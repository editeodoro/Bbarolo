//--------------------------------------------------------------------
// paramguess.cpp: Definitions of functions for the ParamGuess class.
//--------------------------------------------------------------------

/*-----------------------------------------------------------------------
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 BBarolo is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with BBarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning BBarolo may be directed to:
    Internet email: enrico.diteodoro@gmail.com
-----------------------------------------------------------------------*/

#include <iostream>
#include <functional>
#include <cfloat>
#include <Arrays/cube.hh>
#include <Map/detection.hh>
#include <Tasks/galmod.hh>
#include <Tasks/ringmodel.hh>
#include <Tasks/moment.hh>
#include <Tasks/ellprof.hh>
#include <Utilities/paramguess.hh>
#include <Utilities/utils.hh>
#include <Utilities/lsqfit.hh>
#include <Utilities/gnuplot.hh>
#include <Utilities/optimization.hh>



template <class T>
ParamGuess<T>::ParamGuess(Cube<T> *c, Detection<T> *object) {
    
    in  = c;
    obj = object;
    if (!obj->hasParams())
        obj->calcAllParams(in->Array(),in->AxisDim(),in->Head());

    /// Extracting intensity and velocity field
    Vemap = new T[in->DimX()*in->DimY()];
    Intmap = new T[in->DimX()*in->DimY()];
    std::vector<Voxel<T> > voxelList = obj->getPixelSet(in->Array(), in->AxisDim());
    float *fluxint = new float[in->DimX()*in->DimY()];
    float *fluxsum = new float[in->DimX()*in->DimY()];
    for (int i=0; i<in->DimX()*in->DimY();i++) fluxint[i] = fluxsum[i] = Intmap[i]= 0;
    
    for(auto &vox : voxelList) {
        int x = vox.getX();
        int y = vox.getY();
        int z = vox.getZ();
        float flux = in->Array(x,y,z);
        fluxsum[x+y*in->DimX()] += flux;
        fluxint[x+y*in->DimX()] += flux*in->getZphys(z);
        Intmap[x+y*in->DimX()] += flux;
    }
    
    totflux_obs=0;
    for (int i=0; i<in->DimX()*in->DimY();i++) {
        totflux_obs += fluxsum[i];
        Vemap[i] = AlltoVel(fluxint[i]/fluxsum[i], in->Head());
    }
    
    // Initializing radsep to the Beam size
    radsep = in->Head().Bmaj()*arcsconv(in->Head().Cunit(0));
    
    delete [] fluxint;
    delete [] fluxsum;
} 


template <class T>
void ParamGuess<T>::findAll() {
    
    // Front-end function to estimate all geometrical and kinematical
    // parameters needed by BBarolo's 3DFIT task

    // Estimating systemic velocity
    findSystemicVelocity();
    // Estimating centre position
    findCentre();
    // Estimating position angle
    findPositionAngle();
    // Estimating inclination angle
    findInclination();
    // Estimating number of rings and width
    findRings();
    // Estimating rotation velocity
    findRotationVelocity();
}


template <class T>
void ParamGuess<T>::findCentre() {
    
    /// X-Y centres are estimated from the centroids of the
    /// object detected by the source-finding algorithm.
    xcentre = (obj->getXcentre()+obj->getXaverage())/2.;
    ycentre = (obj->getYcentre()+obj->getYaverage())/2.;
}


template <class T>
void ParamGuess<T>::findSystemicVelocity() {
    
    /// Systemic velocity is estimated from the total spectrum of
    /// the object detected by the source-finding algorithm.
    vsystem = obj->getVsys();
}


template <class T>
void ParamGuess<T>::findRotationVelocity() {
    
    /// Rotation velocity is estimated from the W50 of spectrum of
    /// the object detected by the source-finding algorithm.
    vrot=fabs(obj->getW50()/2.);
    if (inclin>40) vrot /= sin(inclin*M_PI/180.);
}


template <class T>
void ParamGuess<T>::findPositionAngle(int algorithm) {
    
    ////////////////////////////////////////////////////////////////////////////
    /// Estimating Position angle using several different algorithms
    /// algorithm = 1: velocity differences from the VSYS
    /// algorithm = 2: regions of highest/lowest velocities
    ///
    /// NB: xcenter, ycenter, vsystem need to be set before calling this function
    /// This function sets posang
    ////////////////////////////////////////////////////////////////////////////
    
    // Getting maximum and minimum velocity in the spectral range of the cube
    double velmin = AlltoVel(in->getZphys(0),in->Head());
    double velmax = AlltoVel(in->getZphys(in->DimZ()-1),in->Head());
    if (velmin>velmax) std::swap(velmin,velmax);
    // Getting maximum and minimum coordinates of detection
    int Xmin=obj->getXmin(), Ymin=obj->getYmin();
    int Xmax=obj->getXmax(), Ymax=obj->getYmax();
    
    if (algorithm == 1) {
        // This algorithm calculates the median difference between a velocity 
        // on the  velocity field and the systemic velocity for any possible 
        // position angle. The PA that returns the highest median value is the 
        // best kinematical PA.
        
        double maxdev=0, bestpa=0, p=0, vl=0, vr=0;
        // Loop over PA with step 0.5 degrees
        while (p<180) {
            std::vector<T> vdev;
            double sumleft=0, sumright=0;
            if (p>45 && p<135) {
                // If 45 < PA < 135 it is better to loop over y
                for (int y=Ymin; y<=Ymax; y++) {
                    // Calculating (x,y) coordinates to sample velocity field
                    int x = lround(1/tan(p*M_PI/180.)*(y-ycentre)+xcentre);
                    if (p==90) x = lround(xcentre);
                    long npix = y*in->DimX()+x; 
                    bool isOK = x>=Xmin && x<=Xmax && !isNaN<T>(Vemap[npix]) && 
                                Vemap[npix]>=velmin && Vemap[npix]<=velmax;
                    if (!isOK) continue;
                    // Collecting the absolute difference from the VSYS
                    vdev.push_back(fabs(Vemap[npix]-vsystem));
                    // Getting info on which side we have the highest velocity
                    if (p<=90) {
                        if (y<ycentre) sumleft += Vemap[npix]-vsystem;
                        else sumright += Vemap[npix]-vsystem;
                    }
                    else {
                        if (x<xcentre) sumleft += Vemap[npix]-vsystem;
                        else sumright += Vemap[npix]-vsystem;
                    }
                }
            }
            else {
                // Loop over x in any other PA case
                for (int x=Xmin; x<=Xmax; x++) {
                    // Calculating (x,y) coordinates to sample velocity field
                    int y = lround(tan(p*M_PI/180.)*(x-xcentre)+ycentre);
                    long npix = y*in->DimX()+x;
                    bool isOK = y>=Ymin && y<=Ymax && !isNaN<T>(Vemap[npix]) && 
                                Vemap[npix]>=velmin && Vemap[npix]<=velmax;
                    if (!isOK) continue;
                    // Collecting the absolute difference from the VSYS
                    vdev.push_back(fabs(Vemap[npix]-vsystem));
                    // Getting info on which side we have the highest velocity
                    if (x<xcentre) sumleft += Vemap[npix]-vsystem;
                    else sumright += Vemap[npix]-vsystem;
                }
            }

            // Calculating the median deviation from VSYS
            T median = findMedian<T>(&vdev[0], vdev.size());
            // If the median is so far the highest, assign best pa
            if (median>maxdev && fabs(median)<1E16) {
                maxdev = median;
                bestpa = p;
                vl = sumleft;
                vr = sumright;
            }

            p+=0.5;
        }

        // Rotate the PA to conform to BBarolo's definition.
        if (vl<vr) {
            if (bestpa<90) posang = 270+bestpa;
            else posang = 90+bestpa;
        }
        else {
            if (bestpa<90) posang = 90+bestpa;
            else posang = bestpa-90;
        }
    }

    else if (algorithm == 2) {
        // This algorithm samples the velocity field in rectangular regions
        // with size equal to the beam size. In each region it calculates the
        // median velocity value and finds the two spots where this median value
        // is the highest and the lowest. Then, a linear regression between
        // these two regions is performed and the position angle estimated.

        float vel_high = vsystem;
        float vel_low  = vsystem;   
        int range = ceil(in->Head().Bmaj()/in->Head().PixScale());
        int coord_high[2]={0,0}, coord_low[2]={0,0};
        int xsize = fabs(obj->getXmax()-obj->getXmin())+1;
        int ysize = fabs(obj->getYmax()-obj->getYmin())+1;
        
        for (int y=range; y<ysize-range; y++) {
            for (int x=range; x<xsize-range; x++) {
                long npix = (y+Ymin)*in->DimX()+x+Xmin;
                if (isNaN<T>(Vemap[npix])) continue;
                std::vector<T> vec;
                for (int yi=y-range; yi<=y+range; yi++) 
                    for (int xi=x-range; xi<=x+range; xi++) 
                        vec.push_back(Vemap[(yi+Ymin)*in->DimX()+xi+Xmin]);
                T median = findMedian<T>(&vec[0], vec.size());
                if (median<vel_low && median>=velmin) {
                    vel_low = median;
                    coord_low[0] = x+Xmin;
                    coord_low[1] = y+Ymin;
                }
                if (median>vel_high && median<=velmax) {
                    vel_high = median;
                    coord_high[0] = x+Xmin;
                    coord_high[1] = y+Ymin;
                }
            }
        }

        std::vector<int> xx(3), yy(3);
        xx[0]=coord_low[0]; xx[1]=coord_high[0]; xx[2]=lround(xcentre);
        yy[0]=coord_low[1]; yy[1]=coord_high[1]; yy[2]=lround(ycentre);
        float rmaj, errmaj[2], pmaj[2];
    
        // Linear regression between the center and the two point found previously.
        // For not including the center, just change the last parameter from
        // 2 to 1 in the function below.
        int a = linear_reg<int>(3, xx, yy, pmaj, errmaj, rmaj, 0, 2);
        
        float ang = atan(pmaj[0]);
        if (coord_high[0]>=xcentre) {   
            if (ang<M_PI_2) posang = (3*M_PI_2+ang)*180/M_PI;
            else posang = (M_PI_2+ang)*180/M_PI;
        }
        else {
            if (ang<M_PI_2) posang = (M_PI_2+ang)*180/M_PI;
            else posang = (ang-M_PI_2)*180/M_PI;    
        }

    }
    else {
        std::cerr << "PARAMGUESS ERROR: unknown algorithm value " << algorithm << std::endl;
        std::terminate();
    }
    
}


template <class T>
void ParamGuess<T>::findInclination(int algorithm) {
    
    ////////////////////////////////////////////////////////////////////////////
    /// Estimating Inclination angle using several different algorithms
    /// algorithm = 1: Length of major/minor axis on the velocity field.
    /// algorithm = 2: Ellipse with largest number of valid pixels.
    /// algorithm = 3: Fit a model intensity map to the observed one.
    ///
    /// NB: xcenter, ycenter, vsystem and posang need to be set before calling
    /// this function.
    /// This function sets inclin, Rmax, nrings, radsep
    ////////////////////////////////////////////////////////////////////////////
    
    // We always start with algorithm==1. This will provide also initial 
    // guesses for other algorithms.
    
    // Estimating axis lengths for major and minor axis
    float pmaj[2], pmin[2];
    setAxesLine(xcentre,ycentre,posang,pmaj,pmin);
    T axmaj = findAxisLength(pmaj, major_max, major_min);
    T axmin = findAxisLength(pmin, minor_max, minor_min);

    if (axmin>axmaj) {
        std::cout << "---------------> WARNING - Finding initial parameters <--------------\n"
                  << " The major axis is shorter than the minor axis. They will be swapped\n"
                  << " for estimating the inclination.\n"
                  << " The galaxy seems to be less elongated in the kinematical axis!!\n\n";
        std::swap(axmin, axmaj);
    }
    
    // Inclination angle and maximum radius
    inclin = acos(axmin/axmaj)*180./M_PI;
    Rmax = axmaj*in->Head().PixScale()*arcsconv(in->Head().Cunit(0));
    
    if (algorithm==1) {
        // We are happy with what we have done above. Just return
        return;
    }
    else if (algorithm==2 || algorithm==3) {
        // Algorithm 2 finds the ellipses that contains the largest number of
        // non-NaN pixels on the velocity field.
        //
        // Algorithm 3 minimizes the difference between a model intensity map
        // and the observed intensity map.
        
        int ndim=2;
    
        std::vector<double> point(ndim), dels(ndim);
        point[0] = Rmax;
        point[1] = inclin;
        //point[2] = posang;
        //point[3] = xcentre;
        //point[4] = ycentre;
    
        // Determine the initial simplex.
        for (int i=0; i<ndim; i++) {
            dels[i]  = 0.1*point[i]; 
            point[i] = point[i]-0.05*point[i];
        }

        // Deciding what function to pass to the optimizer, depending on chosen algorithm.
        auto func2 = std::bind(funcEllipse<T>, std::placeholders::_1, in, posang, xcentre, ycentre, Vemap);
        auto func3 = std::bind(funcIncfromMap<T>, std::placeholders::_1, in, radsep, Rmax, posang, xcentre, ycentre, vsystem, Intmap, totflux_obs);

        NelderMead optimizer;

        double *mymin;
        try {
            if (algorithm==3) mymin = optimizer.minimize(point,dels,func3);
            else mymin = optimizer.minimize(point,dels,func2);
        }
        catch (...) {
            std::cerr << "PARAMGUESS ERROR: Error while estimating inclination." << std::endl;
            std::terminate();
        }

        Rmax   = mymin[0];
        inclin = mymin[1];

        // I am adding a correction for rounding due to the beam
        double Axmin = Rmax*cos(inclin/180.*M_PI);
        double Axmin_corr = Axmin - in->Head().Bmaj()*3600.;
        double Axmaj_corr = Rmax - in->Head().Bmaj()*3600.;
        if (inclin<75) {
            bool isok = Axmin_corr/(in->Head().Bmaj()*3600.) > 1;
            if (Axmin_corr>0 && Axmaj_corr>0 && isok)
                inclin = acos(Axmin_corr/Axmaj_corr)*180/M_PI;
        }
        else {
            bool isok = Axmin_corr/(in->Head().Bmaj()*3600.) > 1 &&  Axmaj_corr/(in->Head().Bmaj()*3600.) > 9;
            if (isok) inclin = acos(Axmin_corr/Axmaj_corr)*180/M_PI;
        }

        if (2*Axmin<in->Head().Bmaj()*3600.) inclin = 80;
        //posang = mymin[2];
        //xcentre= mymin[3];
        //ycentre= mymin[4];
    }
    else {
        std::cerr << "PARAMGUESS ERROR: unknown algorithm value " << algorithm << std::endl;
        std::terminate();
    }
    
    if (inclin<10) inclin=25;

}


template <class T>
void ParamGuess<T>::findRings() {

    nrings = lround(Rmax/radsep);
    if (nrings<20) {
        radsep /= 2.;
        nrings = lround(Rmax/radsep);
    }
    if (nrings>4) nrings -= 1;
}


template <class T>
T ParamGuess<T>::findAxisLength(float *lpar, int *coords_up, int *coords_low) {
    
    // Reset coordinates
    coords_up[0] = coords_up[1] = coords_low[0] = coords_low[1] = 0;
    
    double axis_r_l=0., axis_r_r=0;
    int Xmin=obj->getXmin(), Ymin=obj->getYmin();
    int Xmax=obj->getXmax(), Ymax=obj->getYmax();

    double p = atan(lpar[0])*180./M_PI;
    while (p>=180) p-=180;
    while (p<0) p+=180;
    if (p>45 && p<135) {
        // If 45 < angle < 135 it is better to loop over y
        for (int y=Ymin; y<=Ymax; y++) {
            // Calculating (x,y) coordinates to sample velocity field
            int x = lround((y-lpar[1])/lpar[0]);
            if (fabs(90.-p)<0.1) x = lround(xcentre);
            long npix = y*in->DimX()+x;
            bool isOK = x>=Xmin && x<=Xmax && !isNaN<T>(Vemap[npix]);
            if (!isOK) continue;
            double r = sqrt(pow(double(x-xcentre),2.)+pow(double(y-ycentre),2.));
            if (p-90<0.5) {    // Special case for PA=90
                if (r>axis_r_l && y<=ycentre) {
                    axis_r_l = r;
                    coords_up[0] = x; coords_up[1] = y;
                }
                if (r>axis_r_r && y>ycentre) {
                    axis_r_r = r;
                    coords_low[0] = x; coords_low[1] = y;
                }
            }
            else {
                if (r>axis_r_l && x<=xcentre) {
                    axis_r_l = r;
                    coords_up[0] = x; coords_up[1] = y;
                }
                if (r>axis_r_r && x>xcentre) {
                    axis_r_r = r;
                    coords_low[0] = x; coords_low[1] = y;
                }
            }
        }
    }
    else {
        // Loop over x in any other angle case
        for (int x=Xmin; x<=Xmax; x++) {
            int y = lround(lpar[0]*x+lpar[1]);
            long npix = y*in->DimX()+x; 
            bool isOK = y>=Ymin && y<=Ymax && !isNaN<T>(Vemap[npix]);
            if (!isOK) continue;
            double r = sqrt(pow(double(x-xcentre),2.)+pow(double(y-ycentre),2.));
            if (r>axis_r_l && x<=xcentre) {
                axis_r_l = r;
                coords_up[0] = x; coords_up[1] = y;
            }
            if (r>axis_r_r && x>xcentre) {
                axis_r_r = r;
                coords_low[0] = x; coords_low[1] = y;
            }
        }
    }
        
    return (axis_r_r+axis_r_l)/2.;
    
}


template <class T>
void ParamGuess<T>::setAxesLine(T xcen, T ycen, T pa, float *maj, float *min) {

    /// Calculate angular coefficients and zero-points for major and minor axis.
    float m = pa - 90;
    while (m>=180) m-=180;
    while (m<0) m+=180;
    maj[0] = tan(m*M_PI/180.);
    maj[1] = ycen-maj[0]*xcen;
    min[0] = - 1/maj[0];
    min[1] = ycen-min[0]*xcen;
}


template <class T>
void ParamGuess<T>::tuneWithTiltedRing() {

    /// This function uses a 2D tilted-ring model to get better estimates of
    /// the parameters, in particular the centre and the position angle.
    ///
    /// N.B.: it must be called after all parameters are already set to some
    /// value

    // The tilted ring model will extend to Rmax - 2 pixels
    T rmax = Rmax/(in->Head().PixScale()*arcsconv(in->Head().Cunit(0))) - 2;
    // It uses a fixed ring width of 2 pixels
    T rwidth = 1;
    int nr = rmax/rwidth;
    if (nr<4) return;

    // Initializing rings
    T *radii = new T[nr];
    for (int i=0; i<nr; i++) radii[i] = 2+i*rwidth;

    // Initializing a Ringmodel instance
    Ringmodel<T> tr(nr,radii,rwidth,vsystem,vrot,0,posang,inclin,xcentre,ycentre);
    tr.setfield(Vemap,in->DimX(),in->DimY());
    // Setting free parameters. Order is VSYS, VROT, VEXP, PA, INC, X0, Y0
    // Here I keeping Vsys and inc fixed and fit for the center and pa and vrot
    bool mpar[7] = {true,true,true,true,true,true,true};
    mpar[VSYS] = true;
    mpar[VEXP] = false;
    mpar[INC]  = false;

    tr.setoption(mpar,3,2,15.);
    // Fitting a tilted-ring model
    tr.ringfit(in->pars().getThreads(),false,false);
    //tr.printfinal(std::cout,in->Head());
    std::ofstream fileo(in->pars().getOutfolder()+in->Head().Name()+"_2drings.txt");
    tr.printfinal(fileo,in->Head());

    //tr.set(nr,radii,rwidth,tr.getVsysf(0),tr.getVrotf(0),0,tr.getPosaf(0),tr.getInclf(0),tr.getXposf(0),tr.getYposf(0));
    //mpar[PA] = false;
    //tr.setoption(mpar,3,2,15.);
    //tr.ringfit(in->pars().getThreads(),false,false);
    //tr.printfinal(std::cout,in->Head());

    std::vector<T> xcen,ycen,vsys,posa,incl;
    for (int i=1; i<nr; i++) {  // Skipping first ring
        if (!isNaN(tr.getXposf(i))) xcen.push_back(tr.getXposf(i));
        if (!isNaN(tr.getYposf(i))) ycen.push_back(tr.getYposf(i));
        if (!isNaN(tr.getVsysf(i))) vsys.push_back(tr.getVsysf(i));
        if (!isNaN(tr.getInclf(i))) incl.push_back(tr.getInclf(i));
        if (!isNaN(tr.getPosaf(i))) posa.push_back(tr.getPosaf(i));
    }

    // Setting initial estimates to the median of the tilted-ring model
    if (xcen.size()>1) xcentre = findMedian(&xcen[0],xcen.size());
    if (ycen.size()>1) ycentre = findMedian(&ycen[0],ycen.size());
    if (vsys.size()>1) vsystem = findMedian(&vsys[0],vsys.size());
    if (incl.size()>1) inclin  = findMedian(&incl[0],incl.size());
    if (posa.size()>1) posang  = findMedian(&posa[0],posa.size());

    // Re-estimating inclination angle with the latest parameters
    findInclination(2);
    findRotationVelocity();
    delete [] radii;
}


template <class T>
int ParamGuess<T>::plotGuess(std::string outfile) {
    
    /// Plotting the initial guesses

    std::string outfolder = in->pars().getOutfolder();
    std::vector<T> vec;
    T axmaj_pix = Rmax/(in->Head().PixScale()*arcsconv(in->Head().Cunit(0)));
    T axmin_pix = axmaj_pix*cos(inclin*M_PI/180.);
    float maj[2], min[2];
    int ret = 0;

    // Writing intensity map and velocity field in a text file
    std::ofstream velf((outfolder+"vfield.dat").c_str());
    std::ofstream intf((outfolder+"ifield.dat").c_str());
    short b = 10;
    int xstart = obj->getXmin()-b>=0 ? obj->getXmin()-b : 0;
    int xstop  = obj->getXmax()+b<in->DimX() ? obj->getXmax()+b : in->DimX()-1;
    int ystart = obj->getYmin()-b>=0 ? obj->getYmin()-b : 0;
    int ystop  = obj->getYmax()+b<in->DimY() ? obj->getYmax()+b : in->DimY()-1;

    T *prof = new T[in->DimZ()];
    for (int z=0; z<in->DimZ(); z++) prof[z] = 0;

    for (int x=xstart; x<=xstop; x++) {
        for (int y=ystart; y<=ystop; y++) {
            long npix = x+y*in->DimX(); 
            if (!isNaN(Vemap[npix])) vec.push_back(Vemap[npix]);
            velf << x << "  " << y << "  " << Vemap[npix] << std::endl;
            intf << x << "  " << y << "  " << Intmap[npix] << std::endl;
            for (int z=0; z<in->DimZ(); z++) prof[z] += in->Array(x,y,z);
        }
        velf << std::endl;
    }
    velf.close();
    intf.close();

    // Writing total spectrum to a file
    std::ofstream spf((outfolder+"spec.dat").c_str());
    for (int z=0; z<in->DimZ(); z++)
        spf << AlltoVel(in->getZphys(z),in->Head()) << "  " << prof[z] << std::endl;
    delete [] prof;
    spf.close();

    float maxvel = *max_element(&vec[0], &vec[0]+vec.size());
    float minvel = *min_element(&vec[0], &vec[0]+vec.size());

    double pix[3] = {xcentre, ycentre, 0}, world[3];
    if (pixToWCSSingle(in->Head().WCS(),pix,world)) world[0]=world[1]=0;


#ifdef HAVE_PYTHON
    setAxesLine(xcentre-xstart,ycentre-ystart,posang,maj,min);

    std::ofstream pyf((outfolder+"pyscript_ig.py").c_str());
    pyf << "import numpy as np \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "fsize = 11 \n"
        << "mpl.rc('xtick',direction='in',labelsize=fsize) \n"
        << "mpl.rc('ytick',direction='in',labelsize=fsize) \n"
        << "plt.rc('font',family='sans-serif',serif='Helvetica',size=fsize) \n"
        << "params = {'text.usetex': False, 'mathtext.fontset': 'cm', 'mathtext.default': 'regular'} \n"
        << "plt.rcParams.update(params) \n"
        << std::endl
        << "outdir = '"<< outfolder << "'\n"
        << "xsize, ysize = " << xstop-xstart+1 << " , " << ystop-ystart+1 << "\n"
        << "xpos, ypos, pa, inc = " << xcentre-xstart << " , " << ycentre-ystart
        << " , " << posang << " , " << inclin << "\n"
        << "axmaj, axmin, vsys = " << axmaj_pix << " , " << axmin_pix << " , " << vsystem << "\n"
        << "f = np.genfromtxt('%s/ifield.dat'%outdir,usecols=(2),unpack=True) \n"
        << "v = np.genfromtxt('%s/vfield.dat'%outdir,usecols=(2),unpack=True) \n"
        << "v = v.reshape(xsize,ysize).T \n"
        << "f = f.reshape(xsize,ysize).T \n"
        << "f[f==0] = np.nan \n"
        << "t  = np.linspace(0,2*np.pi,100) \n"
        << "xt = xpos+axmaj*np.cos(np.radians(pa-90))*np.cos(t)-axmin*np.sin(np.radians(pa-90))*np.sin(t) \n"
        << "yt = ypos+axmaj*np.sin(np.radians(pa-90))*np.cos(t)+axmin*np.cos(np.radians(pa-90))*np.sin(t) \n"
        << "fx = lambda x: " << maj[0] << "*x+" << to_string(maj[1]) << std::endl
        << "gx = lambda x: " << min[0] << "*x+" << to_string(min[1]) << std::endl
        << "x_maj, y_maj = ["<< major_min[0]-xstart << "," << major_max[0]-xstart << "], ["
        << major_min[1]-ystart << "," << major_max[1]-ystart << "] \n"
        << "x_min, y_min = ["<< minor_min[0]-xstart << "," << minor_max[0]-xstart << "], ["
        << minor_min[1]-ystart << "," << minor_max[1]-ystart << "] \n"
        << "vel, sp = np.genfromtxt('%s/spec.dat'%outdir,usecols=(0,1),unpack=True) \n"
        << "v20min, v20max = " << obj->getV20Min() << " , " << obj->getV20Max() << "\n"
        << "v50min, v50max = " << obj->getV50Min() << " , " << obj->getV50Max() << "\n"
        << std::endl
        << "fig = plt.figure(figsize=(10,10)) \n"
        << "fig.add_axes([0.00,0.42,0.3,0.3]) \n"
        << "fig.add_axes([0.31,0.42,0.3,0.3]) \n"
        << "fig.add_axes([0.00,0.1,0.4,0.3]) \n"
        << "ax = fig.axes \n"
        << "cax1 = fig.add_axes([ax[0].get_position().x0,ax[0].get_position().y1+0.01,0.3,0.02]) \n"
        << "cax2 = fig.add_axes([ax[1].get_position().x0,ax[1].get_position().y1+0.01,0.3,0.02]) \n"
        << std::endl
        << "im = ax[0].imshow(f,origin='lower',cmap=plt.get_cmap('Spectral_r'),aspect='auto') \n"
        << "cb = plt.colorbar(im,cax=cax1,orientation='horizontal') \n"
        << "cb.set_label('Int flux (std)',labelpad=-50,fontsize=fsize+1) \n"
        << std::endl
        << "im = ax[1].imshow(v,origin='lower',cmap=plt.get_cmap('RdBu_r',25),aspect='auto') \n"
        << "ax[1].contour(v,levels=[vsys],colors='green') \n"
        << "cb = plt.colorbar(im,cax=cax2 ,orientation='horizontal') \n"
        << "cb.set_label(r'V$_\\mathrm{LOS}$ (km/s)',labelpad=-50,fontsize=fsize+1) \n"
        << std::endl
        << "ax[2].plot(vel,sp,'.-') \n"
        << "ax[2].axhline(0,c='k',lw=0.8) \n"
        << "ax[2].axvline(vsys,c='green') \n"
        << "ax[2].axvline(v20min,c='#9932CC',lw=1,ls='--',label='V20') \n"
        << "ax[2].axvline(v20max,c='#9932CC',lw=1,ls='--') \n"
        << "ax[2].axvline(v50min,c='#FF8C00',lw=1,ls='--',label='V50') \n"
        << "ax[2].axvline(v50max,c='#FF8C00',lw=1,ls='--') \n"
        << "ax[2].set_xlabel(r'V$_\\mathrm{LOS}$ (km/s)',fontsize=fsize+1) \n"
        << "ax[2].set_ylabel(r'Flux (std)',fontsize=fsize+1) \n"
        << "ax[2].tick_params(right=True,top=True) \n"
        << std::endl
        << "for a in ax[:2]: \n"
        << "\ta.plot([0,xsize],fx(np.array([0,xsize])),'-',c='#4169E1') \n"
        << "\ta.plot([0,xsize],gx(np.array([0,xsize])),'--',c='#FF7F50') \n"
        << "\ta.plot(xt,yt,'k') \n"
        << "\ta.plot(xpos,ypos,'x',color='#000000',markersize=7,mew=1.5) \n"
        << "\ta.plot(x_maj,y_maj,'o',color='red',markersize=4) \n"
        << "\ta.plot(x_min,y_min,'o',color='red',markersize=4) \n"
        << "\ta.tick_params(right=True,top=True,labelbottom=False,labelleft=False) \n"
        << "\ta.set_xlim(0,xsize) \n"
        << "\ta.set_ylim(0,ysize) \n"
        << std::endl
        << "for a in [cax1,cax2]: \n"
        << "\ta.xaxis.set_ticks_position('top') \n"
        << "\ta.locator_params(nbins=4) \n"
        << std::endl
        << "xpos, ypos = " << xcentre << " , " << ycentre << "\n"
        << "ra, dec    = " << world[0] << " , " << world[1] << "\n"
        << "pstr = [r'X$_\\mathrm{pos}$ = %.2f pix'%xpos, r'Y$_\\mathrm{pos}$ = %.2f pix'%ypos, \n "
        << "        r'$\\alpha$  = %.5f$^\\circ$'%ra, r'$\\delta$ = %.5f$^\\circ$'%dec, \n"
        << "        r'$\\phi$    = %.0f$^\\circ$'%pa,r'inc  = %.1f$^\\circ$'%inc, r'$b/a$ = %.2f'%(axmin/axmaj), \n"
        << "        r'R$_\\mathrm{max}$ = %.2f pix'%axmaj, r'V$_\\mathrm{sys}$ = %.2f km/s'%vsys, \n"
           "        r'W$_{20}$ = %.2f km/s'%(np.fabs(v20max-v20min)),r'W$_{50}$ = %.2f km/s'%(np.fabs(v50max-v50min))] \n"
        << "for i in range (len(pstr)): \n"
        << "\tax[2].text(1.1,0.97-i*0.08,pstr[i],va='top',transform=ax[2].transAxes) \n"
        << std::endl
        << "fig.savefig('%s/"<< outfile << "'%outdir,bbox_inches='tight') \n"
        << std::endl;

    pyf.close();

    std::string cmd = "python \'"+outfolder+"pyscript_ig.py\'";// > /dev/null 2>&1";
    ret = system(cmd.c_str());
    remove((outfolder+"pyscript_ig.py").c_str());
#else
#ifdef HAVE_GNUPLOT
    std::string outf = outfolder+"initial_guesses.eps";
    std::ofstream gnu((outfolder+"gnuscript.gnu").c_str());
    std::string amaj = to_string(axmaj_pix);
    std::string amin = to_string(axmin_pix);
    std::string posa = to_string(posang/180*M_PI-M_PI_2);
    std::string xcen = to_string(xcentre);
    std::string ycen = to_string(ycentre);
    setAxesLine(xcentre,ycentre,posang,maj,min);
    gnu << "unset key\n"
    //  << "set grid\n"
        << "set title 'Axis fitting'\n"
        << "set cbtics scale 0\n"
        << "set palette defined (0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',"
        << " 5 '#ffee00', 6 '#ff7000',7 '#ee0000',8 '#7f0000')\n"
        << "f(x)="<<to_string(maj[0])<<"*x+"<<to_string(maj[1])<<std::endl
        << "g(x)="<<to_string(min[0])<<"*x+"<<to_string(min[1])<<std::endl
        << "set xrange [0:"<<to_string<long>(in->DimX())<<"]\n"
        << "set yrange [0:"<<to_string<long>(in->DimY())<<"]\n"
        << "set cbrange ["<<to_string(minvel)<<":"<<to_string(maxvel)<<"]\n"
        << "set xlabel 'X (pixels)'\n"
        << "set ylabel 'Y (pixels)'\n"
        << "set size square\n"
        << "set parametric\n"
        << "x(t)="+xcen+"+"+amaj+"*cos("+posa+")*cos(t)-"+amin+"*sin("+posa+")*sin(t)\n"
        << "y(t)="+ycen+"+"+amaj+"*sin("+posa+")*cos(t)+"+amin+"*cos("+posa+")*sin(t)\n"
        << "set table '"+outfolder+"ellipse.tab'\n"
        << "plot x(t), y(t)\n"
        << "unset table\n"
        << "unset parametric\n"
        << "set terminal postscript eps enhanced color font 'Helvetica,14'\n"
        << "set output '"<< outf <<"'\n"
        << "plot '"+outfolder+"vfield.dat' w image pixels,"
        << " '"+outfolder+"ellipse.tab' w l ls -1 lw 2, f(x) ls 1 lw 2, g(x) ls 3 lw 2,'-' ls 5, '-' ls 7 \n"
        << xcen+" "+ycen << std::endl
        << "e" << std::endl
        << to_string(major_max[0])+" "+to_string(major_max[1]) << std::endl
        << to_string(major_min[0])+" "+to_string(major_min[1]) << std::endl
        << to_string(minor_max[0])+" "+to_string(minor_max[1]) << std::endl
        << to_string(minor_min[0])+" "+to_string(minor_min[1]) << std::endl
        << "e" << std::endl;

    gnu.close();

    Gnuplot gp;
    gp.begin(); 
    gp.commandln(("load '"+outfolder+"gnuscript.gnu'").c_str());
    gp.end();
    remove((outfolder+"ellipse.tab").c_str());
    remove((outfolder+"gnuscript.gnu").c_str());
#endif
#endif

    remove((outfolder+"vfield.dat").c_str());
    remove((outfolder+"ifield.dat").c_str());
    remove((outfolder+"spec.dat").c_str());

    return ret;
}


template <class T>
double funcEllipse(std::vector<double> &mypar, Cube<T> *c, double pa, double xcen, double ycen, T* Vf) {

    if (mypar[0]<0)  mypar[0]= fabs(mypar[0]);
    if (mypar[1]>90) mypar[1]= 180 - mypar[1];

    T R   = mypar[0]/(c->Head().PixScale()*arcsconv(c->Head().Cunit(0)));
    T inc = mypar[1]*M_PI/180.;
    pa   *= M_PI/180.;//mypar[2];
    //T x0  = mypar[3];
    //T y0  = mypar[4];

    double func = 0;
    for (int x=0; x<c->DimX(); x++) {
        for (int y=0; y<c->DimY(); y++) {
            T xr =  -(x-xcen)*sin(pa)+(y-ycen)*cos(pa);
            T yr = (-(x-xcen)*cos(pa)-(y-ycen)*sin(pa))/cos(inc);
            T r = sqrt(xr*xr+yr*yr);
            bool isIn = r<=R;
            if (!isIn) continue;
            if (isNaN(Vf[x+y*c->DimX()])) func++;
            else func--;
        }
    }
    return func;
}


template <class T>
double funcIncfromMap(std::vector<double> &mypar, Cube<T> *c, double radsep, double Rmax,
                      double pa, double xcen, double ycen, double vsys, T* Imap, double totflux) {


    bool verbosity = c->pars().isVerbose();
    c->pars().setVerbosity(false);

    if (mypar[0]<0) mypar[0]=2*radsep;
    if (mypar[0]>1.5*Rmax) mypar[0]=Rmax;
    if (mypar[1]<0) mypar[1]=1;
    if (mypar[1]>90) mypar[1]=89;

    T RMAX = mypar[0];
    T inc  = mypar[1];

    Rings<T> *rings = new Rings<T>;
    rings->setRings(0,RMAX,radsep,xcen,ycen,vsys,10*DeltaVel(c->Head()),8,0,0,0,0,1E20,0,inc,pa);

    MomentMap<T> *totalmap = new MomentMap<T>;
    totalmap->input(c);
    totalmap->SumMap(true);
    for (int i=0; i<totalmap->NumPix();i++) totalmap->Array()[i] = Imap[i];
    Tasks::Ellprof<T> ell(totalmap,rings);
    ell.RadialProfile();
    delete totalmap;

    double profmin=FLT_MAX;
    for (int i=0; i<rings->nr; i++) {
        double mean = ell.getMean(i);
        if (!isNaN(mean) && profmin>mean && mean>0) profmin = mean;
    }
    float factor = 1;
    while(profmin<0.1) {
        profmin*=10;
        factor *=10;
    }
    while (profmin>10) {
        profmin /= 10;
        factor /= 10;
    }
    for (int i=0; i<rings->nr; i++) {
        rings->dens[i]=factor*fabs(ell.getMean(i))*1E20;
        if (rings->dens[i]==0) rings->dens[i]=profmin*1E20;
    }

    Model::Galmod<T> *mod = new Model::Galmod<T>;
    mod->input(c,rings);
    mod->calculate();
    mod->smooth();

    delete rings;

    T *map_mod = new T[c->DimX()*c->DimY()];

    float totflux_mod=0;
    for (int x=0; x<c->DimX(); x++){
        for (int y=0; y<c->DimY(); y++){
            map_mod[x+y*c->DimX()]=0;
            for (int z=0; z<c->DimZ(); z++)
                map_mod[x+y*c->DimX()]+=mod->Out()->Array(x,y,z);
            totflux_mod += map_mod[x+y*c->DimX()];
        }
    }

    factor = totflux/totflux_mod;

    float res_sum=0;
    for (int i=0; i<c->DimX()*c->DimY();i++) {
        res_sum += fabs(Imap[i]-map_mod[i]*factor);
    }

    delete [] map_mod;

    c->pars().setVerbosity(verbosity);

    return res_sum;
}


template <class T>
ParamGuess<T>* EstimateInitial(Cube<T> *c, GALFIT_PAR *p){
    
    // Running the source finder to detect the source
    if (!c->getIsSearched()) c->search();
    Detection<T> *largest = c->getSources()->LargestDetection();

    bool verb = c->pars().isVerbose();
    c->pars().setVerbosity(false);
    if (verb) std::cout << "\n Estimating initial parameters... " << std::flush;

    if (largest==NULL) {
        std::cout << " 3DFIT ERROR: No sources detected in the datacube. Cannot fit!!! \n";
        std::terminate();
    }

    ParamGuess<T> *ip = new ParamGuess<T>(c,largest);

    // Estimating systemic velocity if not given
    if (p->VSYS!="-1") ip->vsystem = atof(p->VSYS.c_str());
    else ip->findSystemicVelocity();

    // Estimating centre if not given
    string pos[2] = {p->XPOS, p->YPOS};
    double *pixs = getCenterCoordinates(pos, c->Head());
    if (p->XPOS!="-1" && p->YPOS!="-1") {
        ip->xcentre = pixs[0];
        ip->ycentre = pixs[1];
    }
    else ip->findCentre();
    
    // Estimating position angle if not given
    if (p->PHI!="-1") ip->posang = atof(p->PHI.c_str());
    else ip->findPositionAngle(1);

    // Estimating rotation velocity angle if not given
    // In findInclination: 1=axis ratio, 2=ellipse, 3=totalmap
    ip->findInclination(2);
    if (p->INC!="-1") ip->inclin = atof(p->INC.c_str());

    // Estimating rings if not given
    if (p->NRADII!=-1 && p->RADSEP!=-1) {
        ip->radsep = p->RADSEP;
        ip->nrings = p->NRADII;
    }
    else ip->findRings();

    if (p->VROT!="-1") ip->vrot = atof(p->VROT.c_str());
    else ip->findRotationVelocity();

    // This performs an additional step with a 2D tilted ring model
    if (c->pars().getFlagPlots()>=3) ip->plotGuess("initialguesses_"+c->Head().Name()+"_0.pdf");
    if (ip->nrings>3) ip->tuneWithTiltedRing();
    if (c->pars().getFlagPlots()>=2) ip->plotGuess("initialguesses_"+c->Head().Name()+".pdf");

    if (verb) std::cout << "Done." << std::endl;
    c->pars().setVerbosity(verb);
    return ip;
}
template ParamGuess<float>* EstimateInitial(Cube<float> *, GALFIT_PAR *);
template ParamGuess<double>* EstimateInitial(Cube<double> *, GALFIT_PAR *);


// Explicit instantiation of the class
template class ParamGuess<float>;
template class ParamGuess<double>;
