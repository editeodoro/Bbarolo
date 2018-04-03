//--------------------------------------------------------------------
// paramguess.hh: A class to estimate initial parameters for 3D Fit
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
#include <Arrays/cube.hh>
#include <Map/detection.hh>
#include <Tasks/galmod.hh>
#include <Utilities/utils.hh>
#include <Utilities/lsqfit.hh>
#include <Utilities/gnuplot.hh>
#include <Tasks/moment.hh>
#include <Tasks/ellprof.hh>

template <class T>
class ParamGuess 
{
public:
    Cube<T> *in;
    Detection<T> *obj;
    T       xcentre;
    T       ycentre;
    T       vsystem;
    T       inclin;
    T       posang;
    T       Rmax;
    T       vrot;
    int     nrings;
    T       radsep;
    T       axmin;
    T       axmaj;
    T*      Intmap;
    T*      Vemap;
    
    ParamGuess(Cube<T> *c, Detection<T> *object);
    ~ParamGuess();
    
    void setPosang (T v) {posang=v;}
    void setXcentre (T v) {xcentre=v;}
    void setYcentre (T v) {ycentre=v;}


    void findInitial();
    void fitEllipse();
    void fitIncfromMap();
    void plotGuess();
    
private:
    float   pmaj[2];            //Parameter m of linear function y=mx+q
    float   pmin[2];            //Parameter q of linear function y=mx+q
    int     coord_high[2];
    int     coord_low[2];
    int     range;
    int     xsize;
    int     ysize;
    int     zsize;
    int     Xmin;
    int     Ymin;
    int     major_max[2];
    int     minor_max[2];
    int     major_min[2];
    int     minor_min[2];
    float   totflux_obs;

    T     funcEllipse(T *mypar);
    T     funcIncfromMap(T *mypar);
    
    void  findGeometricalParameters();
    bool  fitSimplex(int ndim, T **p);  
    
    typedef T (ParamGuess<T>::*funcPtr) (T *);
    funcPtr func;
};


template <class T>
ParamGuess<T>::ParamGuess(Cube<T> *c, Detection<T> *object) {
    
    in  = c;
    obj = object;
    func = &ParamGuess<T>::funcEllipse;
    
    
    obj->calcFluxes(obj->getPixelSet(in->Array(), in->AxisDim()));
    obj->calcWCSparams(in->Head());
    obj->calcIntegFlux(in->DimZ(), obj->getPixelSet(in->Array(), in->AxisDim()), in->Head());
    
    xsize = fabs(obj->getXmax()-obj->getXmin())+1;
    ysize = fabs(obj->getYmax()-obj->getYmin())+1;
    zsize = fabs(obj->getZmax()-obj->getZmin())+1;
    Xmin = obj->getXmin();
    Ymin = obj->getYmin();
    
        
    /// Extracting intensity and velocity field
    Vemap = new T[in->DimX()*in->DimY()];
    Intmap = new T[in->DimX()*in->DimY()];
    std::vector<Voxel<T> > voxelList = obj->getPixelSet(in->Array(), in->AxisDim());
    float *fluxint = new float[in->DimX()*in->DimY()];
    float *fluxsum = new float[in->DimX()*in->DimY()];
    for (int i=0; i<in->DimX()*in->DimY();i++) fluxint[i] = fluxsum[i] = Intmap[i]= 0;
    
    for(typename std::vector<Voxel<T> >::iterator vox=voxelList.begin();vox<voxelList.end();vox++){
        int x = vox->getX();
        int y = vox->getY();
        int z = vox->getZ();
        float flux = FluxtoJy(in->Array(x,y,z),in->Head());
        fluxsum[x+y*in->DimX()] += flux;
        fluxint[x+y*in->DimX()] += flux*in->getZphys(z);
        Intmap[x+y*in->DimX()] += flux;
    }
    
    totflux_obs=0;
    for (int i=0; i<in->DimX()*in->DimY();i++) {
        totflux_obs += fluxsum[i];
        Vemap[i] = AlltoVel(fluxint[i]/fluxsum[i], in->Head());
    }
    
    delete [] fluxint;
    delete [] fluxsum;
} 

template <class T>
ParamGuess<T>::~ParamGuess() {
    delete [] Vemap;
} 
/*
template <class T>
void ParameterGuess<T>::find_Initial(Detection<T> *obj) {
    
    
    // Initial estimates for fitting parameters

    Rings<T> *init_par = new Rings<T>;
    
    T *array = in->Array(); 
    obj->calcFluxes(obj->getPixelSet(array, in->AxisDim()));
    obj->calcWCSparams(in->Head());
    obj->calcIntegFlux(in->DimZ(), obj->getPixelSet(array, in->AxisDim()), in->Head()); 
    obj->setMass(2.365E5*obj->getIntegFlux()*distance*distance);
    cout << obj->getIntegFlux() << "  " << obj->getMass() << endl << endl;
    
    init_par->xpos.push_back((obj->getXcentre()+obj->getXaverage())/2.);
    init_par->ypos.push_back((obj->getYcentre()+obj->getYaverage())/2.);
    init_par->vsys.push_back(obj->getVsys());
    
    
    int xsize = fabs(obj->getXmax()-obj->getXmin())+1;
    int ysize = fabs(obj->getYmax()-obj->getYmin())+1;
    int zsize = fabs(obj->getZmax()-obj->getZmin())+1;
    
    T *HImap = new T[xsize*ysize];
    T *Vemap = new T[xsize*ysize];
    
    std::vector<bool> isObj(xsize*ysize*zsize,false);
    
    std::vector<Voxel<T> > voxelList = obj->getPixelSet(array, in->AxisDim());
    typename std::vector<Voxel<T> >::iterator vox;
    for(vox=voxelList.begin();vox<voxelList.end();vox++){
        int xm = vox->getX()-obj->getXmin();
        int ym = vox->getY()-obj->getYmin();
        int zm = vox->getZ()-obj->getZmin();
        long pos = xm+ym*xsize+zm*xsize*ysize;
        isObj[pos] = true;
    }
    
    
    for (int x=0; x<xsize; x++) {
        for (int y=0; y<ysize; y++) {
            float fluxint = 0;
            float fluxsum = 0;
            long mappix = x+y*xsize;
            for (int z=0; z<zsize; z++) {
                long cubepix = in->nPix(x+obj->getXmin(), y+obj->getYmin(), z+obj->getZmin());
                HImap[mappix] = 0;
                Vemap[mappix] = 0;
                if (isObj[mappix+z*xsize*ysize]) {
                    float flux = FluxtoJy(array[cubepix],in->Head());
                    fluxsum += flux;
                    float zval = ((z+obj->getZmin()+1-in->Head().Crpix(2))*in->Head().Cdelt(2)+in->Head().Crval(2));
                    fluxint += flux*zval;
                }
            }
            HImap[mappix] = fluxsum*DeltaVel<T>(in->Head());
            Vemap[mappix] = AlltoVel(fluxint/fluxsum, in->Head().Cunit(2), in->Head().Freq0());
        }
    }
    
    
    
    //FitsWrite_2D("DIOCANE1.fits", HImap, xsize, ysize);
    //FitsWrite_2D("DIOCANE2.fits", Vemap, xsize, ysize);

    
    float axmaj=0, axmin=0;
    float posang=0, inclin=0; 
    float pmaj[2] = {0,0}, pmin[2]={0,0};


    float bestmaj = init_par->vsys[0];
    float bestmin = init_par->vsys[0];
    int bestxmaj = 0;
    int bestxmin = 0;
    int bestymaj = 0;
    int bestymin = 0;
    
    for (int x=0; x<xsize; x+=5) {
        for (int y=0; y<ysize; y+=5) {
            if (x<xsize-4 && y<ysize-4) {
                std::vector<float> vec;
                for (int i=x; i<=x+4; i++) {
                    for (int j=y; j<=y+4; j++) {
                        long npix = i+j*xsize;
                        if (isNaN<float>(Vemap[npix])) vec.push_back(0);
                        else vec.push_back(Vemap[npix]);
                    }
                }
                float median = findMedian<float>(vec, vec.size());
                if (median < bestmin && median !=0) {
                    bestmin = median;
                    bestxmin = x+2;
                    bestymin = y+2;
                }
                if (median > bestmaj) {
                    bestmaj = median;
                    bestxmaj = x+2;
                    bestymaj = y+2;
                }
                vec.clear();
            }
        }
    }
    
    
    std::ofstream outmaj, outmin;
    outmaj.open("./output/axmaj.dat");
    outmin.open("./output/axmin.dat");
    for (int x=bestxmin-2; x<=bestxmin+2; x++) {
        for (int y=bestymin-2; y<=bestymin+2; y++) {
            outmaj << x << "   " << y << std::endl;
        }
    }       
    for (int x=bestxmaj-2; x<=bestxmaj+2; x++) {
        for (int y=bestymaj-2; y<=bestymaj+2; y++) {
            outmaj << x << "  " << y << std::endl;
        }
    }       
    outmaj.close();
    
    float vrange = fabs(DeltaVel<float>(in->Head()));

    for (int x=0; x<xsize; x++) {
        for (int y=0; y<ysize; y++) {
            long npix = x+y*xsize;      
            if (Vemap[npix]<= init_par->vsys[0]+vrange && Vemap[npix]>= init_par->vsys[0]-vrange) {
                outmin << x << "   " << y<<std::endl;
            }
        }
    }
    outmin.close();
    
    
    std::vector<int> xx(2);
    xx[0] = bestxmin; xx[1] = bestxmaj;
    std::vector<int> yy(2);
    yy[0] = bestymin; yy[1] = bestymaj;
    float rmaj, errmaj[2];
    int a = linear_reg<int>(2, xx, yy, pmaj, errmaj, rmaj, 0, 1);
    
    if (a==0 && pmaj[0]!=0) {
        pmin[0] = - 1/pmaj[0];
        
        int xlow = 0, ylow = 0;
        int xhigh = xsize-1, yhigh = ysize-1;
        
        float Yxlow = pmaj[0]*xlow + pmaj[1];
        float Yxhigh = pmaj[0]*xhigh + pmaj[1];
        float Xylow = (ylow-pmaj[1])/pmaj[0];
        float Xyhigh = (yhigh-pmaj[1])/pmaj[0];
    
        float x0=0, y0=0, x1=0, y1=0;

        if (Yxlow>=ylow && Yxlow<=yhigh) {
            x0 = xlow;
            y0 = Yxlow;
        }
        else {
            if (Yxlow<ylow && pmaj[0]>=0) {
                x0 = Xylow;
                y0 = ylow;
            }
            else if (Yxlow>yhigh && pmaj[0]<0) {
                x0 = Xyhigh;
                y0 = yhigh;
            }
            else {
                std::cout << "Line is outside the map!!!"<<std::endl;
            }
        }

        if (Yxhigh<=yhigh && Yxhigh>=ylow) {
            x1 = xhigh;
            y1 = Yxhigh;
        }
        else  {
            if (Yxhigh>yhigh && pmaj[0]>=0) {
                x1 = Xyhigh;
                y1 = yhigh;
            }
            else if (Yxhigh<ylow && pmin[0]<0) {
                x1 = Xylow;
                y1 = ylow;
            }
            else {
                std::cout << "Line is outside the map!!!"<<std::endl;
            }
        }
    
        int kmin = y0-pmin[0]*x0;
        int kmax = y1-pmin[0]*x1;
        
        if (kmin>kmax) std::swap(kmin,kmax);

        int best = 0;
        for (int k=kmin; k<kmax; k++) { 
            int count=0;
            for (int x=0; x<xsize; x++) {
                int ymin = pmin[0]*x + k;
                if (ymin>0 && ymin<ysize) {
                    long npix = x+ymin*xsize;
                    if (Vemap[npix]<= init_par->vsys[0]+vrange && Vemap[npix]>= init_par->vsys[0]-vrange) count++;  
                }
            }
            if (count>best) {
                best = count;
                pmin[1] = k;
            }
        }
        
        //float xcenter = (pmin[1]-pmaj[1])/(pmaj[0]-pmin[0]);
        //float ycenter = pmaj[0]*xcenter + pmaj[1];

        int i= lround(init_par->xpos[0]-obj->getXmin());
        int stopmaj=1;
        int stopmin=1;
        float v_mean=0;

        while (i<xsize && (stopmaj>0 || stopmin>0)) {           
            int ymaj = pmaj[0]*i + pmaj[1];
            int ymin = pmin[0]*i + pmin[1];
        
            if (ymaj<ysize && ymaj>0 && stopmaj>0) {
                int npix = i+ymaj*xsize;
                float vmaj = Vemap[npix];
                if (!(vmaj!=vmaj)) {
                    axmaj++;
                    v_mean = (v_mean+vmaj)/2;
                }
                else {
                    int esc=0;
                    for (int k=i+1; k<i+5; k++) {
                        if (k<xsize) {
                            int nymaj = pmaj[0]*k+pmaj[1];
                            if (nymaj<ysize && nymaj>0) {
                                int newpix = k+nymaj*xsize;
                                float nv = Vemap[newpix];
                                if ((nv!=nv)) esc++;
                            }
                            else stopmaj=-1;
                        }
                        else stopmaj=-1;
                    }
                    if (esc>=4) {
                        stopmaj = -1;
                        //std::cout << "STOPMAJ " << i << "  " << ymaj << std::endl;
                    }
                    else axmaj++;
                }
            }
            else stopmaj = -1;
        
            if (ymin<ysize && ymin>0 && stopmin>0) {
                int npix = i+ymin*xsize;
                float vmin = Vemap[npix];
                if (!(vmin!=vmin)) axmin++;
                else {
                    int esc=0;
                    for (int k=i+1; k<i+5; k++) {
                        if (k<xsize) {
                            int nymin = pmin[0]*k+pmin[1];
                            if (nymin<ysize && nymin>0) {
                                int newpix = k+nymin*xsize;
                                float nv = Vemap[newpix];
                                if ((nv!=nv)) esc++;
                            }
                            else stopmin=-1;
                        }
                        else stopmin=-1;
                    }
                    if (esc>=4) {
                        stopmin = -1;
                        //std::cout << "STOPMIN " << i << "  " << ymin << std::endl;
                    }
                    else axmin++;
                }
            }
            else stopmin = -1;
            i++;
        }
        
        if (pmaj[0]>0) {
            float coeffmaj = atan(pmaj[0]);
            axmaj /= cos(coeffmaj);
            axmin /= cos(M_PI_2-coeffmaj);  
        }
        else {
            float coeffmin = atan(pmin[0]);
            axmaj /= cos(M_PI_2-coeffmin);
            axmin /= cos(coeffmin);
            
        }
        
        
        float ang = atan(pmaj[0]);  
    
        if (v_mean>init_par->vsys[0]) { 
            if (ang<M_PI_2) posang = (3*M_PI_2+ang)*180/M_PI;
            else posang = (M_PI_2+ang)*180/M_PI;
        }
        else {
            if (ang<M_PI_2) posang = (M_PI_2+ang)*180/M_PI;
            else posang = (ang-M_PI_2)*180/M_PI;    
        }
        
        if (axmin>axmaj) axmin = axmaj;
        inclin   = in->pars.getParGF().INC!=-1 ? in->pars.getParGF().INC : acos(axmin/axmaj)*180/M_PI;
        
        //double meanw = (fabs(obj->getW20()/2.)+fabs(obj->getW50()/2.))/2.;
        //init_par->vrot.push_back(meanw/sin(inclin*M_PI/180.));
        init_par->vrot.push_back(fabs(obj->getW50()/2.)/sin(inclin*M_PI/180.));
    }
    
    init_par->inc.push_back(inclin);
    init_par->phi.push_back(posang);
    
    double pixScale = ((fabs(in->Head().Cdelt(0))*arcsconv)+
                       (fabs(in->Head().Cdelt(1))*arcsconv))/2.;
    
    float Rmax = axmaj*pixScale;
    init_par->radsep = in->Head().Bmaj()*arcsconv;
    init_par->nr= lround(Rmax/init_par->radsep);
        
    delete [] HImap;
    delete [] Vemap;
    
#ifdef HAVE_GNUPLOT     
    Gnuplot gp;
    
    std::string ouffile = in->pars().getOutfolder()+"axis.eps";
    gp.begin(); 
    gp.commandln("set terminal postscript eps color");
    gp.commandln("unset key");
    gp.commandln("set grid");
    std::string tcmd = "set output '"+ouffile+"'";
    gp.commandln(tcmd.c_str());                   
    gp.commandln("set title 'Axis fitting'");
    
    std::string f = "f(x)="+to_string(pmaj[0])+"*(x-"+to_string(obj->getXmin())+")+"+to_string(pmaj[1]+obj->getYmin());
    std::string g = "g(x)="+to_string(pmin[0])+"*(x-"+to_string(obj->getXmin())+")+"+to_string(pmin[1]+obj->getYmin());
    std::string xr = "set xrange [0:"+to_string<long>(in->DimX())+"]";
    std::string yr = "set yrange [0:"+to_string<long>(in->DimY())+"]";
    std::string xoffset = "($1+"+to_string<int>(obj->getXmin())+")";
    std::string yoffset = "($2+"+to_string<int>(obj->getYmin())+")";
    std::string x_y = "u "+xoffset+":"+yoffset;
    gp.commandln("set xlabel 'X (pixels)'");
    gp.commandln("set ylabel 'Y (pixels)'");
    gp.commandln("set size square");
    gp.commandln(f.c_str());
    gp.commandln(g.c_str());
    gp.commandln(xr.c_str());
    gp.commandln(yr.c_str());
    std::string cmd = "plot './output/axmaj.dat' "+x_y+ ", './output/axmin.dat' "+x_y+ ",f(x), g(x)";
    gp.commandln(cmd.c_str());
    gp.end();
#endif  
    init_par->vdisp.push_back(8);
    init_par->z0.push_back(100);
    showInitial(init_par, std::cout);
    
    remove("./output/axmaj.dat");
    remove("./output/axmin.dat");
    
    return init_par;
    
}
template Rings<float>* Galfit<float>::find_Initial(Detection<float>*);
template Rings<double>* Galfit<double>::find_Initial(Detection<double>*);
*/

template <class T>
void ParamGuess<T>::findInitial() {
    
    // Initial estimates for fitting parameters.
    // This function identify first the major axis and position 
    // angle and then the minor axis and inclination.

    // X-Y centre and systemic velocity are estimated from the 
    // object detected by the source-finding algorithm. 
    xcentre = (obj->getXcentre()+obj->getXaverage())/2.;
    ycentre = (obj->getYcentre()+obj->getYaverage())/2.;
    vsystem = obj->getVsys();
    
    float dist = in->pars().getDistance(); 
    if (dist==-1) dist = VeltoDist(vsystem);
    obj->setMass(2.365E5*obj->getIntegFlux()*dist*dist);
    
    findGeometricalParameters();
    
    radsep = in->Head().Bmaj()*arcsconv(in->Head().Cunit(0));
    nrings = lround(Rmax/radsep);
    vrot=fabs(obj->getW50()/2.)/sin(inclin*M_PI/180.);  
        
}


template <class T>
void ParamGuess<T>::findGeometricalParameters() {
    
    /// Looking for the regions in the velocity field where the characteristic velocities
    /// are the highest and the lowest ones. Then, a linear regression between these regions 
    /// is performed and the position angle estimated.
    T velmin = AlltoVel<T>(in->getZphys(0),in->Head());
    T velmax = AlltoVel<T>(in->getZphys(in->DimZ()-1), in->Head());
    if (velmin>velmax) std::swap(velmin,velmax);
    float vel_high = vsystem;
    float vel_low  = vsystem;   
    range = ceil(in->Head().Bmaj()/in->Head().PixScale());
    
    for (int y=range; y<ysize-range; y++) {
        for (int x=range; x<xsize-range; x++) {
            long npix = (y+Ymin)*in->DimX()+x+Xmin;
            if (isNaN<T>(Vemap[npix])) continue;
            std::vector<T> vec;
            for (int yi=y-range; yi<=y+range; yi++) 
                for (int xi=x-range; xi<=x+range; xi++) 
                    vec.push_back(Vemap[(yi+Ymin)*in->DimX()+xi+Xmin]);
            T median = findMedian<T>(vec, vec.size());
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
    float rmaj, errmaj[2];
    
    /// Linear regression between the center and the two point found previously.
    /// For not including the center, just change the last parameter from 2 to 1 in 
    /// the function below.
    int a = linear_reg<int>(3, xx, yy, pmaj, errmaj, rmaj, 0, 2);
    
    pmin[0] = - 1/pmaj[0];
    pmin[1] = (ycentre)-pmin[0]*xcentre;
    
    float ang = atan(pmaj[0]);  
    if (coord_high[0]>=xcentre) {   
        if (ang<M_PI_2) posang = (3*M_PI_2+ang)*180/M_PI;
        else posang = (M_PI_2+ang)*180/M_PI;
    }
    else {
        if (ang<M_PI_2) posang = (M_PI_2+ang)*180/M_PI;
        else posang = (ang-M_PI_2)*180/M_PI;    
    }
    
    /// Estimating axes lenght, axes ratio and inclination.
    major_max[0] = major_max[1] = minor_max[0] = minor_max[1] = 0;
    major_min[0] = minor_min[0] = in->DimX()-1; 
    major_min[1] = minor_min[1] = in->DimY()-1;
    double major_r_l=0., major_r_r=0;
    double minor_r_l=0., minor_r_r=0;
    
    for (int x=Xmin; x<=obj->getXmax(); x++) {
        int y_maj = lround(pmaj[0]*x+pmaj[1]);
        int y_min = lround(pmin[0]*x+pmin[1]);
        long mpix_maj = x+y_maj*in->DimX(); 
        long mpix_min = x+y_min*in->DimX();
        if (y_maj>=Ymin && y_maj<=obj->getYmax() && !isNaN<T>(Vemap[mpix_maj])) {
            double r = sqrt(pow(double(x-xcentre),2.)+pow(double(y_maj-ycentre),2.));
            if (r>major_r_l && x<=xcentre) {
                major_r_l = r;
                major_max[0] = x;
                major_max[1] = y_maj;
            }
            if (r>major_r_r && x>xcentre) {
                major_r_r = r;
                major_min[0] = x;
                major_min[1] = y_maj;
            }
        }
        if (y_min>=Ymin && y_min<=obj->getYmax() && !isNaN<T>(Vemap[mpix_min])) {
            double r = sqrt(pow(double(x-xcentre),2.)+pow(double(y_min-ycentre),2.));
            if (r>minor_r_l && x<=xcentre) {
                minor_r_l = r;
                minor_max[0] = x;
                minor_max[1] = y_min;
            }
            if (r>minor_r_r && x>xcentre) {
                minor_r_r = r;
                minor_min[0] = x;
                minor_min[1] = y_min;
            }
        }           
    } 
    
    axmaj = (major_r_r+major_r_l)/2.;
    axmin = (minor_r_r+minor_r_l)/2.;
    //axmaj = axmaj-in->Head().Bmaj()/in->Head().PixScale();
    //axmin = axmin-in->Head().Bmaj()/in->Head().PixScale()-atof(in->pars().getZ0().c_str())/3600./in->Head().PixScale();
    //float axmaj = sqrt(pow(double(major_max[0]-major_min[0]),2.)+pow(double(major_max[1]-major_min[1]),2.));
    //float axmin = sqrt(pow(double(minor_max[0]-minor_min[0]),2.)+pow(double(minor_max[1]-minor_min[1]),2.));
    
    if (axmin>axmaj) {
        std::cout << "---------------> WARNING - Finding initial parameters <--------------\n"
                  << " The major axis is shorter than the minor axis. They will be swapped\n"
                  << " for estimating the inclination.\n"
                  << " The galaxy seem to be less elongated in the kynematical axis!!\n\n";
        std::swap(axmin, axmaj);
    }
    
    inclin = in->pars().getParGF().INC!="-1" ? atof(in->pars().getParGF().INC.c_str()) : acos(axmin/axmaj)*180/M_PI;
    
    Rmax = axmaj*in->Head().PixScale()*arcsconv(in->Head().Cunit(0));
    
}

template <class T>
bool ParamGuess<T>::fitSimplex(int ndim, T **p) {
    
    const int NMAX=5000;            
    const double TINY=1.0e-10;
    const double tol =1.E-03;
    
    int mpts=ndim+1;    
    T minimum;
    
    T psum[ndim], x[ndim];
    T *y = new T[mpts];
    
    for (int i=0; i<mpts; i++) {
        for (int j=0; j<ndim; j++) x[j]=p[i][j];
        y[i]=(this->*func)(x); 
    }
    
    int nfunc=0;
    for (int j=0; j<ndim; j++) {
        T sum=0.0;
        for (int i=0; i<mpts; i++) sum += p[i][j];
        psum[j]=sum;
    }
    
    bool Ok = false;
    for (;;) {
        int ihi, inhi;
        int ilo=0;
        ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
        for (int i=0; i<mpts; i++) {
            if (y[i]<=y[ilo]) ilo=i;
            if (y[i]>y[ihi]) {
                inhi=ihi;
                ihi=i;
            } 
            else if (y[i]>y[inhi] && i!=ihi) inhi=i;
        }
        
        double rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        
        if (rtol<tol) {     
            std::swap(y[0],y[ilo]);
            for (int i=0; i<ndim; i++) {
                std::swap(p[0][i],p[ilo][i]);
            }
            minimum=y[0];
            delete [] y;
            Ok = true;
            break;
        }
        
        if (nfunc>=NMAX) {
            delete [] y;
            Ok = false;
            break;
        }
        nfunc += 2;
        
        double fac=-1.0;
        double fac1=(1.0-fac)/ndim;
        double fac2=fac1-fac;
        T ptry[ndim];
        for (int j=0; j<ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
        T ytry=(this->*func)(ptry);         
        if (ytry<y[ihi]) {  
            y[ihi]=ytry;
            for (int j=0; j<ndim; j++) {
                psum[j] += ptry[j]-p[ihi][j];
                p[ihi][j]=ptry[j];
            }
        }
        
        
        if (ytry<=y[ilo]) {
            fac=2.0;
            fac1=(1.0-fac)/ndim;
            fac2=fac1-fac;
            for (int j=0; j<ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
            ytry=(this->*func)(ptry);           
            if (ytry<y[ihi]) {  
                y[ihi]=ytry;
                for (int j=0; j<ndim; j++) {
                    psum[j] += ptry[j]-p[ihi][j];
                    p[ihi][j]=ptry[j];
                }
            }
        }   
        else if (ytry>=y[inhi]) {
            T ysave=y[ihi];
            fac=0.5;
            fac1=(1.0-fac)/ndim;
            fac2=fac1-fac;
            for (int j=0; j<ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
            ytry=(this->*func)(ptry);           
            if (ytry<y[ihi]) {  
                y[ihi]=ytry;
                for (int j=0; j<ndim; j++) {
                    psum[j] += ptry[j]-p[ihi][j];
                    p[ihi][j]=ptry[j];
                }
            }
            if (ytry>=ysave) {      
                for (int i=0; i<mpts; i++) {
                    if (i!=ilo) {
                        for (int j=0; j<ndim; j++)
                            p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                        y[i]=(this->*func)(psum);
                    }
                }
                
                nfunc += ndim;  
                
                for (int j=0;j<ndim;j++) {
                    T sum=0.0;
                    for (int i=0; i<mpts; i++) sum += p[i][j];
                    psum[j]=sum;
                }
            }
        } 
        else --nfunc; 
    }
    
    return Ok;
}


template <class T>
void ParamGuess<T>::fitEllipse() {
    
    int ndim=2; 
    int mpts=ndim+1;
    
    T *point = new T[ndim];
    T *dels  = new T[ndim];
    T **p = allocate_2D<T>(mpts,ndim);
    
    point[0] = Rmax;
    point[1] = inclin;
    //point[2] = posang;
    //point[3] = xcentre;
    //point[4] = ycentre;
    
    /// Determine the initial simplex.
    for (int i=0; i<ndim; i++) {
        dels[i]  = 0.1*point[i]; 
        point[i] = point[i]-0.05*point[i];      
    }
    
    for (int i=0; i<mpts; i++) {
            for (int j=0; j<ndim; j++) p[i][j]=point[j];
        if (i!=0) p[i][i-1] += dels[i-1];
    }
    
    delete [] point;
    delete [] dels;
    
    func = &ParamGuess<T>::funcEllipse;
    
    bool Ok = fitSimplex(ndim, p);
        
    if (Ok) {
        Rmax   = p[0][0];
        inclin = p[0][1];
        //double Axmin = Rmax*cos(inclin/180.*M_PI);
        //double Axmin_corr = Axmin - in->Head().Bmaj()/in->Head().PixScale();
        //double Axmaj_corr = Rmax - in->Head().Bmaj()/in->Head().PixScale();
        //inclin = acos(Axmin_corr/Axmaj_corr)*180/M_PI;
        //posang = parmin[2];
        //xcentre= parmin[3];
        //ycentre= parmin[4];
        deallocate_2D(p,ndim+1);
    }
    else {
        std::cout << "Error while estimating inclination.";
        std::terminate();
    }
     
}


template <class T>
T ParamGuess<T>::funcEllipse(T *mypar) {
    
    double F = M_PI/180.;
    
    T R   = mypar[0]/(in->Head().PixScale()*arcsconv(in->Head().Cunit(0)));
    T inc = mypar[1];
    T phi = posang;//mypar[2];
    T x0  = xcentre;//mypar[3];
    T y0  = ycentre;//mypar[4];

    double func = 0;
    for (int x=0; x<in->DimX(); x++) {
        for (int y=0; y<in->DimY(); y++) {
            T xr =  -(x-x0)*sin(F*phi)+(y-y0)*cos(F*phi);           
            T yr = (-(x-x0)*cos(F*phi)-(y-y0)*sin(F*phi))/cos(F*inc);
            T r = sqrt(xr*xr+yr*yr);
            bool isIn = r<=R;
            if (!isIn) continue;
            if (isNaN(Vemap[x+y*in->DimX()])) func++;
            else func--;
        }
    }
    
    return func;        

}


template <class T>
void ParamGuess<T>::fitIncfromMap() {
    
    int ndim=2;
    int mpts=ndim+1;
    
    T *point = new T[ndim];
    T *dels  = new T[ndim];
    T **p = allocate_2D<T>(mpts,ndim);
    
    point[0] = inclin;
    point[1] = Rmax;

    
    /// Determine the initial simplex.
    for (int i=0; i<ndim; i++) {
        dels[i]  = 0.1*point[i]; 
        point[i] = point[i];        
    }
    
    for (int i=0; i<mpts; i++) {
            for (int j=0; j<ndim; j++) p[i][j]=point[j];
        if (i!=0) p[i][i-1] += dels[i-1];
    }
    
    delete [] point;
    delete [] dels;
    
    func = &ParamGuess<T>::funcIncfromMap;
    
    bool Ok = fitSimplex(ndim, p);
        
    if (Ok) {
        inclin   = p[0][0];
        Rmax     = p[0][1];
        deallocate_2D(p,ndim+1);
    }
    else {
        std::cout << "Error while estimating inclination.";
        std::terminate();
    }
     
}


template <class T> 
T ParamGuess<T>::funcIncfromMap(T *mypar) {
        
    /*
        if (mypar[0]<0) mypar[0]=1;
        if (mypar[0]>90) mypar[0]=89;
        
        T inc = mypar[0];
                
        Rings<T> *rings = new Rings<T>;     
        rings->radsep=radsep;
        rings->nr = int(Rmax/radsep);
        
        T *prof = new T[rings->nr];
        int *counter = new int[rings->nr];

        cout << mypar[0] << " " << radsep << " " << rings->nr << " " << xcentre << " "<< ycentre << " " << posang <<  endl;
        
        for (int i=0; i<rings->nr; i++) {
            rings->radii.push_back(i*rings->radsep);
            float r1 = rings->radii[i]/(in->Head().PixScale()*arcsconv(in->Head().Cunit(0)));
            float r2 = (i+1)*rings->radsep/(in->Head().PixScale()*arcsconv(in->Head().Cunit(0)));
            prof[i]=counter[i]=0.;
        }

        for (int x=0; x<in->DimX();x++) {
            for (int y=0; y<in->DimY(); y++) {
                float xr = -(x-xcentre)*sin(posang*M_PI/360.)+(y-ycentre)*cos(posang*M_PI/360.);
                float yr = (-(x-xcentre)*cos(posang*M_PI/360.)-(y-ycentre)*sin(posang*M_PI/360.))/cos(inc*M_PI/360.);
                float rad = sqrt(xr*xr+yr*yr);
                for (int i=0; i<rings->nr; i++) {
                    float r1 = rings->radii[i]/(in->Head().PixScale()*arcsconv(in->Head().Cunit(0)));
                    float r2 = (i+1)*rings->radsep/(in->Head().PixScale()*arcsconv(in->Head().Cunit(0)));
                    if (rad>=r1 && rad<r2) {
                        //cout << Intmap[x+y*in->DimX()] << endl;
                        prof[i]+=fabs(Intmap[x+y*in->DimX()]);
                        counter[i]++;
                    }
                }
            }
        }       

        for (int i=0; i<rings->nr; i++) prof[i]/=float(counter[i]);
        
        float minval = *min_element(&prof[0],&prof[0]+rings->nr);
        float mulfac = pow(10,-int(log10(minval)));

        for (int i=0; i<rings->nr; i++) {
            if (counter[i]!=0) {
                cout << setprecision(10);
                cout << prof[i] << " " << std::flush;
                prof[i]*= mulfac;
                cout << prof[i] << " " << counter[i] <<endl;
            }
            rings->vrot.push_back(AlltoVel(10*in->Head().Cdelt(2),in->Head()));
            rings->vdisp.push_back(8);
            rings->z0.push_back(0);
            rings->dens.push_back(prof[i]*1E20);
            rings->inc.push_back(inc);
            rings->phi.push_back(posang);
            rings->xpos.push_back(xcentre);
            rings->ypos.push_back(ycentre);
            rings->vsys.push_back(vsystem);
        }

        
        Model::Galmod<T> *mod = new Model::Galmod<T>;

        mod->input(in,rings,1);
        mod->calculate();
        mod->smooth();

        delete rings;
        delete [] prof;
        
        T *map_mod = new T[in->DimX()*in->DimY()];
        
        float totflux_mod=0;
        for (int x=0; x<in->DimX(); x++){
            for (int y=0; y<in->DimY(); y++){
                map_mod[x+y*in->DimX()]=0;
                for (int z=0; z<in->DimZ(); z++)
                    map_mod[x+y*in->DimX()]+=mod->Out()->Array(x,y,z);
                totflux_mod += map_mod[x+y*in->DimX()];
            }           
        }
        
                
        float factor = totflux_obs/totflux_mod;

        float res_sum=0;
        for (int i=0; i<in->DimX()*in->DimY();i++) {
            if (map_mod[i]!=0) res_sum += fabs(Intmap[i]-map_mod[i]*factor);    
        }

        return res_sum;
        */

    bool verbosity = in->pars().isVerbose();
    in->pars().setVerbosity(false);

    if (mypar[0]<0) mypar[0]=1;
    if (mypar[0]>90) mypar[0]=89;
    if (mypar[1]<0) mypar[1]=2*radsep;
    if (mypar[1]>1.5*Rmax) mypar[1]=Rmax;

    T inc  = mypar[0];
    T RMAX = mypar[1];

    Rings<T> *rings = new Rings<T>;
    rings->radsep=radsep/2.;
    rings->nr = int(RMAX/rings->radsep);

    cout << mypar[0] << " " << mypar[1]<<  "  " << radsep << " " << rings->nr << " " << xcentre << " "<< ycentre << " " << posang <<  endl;

    for (int i=0; i<rings->nr; i++) {
        rings->radii.push_back(i*rings->radsep+rings->radsep/2.);
        rings->vrot.push_back(AlltoVel(10*in->Head().Cdelt(2),in->Head()));
        rings->vdisp.push_back(5);
        rings->z0.push_back(0);
        rings->inc.push_back(inc);
        rings->phi.push_back(posang);
        rings->xpos.push_back(xcentre);
        rings->ypos.push_back(ycentre);
        rings->vsys.push_back(vsystem);
        rings->dens.push_back(1E20);
    }

    MomentMap<T> *totalmap = new MomentMap<T>;
    totalmap->input(in);
    totalmap->SumMap(true);
    for (int i=0; i<totalmap->NumPix();i++) totalmap->Array()[i] = Intmap[i];
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
    mod->input(in,rings);
    mod->calculate();
    mod->smooth();

    delete rings;

    T *map_mod = new T[in->DimX()*in->DimY()];

    float totflux_mod=0;
    for (int x=0; x<in->DimX(); x++){
        for (int y=0; y<in->DimY(); y++){
            map_mod[x+y*in->DimX()]=0;
            for (int z=0; z<in->DimZ(); z++)
                map_mod[x+y*in->DimX()]+=mod->Out()->Array(x,y,z);
            totflux_mod += map_mod[x+y*in->DimX()];
        }
    }

    factor = totflux_obs/totflux_mod;

    float res_sum=0;
    for (int i=0; i<in->DimX()*in->DimY();i++) {
        res_sum += fabs(Intmap[i]-map_mod[i]*factor);
    }

    delete [] map_mod;

    in->pars().setVerbosity(verbosity);

    return res_sum;
}

template <class T>
void ParamGuess<T>::plotGuess() {
    

    //cout << posang << " " << xcentre << " " << ycentre << endl;
    /// Plotting axes in a .eps file.
    std::ofstream outmaj1, outmaj2, outmin, velf;
    std::string outfolder = in->pars().getOutfolder();
    outmaj1.open((outfolder+"axmaj1.dat").c_str());
    outmaj2.open((outfolder+"axmaj2.dat").c_str());
    outmin.open((outfolder+"axmin.dat").c_str());
    velf.open((outfolder+"vfield.dat").c_str());
    for (int x=coord_low[0]-range; x<=coord_low[0]+range; x++) 
        for (int y=coord_low[1]-range; y<=coord_low[1]+range; y++) 
            outmaj1 << x << "   " << y << std::endl;
            
    for (int x=coord_high[0]-range; x<=coord_high[0]+range; x++) 
        for (int y=coord_high[1]-range; y<=coord_high[1]+range; y++) 
            outmaj2 << x << "  " << y << std::endl; 
    
    outmaj1.close();
    outmaj2.close();
    
    std::vector<T> vec;
    float vrange = fabs(DeltaVel<float>(in->Head()));
    for (int x=0; x<in->DimX(); x++) {
        for (int y=0; y<in->DimY(); y++) {
            long npix = x+y*in->DimX(); 
            if (!isNaN(Vemap[npix])) vec.push_back(Vemap[npix]); 
            velf<<x<<" "<<y<<" "<<Vemap[npix]<<std::endl;
            if (Vemap[npix]<=vsystem+vrange && Vemap[npix]>=vsystem-vrange) 
                outmin << x << "   " << y << std::endl;
        }
        velf << std::endl;
    }
    outmin.close();
    velf.close();
    
    float maxvel = *max_element(&vec[0], &vec[0]+vec.size());
    float minvel = *min_element(&vec[0], &vec[0]+vec.size());
    vec.clear();
    
    std::string outfile = outfolder+"axis_major.eps";
    std::ofstream gnu((outfolder+"gnuscript.gnu").c_str());
    T Rmaxpix = Rmax/(in->Head().PixScale()*arcsconv(in->Head().Cunit(0)));
    std::string amaj = to_string(Rmaxpix);
    std::string amin = to_string(Rmaxpix*cos(inclin/180*M_PI));
    std::string posa = to_string(posang/180*M_PI-M_PI_2);
    std::string xcenter = to_string(xcentre);
    std::string ycenter = to_string(ycentre); 
    gnu << "unset key\n"
    //  << "set grid\n"
        << "set title 'Axis fitting'\n"
        << "set cbtics scale 0\n"
        << "set palette defined (0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',"
        << " 5 '#ffee00', 6 '#ff7000',7 '#ee0000',8 '#7f0000')\n"
        << "f(x)="<<to_string(pmaj[0])<<"*x+"<<to_string(pmaj[1])<<std::endl
        << "g(x)="<<to_string(pmin[0])<<"*x+"<<to_string(pmin[1])<<std::endl
        << "set xrange [0:"<<to_string<long>(in->DimX())<<"]\n"
        << "set yrange [0:"<<to_string<long>(in->DimY())<<"]\n"
        << "set cbrange ["<<to_string(minvel)<<":"<<to_string(maxvel)<<"]\n"
        << "set xlabel 'X (pixels)'\n"
        << "set ylabel 'Y (pixels)'\n"
        << "set size square\n"
        << "set parametric\n"
        << "x(t)="+xcenter+"+"+amaj+"*cos("+posa+")*cos(t)-"+amin+"*sin("+posa+")*sin(t)\n"
        << "y(t)="+ycenter+"+"+amaj+"*sin("+posa+")*cos(t)+"+amin+"*cos("+posa+")*sin(t)\n"
        << "set table '"+outfolder+"ellipse.tab'\n"
        << "plot x(t), y(t)\n"
        << "unset table\n"
        << "unset parametric\n"
        << "set terminal postscript eps enhanced color font 'Helvetica,14'\n"
        << "set output '"<<outfile<<"'\n"
        << "plot '"+outfolder+"vfield.dat' w image, '"+outfolder+"axmaj1.dat' ls 1 lc 3, "
        << "'"+outfolder+"axmaj2.dat' ls 1 lc 1, '"+outfolder+"axmin.dat' lc 2, "
        << " '"+outfolder+"ellipse.tab' w l ls -1, f(x) ls 1, g(x) ls 3,'-' ls 5, '-' ls 7 \n"
        << to_string(xcentre)+" "+to_string(ycentre) << std::endl
        << "e" << std::endl
        << to_string(major_max[0])+" "+to_string(major_max[1]) << std::endl
        << to_string(major_min[0])+" "+to_string(major_min[1]) << std::endl
        << to_string(minor_max[0])+" "+to_string(minor_max[1]) << std::endl
        << to_string(minor_min[0])+" "+to_string(minor_min[1]) << std::endl
        << "e" << std::endl;
    
#ifdef HAVE_GNUPLOT     
    Gnuplot gp;
    gp.begin(); 
    gp.commandln(("load '"+outfolder+"gnuscript.gnu'").c_str());
    gp.end();
    remove((outfolder+"ellipse.tab").c_str());

#endif  
    remove((outfolder+"axmaj1.dat").c_str());
    remove((outfolder+"axmaj2.dat").c_str());
    remove((outfolder+"axmin.dat").c_str());
    remove((outfolder+"vfield.dat").c_str());
    remove((outfolder+"gnuscript.gnu").c_str());
    
}

