//----------------------------------------------------------
// moment.cpp: Functions of the MomentMap class.
//----------------------------------------------------------

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
#include <fitsio.h>
#include <vector>
#include <Tasks/moment.hh>
#include <Arrays/array.hpp>
#include <Arrays/cube.hh>
#include <Arrays/image.hh>
#include <Utilities/utils.hh>
#include <Utilities/lsqfit.hh>
#include <Utilities/progressbar.hh>


template <class T>
MomentMap<T>::MomentMap() {
    
    this->numPix = 0;
    this->numAxes = 2;
    this->arrayAllocated = false;
    this->headDefined      = false;
    this->statsDefined   = false;
    
}


template <class T> 
void MomentMap<T>::input (Cube<T> *c, int *Blo, int *Bhi, bool *m) {
    
    in = c;
    mask = m;

    for (int i=0; i<3; i++) {
        blo[i] = Blo[i];
        bhi[i] = Bhi[i];
    }
    
    this->axisDim[0] = bhi[0]-blo[0];
    this->axisDim[1] = bhi[1]-blo[1];
    nsubs = bhi[2]-blo[2]; 
    
    this->numPix = this->axisDim[0]*this->axisDim[1];
    this->array = new T[this->numPix];
    this->arrayAllocated = true;
    this->par = c->pars();
}


template <class T> 
void MomentMap<T>::input (Cube<T> *c, bool *m) {
    
    int blo[3] = {0,0,0};
    int bhi[3] = {c->AxesDim(0),c->AxesDim(1),c->AxesDim(2)};
    
    input(c,blo,bhi,m);
    
}


template <class T>
void MomentMap<T>::SumMap (bool msk) {

    // Just sum over che channels
    bool headdef = in->HeadDef();
    in->setHeadDef(false);
    ZeroMoment(msk);
    in->setHeadDef(headdef);
    this->headDefined=setHead(0);
    this->head.setBunit(in->Head().Bunit());

}


template <class T>
void MomentMap<T>::HIMassDensityMap (bool msk) {

    // Write a Mass density map in Msun/pc2. BUNIT must be JY/beam or Kelvin
    
    // Check that input units are jy / beam or Kelvin
    std::string bunit = makelower(in->Head().Bunit());
    bool isJY = bunit.find("/b")!=std::string::npos;
    bool isK  = bunit=="k" || bunit=="kelvin";
    if (isJY || isK) {      
        if(msk && !in->MaskAll()) in->BlankMask();
        bool v = in->pars().isVerbose();
        in->pars().setVerbosity(false);
        if (v) std::cout << " Extracting Mass surfance-density map... " << std::flush;
        SumMap(msk);
        in->pars().setVerbosity(v);
        if (v) std::cout << "Done. " << std::endl;
        in->checkBeam();
        double bmaj = in->Head().Bmaj()*arcsconv(in->Head().Cunit(0));
        double bmin = in->Head().Bmin()*arcsconv(in->Head().Cunit(1));
        double dvel = fabs(DeltaVel<T>(in->Head()));
        for (int i=0; i<this->numPix; i++) {
            T val = this->array[i];
            if (isK) val /= 1360*21.106114*21.106114/(bmaj*bmin);
            this->array[i]=8794*val*dvel/(bmaj*bmin);
        }
        this->head.setBunit("Msun/pc2");
        this->head.setBtype("Mass_surfdens");
        
    }
    else {
        std::cerr << "MOMENT MAPS ERROR: Input datacube must be in JY/BEAM or K.\n";
        std::terminate();
    }
}


template <class T>
void MomentMap<T>::storeMap(bool msk, int whichmap, std::string map_type) {
    
    // This function calculate the requested map (whichmap) and store it 
    // into the Image2D::array.
    
    if (!this->arrayAllocated) {
        std::cout << "MOMENT MAPS error: ";
        std::cout << "Array not allocated. Call 'input' first!!\n";
        std::terminate();
    }
        
    if (!(this->headDefined=setHead(whichmap)) && in->pars().isVerbose()) {
        std::cout<< "MOMENT MAPS warning: cannot create new header.\n";
    }
    
    // Creating mask if it does not exist
    if(msk && mask==nullptr) {
        if (!in->MaskAll()) in->BlankMask();
        mask = in->Mask();
    }

    std::string barstring;
    if (whichmap==0)      barstring = " Extracting intensity map ";
    else if (whichmap==1) barstring = " Extracting velocity map ";
    else if (whichmap==2) barstring = " Extracting velocity dispersion map ";
    
    if (map_type=="GAUSSIAN") {
        map_Type = &MomentMap<T>::fitSpectrum;
        barstring += "(GAUSSIAN)... ";
    }
    else {
        map_Type = &MomentMap<T>::calculateMoments;
        barstring += "(MOMENT)... ";
    }
    
    ProgressBar bar(true,in->pars().isVerbose(),in->pars().getShowbar());

    int nthreads = in->pars().getThreads();    
#pragma omp parallel num_threads(nthreads)
{
    bar.init(barstring,this->axisDim[0]);
#pragma omp for
    for (int x=0; x<this->axisDim[0]; x++) {
        bar.update(x+1);
        for (int y=0; y<this->axisDim[1]; y++) {
            this->array[x+y*this->axisDim[0]]=log(-1);
            double moms[3];
            if((this->*map_Type)(x,y,msk,moms))
                this->array[x+y*this->axisDim[0]] = moms[whichmap];
        }
    }
}
    
    bar.fillSpace("Done.\n");
    
}


template <class T> 
void MomentMap<T>::RMSMap (float level, float sncut) {
    
    // Compute the RMS map, i.e. the RMS in each spectrum.
    // Use an iterative way: calculate rms, cut at 1.5*rms, 
    // start again until convergence at level "level".

    if (!this->arrayAllocated) {
        std::cout << "MOMENT MAPS error: ";
        std::cout << "Array not allocated. Call 'input' first!!\n";
        std::terminate();
    }
        
    if (!(this->headDefined=setHead(3))) {
        std::cout<< "MOMENT MAPS warning: cannot create new header.\n";
    }
    
    int nthreads = in->pars().getThreads();
    bool rob = in->pars().getFlagRobustStats();
    
    // Cube sizes
    size_t xs = in->DimX(), ys = in->DimY(), zs = in->DimZ();
    
    // Progress bar
    ProgressBar bar(true,in->pars().isVerbose(),in->pars().getShowbar());
    
#pragma omp parallel num_threads(nthreads)
{
   bar.init(" Computing RMS map... ",xs*ys);
#pragma omp for
    for (size_t xy=0; xy<xs*ys; xy++) {
        bar.update(xy+1);
            
        // Getting the spectrum at x,y pixel
        std::vector<float> sp(zs);
        for (size_t z=0; z<zs; z++) 
            if (in->Array(xy+z*xs*ys)==in->Array(xy+z*xs*ys)) 
                sp[z] = in->Array(xy+z*xs*ys);
            
        // Start main loop
        float orms = 1E10;
        size_t count = 0;
        while (true) {
            // Calculate median and MADFM
            float rms = 0;
            if (rob) {
                float median = findMedian(&sp[0],sp.size(),false);
                rms = findMADFM(&sp[0],sp.size(),median,false)/0.6745;
            }
            else {
                rms = findStddev(&sp[0],sp.size());
            }
            // Calculate improvement wrt previous step
            float rat = (orms - rms)/orms;
            // If meet criteria, exit 
            if (rat<level || count++>100 || sp.size()<5) break;
            // Else S/N cut the spectrum 
            orms = rms;
            for (int z=sp.size(); z--;) 
                if (sp[z]>sncut*rms) sp.erase(sp.begin()+z);
        }
        this->array[xy] = orms;
    }
}
    bar.fillSpace("Done.\n");

}


template <class T> 
bool MomentMap<T>::setHead(int type) {
    
    if (!in->HeadDef()) this->headDefined = false; 
    else { 
        this->copyHeader(in->Head());
        this->head.setCrpix(0, in->Head().Crpix(0)-blo[0]);
        this->head.setCrpix(1, in->Head().Crpix(1)-blo[1]);
        if (type==0) {              
            this->head.setBtype("intensity");
            std::string bunit;
            if (FluxtoJy(1,in->Head())==1) {
                bunit = in->Head().Bunit() + " * KM/S";
            }
            else bunit = "JY * KM/S";
            this->head.setBunit(bunit);
        }
        else if (type==1 || type==2) {
            if (type==1) this->head.setBtype("velocity");
            else this->head.setBtype("dispersion");
            this->head.setBunit("KM/S");
        }
        else if (type==3) {
            this->head.setBtype("rms");
        }
        this->head.calcArea();
        this->headDefined = true;
    }   
    
    return this->headDefined;
    
}


template <class T>
bool MomentMap<T>::calculateMoments (size_t x, size_t y, bool msk, double *moments) {
    
    // An array to store the spectrum at (x,y) position
    double *spectrum = new double[nsubs]; 
    // Velocities in km/s
    double *vels = new double[nsubs];
    
    T num=0, denom=0;
    for (int z=0; z<nsubs; z++) {
        long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
        vels[z] = z+blo[2];
        if (in->HeadDef()) vels[z] = AlltoVel<T>(in->getZphys(z), in->Head());
                
        if (msk) {
            num   += in->Array(npix)*vels[z]*mask[npix];
            denom += in->Array(npix)*mask[npix];
        }
        else {
            num   += in->Array(npix)*vels[z];
            denom += in->Array(npix);
        }
    }
    
    // If all pixels are masked return all NaNs
    if (denom==0) {
        moments[0] = moments[1] = moments[2] = log(-1);
        delete [] spectrum;
        delete [] vels;
        return true;
    }
    
    // Moment 0th
    moments[0] = denom;
    if (in->HeadDef()) 
        moments[0] = FluxtoJy(denom, in->Head()) * fabs(DeltaVel<T>(in->Head()));
    
    // Moment 1st
    moments[1] = num/denom;

    // Calculating 2nd moment
    num = 0;
    for (int z=0; z<nsubs; z++) {
        long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
        if (msk) num += in->Array(npix)*(vels[z]-moments[1])*(vels[z]-moments[1])*mask[npix];
        else num += in->Array(npix)*(vels[z]-moments[1])*(vels[z]-moments[1]);
    }
    
    // Moment 2nd
    moments[2] = sqrt(num/denom);
    
    delete [] spectrum;
    delete [] vels;
    return true;
}


template <class T>
bool MomentMap<T>::fitSpectrum (size_t x, size_t y, bool msk, double *bestfitpar) {

    // An array to store the spectrum at (x,y) position
    double *spectrum = new double[nsubs]; 
    // Weights and velocities in km/s
    double *ww = new double[nsubs];
    double *vels = new double[nsubs];
    // Parameters of the Gaussian fit and their errors
    double c[3], cerr[3];
    // Paramters to fit
    bool mp[3] = {true,true,true};
    
    float smax = -FLT_MAX, vmax = 0;
    for (int z=0; z<nsubs; z++) {
        ww[z] = 1;
        spectrum[z] = in->Array(x,y,z);
        if (msk) spectrum[z] *= mask[in->nPix(x,y,z)];
        vels[z] = AlltoVel(in->getZphys(z),in->Head());
        // Finding spectrum maximum value and corresponding position
        if (spectrum[z]>smax) {
            smax = spectrum[z];
            vmax = vels[z];
        }
    }
    
    // Initial guesses of the Gaussian fit
    c[0] = smax;
    c[1] = vmax;
    c[2] = 10.;
            
    if (smax==0) return false;
    
    Lsqfit<double> lsq(vels,1,spectrum,ww,nsubs,c,cerr,mp,3,&func_gauss,&derv_gauss);
    int ret;
    ret = lsq.fit();
    if(ret<0) return false;
    
    double integint = sqrt(2)*sqrt(M_PI)*fabs(c[2])*FluxtoJy(c[0], in->Head());
    bestfitpar[0] = integint;        // Integrated intensity
    bestfitpar[1] = c[1];            // Central velocity
    bestfitpar[2] = c[2];            // Velocity dispersion
    
    delete [] spectrum;
    delete [] ww;
    delete [] vels;
    
    return true;
}



template <class T>
std::vector< MomentMap<T> > getAllMoments(Cube<T> *c, bool usemask, bool *mask, string mtype) {

    /// This function computes 0th, 1st and 2nd moment maps in a computationally
    /// efficient way. Maps are stored and returned in a vector of MomentMap
    /// instances.

    std::vector< MomentMap<T> > allmaps(3);

    // Creating mask if it does not exist
    if(mask==nullptr && usemask) {
        if (!c->MaskAll()) c->BlankMask();
        mask = c->Mask();
    }

    // Initiliazing all moment maps
    for (int i=0; i<3; i++) {
        allmaps[i].input(c,mask);
        allmaps[i].setHeadDef(allmaps[i].setHead(i));
    }

    ProgressBar bar(true,c->pars().isVerbose(),c->pars().getShowbar());

    int nthreads = c->pars().getThreads();
#pragma omp parallel num_threads(nthreads)
{
    bar.init(" Deriving kinematic maps... ",c->DimY());
#pragma omp for
    for (int y=0; y<c->DimY(); y++) {
        bar.update(y+1);
        for (int x=0; x<c->DimX(); x++) {
            double moms[3];
            if (mtype=="GAUSSIAN") allmaps[0].fitSpectrum(x,y,usemask,moms);
            else allmaps[0].calculateMoments(x,y,usemask,moms);
            allmaps[0].Array(x,y) = moms[0];
            allmaps[1].Array(x,y) = moms[1];
            allmaps[2].Array(x,y) = moms[2];
        }
    }
}

    bar.fillSpace("Done.\n");

    return allmaps;
}
template std::vector< MomentMap<short> > getAllMoments(Cube<short>*,bool,bool*,std::string);
template std::vector< MomentMap<int> > getAllMoments(Cube<int>*,bool,bool*,std::string);
template std::vector< MomentMap<long> > getAllMoments(Cube<long>*,bool,bool*,std::string);
template std::vector< MomentMap<float> > getAllMoments(Cube<float>*,bool,bool*,std::string);
template std::vector< MomentMap<double> > getAllMoments(Cube<double>*,bool,bool*,std::string);

////////////////////////////////////////////////////////////////////////////////////////
// Functions for PvSlice class
////////////////////////////////////////////////////////////////////////////////////////
template <class T>
PvSlice<T>::PvSlice(Cube<T> *c) {
    
    in = c;
    Param &p = c->pars();
    
    bool valid2points = p.getP1_PV(0)>=0 && p.getP1_PV(1)>=0 &&
                        p.getP2_PV(0)>=0 && p.getP2_PV(1)>=0;
        
    if (valid2points) {
        // The slice is defined with two points
        x1 = p.getP1_PV(0); y1 = p.getP1_PV(1);
        x2 = p.getP2_PV(0); y2 = p.getP2_PV(1);
        isAngle = false;
    }
    else if (p.getXPOS_PV()!="-1" && p.getYPOS_PV()!="-1") {
        // The slice is defined with 1 point and a angle
        std::string pos[2] = {p.getXPOS_PV(), p.getYPOS_PV()};
        double *pixs  = getCenterCoordinates(pos, c->Head());
        x0 = pixs[0];
        y0 = pixs[1];
        angle = c->pars().getPA_PV();
        isAngle = true;
    }
    else {
        throw("PVSLICE ERROR: no slice has been defined!");
    }

}


template <class T>
bool PvSlice<T>::slice() {
    
    // Front-end function to slice the cube and extract the position-velocity.
    // Slice is made through the data provided in the constructor.
    
    xpix = in->DimX();
    ypix = in->DimY();
    zpix = in->DimZ();
    
    
    if (!(xpix*ypix*zpix > 0)) {
        fprintf (stderr, "PvSlice ERROR: input cube dimensions are wrong.\n");
        return false;
    }
    
    // If the slice is defined through a point and a angle, the slice is always
    // taken across the entire cube, so calculate intercepts at 0 and xsize-1
    if (isAngle) {
        if (x0<0 || x0>=xpix || y0<0 || y0>=ypix) {
            fprintf (stderr, "PvSlice ERROR: center (x0,y0) is outside the cube.\n");
            return false;
        }
        while (angle>=360) angle-=360;
        while (angle<=-360) angle+=360;
        if (angle==0 || angle==180) {
            x1 = x2 = x0;
            y1 = 0; y2 = ypix-1;
        }
        else {
            float theta = (angle+90)*M_PI/180.;
            x1 = 0;      y1 = tan(theta)*(x1-x0)+y0;
            x2 = xpix-1; y2 = tan(theta)*(x2-x0)+y0;
        }
    }
    
    // Get the slice locus 
    if(!define_slice(x1,y1,x2,y2)) return false;
    
    if (num_points>0) {
        // Setting the PV array
        int dimen[2] = {num_points,zpix};
        this->setImage(dimen);
        
        // Extract slice from cube 
        if (!pvslice()) return false;
        define_header();
    }
    
    return true;
}


template <class T>
bool PvSlice<T>::define_slice(int x1, int y1, int x2, int y2) {

    // Determine a locus of pixels in a slice line, from two endpoints 
    
    int    blx, bly, Trx, Try;
    float  theta, ctheta, stheta;

    if (!check_bounds(&blx,&bly,&Trx,&Try) ) return false;

    // Now calculate slice length 
    float dx    = (float) (Trx-blx+1);
    float dy    = (float) (Try-bly+1);
    if (blx==Trx) {
        theta = M_PI/2.0; ctheta=0.0; stheta=1.0;
        if (bly>Try) {theta = 3.0*M_PI/2.0; stheta=-1.0;}
    }
    else {theta = atan2(dy,dx); ctheta = cos(theta); stheta = sin(theta);}

    float slicelength = sqrt(dx*dx+dy*dy);
    num_points = lround(slicelength);

    if (locusAllocated) delete [] x_locus;
    if (locusAllocated) delete [] y_locus;
    x_locus = new float[num_points];
    y_locus = new float[num_points];
    locusAllocated = true;
    
    for (int i=0; i<num_points; i++) {
        x_locus[i] = i * ctheta + blx;
        y_locus[i] = i * stheta + bly;
    }

    return true;
} 


template <class T>
bool PvSlice<T>::check_bounds (int *blx, int *bly, int *Trx, int *Try) {
    
    // Checks bounds & if endpoints are outside the area of cube's front face,
    // computes overlapping section of slice 

    float f, x, y;
    int swapped=0, v[5], w[5];
    float ix1 = x1, iy1 = y1,  ix2 = x2, iy2 = y2;

    // Do everything assuming first point is left of second, swap back later 
    if (ix2 < ix1) {
        ix1 = x2; ix2 = x1; iy1 = y2; iy2 = y1; swapped=1;
    } 
    else {
        if (ix1 == ix2 && iy1 > iy2) {
            iy1 = y2; iy2 = y1; swapped=2;
        }
    }
   
    *blx = ix1; *bly = iy1; *Trx = ix2; *Try = iy2;

    // Handle no-overlap cases: the extended slice line may overlap eventually,
    // but the piece within the endpoints does not. So stop here.
    if ( (ix1 < 0 && ix2 < 0) || (ix1 >= xpix && ix2 >= xpix) ) {
        fprintf(stderr, "PvSlice ERROR: Slice does not intercept cube.\n"); return false;
    }
    if ( (iy1 < 0 && iy2 < 0) || (iy1 >= ypix && iy2 >= ypix) ) {
        fprintf(stderr, "PvSlice ERROR: Slice does not intercept cube.\n"); return false;
    }
    // Single pixel case
    if ( ix1 == ix2 && iy1 == iy2 ) {
        fprintf(stderr, "PvSlice ERROR: Line segment too short\n"); return false;
    }

    // Handle the case of a vertical slice 
    if ( ix1 == ix2 ) {

        if (ix1 < 0 || ix1 >= xpix) {
            fprintf(stderr, "PvSlice ERROR: Slice does not intercept cube.\n"); return false;
        }
        else {
            *blx = ix1; *Trx = ix2;
            // Now check if either of the y coords are good 
            if (iy1 < 0 || iy1 >= ypix) { *bly = 0;}
            if (iy2 < 0 || iy2 >= ypix) { *Try = ypix-1;}
        }
    }
    else {
        // Handle the case of a horizontal slice
        if ( iy1 == iy2 ) {

            if (iy1 < 0 || iy1 >= ypix) {
                fprintf(stderr, "PvSlice ERROR: Slice does not intercept cube.\n"); return false;
            }
            else {
                *bly = iy1; *Try = iy2;
                // Now check if either of the x coords are good 
                if (ix1 < 0 || ix1 >= xpix) { *blx = 0;}
                if (ix2 < 0 || ix2 >= xpix) { *Trx = xpix-1;}
            }
        }
        else {
            // Calculate intercepts with x=0, y=0, x=xpix-1, y=ypix-1
            // store x-intercepts in v[] and y-intercepts in w[]       
            f = ( (float)(iy2) - (float)(iy1) ) / ( (float)(ix2) - (float)(ix1) );

            v[1]=0; w[2]=0; v[3]=xpix-1; w[4]=ypix-1;
            x=0.0;    y=(float)(iy1) + (x-(float)(ix1))*f;          w[1]=(int)(y+0.5);
            y=0.0;    x=(float)(ix1) + (y-(float)(iy1))/f;          v[2]=(int)(x+0.5);
            x=(float)(xpix-1); y=(float)(iy1) + (x-(float)(ix1))*f; w[3]=(int)(y+0.5);
            y=(float)(ypix-1); x=(float)(ix1) + (y-(float)(iy1))/f; v[4]=(int)(x+0.5);


            // Things are different depending on if slope is +ve or -ve 
            // For either, there are six cases of where the intercepts of the
            // slice line (extended to infinity in both directions) lie w.r.t. the
            // four lines noted above
            if ( f >= 0.0 ) {
                if (v[2] < 0 ) {
                    if (w[1] < ypix) {
                        *blx = 0; *bly = w[1];
                        if ( w[3] < ypix ) {*Trx = xpix-1; *Try = w[3];}
                        else{*Trx = v[4]; *Try = ypix-1;}
                    }
                    else {fprintf(stderr, "PvSlice ERROR: Out of bounds\n"); return false;}
                }
                else{
                    if (v[2] < xpix) {
                        *bly = 0; *blx = v[2];
                        if (v[4] < xpix) {*Try = ypix-1; *Trx = v[4];}
                        else {*Try = w[3]; *Trx = xpix-1;}
                    }
                    else{fprintf(stderr, "PvSlice ERROR: Out of bounds\n"); return false;}
                }
            }
            else {
                if (v[4] < 0) {
                    if ( w[1] > 0 ) {
                        *blx = 0; *bly = w[1];
                        if (w[3] > 0 ) {*Trx = xpix-1; *Try = w[3];}
                        else {*Trx = v[2]; *Try = 0;}
                    }
                    else {fprintf(stderr, "PvSlice ERROR: Out of bounds\n"); return false;}
                }
                else {
                    if (v[4] < xpix) {
                        *blx = v[4]; *bly = ypix-1;
                        if ( v[2] < xpix ) {*Trx = v[2]; *Try = 0;}
                        else {*Trx = xpix-1; *Try = w[3];}
                    }
                    else {fprintf(stderr, "PvSlice ERROR: Out of bounds\n"); return false;}
                }
            } 
        } 
    } 
    
    // now check to see if either endpoint NEEDS to be replaced:
    // if they are legal return the input value 
    if ( (ix1 >= 0 && ix1 < xpix) && (iy1 >= 0 && iy1 < ypix) ) {
        *blx = ix1; *bly = iy1;
    }
    if ( (ix2 >= 0 && ix2 < xpix) && (iy2 >= 0 && iy2 < ypix) ) {
        *Trx = ix2; *Try = iy2;
    }

    // Swap ends back if necessary
    if (swapped > 0) {
        ix1 = *blx; *blx = *Trx; *Trx = ix1;
        iy1 = *bly; *bly = *Try; *Try = iy1;
    }

    return true;

} 


template <class T>
bool PvSlice<T>::pvslice () {

    // Extract pixels from a cube, along an input locus and write the PV in 
    // the main array. Result is antialiased.
    // 
    // The cube is assumed to have axes in x,y,v order.
    // The output array has the same spatial scale as input cube, ie
    // the scale is assumed to be the same for both spatial axes,
    // and velocity pixels are given the same width as the channel spacing.
    // Output is a weighted sum over all pixels nearby the point where the
    // slice locus passes, to reduce aliasing effects.
    //
    // There are four cases to consider for the neighbour pixels, depending
    // where within a pixel the hit occurs:
    //
    //         |-----------|-----------|     |-----------|-----------|
    //         |           |           |	   |           |           |
    //    2.0  -     3     |     2     | 	   |     3     |     2     |
    //         |           |           |	   |       xy  |           |
    //         |-----------|-----------|	   |-----------|-----------|
    //         |           | xy        |	   |           |           |
    //    1.0  -     4     |     1     |	   |     4     |     1     |
    //         |           |           |	   |           |           |
    //         |-----|-----|-----|-----|	   |-----------|-----------|
    //              1.0         2.0
    //
    //         |-----------|-----------|     |-----------|-----------|
    //         |           |           |	   |           |           |
    //         |     3     |     2     | 	   |     3     |     2     |
    //         |           | xy        |	   |           |           |
    //         |-----------|-----------|	   |-----------|-----------|
    //         |           |           |	   |       xy  |           |
    //         |     4     |     1     |	   |     4     |     1     |
    //         |           |           |	   |           |           |
    //         |-----------|-----------|	   |-----------|-----------|
    //
    // All these can be handled by int(x+/-0.5), int(y+/-0.5)
    // Pixel coordinates are associated with centres of pixels. In pixels that 
    // don't have 4 neighbours, missing neighbours are aasigned a zero weight.

    
    int    xp, yp, xn, yn;
    float  xc, yc, xx, yy, sw, f;
    size_t MAXNB = 5;    // Number of neighbours for antialiasing
    

    if (num_points < 2) return false;
    
    float2D wt (MAXNB,size_t(num_points));
    int3D nb(size_t(2),MAXNB,size_t(num_points));

    // Find the neighbour pixels & antialiasing weights.
    // this is done for one channel only, then list applied to all channels
    for (int i = 0; i < num_points; i++) {
        xc = x_locus[i];       yc = y_locus[i];
        xp = (int)( xc+0.5 );  yp = (int)( yc+0.5 );
        // here I only do the 4 nearest pixels
        nb(0,0,i) = xp; nb(1,0,i) = yp;
        nb(0,1,i) = (int)(xc+0.5); nb(1,1,i) = (int)(yc-0.5); // neigbours
        nb(0,2,i) = (int)(xc+0.5); nb(1,2,i) = (int)(yc+0.5);
        nb(0,3,i) = (int)(xc-0.5); nb(1,3,i) = (int)(yc+0.5);
        nb(0,4,i) = (int)(xc-0.5); nb(1,4,i) = (int)(yc-0.5);

        // calculate the weight for each neighbour
        for (int j = 1; j < MAXNB; j++) {
            xp = nb(0,j,i); yp = nb(1,j,i);
            xx = (float) xp;   yy = (float) yp;
            if ( (xp >= 0) && (xp < xpix) && (yp >= 0) && (yp < ypix) )
                wt(j, i) = weight (xx, yy, xc, yc);
            else wt(j, i) = 0.0;
        }   
    }   

    
    for (int z = 0; z < zpix; ++z) {
        for (int i=0; i<num_points; i++) { // Start slice loop 
            f = 0.0; sw = 0.0;
            xp = nb(0,0,i);
            yp = nb(1,0,i);

            for (int j=1; j<MAXNB; j++) {
                xn = nb(0,j,i);
                yn = nb(1,j,i);
                // The first test protects the second; it is possible for xn and yn
                // to go out of range (neighbours of a pixel at edge of cube) 
                if ( wt(j,i) >0.0 && in->Array(xn,yn,z)!=1.0e30) {
                    sw += wt(j,i);
                    f  += in->Array(xn,yn,z) * wt(j,i);
                }
            } 
            
            if (sw>0.0) this->array[i+z*num_points] = f/sw;
            else {
                if ( xp >= 0 && xp<xpix && yp>=0 && yp<ypix ) this->array[i+z*num_points] = in->Array(xn,yn,z);
                else this->array[i+z*num_points] = 0.0;
            }
            
        } 
    } 
    
    return true;
} 


template <class T>
bool PvSlice<T>::slice_old() {
    
    // This is my first function to extract PV. It is less polished than the newer function
    // I am still using this in 3DFIT, but I HAVE TO CHECK!
    
    float phi = angle;
    while(phi>=180) phi -= 180;
    while(phi<0) phi += 180;

    double P = phi*M_PI/180.;
    xpix = in->DimX();
    ypix = in->DimY();
    zpix = in->DimZ();
    
    int dim[2] = {0, zpix}; 
    int xdim=xpix, ydim=ypix;
    int xmax=xpix, ymax=ypix;
    int xmin=0, ymin=0;

    Header &h = in->Head();
    float cdelt0;
    
    if (phi==90) {
        dim[0] = xpix;
        this->setImage(dim);
        for (int x=0; x<dim[0]; x++)
            for (int z=0; z<dim[1]; z++)
                this->array[x+z*dim[0]] = in->Array(x,y0,z);
        cdelt0 = fabs(h.Cdelt(0));
    }
    else {
        std::vector<int> xx, yy; 
        double mx=0, my=0;
    
        if (phi<90) my = tan(P+M_PI_2);
        else my = tan(P-M_PI_2);

        mx = 1./my;
    
        int x_1 = lround(x0+(ypix-y0)/my); 
        int x_0 = lround(x0-y0/my);
        xmax = my>0 ? x_1 : x_0;
        xmin = my>0 ? x_0 : x_1;
        if (xmax<xmin) std::swap(xmax,xmin);
        if (xmax>xpix) xmax=xpix;
        if (xmin<0) xmin=0;
        xdim = fabs(xmax-xmin);
            
        int y_1 = lround(y0+(xpix-x0)/mx); 
        int y_0 = lround(y0-x0/mx);
        ymax = mx>0 ? y_1 : y_0;
        ymin = mx>0 ? y_0 : y_1;
        if (ymax<ymin) std::swap(ymax,ymin);
        if (ymax>ypix) ymax=ypix;
        if (ymin<0) ymin=0;
        ydim = fabs(ymax-ymin);
    
        if (xdim>=ydim) {
            int nxdim=0;
            for (int x=0; x<xpix; x++) { 
                int y1 = lround(my*(x-x0)+y0);
                bool isin = y1>=0 && y1<ypix;  
                if (isin) {
                    xx.push_back(x);
                    yy.push_back(y1);
                    nxdim++;
                }
            }
            dim[0]=nxdim;
        }
        else {
            int nxdim=0;
            for (int y=0; y<ypix; y++) { 
                int x1 = lround(mx*(y-y0)+x0);
                bool isin = x1>=0 && x1<xpix;  
                if (isin) {
                    xx.push_back(x1);
                    yy.push_back(y);
                    nxdim++;
                }
            }
            dim[0]=nxdim;
        }

        this->setImage(dim);
        for (int i=0; i<dim[0]; i++) { 
            for (int z=0; z<zpix; z++) {
                this->array[i+z*dim[0]] = in->Array(xx[i],yy[i],z);
            } 
        }

        float xdom = xdim*h.Cdelt(0);
        float ydom = ydim*h.Cdelt(1);
        cdelt0 = sqrt(xdom*xdom+ydom*ydom)/dim[0];
    }

    this->copyHeader(h);
    
    float crpix0 = xdim>ydim ? x0-xmin : y0-ymin;
    if (phi==90) crpix0 = x0;
    this->Head().setCrpix(0, crpix0+1);
    
    this->Head().setCrval(0, 0);
    this->Head().setCdelt(0, cdelt0);
    this->Head().setCtype(0, "Offset");
    this->Head().setCrpix(1, h.Crpix(2));
    this->Head().setCdelt(1, h.Cdelt(2));
    this->Head().setCrval(1, h.Crval(2));
    this->Head().setCunit(1, h.Cunit(2));
    this->Head().setCtype(1, h.Ctype(2));
    this->Head().setMinMax(0.,0.);

    std::string name = h.Name()+"_pv"+to_string(phi,0);
    this->Head().setName(name);
    this->setHeadDef(true);
    
    return true;
    
}

template <class T>
void PvSlice<T>::define_header () {
    
    // Defining the Header for the PV slice.
    
    ///* Calculating x-axis WCS
    double *coord_start = in->getXYphys(x_locus[0],y_locus[0]);
    double *coord_end   = in->getXYphys(x_locus[num_points-1],y_locus[num_points-1]);
    double coord_center[2] = {0.5*(coord_start[0]+coord_end[0]),0.5*(coord_start[1] + coord_end[1])};
    
    double ra_off  = (coord_start[0]-coord_center[0])*cos(coord_center[1]*M_PI/180.);
    double dec_off = coord_end[1] - coord_center[1];
    double position_angle = atan2(-ra_off, dec_off)*180./M_PI - 90.0;
    if (position_angle < 0) position_angle += 360.0;
    double first_coord  = -sqrt( ra_off*ra_off + dec_off*dec_off );
    
    double cdelt1 = fabs(2*first_coord/(num_points-1));
    double pcent = 0.5*num_points;
    
    if (isAngle) {
        // If the center is given, let's center the WCS on it.
        double dx_0 = sqrt( (x0-x_locus[0])*(x0-x_locus[0]) + (y0-y_locus[0])*(y0-y_locus[0]) ) ;
        double dx_t = sqrt( (x_locus[0]-x_locus[num_points-1])*(x_locus[0]-x_locus[num_points-1]) 
                          + (y_locus[0]-y_locus[num_points-1])*(y_locus[0]-y_locus[num_points-1]) );
        pcent =  dx_0/dx_t*num_points;
    }
    
    
    // Setting spatial coordinates
    this->Head().setCtype(0,"OFFSET");
    this->Head().setCrpix(0,pcent);
    this->Head().setCrval(0,0);
    this->Head().setCdelt(0,cdelt1);
    this->Head().setCunit(0,in->Head().Cunit(0));
    
    // Setting spectral coordinates
    this->Head().setCtype(1,in->Head().Ctype(2));
    this->Head().setCrpix(1,in->Head().Crpix(2));
    this->Head().setCrval(1,in->Head().Crval(2));
    this->Head().setCdelt(1,in->Head().Cdelt(2));
    this->Head().setCunit(1,in->Head().Cunit(2));
    
    // Setting other properties
    this->Head().setBmaj(in->Head().Bmaj());
    this->Head().setBmin(in->Head().Bmin());
    this->Head().setBpa(in->Head().Bpa());
    this->Head().setBunit(in->Head().Bunit());
    this->Head().setBtype(in->Head().Btype());
    this->Head().setEpoch(in->Head().Epoch());
    this->Head().setTelesc(in->Head().Telesc());
    std::string name = in->Head().Name()+"_pv"+to_string(angle,0);
    this->Head().setName(name);
    
    std::string str;
    if (isAngle) {
        double *crt   = in->getXYphys(x0,y0);
        str = "RA="+to_string(crt[0])+" deg, DEC="+to_string(crt[1])+" deg, PA="+to_string(angle,1)+" deg";
    }
    else 
        str = "P1=("+to_string(x1,1)+","+to_string(y1,1)+") P2=("+to_string(x2,1)+","+to_string(y2,1)+")";
    
    this->Head().Keys().push_back("HISTORY BBAROLO PVSLICE: "+str);
    this->setHeadDef(true);

}



template <class T>
PvSlice<T>* PositionVelocity (Cube<T> *c, float x0, float y0, float Phi, bool oldmethod) {
    
    PvSlice<T> *pv = new PvSlice<T>(c,x0,y0,Phi);
    if (oldmethod) pv->slice_old();
    else pv->slice();
    return pv;
    
}
template PvSlice<short>* PositionVelocity (Cube<short>*,float,float,float,bool);
template PvSlice<int>* PositionVelocity (Cube<int>*,float,float,float,bool);
template PvSlice<long>* PositionVelocity (Cube<long>*,float,float,float,bool);
template PvSlice<float>* PositionVelocity (Cube<float>*,float,float,float,bool);
template PvSlice<double>* PositionVelocity (Cube<double>*,float,float,float,bool);



// Explicit instantiation of the classes
template class MomentMap<short>;
template class MomentMap<int>;
template class MomentMap<long>;
template class MomentMap<float>;
template class MomentMap<double>;

template class PvSlice<short>;
template class PvSlice<int>;
template class PvSlice<long>;
template class PvSlice<float>;
template class PvSlice<double>;
