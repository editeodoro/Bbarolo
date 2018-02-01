//----------------------------------------------------------
// mmaps.hh: Definition of the MomentMap class.
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

#ifndef MOMENTMAP_HH_
#define MOMENTMAP_HH_

#include <iostream>
#include <fitsio.h>
#include <vector>

#include <Arrays/cube.hh>
#include <Arrays/image.hh>
#include <Utilities/utils.hh>
#include <Utilities/progressbar.hh>


template <class T> 
class MomentMap : public Image2D <T> 
{
public:
    MomentMap();
    ~MomentMap() {}
    MomentMap(const MomentMap &i);      
    MomentMap& operator=(const MomentMap &i);

    void input (Cube<T> *c, int *Blo, int *Bhi);
    void input (Cube<T> *c);
    void SumMap (bool msk);
    void HIMassDensityMap (bool msk);
    void ZeroMoment (bool msk);
    void FirstMoment (bool msk);
    void SecondMoment (bool msk);
    void RMSMap (float level=0.1, float sncut = 1.5);

    bool setHead(int type); 
    
private:
    Cube<T> *in;
    int blo[3],bhi[3];
    int nsubs;
    
};


template <class T>
MomentMap<T>::MomentMap() {
    
    this->numPix = 0;
    this->numAxes = 2;
    this->arrayAllocated = false;
    this->headDefined      = false;
    this->statsDefined   = false;
    
}


template <class T> 
void MomentMap<T>::input (Cube<T> *c, int *Blo, int *Bhi) {
    
    in = c;
    
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
    
}


template <class T> 
void MomentMap<T>::input (Cube<T> *c) {
    
    in = c;
    
    for (int i=0; i<3; i++) {
        blo[i] = 0;
        bhi[i] = in->AxesDim(i);
    }
    
    this->axisDim[0] = bhi[0]-blo[0];
    this->axisDim[1] = bhi[1]-blo[1];
    nsubs = bhi[2]-blo[2]; 
    
    this->numPix = this->axisDim[0]*this->axisDim[1];
    this->array = new T[this->numPix];
    this->arrayAllocated = true;
    
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

    // Write a Mass density map in Msun/pc2. BUNIT must be JY/beam
    
    // Check that input units are jy / beam
    std::string bunit = makelower(in->Head().Bunit());
    bool isOK = bunit.find("jy/b")>=0 || bunit.find("j/b")>=0;
    if (isOK) {      
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
        for (int i=0; i<this->numPix; i++)
            this->array[i]=8794*this->array[i]*dvel/(bmaj*bmin);

        this->head.setBunit("Msun/pc2");
        this->head.setBtype("Mass_surfdens");
        
    }
    else {
        std::cerr << "MOMENT MAPS ERROR: Input datacube must be in JY/BEAM.\n";
        std::terminate();
    }
}

template <class T>
void MomentMap<T>::ZeroMoment (bool msk) {
    
    if (!this->arrayAllocated) {
        std::cout << "MOMENT MAPS error: ";
        std::cout << "Array not allocated. Call 'input' first!!\n";
        std::terminate();
    }
        
    if (!(this->headDefined=setHead(0)) && in->pars().isVerbose()) {
        std::cout<< "MOMENT MAPS warning: cannot create new header.\n";
    }
    
    bool isVerbose = in->pars().isVerbose();
    
    if(msk && !in->MaskAll()) in->BlankMask();
    
    float deltaV = 1;
    if (in->HeadDef()) deltaV = fabs(DeltaVel<T>(in->Head()));
    
    ProgressBar bar(" Extracting 0th moment map... ", true);
    bar.setShowbar(in->pars().getShowbar());
    if (isVerbose) bar.init(this->axisDim[0]);

    for (int x=0; x<this->axisDim[0]; x++) {
        if (isVerbose) bar.update(x+1);
        for (int y=0; y<this->axisDim[1]; y++) {
            float fluxsum = 0;
            this->array[x+y*this->axisDim[0]]=0;
            for (int z=0; z<nsubs; z++) {
                long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
                if (msk)fluxsum += in->Array(npix)*in->Mask(npix);
                else fluxsum += in->Array(npix);
            }
            if (fluxsum==0) fluxsum /= fluxsum;
            if (in->HeadDef())
                this->array[x+y*this->axisDim[0]] = FluxtoJy(fluxsum, in->Head())*deltaV;
            else this->array[x+y*this->axisDim[0]] = fluxsum;
        }
    }

    if (isVerbose) bar.fillSpace("Done.\n");
    
}


template <class T>
void MomentMap<T>::FirstMoment (bool msk) {
    
    
    if (!this->arrayAllocated) {
        std::cout << "MOMENT MAPS error: ";
        std::cout << "Array not allocated. Call 'input' first!!\n";
        std::terminate();
    }
        
    if (!(this->headDefined=setHead(1))) {
        std::cout<< "MOMENT MAPS warning: cannot create new header.\n";
    }
    
    bool isVerbose = in->pars().isVerbose();
    
    if(msk && !in->MaskAll()) in->BlankMask();
            
    ProgressBar bar(" Extracting 1st moment map... ", true);
    bar.setShowbar(in->pars().getShowbar());
    if (isVerbose) bar.init(this->axisDim[0]);
    
    
    for (int x=0; x<this->axisDim[0]; x++) {        
        if (isVerbose) bar.update(x+1);
        for (int y=0; y<this->axisDim[1]; y++) {
            T num = 0;
            T denom = 0;
            this->array[x+y*this->axisDim[0]]=0;
            for (int z=0; z<nsubs; z++) {
                long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
                T VEL;
                if (in->HeadDef()) {
                    double crpix2 = in->Head().Crpix(2);
                    double cdelt2 = in->Head().Cdelt(2);
                    double crval2 = in->Head().Crval(2);
                    T zval = (z+1+blo[2]-crpix2)*cdelt2+crval2;
                    VEL = AlltoVel<T>(zval, in->Head());
                }
                else VEL = z+blo[2];
                if (msk) {
                    num += in->Array(npix)*VEL*in->Mask(npix);
                    denom += in->Array(npix)*in->Mask(npix);
                }
                else {
                    num += in->Array(npix)*VEL;
                    denom += in->Array(npix);   
                }
            }
            this->array[x+y*this->axisDim[0]]=num/denom;    
        }
    }
    
    if (isVerbose) bar.fillSpace("Done.\n");

}


template <class T>
void MomentMap<T>::SecondMoment (bool msk) {
    
    if (!this->arrayAllocated) {
        std::cout << "MOMENT MAPS error: ";
        std::cout << "Array not allocated. Call 'input' first!!\n";
        std::terminate();
    }
        
    if (!(this->headDefined=setHead(2))) {
        std::cout<< "MOMENT MAPS warning: cannot create new header.\n";
    }
    
    bool isVerbose = in->pars().isVerbose();
    
    if(msk && !in->MaskAll()) in->BlankMask();

    ProgressBar bar(" Extracting 2nd moment map... ", true);
    bar.setShowbar(in->pars().getShowbar());
    if (isVerbose) bar.init(this->axisDim[0]);
    
    for (int x=0; x<this->axisDim[0]; x++) {        
        if (isVerbose) bar.update(x+1);
        for (int y=0; y<this->axisDim[1]; y++) {
            T num=0, denom=0;
            T numfrst=0, denomfrst=0;
            T firstmoment=0;
            this->array[x+y*this->axisDim[0]]=0;
            for (int z=0; z<nsubs; z++) {
                long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
                T VEL;
                if (in->HeadDef()) {
                    double crpix2 = in->Head().Crpix(2);
                    double cdelt2 = in->Head().Cdelt(2);
                    double crval2 = in->Head().Crval(2);
                    T zval = (z+1+blo[2]-crpix2)*cdelt2+crval2;
                    VEL = AlltoVel<T>(zval, in->Head());
                }
                else VEL = z+blo[2];
                if (msk) {  
                    if (in->Mask(npix)) {
                        numfrst += in->Array(npix)*VEL;
                        denomfrst += in->Array(npix);               
                    }
                }
                else {
                    numfrst += in->Array(npix)*VEL;
                    denomfrst += in->Array(npix);   
                }
            }
            if(denomfrst!=0) firstmoment = numfrst/denomfrst;
            else firstmoment = 0;
            for (int z=0; z<nsubs; z++) {
                long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
                T VEL;
                if (in->HeadDef()) {
                    double crpix2 = in->Head().Crpix(2);
                    double cdelt2 = in->Head().Cdelt(2);
                    double crval2 = in->Head().Crval(2);
                    T zval = (z+1+blo[2]-crpix2)*cdelt2+crval2;
                    VEL = AlltoVel<T>(zval, in->Head());
                }
                else VEL = z+blo[2];
                if (msk) {
                    if (in->Mask(npix)) { 
                        num += in->Array(npix)*(VEL-firstmoment)*(VEL-firstmoment);
                        denom += in->Array(npix);
                    }
                }
                else {
                    num += in->Array(npix)*(VEL-firstmoment)*(VEL-firstmoment);
                    denom += in->Array(npix);
                }
            }
            this->array[x+y*this->axisDim[0]]=sqrt(num/denom);          
        }
    }
    if (isVerbose) bar.fillSpace("Done.\n");
    
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
    
    bool isVerbose = in->pars().isVerbose();
    int nthreads = in->pars().getThreads();
    bool rob = in->pars().getFlagRobustStats();
    
    // Cube sizes
    size_t xs = in->DimX(), ys = in->DimY(), zs = in->DimZ();
    
    // Progress bar
    ProgressBar bar(" Computing RMS map... ", true);
    bar.setShowbar(in->pars().getShowbar());
    if (nthreads>1) bar.setShowbar(false);
    if (isVerbose) bar.init(xs*ys);
    
#pragma omp parallel for num_threads(nthreads)   
    for (int xy=0; xy<xs*ys; xy++) {
        if (isVerbose) bar.update(xy+1);
            
        // Getting the spectrum at x,y pixel
        std::vector<float> sp(zs);
        for (int z=0; z<zs; z++) 
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
    
    if (isVerbose) bar.fillSpace("Done.\n");


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
            if (in->Head().BeamArea()!=0) 
                this->head.setBunit("JY * KM/S");
            else this->head.setBunit("JY/BEAM * KM/S");
        }
        else if (type==1 || type==2) {
            if (type==1) this->head.setBtype("velocity");
            else this->head.setBtype("dispersion");
            this->head.setBunit("KM/S");
        }
        else if (type==3) {
            this->head.setBtype("rms");
        }
        this->headDefined = true;
    }   
    
    return this->headDefined;
    
}


template <class T>
Image2D<T>* PositionVelocity (Cube<T> *c, T x0, T y0, T Phi) {
    
    T phi = Phi;
    while(phi>=180) phi -= 180;
    while(phi<0) phi += 180;
    
    double P = phi*M_PI/180.;
    int dim[2] = {0, c->DimZ()}; 
    Image2D<T> *pv;
    int xdim=c->DimX(), ydim=c->DimY();
    int xmax=xdim, ymax=ydim;
    int xmin=0, ymin=0; 
    
    Header &h = c->Head();
    float cdelt0;

    if (phi==90) {
        dim[0] = xdim;
        pv  = new Image2D<T>(dim);
        for (int x=0; x<dim[0]; x++)
            for (int z=0; z<dim[1]; z++)
                pv->Array()[x+z*dim[0]] = c->Array(c->nPix(x,y0,z));
        cdelt0 = fabs(h.Cdelt(0));
    }
    else {
        std::vector<int> xx, yy; 
        double mx=0, my=0;
        
        if (phi<90) my = tan(P+M_PI_2);
        else my = tan(P-M_PI_2);

        mx = 1./my;
        
        int x_1 = lround(x0+(c->DimY()-y0)/my); 
        int x_0 = lround(x0-y0/my);
        xmax = my>0 ? x_1 : x_0;
        xmin = my>0 ? x_0 : x_1;
        if (xmax<xmin) std::swap(xmax,xmin);
        if (xmax>c->DimX()) xmax=c->DimX();
        if (xmin<0) xmin=0;
        xdim = fabs(xmax-xmin);
                
        int y_1 = lround(y0+(c->DimX()-x0)/mx); 
        int y_0 = lround(y0-x0/mx);
        ymax = mx>0 ? y_1 : y_0;
        ymin = mx>0 ? y_0 : y_1;
        if (ymax<ymin) std::swap(ymax,ymin);
        if (ymax>c->DimY()) ymax=c->DimY();
        if (ymin<0) ymin=0;
        ydim = fabs(ymax-ymin);
        
        if (xdim>=ydim) {
            int nxdim=0;
            for (int x=0; x<c->DimX(); x++) { 
                int y1 = lround(my*(x-x0)+y0);
                bool isin = y1>=0 && y1<c->DimY();  
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
            for (int y=0; y<c->DimY(); y++) { 
                int x1 = lround(mx*(y-y0)+x0);
                bool isin = x1>=0 && x1<c->DimX();  
                if (isin) {
                    xx.push_back(x1);
                    yy.push_back(y);
                    nxdim++;
                }
            }
            dim[0]=nxdim;
        }

        pv  = new Image2D<T>(dim);  
        for (int i=0; i<dim[0]; i++) { 
            for (int z=0; z<c->DimZ(); z++) {
                pv->Array()[i+z*dim[0]] = c->Array(c->nPix(xx[i],yy[i],z));     
            } 
        }

        float xdom = xdim*h.Cdelt(0);
        float ydom = ydim*h.Cdelt(1);
        cdelt0 = sqrt(xdom*xdom+ydom*ydom)/pv->DimX();
    }
    
    pv->copyHeader(h);
    float crpix0 = xdim>ydim ? x0-xmin : y0-ymin;
    if (phi==90) crpix0 = x0;
    pv->Head().setCrpix(0, crpix0+1);
    pv->Head().setCrval(0, 0);
    pv->Head().setCdelt(0, cdelt0);
    pv->Head().setCtype(0, "Offset");
    pv->Head().setCrpix(1, h.Crpix(2));
    pv->Head().setCdelt(1, h.Cdelt(2));
    pv->Head().setCrval(1, h.Crval(2));
    pv->Head().setCunit(1, h.Cunit(2));
    pv->Head().setCtype(1, h.Ctype(2));
    pv->Head().setMinMax(0.,0.);
    
    std::string name = c->Head().Name()+"_pv"+to_string(Phi);
    pv->Head().setName(name);

    return pv;
    
}


#endif
