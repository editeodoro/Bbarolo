//--------------------------------------------------------------------
// galmod.cpp: Member functions for the Galmod class
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
#include <cmath>
#include <vector>
#include <Tasks/galmod.hh>
#include <Arrays/cube.hh>
#include <Tasks/smooth3D.hh>
#include <Tasks/moment.hh>
#include <Utilities/utils.hh>


#define C  2.99792458E08            // Speed of light in M/S
#define nran 16777216

namespace Model {

template <class T>
void Galmod<T>::defaults() {
    
    subAllocated    = false;
    outDefined      = false;
    readytomod      = false;
    ringDefined     = false;
    modCalculated   = false;
    ltype           = 1;
    cmode           = 1;    
    crota2          = 0; 
    cdens           = 1.0;
    iseed           = -1;   
    freq0           = 0.1420405751786E10;
    velsys          = 2;    //Optical or radio definition of the velocities resp. 1 2.      
    
}
template void Galmod<float>::defaults();
template void Galmod<double>::defaults();


template <class T>
Galmod<T>::Galmod() {
    
    defaults();
}
template Galmod<float>::Galmod();
template Galmod<double>::Galmod();


template <class T>
Galmod<T>::~Galmod() {
    
    if (subAllocated) delete [] cd2i;
    if (outDefined) delete out;
    if (ringDefined) delete r;
    
}
template Galmod<float>::~Galmod();
template Galmod<double>::~Galmod();


template <class T>
Galmod<T>::Galmod(const Galmod<T> &g) {

    operator=(g);
}
template Galmod<float>::Galmod(const Galmod<float>&);
template Galmod<double>::Galmod(const Galmod<double>&);


template <class T>
Galmod<T>& Galmod<T>::operator=(const Galmod &g) {
    
    if(this==&g) return *this;
    
    this->in        = g.in;
    if (this->outDefined) delete out;
    this->outDefined = g.outDefined;
    if (this->outDefined) *this->out = *g.out;
    this->crpix3    = g.crpix3;
    this->crval3    = g.crval3;
    this->drval3    = g.drval3;
    this->cdelt3    = g.cdelt3;
    this->cunit3    = g.cunit3;
    this->ctype3    = g.ctype3;
    
    for(int i=0; i<2; i++) {
        this->ctype[i]  = g.ctype[i];
        this->cunit[i]  = g.cunit[i];
        this->cdelt[i]  = g.cdelt[i];
        this->blo[i]    = g.blo[i];
        this->bhi[i]    = g.bhi[i];
        this->bsize[i]  = g.bsize[i];
        this->pixsize[i]    = g.pixsize[i];
    } 
    
    this->pixarea   = g.pixarea;
    this->velsys    = g.velsys;
    this->axtyp     = g.axtyp;
    this->freq0     = g.freq0;
    this->crota2    = g.crota2;
    this->readytomod= g.readytomod;
    
    this->nsubs     = g.nsubs;
    if(this->subAllocated) delete [] cd2i;
    this->subAllocated = g.subAllocated;
    if (this->subAllocated) {
        this->cd2i  = new float[nsubs];
        for (int i=0; i<nsubs; i++) {
            this->cd2i[i]   = g.cd2i[i];
        }
    }
        
    this->ringDefined = g.ringDefined;
    if (ringDefined) this->r = g.r;
        
    this->ltype = g.ltype;
    this->cmode = g.cmode;
    this->cdens = g.cdens;
    this->iseed = g.iseed;
    this->arcmconv = g.arcmconv;
    this->modCalculated = g.modCalculated;

    
    return *this;
}
template Galmod<float>& Galmod<float>::operator=(const Galmod<float>&);
template Galmod<double>& Galmod<double>::operator=(const Galmod<double>&);


template <class T>
void Galmod<T>::input(Cube<T> *c, int *Boxup, int *Boxlow, Rings<T> *rings, 
                         int NV, int LTYPE, int CMODE, float CDENS, int ISEED) {
    
    /// This function sets all parameters needed to use the Galmod Object. A description
    /// of input is given in the Galmod class definition file (galmod.hh). 
    
        
    initialize(c, Boxup, Boxlow);
    
    ringIO(rings);
    
    ltype=LTYPE;
    cmode=CMODE;
    if (cmode<0 || cmode>2) {
        std::cout << "GALMOD warning: CMODE must be 0,1 or 2. Assuming 1.\n";
        cmode=1;
    }
    
    cdens = CDENS;
    if (cdens<0) {
        std::cout << "GALMOD warning: CDENS must greater than 0. Assuming 1.0.\n";
        cdens=1.;
    }
    
    int nvtmp = NV;
    if (nvtmp<1) {
        std::cout << "GALMOD warning: NV must greater than 0. Assuming "<< nsubs;
        std::cout << std::endl;
        nvtmp=nsubs;
    }
    for (int i=0; i<r->nr; i++) {
        if (r->vdisp[i]==0) r->nv.push_back(1);
        else r->nv.push_back(nvtmp);
    }
    
    iseed = ISEED;
    if (iseed>=0) {
        std::cout << "GALMOD warning: ISEED must be negative. Assuming -1.\n";
        iseed=-1;
    }
    
    readytomod=true;
}
template void Galmod<float>::input(Cube<float>*,int*,int*,Rings<float>*,int,int,int,float,int);
template void Galmod<double>::input(Cube<double>*,int*,int*,Rings<double>*,int,int,int,float,int);


template <class T>
void Galmod<T>::input(Cube<T> *c, Rings<T> *rings, int LTYPE) {
    
    /// Alternative semplified user interface function. 
    /// It sets the box to all the Cube and uses the default 
    /// values for cmode, cdens, nv and iseed.
    
    int Blo[2] = {0,0};
    int Bup[2] = {c->DimX(),c->DimY()};
    
    initialize(c, Bup, Blo);
    
    ringIO(rings);
    
    ltype = LTYPE;
    cmode = 1;
    cdens = 1;
    
    int nvtmp = nsubs;
    for (int i=0; i<r->nr; i++) {
        if (r->vdisp[i]==0) r->nv.push_back(1);
        else r->nv.push_back(nvtmp);
    }
    
    iseed = -1;
    
    readytomod=true;
    
}
template void Galmod<float>::input(Cube<float>*,Rings<float>*,int);
template void Galmod<double>::input(Cube<double>*,Rings<double>*,int);


template <class T>
void Galmod<T>::calculate() {
    
    /// Front end function to calculate the model.
    
    if (readytomod) {
        NHItoRAD();
        galmod();
    }
    else {
        std::cout<< "GALMOD error: wrong or unknown input parameter.\n";
        std::terminate();
    }

}
template void Galmod<float>::calculate();
template void Galmod<double>::calculate();



template <class T>
bool Galmod<T>::smooth(bool usescalefac) {
        
    if (!modCalculated) {
        std::cout << "GALMOD smoothing error: the model has not been ";
        std::cout << "calculated yet.\n";
        return false; 
    }
    
    //Beam oldbeam = {pixsize[0]*60, pixsize[1]*60, 0};
    Beam oldbeam = {0., 0., 0};

    if (oldbeam.bmaj<oldbeam.bmin) 
        std::swap(oldbeam.bmaj,oldbeam.bmin);
    
    Beam newbeam = {in->Head().Bmaj()*arcmconv*60,
                    in->Head().Bmin()*arcmconv*60,
                    in->Head().Bpa()};
    
    Smooth3D<T> *smoothed = new Smooth3D<T>;    
    smoothed->setUseScalefac(usescalefac);
    smoothed->setUseBlanks(false);
    if(!smoothed->smooth(out, oldbeam, newbeam, out->Array(), out->Array()))
        return false;   
    
    for (int i=0; i<out->NumPix(); i++)
        if (out->Array(i)<1.E-12) out->Array()[i] = 0;
    
    delete smoothed;
    
    return true;
}
template bool Galmod<float>::smooth(bool);
template bool Galmod<double>::smooth(bool);


template <class T>
void Galmod<T>::normalize() {
    
    if (!modCalculated) {
        std::cout << "GALMOD ERROR: the model has not been calculated yet.\n";
        return; 
    }
    
    bool verb = in->pars().isVerbose();
    if (verb) {
        in->pars().setVerbosity(false);
        out->pars().setVerbosity(false);
    }
    
    if (!in->StatsDef()) in->setCubeStats();
    in->stat().setThresholdSNR(3);
    T thres = in->stat().getThreshold();
    
    if (verb) 
        std::cout << " Normalizing model to data...";
        
    for (int x=0; x<bsize[0]; x++) { 
        for (int y=0; y<bsize[1]; y++) {
            T modSum = 0;
            T obsSum = 0;
            T factor = 0;   
            for (int z=0; z<nsubs; z++) {
                long modPix = out->nPix(x,y,z);
                long obsPix = in->nPix(x+blo[0],y+blo[1],z);
                modSum += out->Array(modPix);
                if (in->Array(obsPix)>thres) 
                    obsSum += in->Array(obsPix);
            }
            if (modSum!=0) factor = obsSum/modSum;
            else continue;
            for (int z=0; z<nsubs; z++) {
                long modPix = out->nPix(x,y,z);
                out->Array()[modPix] *= factor;
            }
        }
    }
        
    
    if (verb) {
        in->pars().setVerbosity(true);
        out->pars().setVerbosity(true);
    }
    
}
template void Galmod<float>::normalize();
template void Galmod<double>::normalize();


template <class T>
void Galmod<T>::initialize(Cube<T> *c, int *Boxup, int *Boxlow) {

    in = c;
    nsubs  = c->DimZ();
    
    if (!subAllocated) cd2i = new float[nsubs];
    subAllocated=true;
    
    /// Information about RA-DEC and conversions.
    arcmconv=0;
    for (int i=0; i<2; i++) {                        
        bhi[i] = Boxup[i];                  
        blo[i] = Boxlow[i];
        ctype[i] = c->Head().Ctype(i);
        cunit[i] = c->Head().Cunit(i);
        cdelt[i] = c->Head().Cdelt(i);

        if (cunit[i]=="DEGREE" || cunit[i]=="DEGREES" || cunit[i]=="DEG" ||
            cunit[i]=="degree" || cunit[i]=="degrees" || cunit[i]=="deg") 
                arcmconv = 60;
        else if (cunit[i]=="ARCSEC" || cunit[i]=="ARCS" ||
                 cunit[i]=="arcsec" || cunit[i]=="arcs") 
                arcmconv = 1/60;
        else if (cunit[i]=="ARCMIN" || cunit[i]=="ARCM" ||
                 cunit[i]=="arcmin" || cunit[i]=="arcm") 
                arcmconv = 1;
        else {
            std::cout << "GALMOD error (unknown CUNIT for RA-DEC): ";
            std::cout << "cannot convert to ARCMIN.\n";
            std::cout << cunit[i];
            std::terminate(); 
        }
        cdelt[i]   = cdelt[i]*arcmconv;
        pixsize[i] = fabs(cdelt[i]); 
        bsize[i]   = bhi[i]-blo[i]; 
    }
    pixarea=pixsize[0]*pixsize[1];
    
    crota2=in->Head().Crota();
    crota2=crota2*M_PI/180.;
    
    int ax[3]={bsize[0],bsize[1],nsubs};
    out = new Cube<T>(ax);
    out->saveHead(c->Head());
    out->saveParam(c->pars());
    out->Head().setDimAx(0, bsize[0]);
    out->Head().setDimAx(1, bsize[1]);
    out->Head().setDimAx(2, nsubs);
    out->Head().setCrpix(0, c->Head().Crpix(0)-blo[0]);
    out->Head().setCrpix(1, c->Head().Crpix(1)-blo[1]);
    outDefined = true;
    
    /// Information about frequency/velocity axis and conversions.
    freq0 = c->Head().Freq0();      
    if (freq0==0) {
        freq0 = 0.1420405751786E10;
        std::cout << "Header item FREQ0 not found. Assuming " << freq0;
        std::cout << std::endl; 
    }                                                           
    ctype3 = makelower(c->Head().Ctype(2));
    cunit3 = makelower(c->Head().Cunit(2));
    cdelt3 = c->Head().Cdelt(2);
    crval3 = c->Head().Crval(2);
    crpix3 = c->Head().Crpix(2);
    drval3 = c->Head().Drval3();
    
    if (ctype3=="wave" || ctype3=="awav" || cunit3=="um" || cunit3=="nm" || cunit3=="ang" ) axtyp =2;
    else if (ctype3=="freq" || cunit3=="hz" || cunit3=="mhz") axtyp = 3;
    else if (ctype3=="velo" || ctype3=="velo-helo" || ctype3=="velo-hel" ||
             cunit3=="m/s" || cunit3=="km/s") {
        axtyp = 4;
    }
    else {
        std::cout << "GALMOD error (unknown CUNIT for spectral axis): cannot convert.";
        std::terminate(); 
    }
    
    if (axtyp==2) {
        float mconv=0;
        if (cunit3=="um"||cunit3=="mum"||cunit3=="micron") mconv = 1.E-06;
        else if (cunit3=="nm"||cunit3=="nanom") mconv = 1.0E-09;
        else if (cunit3=="a" ||cunit3=="ang"||cunit3=="angstrom") mconv = 1.0E-10;
        else {
            std::cout << "GALMOD error (unknown CUNIT3): cannot convert to M.\n";
            std::cout << cunit3;
            std::terminate();
        }

        crval3=crval3*mconv;
        cdelt3=cdelt3*mconv;

        double crvalfreq = C/crval3;

        // If redshift and wavelength parameters are set, take them for freq0, otherwise central channel
        double restw = in->pars().getRestwave(), reds = in->pars().getRedshift();
        if (restw!=-1 && reds!=-1) freq0 = C/(restw*(1+reds)*mconv);
        else freq0 = crvalfreq;                      // Velocity is 0 always at the reference channel

        drval3 = C*(freq0*freq0-crvalfreq*crvalfreq)/(freq0*freq0+crvalfreq*crvalfreq);

    }
    else if (axtyp==3) {
        
        float hzconv=0;
        if (cunit3=="hz") hzconv = 1;
        else if (cunit3=="khz") hzconv = 1.0E03;
        else if (cunit3=="mhz") hzconv = 1.0E06;
        else if (cunit3=="ghz") hzconv = 1.0E09;
        else {
            std::cout << "GALMOD error (unknown CUNIT3): cannot convert to Hz.\n";
            std::cout << cunit3;
            std::terminate(); 
        }
        
        crval3=crval3*hzconv;
        cdelt3=cdelt3*hzconv;  
        
        double crvalfreq = c->Head().Crval(2)*hzconv;
        drval3 = C*(freq0*freq0-crvalfreq*crvalfreq)/(freq0*freq0+crvalfreq*crvalfreq);
    }
    else if (axtyp==4) {
        
        float msconv=0;
        if (cunit3=="m/s" || cunit3=="ms") msconv = 1;
        else if (cunit3=="km/s") msconv = 1.0E03;
        else if ("cm/s") msconv = 1.0E-03;
        else {
            std::cout << "GALMOD error (unknown CUNIT3): cannot convert to M/S.\n";
            std::cout << cunit3;
            std::terminate(); 
        }
        
        crval3=crval3*msconv;
        cdelt3=cdelt3*msconv;
        double crvalvel = c->Head().Crval(2)*msconv;
        drval3=freq0*sqrt((C-crvalvel)/(C+crvalvel));
        //std::cout << "Give frequency at reference grid in HZ. DRVAL3= ";
        //std::cin >> drval3;
    }
    else { 
        std::cout << "Unknown axis type: no velocities along spectral axis.\n";
        std::terminate();
    }
 
    //std::string bunit = c->Head().Bunit();
    //if (bunit=="NONE") 
        //std::cout << "GALMOD warning: No units in the maps found.\n";
    //if (bunit!="W.U." && bunit!="w.u.")
    //  std::cout << "GALMOD warning: Units in the maps should be W.U..\n";
    
    // Get the instrumental broadnening: when Hanning smoothing has been applied, 
    // the instrumental profile can be approximated by a gaussian with a FWHM of 
    // twice the channnel separation. Thus, the sig_instr is FWHM/2.355.
    float nch = in->pars().getLinear();
    if (nch==-1) nch=2./(2*sqrt(2*log(2)));
    chwidth = fabs(DeltaVel<double>(in->Head()))*1000;
    sig_instr = nch*chwidth;
    
}
template void Galmod<float>::initialize(Cube<float>*,int*,int*);
template void Galmod<double>::initialize(Cube<double>*,int*,int*);


template <class T>
void Galmod<T>::ringIO(Rings<T> *rings) {
    
    /// This function creates the set of rings to be used in the model.
    /// - Velocities are given in km/s and then converted in m/s.
    /// - Radii are given in arcsec and then converted in arcmin.
    /// - Column densities are given in HI atoms/cm^2 and the converted 
    ///   in units of 1E20 atoms/cm^2.
    /// - Angles are given in grades and then converted in radians.
    
    r = new Rings<T>;

    int nur = rings->nr;
    T *uradii = new T[nur];
    T *uvrot  = new T[nur];
    T *uvdisp = new T[nur];
    T *uvrad  = new T[nur];
    T *uvvert = new T[nur];
    T *udvdz  = new T[nur];
    T *uzcyl  = new T[nur];
    T *udens  = new T[nur];
    T *uz0    = new T[nur];
    T *uinc   = new T[nur];
    T *uphi   = new T[nur];
    T *uxpos  = new T[nur];
    T *uypos  = new T[nur];
    T *uvsys  = new T[nur];
    
    r->radsep=0.75*min(pixsize[0],pixsize[1]);
    uradii[0]=rings->radii[0]/60.;
    for (int i=1; i<nur; i++) {
        uradii[i]=rings->radii[i]/60.;
        if (uradii[i-1]+r->radsep>=uradii[i]) {
            if (uradii[i-1]>uradii[i]) {
                std::cout << "GALMOD error: Radii not in increasing order.\n";
                std::terminate();
            }
            else {
                std::cout << "GALMOD error: Radius separation too small.\n";
                std::terminate();
            }
        }
    }
    
    bool empty = true;
    /*
    if (uradii[0]!=0) {
        std::cout << "Leave innerpart of first ring empty? [Y,N]";
        std::string a;
        cin >> a;
        if (a=="N" || a=="n" || a=="NO" || a=="no") empty=true;
    }
    */

    for (int i=0; i<nur; i++) {
        uvrot[i]  = rings->vrot[i]*1000;   
        uvrad[i]  = rings->vrad[i]*1000;   
        uvvert[i] = rings->vvert[i]*1000;
        
        uvdisp[i] = rings->vdisp[i]*1000;
        if (uvdisp[i]<0) {
            std::cout << "GALMOD error: Negative velocity dispersion not allowed.\n";
            std::terminate();
        }
        
        udens[i]=rings->dens[i]/1.0E20;

        if (udens[i]<0) {
            std::cout << "GALMOD error: Negative column-density not allowed.\n";
            std::terminate(); 
        }
        
        uz0[i]=rings->z0[i]/60.;
        if (uz0[i]<0) {
            std::cout << "GALMOD error: Negative scale height not allowed.\n";
            std::terminate(); 
        }
        
        udvdz[i]  = rings->dvdz[i]*1000*60; // m/s/arcmin
        uzcyl[i]  = rings->zcyl[i]/60.;
        uinc[i]   = rings->inc[i]*M_PI/180.;
        uphi[i]   = rings->phi[i]*M_PI/180.;
        uxpos[i]  = rings->xpos[i]+1;
        uypos[i]  = rings->ypos[i]+1;
        uvsys[i]  = rings->vsys[i]*1000;
    }
    
    if (uradii[0]!=0 && empty) r->radii.push_back(uradii[0]+r->radsep/2.0);
    else {
        r->radii.push_back(r->radsep/2.0);
        while (r->radii.back()<uradii[0]) {
            r->vrot.push_back(uvrot[0]);
            r->vrad.push_back(uvrad[0]);
            r->vvert.push_back(uvvert[0]);
            r->vdisp.push_back(uvdisp[0]);
            r->dens.push_back(udens[0]);
            r->z0.push_back(uz0[0]);
            r->dvdz.push_back(udvdz[0]);
            r->zcyl.push_back(uzcyl[0]);
            r->inc.push_back(uinc[0]);
            r->phi.push_back(uphi[0]);
            r->xpos.push_back(uxpos[0]);
            r->ypos.push_back(uypos[0]);
            r->vsys.push_back(uvsys[0]);
            r->radii.push_back(r->radii.back()+r->radsep);
        }
    }
    
    for (int i=1; i<nur; i++) {
        T dur       = uradii[i]-uradii[i-1];
        T dvrotdr   = (uvrot[i]-uvrot[i-1])/dur;
        T dvraddr   = (uvrad[i]-uvrad[i-1])/dur;
        T dvvertdr  = (uvvert[i]-uvvert[i-1])/dur;
        T dvdispdr  = (uvdisp[i]-uvdisp[i-1])/dur;
        T ddensdr   = (udens[i]-udens[i-1])/dur;
        T dz0dr     = (uz0[i]-uz0[i-1])/dur;
        T ddvdzdr   = (udvdz[i]-udvdz[i-1])/dur;
        T dzcyldr   = (uzcyl[i]-uzcyl[i-1])/dur;
        T dincdr    = (uinc[i]-uinc[i-1])/dur;
        T dphidr    = (uphi[i]-uphi[i-1])/dur;
        T dxposdr   = (uxpos[i]-uxpos[i-1])/dur;
        T dyposdr   = (uypos[i]-uypos[i-1])/dur;
        T dvsysdr   = (uvsys[i]-uvsys[i-1])/dur;
        while (r->radii.back()<uradii[i]) {
            T dr = r->radii.back()-uradii[i-1];
            r->vrot.push_back(uvrot[i-1]+dvrotdr*dr);
            r->vrad.push_back(uvrad[i-1]+dvraddr*dr);
            r->vvert.push_back(uvvert[i-1]+dvvertdr*dr);
            r->vdisp.push_back(uvdisp[i-1]+dvdispdr*dr);
            r->dens.push_back(udens[i-1]+ddensdr*dr);
            r->z0.push_back(uz0[i-1]+dz0dr*dr);
            r->dvdz.push_back(udvdz[i-1]+ddvdzdr*dr);
            r->zcyl.push_back(uzcyl[i-1]+dzcyldr*dr);
            r->inc.push_back(uinc[i-1]+dincdr*dr);
            r->phi.push_back(uphi[i-1]+dphidr*dr);
            r->pa.push_back(r->phi.back()+crota2);
            r->xpos.push_back(uxpos[i-1]+dxposdr*dr);
            r->ypos.push_back(uypos[i-1]+dyposdr*dr);
            r->vsys.push_back(uvsys[i-1]+dvsysdr*dr);
            r->radii.push_back(r->radii.back()+r->radsep);
        }
    }
    
    r->radii.pop_back();
    r->nr=r->radii.size();      
    
    ringDefined = true;
        
    delete [] uradii;
    delete [] uvrot;
    delete [] uvrad;
    delete [] uvvert;
    delete [] uvdisp;
    delete [] udens;
    delete [] uz0;
    delete [] udvdz;
    delete [] uzcyl;
    delete [] uinc;
    delete [] uphi;
    delete [] uxpos;
    delete [] uypos;
    delete [] uvsys;
}
template void Galmod<float>::ringIO(Rings<float>*);
template void Galmod<double>::ringIO(Rings<double>*);


template <class T>
void Galmod<T>::galmod() {
    
    const double twopi = 2*M_PI;
    const int buflen=bsize[0]*bsize[1]*nsubs+1;
    T *datbuf = new T [buflen];
            
    bool verb = in->pars().isVerbose();
    ProgressBar bar(" Modeling... ",false);;
    bar.setShowbar(in->pars().getShowbar());
    if (verb) bar.init(r->nr);
    
//  Reinitialize random number generator for each new model.
    int isd = iseed;
//  Get number of velocity profiles that will be done.
    int nprof = bsize[0]*bsize[1];
//  Initialize data buffer on zero.
    for (int i=0; i<buflen; i++) datbuf[i]=0.0;

    // ==>> Loop over standard rings.
    for (int ir=0; ir<r->nr; ir++) {
        if (verb) bar.update(ir+1);
//      Get radius
        double rtmp = r->radii[ir];
//      Get number of clouds inside ring.
        int nc = lround(cdens*pow(r->dens[ir],cmode)*twopi*rtmp*r->radsep/pixarea); 
        if (nc==0) {
            std::cerr << " GALMOD ERROR: No clouds used. Choose higher CDENS " << std::endl;
            std::terminate();
//          Do next ring, jump to end of loop for rings.
            continue;
        }
//      Get values of ir-radius.            
        double vrottmp = r->vrot[ir];
        double vradtmp = r->vrad[ir];
        double vverttmp= r->vvert[ir];
        //  The VDISP should be such that VDISP^2 + chwidth^2 / 12 = sig_instr^2 + sig_v^2
        //double vdisptmp= sqrt(r->vdisp[ir]*r->vdisp[ir]+sig_instr*sig_instr-(chwidth*chwidth)/12.);
        double vdisptmp= sqrt(r->vdisp[ir]*r->vdisp[ir]+sig_instr*sig_instr);
        double vsystmp = r->vsys[ir];
        double z0tmp   = r->z0[ir];
        double dvdztmp = r->dvdz[ir];
        double zcyltmp = r->zcyl[ir];
        double sinc    = sin(r->inc[ir]);
        double cinc    = cos(r->inc[ir]);
        double spa     = sin(r->phi[ir])*cos(crota2)+cos(r->phi[ir])*sin(crota2);
        double cpa     = cos(r->phi[ir])*cos(crota2)-sin(r->phi[ir])*sin(crota2);
        double nvtmp   = r->nv[ir];
        float  fluxsc  = r->dens[ir]*twopi*rtmp*r->radsep/(nc*nvtmp);
            
// ==>> Loop over clouds inside each ring.
        for (int ic=0; ic<nc; ic++) {
//          Get radius inside ring. The range includes the inner boundary,
//          excludes the outer boundary. The probability of a radius inside
//          a ring is proportional to the total radius and thus the 
//          surface density of the clouds is constant over the area of the ring.
            double ddum = double(iran(isd))/double(nran); 
            double R    = sqrt(pow((rtmp-0.5*r->radsep),2)+2*r->radsep*rtmp*ddum);
//          Get azimuth and its sine and cosine.
            double az   = twopi*double(iran(isd))/double(nran); 
            double saz  = sin(az); 
            double caz  = cos(az); 
//          Get height above the plane of the ring using a random deviate
//          drawn from density profile of the laye.
            double z    = fdev(isd)*z0tmp;
//          Get position in the plane of the sky with respect to the major
//          and minor axes of the spiral galaxy.
            double x    = R*caz;
            double y    = R*saz*cinc-z*sinc;                                                 
//          Get grid of this position, check if it is inside area of box.
            long grid[2] = {lround(r->xpos[ir]+(x*spa-y*cpa)/cdelt[0]),
                            lround(r->ypos[ir]+(x*cpa+y*spa)/cdelt[1])};             
            if (grid[0]<=blo[0] || grid[0]>bhi[0]) continue;
            if (grid[1]<=blo[1] || grid[1]>bhi[1]) continue;

/*
            /////////// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            float rr, theta;
            float xr = (-(grid[0]-r->pos[0])*spa+(grid[1]-r->pos[1])*cpa);
            float yr = (-(grid[0]-r->pos[0])*cpa-(grid[1]-r->pos[1])*spa)/cinc;
            rr = sqrt(xr*xr+yr*yr);
            if (rr<0.1) theta = 0.0;            
            else theta = atan2(yr, xr)/M_PI*180;    
                
            int side = 1;
            bool use;
            switch (side) {                     // Which side of galaxy.                    
                case 1:                         //< Receding half.                              
                    use = (fabs(theta)<=90.0);      
                    break;
                case 2:                         //< Approaching half. 
                    use = (fabs(theta)>=90.0);
                    break;
                case 3:                         //< Both halves.
                    use = 1;
                    break;
                default: 
                    break;  
            }
                
            if (!use) continue;
*/                  
            //////////// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
//              
//          Get profile number of current pixel and check if position is in
//          range of positions of profiles that are currently being done.
//          If outside any of the ranges, jump to next cloud.
            int iprof = (grid[1]-blo[1])*bsize[0]+grid[0]-blo[0]-bsize[0];  
//          Get systematic velocity of cloud.
            double vsys = vsystmp+(vrottmp*caz+vradtmp*saz)*sinc;
//          Adding vertical rotational gradient after zcyl
            if (abs(z)>zcyltmp) vsys = vsystmp+((vrottmp-dvdztmp*(abs(z)-zcyltmp))*caz+vradtmp*saz)*sinc;
//          Adding vertical velocity
            if (z>0.) vsys += vverttmp*cinc;
            else vsys += vverttmp*cinc;     //// <--- Check this. In galmod is a +
                
                    
/* ORIGINAL BUILDING PROFILES
// ==>>     Build velocity profile.
            for (int iv=0; iv<nvtmp; iv++) {
//              Get deviate drawn from gaussian velocity profile and add
//              to the systematic velocity.
                double v     = vsys+gasdev(isd)*vdisptmp;
//              Get grid of velocity along FREQ-OHEL or VELO axis.
//              If a grid is not in the range, jump to next velocity profile.
                int isubs = lround(velgrid(v)+crpix3-1);
                if (isubs<0 || isubs>=nsubs) continue;                  
                int idat  = iprof+isubs*nprof;
//              Convert HI atom flux per pixel to flux per pixel of 21cm
//              radiation expressed in W.U. and add subcloud to the data
//              buffer.
                datbuf[idat] = datbuf[idat]+fluxsc*cd2i[isubs];
            }

*/

//          PARTE PER I DOPPIETTI CHE SOSTITUISCE IL BUILDING PROFILES DI SOPRA
            uint nlines = in->pars().getParGM().NLINES;
            float relvel_lines[2] = {0,220000};
            float relint_lines[2] = {1,1.70};

            for (int iv=0; iv<nvtmp; iv++) {
                double vdev = gasdev(isd)*vdisptmp;
                for (int nl=0; nl<nlines; nl++) {
                    double v     = vsys+vdev+relvel_lines[nl];
                    int isubs = lround(velgrid(v)+crpix3-1);
                    if (isubs<0 || isubs>=nsubs) continue;
                    int idat  = iprof+isubs*nprof;
                    datbuf[idat] = datbuf[idat]+relint_lines[nl]*fluxsc*cd2i[isubs];
                }
            }
        }
    }
//  Write data to output Cube.      
    for (int isubs=0; isubs<nsubs; isubs++) {
        int pixstart=isubs*bsize[0]*bsize[1];
        int idat=1+isubs*nprof;
        for (int i=0; i<nprof; i++)
            out->Array(pixstart+i)=datbuf[idat+i];
    }
    
    if (verb) bar.fillSpace("OK.\n");
    
    delete [] datbuf;

    modCalculated=true;

}
template void Galmod<float>::galmod();
template void Galmod<double>::galmod();
    

template <class T>
void Galmod<T>::NHItoRAD(){
    
/// Get conversion from column density of the HI and intensity of 21cm
/// radiation. (See Mihalas and Binney 'Galactic Astronomy'. p.489 or
/// Binney & Merrifield Chapter 8.) 
/// The relation between column density and brightness temperature 
/// in cgs is (see Roberts75):
/// 
///   N_HI = 1.823E13 * integ{Tb(v)*dv}      integral in velocity
///   N_HI = 3.848E14 * integ{Tb(nu)*dnu}    integral in frequency
/// 
/// where the brightness temperature is proportional to specific intensity:
///
///   I = (2*k/lamda**2)*Tb = (2*k*nu**2/c**2)*Tb
///    
/// where k is the Boltzmann constant and c the speed of light.
/// The intensity profile is assumed to be a delta peak since all atoms 
/// are on one frequency range inside the subset. Thus:
///    
///  I = (2*k/lamda**2) / (cdelt3*fac) * N_HI
///
/// where fac = 1.823E13 or 3.848E14.
/// N.B.:
/// The column density of the HI is in units of 1E20 atoms per cm2.
/// Velocities are in M/S and frequencies in HZ. Solid angle is in steradians.
/// Intensity is flux per one steradian solid angle. It is defined
/// as the flux per beam of one square arcminute in W.U..

    double const wu = 5E-29;             // Conversion factor to W.U.
    double const K  = 1.3806488E-23;     // Boltzmann constant in Joule/Kelvin 
    
    for (int isubs=0; isubs<nsubs; isubs++) {
        double labsubs=0, fac=0;

        if (axtyp==2) {    // Z-axis is wavelength
            labsubs = crval3+cdelt3*(isubs+1-crpix3);
            fac = 3.8475E-06;     // This is wrong for wavelength
        }
        else if (axtyp==3) {    // Z-axis is frequency
            double fsubs = crval3+cdelt3*(isubs+1-crpix3);
            labsubs = C/fsubs;
            fac = 3.8475E-06;     // This is 3.848E14 / 1E20
        }
        else if (axtyp==4) {     // Z-axis is velocity
            double vsubs=crval3+cdelt3*(isubs+1-crpix3);
            double fsubs=freq0*C*drval3/(drval3*(vsubs-crval3)+freq0*C);
            labsubs=C/fsubs;
            fac = 1.823E-05;       // This is 1.823E13 / 1E20 * 1E02
        } 

        double ddum = 2*(K/wu)/(fac*(labsubs*labsubs))/fabs(cdelt3);
        // Following line removes dependency from lambda
        //double ddum=2*(K/wu)/fac/fabs(cdelt3);
        // The solid angle in the definition of the intensity is 
        // in steradians, thus convert it to square arcminutes.
        cd2i[isubs]=ddum*pow((M_PI/(180.0*60.0)),2);
    }

}
template void Galmod<float>::NHItoRAD();
template void Galmod<double>::NHItoRAD();


template <class T>  
double Galmod<T>::velgrid(double v) {
    
    /// Function to transform a velocity to a grid.
        
    double velg=0, fdv;
    
    if (axtyp==2) {                     //< Wavelength axis
        double f = freq0*sqrt((C-v)/(C+v));
        double l = C/f;
        velg = (l-crval3)/cdelt3;
    }
    else if (axtyp==3) {                        //< Frequency axis.
        if (velsys==1) {                //< Optical definition of velocities.
            fdv  = (drval3-v)*crval3;
            velg = crval3*fdv/((freq0*C-fdv)*cdelt3);
        }
        else if (velsys==2) velg = (drval3-v)*freq0/(C*cdelt3); //< Radio definition
    }
    else if (axtyp==4) velg = (v-crval3)/cdelt3;
    else {
        std::cout << "Invalid frequency/velocity system." << std::endl;
    }
   
    return velg;

}
template double Galmod<float>::velgrid(double);
template double Galmod<double>::velgrid(double);


template <class T>
double Galmod<T>::fdev(int& idum){
    
    /// Function to get random deviates for various functions.
    /// The double precision variable Fdev contains the random deviate.
    /// If nran is within a factor two of the largest possible integer
    /// then integer overflow could occur.

    double Fdev=0, x=0;
    
    if (ltype==1) {                                 /// Gaussian function: exp(-0.5*x^2) 
        Fdev = gasdev(idum);
    }
    else {
        x = double(1+2*iran(idum)-nran)/double(nran);
        if (ltype==2) Fdev = atanh(x);              /// Sech2 function: sech2(x)
        else if (ltype==3) {                        /// Exponential function: exp(-|x|)
           if (x>=0.0) Fdev = - log(x);
           if (x<0.0)  Fdev = log(-x);
        }
        else if (ltype==4) Fdev = tan(M_PI_2*x);    /// Lorentzian function: 1/(1+x**2)
        else if (ltype==5) Fdev = x;                /// Box function.
        else cout << "Unknown function." << endl;
    }

    return Fdev;
}
template double Galmod<float>::fdev(int&);
template double Galmod<double>::fdev(int&);
    

template <class T>
double Galmod<T>::gasdev(int &idum) {

    /// Function to get random deviates from a gaussian distribution.
    /// Drawn from 'Numerical Recipes' but modified a bit.

    static double v1, v2, r;
    static double fac, gset;
    static int iset=0;
    double Gasdev=0;
    
    if (iset==0) {
        int ctrl = 0;
        while (ctrl==0) {
            v1 = double(1+2*iran(idum)-nran)/double(nran);
            v2 = double(1+2*iran(idum)-nran)/double(nran);
            r=v1*v1+v2*v2;
            if (r>=1 || r==0) ctrl=0;
            else ctrl=1;
        }
        fac  = sqrt(-2.0*log(r)/r);
        gset = v1*fac;
        Gasdev = v2*fac;
        iset = 1;
    }
    else { 
        Gasdev=gset;
        iset=0;
    }
    
    return Gasdev;

}
template double Galmod<float>::gasdev(int&);
template double Galmod<double>::gasdev(int&);

///*
template <class T>
int Galmod<T>::iran(int &idum) {
    
    /// Random number generator due to Knuth.
    /// The choices for MBIG and MSEED are not particularly important 
    /// since the they are only used 55 times in a linear congruential 
    /// algorithm to initialize the array ma.
    
        
    const int MBIG = 16777216;
    const int MSEED= 1618033;

    static int ma[56];
    static int inext, inextp, mk;
    int Iran=0; 

    if (idum<0) {
        Iran = MSEED-fabs(idum);
        Iran = Iran % MBIG;
        ma[55]=Iran;
        mk=1;
        for (int i=1; i<=54; i++) {
            int ii = (21*i) % 55;
            ma[ii] = mk;
            mk = Iran-mk;
            if (mk<0) mk = mk+MBIG;
             Iran = ma[ii];
        }
        for (int kk=1; kk<=4; kk++) {
            for (int i=1; i<=55; i++) {
                int j= (i+30) % 55;             
                ma[i]=ma[i]-ma[1+j];
                if (ma[i]<0) ma[i] = ma[i]+MBIG;
            }
        }

        inext=0;
        inextp=31;
        idum=1;
    }

    inext=inext+1;
    if (inext==56) inext=1;
    inextp=inextp+1;
    if (inextp==56) inextp=1;
    Iran = ma[inext]-ma[inextp];
    if (Iran<0) Iran=Iran+MBIG;
    ma[inext] = Iran;
    
    return Iran;
}
template int Galmod<float>::iran(int&);
template int Galmod<double>::iran(int&);
//*/

/*
template <class T>
int Galmod<T>::iran(int &idum) {
    
    /// Random number generator due to Knuth.
    /// The choices for MBIG and MSEED are not particularly important 
    /// since the they are only used 55 times in a linear congruential 
    /// algorithm to initialize the array ma.
    
        
    const int MBIG = 16777216;
    const int MSEED= 1618033;

    static int ma[55];
    static int inext, inextp, mk;
    int Iran=0; 

    if (idum<0) {
        Iran = MSEED-fabs(idum);
        Iran = Iran % MBIG;
        ma[54]=Iran;
        mk=1;
        for (int i=0; i<54; i++) {
            int ii = (21*(i+1)) % 55;
            ma[ii] = mk;
            mk = Iran-mk;
            if (mk<0) mk = mk+MBIG;
             Iran = ma[ii];
        }
        for (int kk=1; kk<=4; kk++) {
            for (int i=0; i<55; i++) {
                int j= (i+1+30) % 55;             
                ma[i]=ma[i]-ma[1+j];
                if (ma[i]<0) ma[i] = ma[i]+MBIG;
            }
        }

        inext=0;
        inextp=31;
        idum=1;
    }


       
    inext=inext+1;
    if (inext==55) inext=1;
    inextp=inextp+1;
    if (inextp==55) inextp=1;
    Iran = ma[inext]-ma[inextp];
    if (Iran<0) Iran=Iran+MBIG;
    ma[inext] = Iran;
    
    return Iran;
}
template int Galmod<float>::iran(int&);
template int Galmod<double>::iran(int&);
*/

}

#undef C  
#undef nran 
