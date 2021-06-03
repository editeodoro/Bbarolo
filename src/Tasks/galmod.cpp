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

// NB: Random number generators supported:
// 1) GALMOD classic: iran(), gasdev() -> Uncomment lines with "Classic galmod"
// 2) STD c++ library -> Uncomment lines with "STD library"
// There are 5 lines to change: 3 in galmod() and 2 in fdev()

#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <functional>
#include <Tasks/galmod.hh>
#include <Arrays/cube.hh>
#include <Tasks/smooth3D.hh>
#include <Tasks/moment.hh>
#include <Utilities/utils.hh>
#include <Utilities/progressbar.hh>

#include <sys/socket.h>
#include <math.h>
#include <stdlib.h>
#define C  2.99792458E08            // Speed of light in M/S
#define nran 16777216

int iran(int &idum) {
    
    /// Random number generator due to Knuth.
    /// The choices for MBIG and MSEED are not particularly important 
    /// since the they are only used 55 times in a linear congruential 
    /// algorithm to initialize the array ma.
    
        
    const int MBIG = 16777216;
    const int MSEED= 1618033;

    static int ma[56];
    static int inext, inextp, mk;
#pragma omp threadprivate(ma,inext,inextp,mk)

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
            if (mk<0) mk += MBIG;
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


double gasdev(int &idum) {

    /// Function to get random deviates from a gaussian distribution.
    /// Drawn from 'Numerical Recipes' but modified a bit.

    static double v1, v2, r;
    static double fac, gset;
    static int iset=0;
#pragma omp threadprivate(v1,v2,r,fac,gset,iset)
    
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


namespace Model {

template <class T>
void Galmod<T>::defaults() {
    
    subAllocated  = false;
    outDefined    = false;
    readytomod    = false;
    ringDefined   = false;
    modCalculated = false;
    ltype         = 1;
    cmode         = 1; 
    iseed           = -1;      
    crota2        = 0; 
    cdens         = 1.0;
    freq0         = 0.1420405751786E10;
    velsys        = 2;    //Optical or radio definition of the velocities resp. 1 2.
    nlines        = 1;
    relvel.push_back(0);
    relint.push_back(1);
           
}


template <class T>
Galmod<T>::Galmod() {
    
    defaults();
}


template <class T>
Galmod<T>::~Galmod() {
    
    if (subAllocated) delete [] cd2i;
    if (outDefined) delete out;
    if (ringDefined) delete r;
}


template <class T>
Galmod<T>::Galmod(const Galmod<T> &g) {

    operator=(g);
}


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
    this->nlines    = g.nlines;
    this->relvel    = g.relvel;
    this->relint    = g.relint;
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


template <class T>
void Galmod<T>::input(Cube<T> *c, int *Boxup, int *Boxlow, Rings<T> *rings, 
                         int NV, int LTYPE, int CMODE, float CDENS, int ISEED) {
    
    /// This function sets all parameters needed to use the Galmod Object. A description
    /// of input is given in the Galmod class definition file (galmod.hh). 

    initialize(c, Boxup, Boxlow);

    ringIO(rings);

    setOptions(LTYPE, CMODE, CDENS, ISEED);
    
    int nvtmp = NV;
    if (nvtmp<1) nvtmp=nsubs;
    for (int i=0; i<r->nr; i++) {
        if (r->vdisp[i]==0) nv.push_back(1);
        else nv.push_back(nvtmp);
    }   
    
    NHItoRAD();
    
    readytomod=true;
    
}


template <class T>
void Galmod<T>::input(Cube<T> *c, Rings<T> *rings, int NV, int LTYPE, int CMODE, float CDENS, int ISEED) {
    
    /// Alternative semplified user interface function. 
    
    int Blo[2] = {0,0};
    int Bup[2] = {c->DimX(),c->DimY()};
    
    input(c,Bup,Blo,rings,NV,LTYPE,CMODE,CDENS,ISEED);
}


template <class T>
void Galmod<T>::setOptions(int LTYPE, int CMODE, float CDENS, int ISEED) {
    
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
    
    iseed = ISEED;
    if (iseed>=0) {
        std::cout << "GALMOD warning: ISEED must be negative. Assuming -1.\n";
        iseed=-1;
    }
    
    // Initializing random number engines
    generator.seed(iseed);
    uniform   = std::uniform_real_distribution<float>(-1.,1.);
    gaussia   = std::normal_distribution<float>(0.,1.);
        
}


template <class T>
bool Galmod<T>::calculate() {
    
    /// Front end function to calculate the model.
    
    if (readytomod) {
        galmod();
    }
    else {
        std::cout<< "GALMOD error: wrong or unknown input parameters.\n";
        return false;
    }
    return true;
}


template <class T>
bool Galmod<T>::smooth(bool usescalefac) {
    
    if (!modCalculated) {
        std::cout << "GALMOD smoothing error: the model has not been ";
        std::cout << "calculated yet.\n";
        return false; 
    }
    
    //Beam oldbeam = {pixsize[0]*60, pixsize[1]*60, 0};
    Beam oldbeam = {0., 0., 0};    
    Beam newbeam = {in->Head().Bmaj()*3600.,
                    in->Head().Bmin()*3600.,
                    in->Head().Bpa()};

    Smooth3D<T> *smoothed = new Smooth3D<T>;
    smoothed->setUseScalefac(usescalefac);
    smoothed->setUseBlanks(false);
    if(!smoothed->smooth(out, oldbeam, newbeam, out->Array(), out->Array()))
        return false;   
    
    for (size_t i=0; i<out->NumPix(); i++)
        if (out->Array(i)<1.E-12) out->Array()[i] = 0;
    
    delete smoothed;
    
    return true;
}


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
                arcmconv = 60.;
        else if (cunit[i]=="ARCSEC" || cunit[i]=="ARCS" ||
                 cunit[i]=="arcsec" || cunit[i]=="arcs")
                arcmconv = 1/60.;
        else if (cunit[i]=="ARCMIN" || cunit[i]=="ARCM" ||
                 cunit[i]=="arcmin" || cunit[i]=="arcm") 
                arcmconv = 1.;
        else {
            std::cerr << "GALMOD error (unknown CUNIT for RA-DEC): ";
            std::cerr << "cannot convert to ARCMIN.\n";
            std::cerr << cunit[i];
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
    out->Head().setCrpix(0, c->Head().Crpix(0)-blo[0]);
    out->Head().setCrpix(1, c->Head().Crpix(1)-blo[1]);
    outDefined = true;
    
    double reds = in->pars().getRedshift();
    ctype3 = makelower(c->Head().Ctype(2));
    cunit3 = makelower(c->Head().Cunit(2));
    cdelt3 = c->Head().Cdelt(2);
    crval3 = c->Head().Crval(2);
    crpix3 = c->Head().Crpix(2);
    drval3 = c->Head().Drval3();
    
    if (ctype3=="wave" || ctype3=="awav" || ctype3=="wavelength" ||
        cunit3=="um" || cunit3=="nm" || cunit3=="ang" ) axtyp =2;
    else if (ctype3=="freq" || cunit3=="hz" || cunit3=="mhz") axtyp = 3;
    else if (ctype3=="velo" || ctype3=="velo-helo" || ctype3=="velo-hel" ||
             cunit3=="m/s" || cunit3=="km/s") {
        axtyp = 4;
    }
    else {
        std::cout << "GALMOD error (unknown CUNIT for spectral axis): cannot convert." << std::endl;
        std::terminate(); 
    }
    
    if (axtyp==2) {                 // Wavelength axis
        float mconv=0;
        if (cunit3=="um"||cunit3=="mum"||cunit3=="micron") mconv = 1.0E-06;
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

        // If wavelength parameter is set, take them for freq0, otherwise central channel
        double restw = in->pars().getRestwave(); 
        if (restw!=-1) freq0 = C/(restw*(1+reds)*mconv);
        else freq0 = crvalfreq; // Velocity is 0 always at the reference channel

        drval3 = C*(freq0*freq0-crvalfreq*crvalfreq)/(freq0*freq0+crvalfreq*crvalfreq);

        // Set number of emission lines
        nlines = in->pars().getParGM().RESTWAVE.size();
        if (nlines>1) {
            relint.clear();
            relvel.clear();
            for (auto i=0; i<nlines; i++) {
                double wred = in->pars().getParGM().RESTWAVE[i]*(1+reds);
                relvel.push_back(AlltoVel(wred,in->Head())*1000.);
                relint.push_back(in->pars().getParGM().RELINT[i]);
            }
        }
    }
    else if (axtyp==3) {                // Frequency axis
        
        float hzconv=0;
        if (cunit3=="hz") hzconv = 1;
        else if (cunit3=="khz") hzconv = 1.0E03;
        else if (cunit3=="mhz") hzconv = 1.0E06;
        else if (cunit3=="ghz") hzconv = 1.0E09;
        else {
            std::cerr << "GALMOD error (unknown CUNIT3): cannot convert to Hz.\n";
            std::cerr << cunit3;
            std::terminate(); 
        }
        
        crval3=crval3*hzconv;
        cdelt3=cdelt3*hzconv;  
        
        /// Information about frequency/velocity axis and conversions.
        freq0 = c->Head().Freq0()/(1+c->Head().Redshift());
        if (freq0==0) {
            freq0 = 0.1420405751786E10;
            std::cerr << "Header item FREQ0 not found. Assuming " << freq0;
            std::cerr << std::endl;
        }
        
        double crvalfreq = c->Head().Crval(2)*hzconv;
        drval3 = C*(freq0*freq0-crvalfreq*crvalfreq)/(freq0*freq0+crvalfreq*crvalfreq);
        
        // Set number of emission lines
        nlines = in->pars().getParGM().RESTFREQ.size();
        if (nlines>1) {
            relint.clear();
            relvel.clear();
            for (auto i=0; i<nlines; i++) {
                double fred = in->pars().getParGM().RESTFREQ[i]*(1+reds);
                relvel.push_back(AlltoVel(fred,in->Head())*1000.);
                relint.push_back(in->pars().getParGM().RELINT[i]);                
            }
        }
    }
    else if (axtyp==4) {                // Velocity axis
        
        float msconv=0;
        if (cunit3=="m/s" || cunit3=="ms") msconv = 1;
        else if (cunit3=="km/s") msconv = 1.0E03;
        else if ("cm/s") msconv = 1.0E-03;
        else {
            std::cerr << "GALMOD error (unknown CUNIT3): cannot convert to M/S.\n";
            std::cerr << cunit3;
            std::terminate(); 
        }
        
        /// Information about frequency/velocity axis and conversions.
        freq0 = c->Head().Freq0()/(1+c->Head().Redshift());      
        if (freq0==0) {
            freq0 = 0.1420405751786E10;
            std::cerr << "Header item FREQ0 not found. Assuming " << freq0;
            std::cerr << std::endl;
        }
        
        crval3=crval3*msconv;
        cdelt3=cdelt3*msconv;
        double crvalvel = c->Head().Crval(2)*msconv;
        drval3=freq0*sqrt((C-crvalvel)/(C+crvalvel));
        
    }
    else { 
        std::cerr << "Unknown axis type: no velocities along spectral axis.\n";
        std::terminate();
    }
    
    // Get the instrumental broadnening: when Hanning smoothing has been applied, 
    // the instrumental profile can be approximated by a gaussian with a FWHM of 
    // twice the channnel separation. Thus, the sig_instr is FWHM/2.355.
    float nch = in->pars().getLinear();
    if (nch==-1) nch=2./(2*sqrt(2*log(2)));
    chwidth = fabs(DeltaVel<double>(in->Head()))*1000;
    sig_instr = nch*chwidth;
    
}


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
                std::cerr << "GALMOD error: Radii not in increasing order.\n";
                std::terminate();
            }
            else {
                std::cerr << "GALMOD error: Radius separation too small.\n";
                std::terminate();
            }
        }
    }
    
    // This is better as true when 3DFIT (otherwise it lowers the VROT and increase VDISP)
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
            std::cerr << "GALMOD error: Negative velocity dispersion not allowed.\n";
            std::terminate();
        }
        
        udens[i]=rings->dens[i]/1.0E20;

        if (udens[i]<0) {
            std::cerr << "GALMOD error: Negative column-density not allowed.\n";
            std::terminate(); 
        }
        
        uz0[i]=rings->z0[i]/60.;
        if (uz0[i]<0) {
            std::cerr << "GALMOD error: Negative scale height not allowed.\n";
            std::terminate(); 
        }
        
        udvdz[i]  = rings->dvdz[i]*1000*60; // m/s/arcmin
        uzcyl[i]  = rings->zcyl[i]/60.;
        uinc[i]   = rings->inc[i]*M_PI/180.;
        uphi[i]   = rings->phi[i]*M_PI/180.;
        uxpos[i]  = rings->xpos[i];
        uypos[i]  = rings->ypos[i];
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

/*
template <class T>
void Galmod<T>::galmod() {

    //////////
    ////////// NEED TO UPDATE LIKE NEW GALMOD BELOW 
    //////////

    std::cout << "GALMOD PROTO" << std::endl;

    const double twopi = 2*M_PI;
    const int buflen=bsize[0]*bsize[1]*nsubs+1;
    T *datbuf = new T [buflen];
            
    ProgressBar bar(false,in->pars().isVerbose(),in->pars().getShowbar());
    bar.init(" Modeling... ",r->nr);
    
    int isd = iseed;
//  Get number of velocity profiles that will be done.
    int nprof = bsize[0]*bsize[1];
//  Initialize data buffer on zero.
    for (int i=0; i<buflen; i++) datbuf[i]=0.0;
//  Convenient random generator functions
    auto fran = std::bind(uniform, generator);
    
    // ==>> Loop over standard rings.
    for (int ir=0; ir<r->nr; ir++) {
        bar.update(ir+1);
        if (r->dens[ir]==0) continue;
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
        double z0tmp   = 0;//r->z0[ir];
        double dvdztmp = r->dvdz[ir];
        double zcyltmp = r->zcyl[ir];
        double sinc    = sin(r->inc[ir]);
        double cinc    = cos(r->inc[ir]);
        double spa     = sin(r->phi[ir])*cos(crota2)+cos(r->phi[ir])*sin(crota2);
        double cpa     = cos(r->phi[ir])*cos(crota2)-sin(r->phi[ir])*sin(crota2);
        double nvtmp   = nv[ir];
        float  fluxsc  = r->dens[ir]*twopi*rtmp*r->radsep/(nc*nvtmp);


//        double zs = 0.14*rtmp;
        double zs = 2*pow(rtmp/20,0.5);

        double zstart[2] = {r->z0[ir],-r->z0[ir]};
        //double zstart[2] = {zs,-zs};

        //std::cout << rtmp << " "  << zstart[0] << "  " << zstart[1] << std::endl;

        for (int kk=0; kk<2; kk++) {
            if (kk==1) fluxsc /= 2;

            // ==>> Loop over clouds inside each ring.
            for (int ic=0; ic<nc; ic++) {
                //          Get radius inside ring. The range includes the inner boundary,
                //          excludes the outer boundary. The probability of a radius inside
                //          a ring is proportional to the total radius and thus the
                //          surface density of the clouds is constant over the area of the ring.
                double ddum = fabs(fran());                        // STD library
                //double ddum = double(iran(isd))/double(nran);       // Classic galmod
                double R    = sqrt(pow((rtmp-0.5*r->radsep),2)+2*r->radsep*rtmp*ddum);
                //          Get azimuth and its sine and cosine.
                double az   = twopi*fabs(fran());                 // STD library
                //double az   = twopi*double(iran(isd))/double(nran); // Classic galmod
                double saz  = sin(az);
                double caz  = cos(az);
                //          Get height above the plane of the ring using a random deviate
                //          drawn from density profile of the layer.
                double z    = zstart[kk]+fdev(isd)*z0tmp;
                //          Get position in the plane of the sky with respect to the major
                //          and minor axes of the spiral galaxy.
                double x    = R*caz;
                double y    = R*saz*cinc-z*sinc;
                //          Get grid of this position, check if it is inside area of box.
                long grid[2] = {lround(r->xpos[ir]+(x*spa-y*cpa)/cdelt[0]),
                                lround(r->ypos[ir]+(x*cpa+y*spa)/cdelt[1])};
                if (grid[0]<=blo[0] || grid[0]>bhi[0]) continue;
                if (grid[1]<=blo[1] || grid[1]>bhi[1]) continue;

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
                else vsys -= vverttmp*cinc;     // <--- The sign here depends on the meaning of vvert
                

                //          MODIFIED BUILDING PROFILE FOR MULTIPLE LINES
                for (int iv=0; iv<nvtmp; iv++) {
                    double vdev = gaussia(generator)*vdisptmp;        // STD library
                    //double vdev = gasdev(isd)*vdisptmp;                 // Classic galmod
                    for (int nl=0; nl<nlines; nl++) {
                        double v     = vsys+vdev+relvel[nl];
                        int isubs = lround(velgrid(v)+crpix3-1);
                        if (isubs<0 || isubs>=nsubs) continue;
                        int idat  = iprof+isubs*nprof;
                        datbuf[idat] = datbuf[idat]+relint[nl]*fluxsc*cd2i[isubs];
                    }
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
    
    bar.fillSpace("OK.\n");
    
    delete [] datbuf;

    modCalculated=true;

}
*/

///*
template <class T>
void Galmod<T>::galmod() {

    const double twopi = 2*M_PI;
    T *array = out->Array();
//  Initialize output array on zero.
    for (int i=0; i<out->NumPix(); i++) array[i] = 0.;

    ProgressBar bar(false,in->pars().isVerbose(),in->pars().getShowbar());
    bar.init(" Modeling... ",r->nr);

    int isd = iseed;
//  Get number of velocity profiles that will be done.
    int nprof = bsize[0]*bsize[1];

//  Convenient random generator functions
    auto fran = std::bind(uniform, generator);

    // ==>> Loop over standard rings.
    for (int ir=0; ir<r->nr; ir++) {
        bar.update(ir+1);
        if (r->dens[ir]==0) continue;
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
        double nvtmp   = nv[ir];
        float  fluxsc  = r->dens[ir]*twopi*rtmp*r->radsep/(nc*nvtmp);

// ==>> Loop over clouds inside each ring.
        for (int ic=0; ic<nc; ic++) {
//          Get radius inside ring. The range includes the inner boundary,
//          excludes the outer boundary. The probability of a radius inside
//          a ring is proportional to the total radius and thus the
//          surface density of the clouds is constant over the area of the ring.
            double ddum = fabs(fran());                        // STD library
            //double ddum = double(iran(isd))/double(nran);       // Classic galmod
            double R    = sqrt(pow((rtmp-0.5*r->radsep),2)+2*r->radsep*rtmp*ddum);
//          Get azimuth and its sine and cosine.
            double az   = twopi*fabs(fran());                 // STD library
            //double az   = twopi*double(iran(isd))/double(nran); // Classic galmod
            double saz  = sin(az);
            double caz  = cos(az);
//          Get height above the plane of the ring using a random deviate
//          drawn from density profile of the layer.
            double z    = fdev(isd)*z0tmp;
//          Get position in the plane of the sky with respect to the major
//          and minor axes of the spiral galaxy.
            double x    = R*caz;
            double y    = R*saz*cinc-z*sinc;
//          Get grid of this position, check if it is inside area of box.
            long grid[2] = {lround(r->xpos[ir]+(x*spa-y*cpa)/cdelt[0]),
                            lround(r->ypos[ir]+(x*cpa+y*spa)/cdelt[1])};
//          If outside any of the ranges, jump to next cloud.
            if (grid[0]<blo[0] || grid[0]>=bhi[0]) continue;
            if (grid[1]<blo[1] || grid[1]>=bhi[1]) continue;

//            /////////// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//            float rr, theta;
//            float xr = (-(grid[0]-r->pos[0])*spa+(grid[1]-r->pos[1])*cpa);
//            float yr = (-(grid[0]-r->pos[0])*cpa-(grid[1]-r->pos[1])*spa)/cinc;
//            rr = sqrt(xr*xr+yr*yr);
//            if (rr<0.1) theta = 0.0;
//            else theta = atan2(yr, xr)/M_PI*180;
//
//            int side = 1;
//            bool use;
//            switch (side) {                     // Which side of galaxy.
//                case 1:                         //< Receding half.
//                    use = (fabs(theta)<=90.0);
//                    break;
//                case 2:                         //< Approaching half.
//                    use = (fabs(theta)>=90.0);
//                    break;
//                case 3:                         //< Both halves.
//                    use = 1;
//                    break;
//                default:
//                    break;
//            }
//
//            if (!use) continue;

//          //////////// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
//          Get profile number of current pixel and check if position is in
//          range of positions of profiles that are currently being done.
            int iprof = (grid[1]-blo[1])*bsize[0]+grid[0]-blo[0];
//          Get systematic velocity of cloud.
            double vsys = vsystmp+(vrottmp*caz+vradtmp*saz)*sinc;
//          Adding vertical rotational gradient after zcyl
            if (abs(z)>zcyltmp) vsys = vsystmp+((vrottmp-dvdztmp*(abs(z)-zcyltmp))*caz+vradtmp*saz)*sinc;
//          Adding vertical velocity
            if (z>0.) vsys += vverttmp*cinc;
            else vsys += vverttmp*cinc;     // <--- The sign here depends on the meaning of vvert. If different from the one above, GALWIND will not work (should have a flag for GALWIND execution)


//          ORIGINAL GALMOD BUILDING PROFILES
// ==>>     Build velocity profile.
//            for (int iv=0; iv<nvtmp; iv++) {
//              Get deviate drawn from gaussian velocity profile and add
//              to the systematic velocity.
//                double v     = vsys+gasdev(isd)*vdisptmp;
//              Get grid of velocity along FREQ-OHEL or VELO axis.
//              If a grid is not in the range, jump to next velocity profile.
//                int isubs = lround(velgrid(v)+crpix3-1);
//                if (isubs<0 || isubs>=nsubs) continue;
//               int idat  = iprof+isubs*nprof;
//              Convert HI atom flux per pixel to flux per pixel of 21cm
//              radiation expressed in Jy and add subcloud to the data
//              buffer.
//                datbuf[idat] = datbuf[idat]+fluxsc*cd2i[isubs];
//            }

//          MODIFIED BUILDING PROFILE FOR MULTIPLE LINES
            for (int iv=0; iv<nvtmp; iv++) {
                double vdev = gaussia(generator)*vdisptmp;        // STD library
                //double vdev = gasdev(isd)*vdisptmp;                 // Classic galmod
                for (int nl=0; nl<nlines; nl++) {
                    double v     = vsys+vdev+relvel[nl];
                    int isubs = lround(velgrid(v)+crpix3-1);
                    if (isubs<0 || isubs>=nsubs) continue;
                    size_t idat  = iprof+isubs*nprof;
                    array[idat] += relint[nl]*fluxsc*cd2i[isubs];
                }
            }
        }
    }

    bar.fillSpace("OK.\n");

    modCalculated=true;

}
//*/

/*
// GALMOD FOR SMC
template <class T>
void Galmod<T>::galmod() {

    const double twopi = 2*M_PI;
    const int buflen=bsize[0]*bsize[1]*nsubs+1;
    T *datbuf = new T [buflen];
            
    bool verb = in->pars().isVerbose();
    ProgressBar bar(" Modeling... ",false);;
    bar.setShowbar(in->pars().getShowbar());
    if (verb) bar.init(r->nr);
    
    double *ras  = new double[in->DimX()*in->DimY()];
    double *decs = new double[in->DimX()*in->DimY()];
    double *cosphi = new double[in->DimX()*in->DimY()];
    double *cosrho = new double[in->DimX()*in->DimY()];
    double *sinphi = new double[in->DimX()*in->DimY()];
    double *sinrho = new double[in->DimX()*in->DimY()];  
    
    double pc[3] = {double(r->xpos[0]),double(r->ypos[0]),0}; 
    double *wc = new double[3]; 
    pixToWCSSingle(in->Head().WCS(),pc,wc);
    double xcc = wc[0]*M_PI/180.;       // RA of galaxy centre 
    double ycc = wc[1]*M_PI/180.;       // DEC of galaxy centre
    
    for (auto x=0; x<in->DimX(); x++) {
        for (auto y=0; y<in->DimY(); y++) {
            double p[3] = {double(x),double(y),0}; 
            double *w = new double[3]; 
            pixToWCSSingle(in->Head().WCS(),p,w);
            size_t a  = x+y*in->DimX();
            ras[a]    = w[0]*M_PI/180.;
            decs[a]   = w[1]*M_PI/180.;
            cosrho[a] = cos(decs[a])*cos(ycc)*cos(ras[a]-xcc)+sin(decs[a])*sin(ycc);
            sinrho[a] = sin(acos(cosrho[a]));
            cosphi[a] = -cos(decs[a])*sin(ras[a]-xcc)/sinrho[a];
            sinphi[a] = (sin(decs[a])*cos(ycc)-cos(decs[a])*sin(ycc)*cos(ras[a]-xcc))/sinrho[a];
        }
    }
    
    //FitsWrite_2D("diorho.fits",rhos,in->DimX(),in->DimY());
    //FitsWrite_2D("diophi.fits",phis,in->DimX(),in->DimY());
    
    // Kallivayalil et al. (2013)
    const double mu_w = -0.772;
    const double mu_n = -1.117;
    const double D0   = 63.;
    const double didt = -281;
    const double dthdt= 0.;
    const double c1   = 3.08568E19/3.15576E07/3.6E06*M_PI/180.; // Conversion kpc * mas/yr  -> m/s
    const double c2   = 3.08568E19/3.15576E16*M_PI/180.;        // Conversion kpc * deg/Gyr -> m/s 
    
    // Transverse velocity and direction
    const double V_t     = c1*D0*sqrt(mu_w*mu_w+mu_n*mu_n);
    const double Theta_t = -atan2(mu_w,mu_n);
    
    int isd = iseed;
//  Get number of velocity profiles that will be done.
    int nprof = bsize[0]*bsize[1];
//  Initialize data buffer on zero.
    for (int i=0; i<buflen; i++) datbuf[i]=0.0;
//  Convenient random generator functions
    auto fran = std::bind(uniform, generator);
    
    // ==>> Loop over standard rings.
    for (int ir=0; ir<r->nr; ir++) {
        //std::cout << "WORKING ON RING #" << ir+1 <<  "/" << r->nr << std::endl;
        
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
        double nvtmp   = nv[ir];
        float  fluxsc  = r->dens[ir]*twopi*rtmp*r->radsep/(nc*nvtmp);

// ==>> Loop over clouds inside each ring.
        for (int ic=0; ic<nc; ic++) {
//          Get radius inside ring. The range includes the inner boundary,
//          excludes the outer boundary. The probability of a radius inside
//          a ring is proportional to the total radius and thus the 
//          surface density of the clouds is constant over the area of the ring.
            double ddum = fabs(fran());                        // STD library
            //double ddum = double(iran(isd))/double(nran);       // Classic galmod
            double R    = sqrt(pow((rtmp-0.5*r->radsep),2)+2*r->radsep*rtmp*ddum);
//          Get azimuth and its sine and cosine.
            double az   = twopi*fabs(fran());                 // STD library
            //double az   = twopi*double(iran(isd))/double(nran); // Classic galmod
            double saz  = sin(az); 
            double caz  = cos(az); 
//          Get height above the plane of the ring using a random deviate
//          drawn from density profile of the layer.
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

//          Get profile number of current pixel and check if position is in
//          range of positions of profiles that are currently being done.
//          If outside any of the ranges, jump to next cloud.
            int iprof = (grid[1]-blo[1])*bsize[0]+grid[0]-blo[0]-bsize[0];  

            // PART FOR SMC
            int np = grid[0]+grid[1]*bsize[0];
                
//          Getting vtc, vts and wts    
            double mu_c = -mu_w*spa + mu_n*cpa;
            double mu_s = -mu_w*cpa - mu_n*spa;
            double vtc  = c1*D0*mu_c;    
            double vts  = c1*D0*mu_s;
            double wts  = vts + c2*D0*didt;
            double theta   = r->phi[ir] + M_PI_2;
            double theta_t = Theta_t + M_PI_2;
            double cospt  = cosphi[np]*cos(theta)+sinphi[np]*sin(theta); // = cos(phi-theta)
            double sinpt  = sinphi[np]*cos(theta)-cosphi[np]*sin(theta); // = sin(phi-theta)
            double cosptt = cosphi[np]*cos(theta_t)+sinphi[np]*sin(theta_t); // = cos(phi-theta_t)
            double sinptt = sinphi[np]*cos(theta_t)-cosphi[np]*sin(theta_t); // = sin(phi-theta)

            //##############################################################
            // # Center of mass velocity components (eq. 13 vdM+2002)
            //##############################################################
            double v1_CM =  V_t*sinrho[np]*cosptt + vsystmp*cosrho[np];
            double v2_CM =  V_t*cosrho[np]*cosptt - vsystmp*sinrho[np];
            double v3_CM = -V_t           *sinptt;

            //##############################################################
            //# Precession/nutation velocity components (eq. 16 vdM+2002)
            //##############################################################
            //# Constant factor 
            double fac   = c2*D0*sinrho[np]/(cinc*cosrho[np] + sinc*sinrho[np]*sinpt);
            //# Components
            double v1_PN = fac*(didt*sinpt*( cinc*cosrho[np] + sinc*sinrho[np]*sinpt));
            double v2_PN = fac*(didt*sinpt*(-cinc*sinrho[np] + sinc*cosrho[np]*sinpt));
            double v3_PN = fac*(didt*sinpt*(                 + sinc           *cospt) + dthdt*cinc);

            //##############################################################
            //# Internal motion velocity components (eq. 21 vdM+2002)
            //##############################################################
            //# Spin factor (Negative if PA is defined from the receding side)
            double s     = 1.;
            //# Constant factor 
            double fac2   = s*vrottmp/sqrt(cinc*cinc*cospt*cospt+sinpt*sinpt);
            //# Components
            double v1_IN = fac2*( sinc*cospt*(cinc*cosrho[np] + sinc*sinrho[np]*sinpt));
            double v2_IN = fac2*(-sinc*cospt*(cinc*sinrho[np] - sinc*cosrho[np]*sinpt));
            double v3_IN = fac2*(-(cinc*cinc*cospt*cospt+sinpt*sinpt));

            //##############################################################
            //# Observed velocity components (eq. 11 vdM+2002)
            //##############################################################
            double v1 = v1_CM + v1_PN + v1_IN;
            double v2 = v2_CM + v2_PN + v2_IN;
            double v3 = v3_CM + v3_PN + v3_IN;

//          Get systematic velocity of cloud.
            double vlos  = v1;
            
            for (int iv=0; iv<nvtmp; iv++) {
                double vdev = gaussia(generator)*vdisptmp;        // STD library
                //double vdev = gasdev(isd)*vdisptmp;                 // Classic galmod
                for (int nl=0; nl<nlines; nl++) {
                    double v     = vlos+vdev+relvel[nl];
                    int isubs = lround(velgrid(v)+crpix3-1);
                    if (isubs<0 || isubs>=nsubs) continue;
                    int idat  = iprof+isubs*nprof;
                    datbuf[idat] = datbuf[idat]+relint[nl]*fluxsc*cd2i[isubs];
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
    delete [] ras;
    delete [] decs;
    delete [] cosphi;
    delete [] cosrho;
    delete [] sinphi;
    delete [] sinrho;  

    modCalculated=true;
    
}
*/



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
/// as the flux per beam of one square arcminute in Jy.

    double const jy = 1E-26;             // Conversion factor to Jy
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

        double ddum = 2*(K/jy)/(fac*(labsubs*labsubs))/fabs(cdelt3);
        // Following line removes dependency from lambda
        //double ddum=2*(K/jy)/fac/fabs(cdelt3);
        // The solid angle in the definition of the intensity is 
        // in steradians, thus convert it to square arcminutes.
        cd2i[isubs]= ddum*pow((M_PI/(180.0*60.0)),2);
    }

}


template <class T>  
double Galmod<T>::velgrid(double v) {
    
    /// Function to transform a velocity to a grid.
        
    double velg=0, fdv;
    
    if (axtyp==2) {                     //< Wavelength axis
        double f = freq0*sqrt((C-v)/(C+v));
        double l = C/f;
        velg = (l-crval3)/cdelt3;
    }
    else if (axtyp==3) {                //< Frequency axis.
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


template <class T>
double Galmod<T>::fdev(int &idum){
    
    /// Function to get random deviates for various functions.
    /// The double precision variable Fdev contains the random deviate.
    /// If nran is within a factor two of the largest possible integer
    /// then integer overflow could occur.

    double Fdev=0, x=0;
    
    if (ltype==1) {                                 /// Gaussian function: exp(-0.5*x^2) 
        Fdev = gaussia(generator);                // STD library
        //Fdev = gasdev(idum);                        // Classic galmod
    }
    else {
        x = uniform(generator);                    // STD library
        //x = double(1+2*iran(idum)-nran)/double(nran);// Classic galmod
        if (ltype==2) Fdev = atanh(x);               // Sech2 function: sech2(x)
        else if (ltype==3) {                         // Exponential function: exp(-|x|)
           if (x>=0.0) Fdev = -log(x);
           if (x<0.0)  Fdev = log(-x);
        }
        else if (ltype==4) Fdev = tan(M_PI_2*x);     // Lorentzian function: 1/(1+x**2)
        else if (ltype==5) Fdev = x;                 // Box function.
        else cout << "Unknown function." << endl;
    }

    return Fdev;
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// Members of Galmod_wind class
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
template <class T>
void Galmod_wind<T>::input(Cube<T> *c, int *Boxup, int *Boxlow, Shells<T> *shells, 
                      int NV, int LTYPE, int CMODE, float CDENS, int ISEED) {
    
    /// This function sets all parameters needed to use the Galmod Object. A description
    /// of input is given in the Galmod class definition file (galmod.hh). 
        
    this->initialize(c, Boxup, Boxlow);
    shellIO(shells);
    this->setOptions(LTYPE, CMODE, CDENS, ISEED);
    
    int nvtmp = NV;
    if (nvtmp<1) nvtmp=this->nsubs;
    for (int i=0; i<s->ns; i++) {
        if (s->vdisp[i]==0) this->nv.push_back(1);
        else this->nv.push_back(nvtmp);
    }   

    this->NHItoRAD();
    this->readytomod=true;
    
}


template <class T>
void Galmod_wind<T>::shellIO(Shells<T> *shells) {
    
    /// This function creates the set of shell to be used in the model.
    /// - Velocities are given in km/s and then converted in m/s.
    /// - Radii are given in arcsec and then converted in arcmin.
    /// - Column densities are given in HI atoms/cm^2 and the converted 
    ///   in units of 1E20 atoms/cm^2.
    /// - Angles are given in grades and then converted in radians.
    
    s = new Shells<T>;

    int nur = shells->ns;
    T *uradii = new T[nur];
    T *uvrot  = new T[nur];
    T *uvdisp = new T[nur];
    T *uvsph  = new T[nur];
    T *udens  = new T[nur];
    T *uinc   = new T[nur];
    T *uphi   = new T[nur];
    T *uxpos  = new T[nur];
    T *uypos  = new T[nur];
    T *uvsys  = new T[nur];
    T *uvoan  = new T[nur];
    
    s->sep=0.75*min(this->pixsize[0],this->pixsize[1]);
    uradii[0]=shells->radii[0]/60.;
    for (int i=1; i<nur; i++) {
        uradii[i]=shells->radii[i]/60.;
        if (uradii[i-1]+s->sep>=uradii[i]) {
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

    for (int i=0; i<nur; i++) {
        uvrot[i]  = shells->vrot[i]*1000;   
        uvsph[i]  = shells->vsph[i]*1000;   
        uvdisp[i] = shells->vdisp[i]*1000;
        if (uvdisp[i]<0) {
            std::cout << "GALMOD error: Negative velocity dispersion not allowed.\n";
            std::terminate();
        }
        
        udens[i]=shells->dens[i]/1.0E20;

        if (udens[i]<0) {
            std::cout << "GALMOD error: Negative column-density not allowed.\n";
            std::terminate(); 
        }
        
        uinc[i]   = shells->inc[i]*M_PI/180.;
        uphi[i]   = shells->pa[i]*M_PI/180.;
        uvoan[i]  = shells->openang[i]*M_PI/180.;
        uxpos[i]  = shells->xpos[i]+1;
        uypos[i]  = shells->ypos[i]+1;
        uvsys[i]  = shells->vsys[i]*1000;
    }
    
    if (uradii[0]!=0 && empty) s->radii.push_back(uradii[0]+s->sep/2.0);
    else {
        s->radii.push_back(s->sep/2.0);
        while (s->radii.back()<uradii[0]) {
            s->vrot.push_back(uvrot[0]);
            s->vsph.push_back(uvsph[0]);
            s->vdisp.push_back(uvdisp[0]);
            s->dens.push_back(udens[0]);
            s->inc.push_back(uinc[0]);
            s->pa.push_back(uphi[0]);
            s->openang.push_back(uvoan[0]);
            s->xpos.push_back(uxpos[0]);
            s->ypos.push_back(uypos[0]);
            s->vsys.push_back(uvsys[0]);
            s->radii.push_back(s->radii.back()+s->sep);
        }
    }
    
    for (int i=1; i<nur; i++) {
        T dur       = uradii[i]-uradii[i-1];
        T dvrotdr   = (uvrot[i]-uvrot[i-1])/dur;
        T dvsphdr   = (uvsph[i]-uvsph[i-1])/dur;
        T dvdispdr  = (uvdisp[i]-uvdisp[i-1])/dur;
        T ddensdr   = (udens[i]-udens[i-1])/dur;
        T dincdr    = (uinc[i]-uinc[i-1])/dur;
        T dphidr    = (uphi[i]-uphi[i-1])/dur;
        T dopendr   = (uvoan[i]-uvoan[i-1])/dur;
        T dxposdr   = (uxpos[i]-uxpos[i-1])/dur;
        T dyposdr   = (uypos[i]-uypos[i-1])/dur;
        T dvsysdr   = (uvsys[i]-uvsys[i-1])/dur;
        while (s->radii.back()<uradii[i]) {
            T dr = s->radii.back()-uradii[i-1];
            s->vrot.push_back(uvrot[i-1]+dvrotdr*dr);
            s->vsph.push_back(uvsph[i-1]+dvsphdr*dr);
            s->vdisp.push_back(uvdisp[i-1]+dvdispdr*dr);
            s->dens.push_back(udens[i-1]+ddensdr*dr);
            s->inc.push_back(uinc[i-1]+dincdr*dr);
            s->pa.push_back(uphi[i-1]+dphidr*dr);
            s->openang.push_back(uvoan[i-1]+dopendr*dr);
            s->xpos.push_back(uxpos[i-1]+dxposdr*dr);
            s->ypos.push_back(uypos[i-1]+dyposdr*dr);
            s->vsys.push_back(uvsys[i-1]+dvsysdr*dr);
            s->radii.push_back(s->radii.back()+s->sep);
        }
    }
    
    s->radii.pop_back();
    s->ns=s->radii.size();      
    
    shellDefined = true;
        
    delete [] uradii;
    delete [] uvrot;
    delete [] uvsph;
    delete [] uvdisp;
    delete [] udens;
    delete [] uinc;
    delete [] uphi;
    delete [] uxpos;
    delete [] uypos;
    delete [] uvsys;
    delete [] uvoan;
}


template <class T>
bool Galmod_wind<T>::calculate() {
    
    /// Front end function to calculate the model.
    
    if (this->readytomod) galmod_wind();
    else {
        std::cout<< "GALMOD_WIND error: wrong or unknown input parameters.\n";
        return false;
    }
    return true;
}


template <class T>
void Galmod_wind<T>::galmod_wind() {

    // GALMOD FOR SPHERICAL WIND
    const int buflen=this->bsize[0]*this->bsize[1]*this->nsubs+1;
    T *datbuf = new T [buflen];
//  Initialize data buffer on zero.
    for (int i=0; i<buflen; i++) datbuf[i]=0.0;
    
    int isd = this->iseed;
//  Get number of velocity profiles that will be done.
    int nprof = this->bsize[0]*this->bsize[1];
//  Convenient random generator functions
    auto fran = std::bind(this->uniform, this->generator);

    ProgressBar bar(false,this->in->pars().isVerbose(),this->in->pars().getShowbar());

    int nthreads = this->in->pars().getThreads();
#pragma omp parallel num_threads(nthreads)
{
    bar.init(" Modeling outflow... ",s->ns);
    
#pragma omp for schedule (dynamic)
    // ==>> Loop over standard spherical shells.
    for (int ir=0; ir<s->ns; ir++) {
        bar.update(ir+1);
//      Get radius
        double rtmp = s->radii[ir];
//      Get number of clouds inside ring.
        int nc = lround(this->cdens*pow(s->dens[ir],this->cmode)*2*M_PI*rtmp*s->sep/this->pixarea); 
        if (nc==0) {
            std::cerr << " GALMOD ERROR: No clouds used. Choose higher CDENS " << std::endl;
            std::terminate();
        }
//      Get values of ir-radius.
        double vrottmp = s->vrot[ir];
        double vsphtmp = s->vsph[ir];
        double opentmp = s->openang[ir]/2.;
        double vdisptmp= sqrt(s->vdisp[ir]*s->vdisp[ir]+this->sig_instr*this->sig_instr);
        double vsystmp = s->vsys[ir];
        double sinc    = sin(s->inc[ir]);
        double cinc    = cos(s->inc[ir]);
        double spa     = sin(s->pa[ir])*cos(this->crota2)+cos(s->pa[ir])*sin(this->crota2);
        double cpa     = cos(s->pa[ir])*cos(this->crota2)-sin(s->pa[ir])*sin(this->crota2);
        double nvtmp   = this->nv[ir];
        double fluxsc  = s->dens[ir]*2*M_PI*rtmp*s->sep/(nc*nvtmp);
        
// ==>> Loop over clouds inside each ring.
        for (int ic=0; ic<nc; ic++) {
//          Get radius inside shell. The range includes the inner boundary,
//          excludes the outer boundary.
            double ddum = fabs(fran());
            double rad  = sqrt(pow((rtmp-0.5*s->sep),2)+2*s->sep*rtmp*ddum);
//          Get azimuth and its sine and cosine.
            double az   = 2*M_PI*fabs(fran());
            double saz  = sin(az); 
            double caz  = cos(az);
            // Get polar angle
            double phi  = opentmp*fabs(fran());
            if (ic%2==0) phi = M_PI - phi;
            double sphi = sin(phi);
            double cphi = cos(phi);
//          Get cylindrical coordinates
            double R    = rad*sphi;
            double z    = rad*cphi;
//          Get position in the plane of the sky with respect to the major
//          and minor axes of the spiral galaxy.
            double x    = R*caz;
            double y    = R*saz*cinc-z*sinc;                                                 
//          Get grid of this position, check if it is inside area of box.
            long grid[2] = {lround(s->xpos[ir]+(x*spa-y*cpa)/this->cdelt[0]),
                            lround(s->ypos[ir]+(x*cpa+y*spa)/this->cdelt[1])};             
            if (grid[0]<=this->blo[0] || grid[0]>this->bhi[0]) continue;
            if (grid[1]<=this->blo[1] || grid[1]>this->bhi[1]) continue;
  
//          Get profile number of current pixel and check if position is in
//          range of positions of profiles that are currently being done.
            int iprof = (grid[1]-this->blo[1])*this->bsize[0]+grid[0]-this->blo[0]-this->bsize[0];
//          Get LOS velocity of cloud.
            double vsys = vsystmp+vsphtmp*(cphi*cinc+sphi*saz*sinc)+(vrottmp*caz*sinc);
                
            for (int iv=0; iv<nvtmp; iv++) {
                double vdev = this->gaussia(this->generator)*vdisptmp;
                for (int nl=0; nl<this->nlines; nl++) {
                    double v  = vsys+vdev+this->relvel[nl];
                    int isubs = lround(this->velgrid(v)+this->crpix3-1);
                    if (isubs<0 || isubs>=this->nsubs) continue;
                    int idat  = iprof+isubs*nprof;
                    datbuf[idat] = datbuf[idat]+this->relint[nl]*fluxsc*this->cd2i[isubs];
                }
            }
        }
        
    }
}

//  Write data to output Cube.      
    for (int isubs=0; isubs<this->nsubs; isubs++) {
        int pixstart=isubs*this->bsize[0]*this->bsize[1];
        int idat=1+isubs*nprof;
        for (int i=0; i<nprof; i++)
            this->out->Array(pixstart+i)=datbuf[idat+i];
    }
    
    bar.fillSpace("OK.\n");
    delete [] datbuf;
    
    this->modCalculated=true;

}


// Explicit instantiation of the classes
template class Galmod<float>;
template class Galmod<double>;
template class Galmod_wind<float>;
template class Galmod_wind<double>;

}

#undef C  
#undef nran 
