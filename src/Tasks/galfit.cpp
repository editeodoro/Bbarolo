//--------------------------------------------------------------------
// galfit.cpp: Members functions of the Galfit class.
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
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <ctime>
#include <string>
#include <Arrays/param.hh>
#include <Arrays/cube.hh>
#include <Map/detection.hh>
#include <Tasks/galfit.hh>
#include <Tasks/galmod.hh>
#include <Tasks/smooth3D.hh>
#include <Utilities/utils.hh>
#include <Utilities/lsqfit.hh>
#include <Utilities/progressbar.hh>
#include <Utilities/paramguess.hh>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Model {
    
template <class T>
void Galfit<T>::defaults() {
    
    inDefined   = false;
    outDefined  = false;
    line_imDefined = false;
    cfieldAllocated = false;
    chan_noiseAllocated = false;
    wpow = 1;
    anglepar = 3;
    second = false;
    global=false;
    verb = true;
    func_norm = &Model::Galfit<T>::norm_local;
    mins = new T[MAXPAR];
    maxs = new T[MAXPAR];
    mpar = new bool[MAXPAR];
}
template void Galfit<float>::defaults();
template void Galfit<double>::defaults();


template <class T>
Galfit<T>::~Galfit () {
    if (outDefined) delete outr;
    if (inDefined) delete inr;
    if (line_imDefined) delete line_im;
    if (cfieldAllocated) delete [] cfield; 
    if (chan_noiseAllocated) delete [] chan_noise;
    delete [] mins;
    delete [] maxs;
    delete [] mpar;
}
template Galfit<float>::~Galfit();
template Galfit<double>::~Galfit();


template <class T>
Galfit<T>& Galfit<T>::operator=(const Galfit &g) {
    
    if(this==&g) return *this;
    
    this->in    = g.in;
    
    if (inDefined) delete inr;
    this->inDefined = g.inDefined;
    if (inDefined) *inr = *g.inr;
    else this->inr = g.inr;
    
    if (outDefined) delete outr;
    this->outDefined = g.outDefined;
    if (outDefined) {
        outr = new Rings<T>;
        *outr = *g.outr;
    }
    
    this->mask = g.mask;
    
    for (int i=0; i<MAXPAR; i++) {
        this->mpar[i] = g.mpar[i];
        this->maxs[i] = g.maxs[i];
        this->mins[i] = g.mins[i];
    }
    
    this->par       = g.par;
    this->nfree     = g.nfree;
    this->arcconv   = g.arcconv;                            
    this->distance  = g.distance;   
    this->NconX     = g.NconX;
    this->NconY     = g.NconY;
    this->wpow      = g.wpow;
    this->global    = g.global;

    if (cfieldAllocated) delete [] cfield;
    this->cfieldAllocated = g.cfieldAllocated;
    if (cfieldAllocated) {
        this->cfield = new double[NconX*NconY];
        for (int i=0; i<this->NconX*this->NconY; i++)
            this->cfield[i] = g.cfield[i];
    }

    if (chan_noiseAllocated) delete [] chan_noise;
    this->chan_noiseAllocated = g.chan_noiseAllocated;
    if (chan_noiseAllocated) {
        this->chan_noise = new float[this->in->DimZ()];
        for (int i=0; i<this->in->DimZ(); i++)
            this->chan_noise[i] = g.chan_noise[i];
    }
    
    func_norm = g.func_norm;

    if (line_imDefined) delete line_im;
    this->line_imDefined = g.line_imDefined;
    if (this->line_imDefined) {
        this->line_im = g.line_im;
    }

    return *this;
}
template Galfit<float>& Galfit<float>::operator=(const Galfit<float>&);
template Galfit<double>& Galfit<double>::operator=(const Galfit<double>&);


template <class T>
Galfit<T>::Galfit(Cube<T> *c) {
    
    // This constructor reads all needed parameters from the Cube object.
    // Cube object must contain a Param method with all information.
    
    par = c->pars().getParGF();
    
    c->checkBeam();

    // Building Rings object 
    // Try to read ring information from an input file
    Rings<T> file_rings;
    bool radii_b = getDataColumn(file_rings.radii,par.RADII);
    bool xpos_b  = getDataColumn(file_rings.xpos,par.XPOS);
    bool ypos_b  = getDataColumn(file_rings.ypos,par.YPOS);
    bool vsys_b  = getDataColumn(file_rings.vsys,par.VSYS);
    bool vrot_b  = getDataColumn(file_rings.vrot,par.VROT);
    bool vrad_b  = getDataColumn(file_rings.vrad,par.VRAD);
    bool vdisp_b = getDataColumn(file_rings.vdisp,par.VDISP);
    bool z0_b    = getDataColumn(file_rings.z0,par.Z0);
    bool dens_b  = getDataColumn(file_rings.dens,par.DENS);
    bool inc_b   = getDataColumn(file_rings.inc,par.INC);
    bool pa_b    = getDataColumn(file_rings.phi,par.PHI);
    bool vvert_b = getDataColumn(file_rings.vvert,par.VVERT);
    bool dvdz_b  = getDataColumn(file_rings.dvdz,par.DVDZ);
    bool zcyl_b  = getDataColumn(file_rings.zcyl,par.ZCYL);
    bool onefile = radii_b||xpos_b||ypos_b||vsys_b||vrot_b||vdisp_b||z0_b||dens_b||inc_b||pa_b||vrad_b||vvert_b||dvdz_b||zcyl_b;

    size_t size[MAXPAR+4] = {file_rings.radii.size(),file_rings.xpos.size(), file_rings.ypos.size(), 
                             file_rings.vsys.size(),file_rings.vrot.size(),file_rings.vdisp.size(),
                             file_rings.z0.size(),file_rings.dens.size(),file_rings.inc.size(),file_rings.phi.size(),
                             file_rings.vrad.size(),file_rings.vvert.size(),file_rings.dvdz.size(),file_rings.zcyl.size()};

    size_t max_size=UINT_MAX;
    for (int i=0; i<MAXPAR+4; i++) if (size[i]!=0 && size[i]<max_size) max_size=size[i];
    
    int nr=0;
    T radsep, xpos, ypos, vsys, vrot, vdisp, z0, dens, inc, pa, vrad, vvert, zcyl, dvdz;

    bool toEstimate =  (par.RADII=="-1" && (par.NRADII==-1 || par.RADSEP==-1)) ||
                        par.XPOS=="-1" || par.YPOS=="-1" || par.VSYS=="-1" ||
                        par.VROT=="-1" || par.PHI=="-1"  || par.INC=="-1";

    if (toEstimate) {

        if (verb) std::cout << "\n Estimating initial parameters... " << std::flush;
        T *init_par = EstimateInitial(c,&par);
        if (verb) std::cout << "Done." << std::endl;
        
        string pos[2] = {par.XPOS, par.YPOS};
        double *pixs = getCenterCoordinates(pos, c->Head());
        
        nr    = par.NRADII!=-1 ? par.NRADII : init_par[0];
        radsep= par.RADSEP!=-1 ? par.RADSEP : init_par[1];
        xpos  = par.XPOS!="-1" ? pixs[0] : init_par[2];
        ypos  = par.YPOS!="-1" ? pixs[1] : init_par[3];
        vsys  = par.VSYS!="-1" ? atof(par.VSYS.c_str()) : init_par[4];
        if (distance==-1) distance = VeltoDist(fabs(vsys));
        vrot  = par.VROT!="-1" ? atof(par.VROT.c_str()) : init_par[5];
        vdisp = par.VDISP!="-1" ? atof(par.VDISP.c_str()): 8.;// default is 8 km/s
        z0    = par.Z0!="-1" ? atof(par.Z0.c_str()) : 0.; // default is infinitely thin disk
        dens  = par.DENS!="-1" ? atof(par.DENS.c_str()) : 1.;
        inc   = par.INC!="-1" ? atof(par.INC.c_str()) : init_par[6];
        pa    = par.PHI!="-1" ? atof(par.PHI.c_str()) : init_par[7];
        vrad  = par.VRAD!="-1" ? atof(par.VRAD.c_str()) : 0.;
        delete [] init_par;
    }
    else {
        nr    = par.NRADII;
        radsep= par.RADSEP;
        string pos[2] = {par.XPOS, par.YPOS};
        double *pixs = getCenterCoordinates(pos, c->Head());
        xpos  = pixs[0];
        ypos  = pixs[1];
        vsys  = atof(par.VSYS.c_str());
        if (par.DISTANCE==-1) par.DISTANCE = VeltoDist(fabs(vsys));
        vrot  = atof(par.VROT.c_str());
        vdisp = par.VDISP!="-1" ? atof(par.VDISP.c_str()): 8.;                  // default is 8 km/s
        z0    = par.Z0!="-1" ? atof(par.Z0.c_str()) : 0.; // default is infinitely thin disk
        vrad  = par.VRAD!="-1" ? atof(par.VRAD.c_str()) : 0.;
        dens  = par.DENS!="-1" ? atof(par.DENS.c_str()) : 1.;
        inc   = atof(par.INC.c_str());
        pa    = atof(par.PHI.c_str());
    }
    
    vvert = par.VVERT!="-1" ? atof(par.VVERT.c_str()) : 0.;
    dvdz  = par.DVDZ!="-1" ? atof(par.DVDZ.c_str()) : 0.;
    zcyl  = par.ZCYL!="-1" ? atof(par.ZCYL.c_str()) : 0.;
    
    if (nr==0) {
        std::cout << "\n 3DFIT ERROR: The number of radii must be > 0! " << std::endl;
        std::terminate();
    }

    nr = nr>0 && nr<max_size ? nr : max_size;
    if (radii_b) {
        radsep = 0;
        for (uint i=1; i<file_rings.radii.size()-1; i++)
            radsep += file_rings.radii[i+1]-file_rings.radii[i];
        radsep/=(file_rings.radii.size()-2);
    }

    Rings<T> *inR = new Rings<T>;
    inR->nr     = nr;
    inR->radsep = radsep;
    for (int i=0; i<inR->nr; i++) {
        if (radii_b) inR->radii.push_back(file_rings.radii[i]);
        else inR->radii.push_back(i*radsep+radsep/2.);
        if (vrot_b) inR->vrot.push_back(file_rings.vrot[i]);
        else inR->vrot.push_back(vrot);
        if (vdisp_b) inR->vdisp.push_back(file_rings.vdisp[i]);
        else inR->vdisp.push_back(vdisp);
        if (z0_b) inR->z0.push_back(file_rings.z0[i]);
        else inR->z0.push_back(z0);
        if (dens_b) inR->dens.push_back(file_rings.dens[i]*1.E20);
        else inR->dens.push_back(dens*1.E20);
        if (inc_b) inR->inc.push_back(file_rings.inc[i]);
        else inR->inc.push_back(inc);
        if (pa_b) inR->phi.push_back(file_rings.phi[i]);
        else inR->phi.push_back(pa);
        if (xpos_b) inR->xpos.push_back(file_rings.xpos[i]);
        else inR->xpos.push_back(xpos);
        if (ypos_b) inR->ypos.push_back(file_rings.ypos[i]);
        else inR->ypos.push_back(ypos);
        if (vsys_b) inR->vsys.push_back(file_rings.vsys[i]);
        else inR->vsys.push_back(vsys);
        if (vrad_b) inR->vrad.push_back(file_rings.vrad[i]);
        else inR->vrad.push_back(vrad);
        
        // In the current version, vertical motions, and gradients are not fitted
        if (vvert_b) inR->vvert.push_back(file_rings.vvert[i]);
        else inR->vvert.push_back(vvert);
        if (dvdz_b) inR->vvert.push_back(file_rings.dvdz[i]);
        else inR->dvdz.push_back(dvdz);
        if (zcyl_b) inR->vvert.push_back(file_rings.zcyl[i]);
        else inR->zcyl.push_back(zcyl);
    }
    
    if (!c->pars().getflagGalMod()) {
        if (!onefile) showInitial(inR, std::cout);
        else printInitial(inR, c->pars().getOutfolder()+"initial_rings.txt");
    }


    // Setup all needed parameters
    setup(c,inR,&par);
  
}
template Galfit<float>::Galfit(Cube<float>*);
template Galfit<double>::Galfit(Cube<double>*);


template <class T>
Galfit<T>::Galfit (Cube<T> *c, Rings<T> *inrings, float DELTAINC, float DELTAPHI, int LTYPE, int FTYPE, 
                   int WFUNC, int BWEIGHT, int NV, double TOL, int CDENS, int STARTRAD, string MASK, 
                   string NORM, string FREE, string SIDE, bool TWOSTAGE, string POLYN, bool ERRORS, 
                   bool SMOOTH, float DISTANCE, double REDSHIFT, double RESTWAVE, string OUTFOLD, int NTHREADS) {
                       
    // Setting all parameters in the GALFIT_PAR container
    par.DELTAINC = DELTAINC;
    par.DELTAPHI = DELTAPHI;
    par.LTYPE    = LTYPE;
    par.FTYPE    = FTYPE;
    par.WFUNC    = WFUNC;
    par.BWEIGHT  = BWEIGHT;
    par.NV       = NV;
    par.TOL      = TOL;
    par.CDENS    = CDENS;
    par.STARTRAD = STARTRAD;
    par.NORM     = NORM;
    par.FREE     = FREE;
    par.SIDE     = SIDE;
    par.TWOSTAGE = TWOSTAGE;
    par.POLYN    = POLYN;
    par.SM       = SMOOTH;
    par.DISTANCE = DISTANCE;
    par.REDSHIFT = REDSHIFT;
    par.RESTWAVE[0] = RESTWAVE;   //@TODO Update for doublets from pyBBarolo
    par.flagERRORS = ERRORS;
        
    c->pars().getParGF() = par;
    
    // Create directory tree if it does not exist
    checkHome(OUTFOLD);    
    c->pars().setOutfolder(OUTFOLD);
    c->pars().setThreads(NTHREADS);
    c->pars().setMASK(MASK);
    
    mkdirp(OUTFOLD.c_str());
    
    setup(c,inrings,&par);
    
}
template Galfit<float>::Galfit(Cube<float>*,Rings<float> *,float,float,int,int,int,int,int,double,int,int,
                               string,string,string,string,bool,string,bool,bool,float,double,double,string,int);
template Galfit<double>::Galfit(Cube<double>*,Rings<double> *,float,float,int,int,int,int,int,double,int,int,
                                string,string,string,string,bool,string,bool,bool,float,double,double,string,int);


template <class T>
void Galfit<T>::setup (Cube<T> *c, Rings<T> *inrings, GALFIT_PAR *p) {
    
    defaults();    
    in = c;
    par = *p;
    inr = new Rings<T>;
    *inr = *inrings;
    inDefined = true;

    // Check that radii are ok.
    for (int ir=0; ir<inr->nr-1; ir++) {
        if (ir!=inr->nr-1) {
            if (inr->radii[ir+1]<=inr->radii[ir]) {
                cout << "3DFIT WARNING: Radii not in increasing order.\n";
            }
        }
        if (inr->radii[ir]<0) {
            cout << "3DFIT ERROR: Negative radius!!!\n";
            std::terminate();
        }
    }
    
    // Checking that the beam has all information
    in->checkBeam();

    // Setting other GALFIT variables
    verb = in->pars().isVerbose();
    arcconv = arcsconv(in->Head().Cunit(0));
    distance = par.DISTANCE==-1 ? VeltoDist(fabs(inr->vsys[0])) : par.DISTANCE;
    chan_noise = new float[in->DimZ()];
    chan_noiseAllocated = true;
    for (int z=0; z< in->DimZ(); z++) chan_noise[z]=1;
    
    wpow = par.WFUNC;
    string polyn = makelower(par.POLYN);
    if (polyn=="bezier") anglepar=-1;
    else anglepar = 1+atoi(polyn.c_str());
    
    // Read par.FREE and set free parameters
    setFree();
    
    // Choose right function for normalization 
    if (par.NORM=="NONE") func_norm = &Model::Galfit<T>::norm_none;
    else if (par.NORM=="AZIM") func_norm = &Model::Galfit<T>::norm_azim;
    else func_norm = &Model::Galfit<T>::norm_local;

    // Creating mask if does not exist and write it in a fitsfile.
    if (!in->MaskAll() || in->pars().getMASK()=="NEGATIVE") in->BlankMask(chan_noise);
    mask = in->Mask();
    
    // Setting limits for fitting parameters
    double kpcperarc = KpcPerArc(distance);
    maxs[VROT]  = 1000;
    mins[VROT]  = 0;
    maxs[VDISP] = 500;
    mins[VDISP] = 0.01;
    maxs[Z0] = 5/kpcperarc;         // Max scaleheight allowed is 5 Kpc.  
    mins[Z0] = 0.;                  // Min scaleheight allowed is 0 pc.
    maxs[INC] = *max_element(inr->inc.begin(),inr->inc.end())+par.DELTAINC;
    mins[INC] = *min_element(inr->inc.begin(),inr->inc.end())-par.DELTAINC;
    maxs[PA]  = *max_element(inr->phi.begin(),inr->phi.end())+par.DELTAPHI;
    mins[PA]  = *min_element(inr->phi.begin(),inr->phi.end())-par.DELTAPHI; 
    if (maxs[INC]>90) maxs[INC] = 90;
    if (mins[INC]<0)  mins[INC] = 0;
    if (maxs[PA]>360) maxs[PA]  = 360;
    if (mins[PA]<-360) mins[PA]  = -360;
    maxs[XPOS] = *max_element(inr->xpos.begin(),inr->xpos.end())+10;
    mins[XPOS] = *min_element(inr->xpos.begin(),inr->xpos.end())-10;
    maxs[YPOS] = *max_element(inr->ypos.begin(),inr->ypos.end())+10;
    mins[YPOS] = *min_element(inr->ypos.begin(),inr->ypos.end())-10;
    maxs[VSYS] = AlltoVel(in->getZphys(in->DimZ()-1), in->Head());
    mins[VSYS] = AlltoVel(in->getZphys(0), in->Head());
    if (maxs[XPOS]>in->DimX()) maxs[XPOS] = in->DimX();
    if (maxs[YPOS]>in->DimY()) maxs[YPOS] = in->DimY();
    if (mins[XPOS]<0)  mins[XPOS] = 0;
    if (mins[YPOS]<0)  mins[YPOS] = 0;
    if (maxs[VSYS]<mins[VSYS]) std::swap(maxs[VSYS],mins[VSYS]);
    maxs[VRAD]  = 100;
    mins[VRAD]  = -100;
    
    // Setting the convolution field
    if (par.SM) {
        if (!setCfield()) {
            std::cout << "3DFIT WARNING: can not set an appropriate convolution "
                      << "field. Turning off the convolution step.\n";
            par.SM = false;
        }
    }
    
    // Allocate output Rings
    outr = new Rings<T>;
    *outr = *inr;
    outDefined = true;
}
template void Galfit<float>::setup (Cube<float>*, Rings<float>*, GALFIT_PAR*);
template void Galfit<double>::setup (Cube<double>*, Rings<double>*, GALFIT_PAR*);
    

template <class T>
void Galfit<T>::galfit() {

    using namespace std;
    verb = in->pars().isVerbose();
    
    static int n=0;
    n = n==1 ? 2 : 1;
    std::string fileo = in->pars().getOutfolder()+"ringlog"+to_string(n)+".txt";
    remove(fileo.c_str());

    std::ofstream fileout;
    fileout.open(fileo.c_str(), std::ios_base::app);
    
    int m=10;
    fileout << left << setw(m) << "#RAD(Kpc)"
            << setw(m+1) << "RAD(arcs)"
            << setw(m+1) << "VROT(km/s)"
            << setw(m+1) << "DISP(km/s)"
            << setw(m)   << "INC(deg)" 
            << setw(m)   << "P.A.(deg)" 
            << setw(m)   << "Z0(pc)"
            << setw(m)   << "Z0(arcs)"
            << setw(m)   << "SIG(E20)"
            << setw(m)   << "XPOS(pix)"
            << setw(m)   << "YPOS(pix)"
            << setw(m+1) << "VSYS(km/s)"
            << setw(m+1) << "VRAD(km/s)";
            
    if (par.flagERRORS) {       
        if (mpar[VROT])  fileout << setw(m) << "E_VROT1" << setw(m) << "E_VROT2";
        if (mpar[VDISP]) fileout << setw(m) << "E_DISP1" << setw(m) << "E_DISP2";
        if (mpar[DENS])  fileout << setw(m) << "E_DENS1" << setw(m) << "E_DENS2";
        if (mpar[Z0])    fileout << setw(m) << "E_Z01"   << setw(m) << "E_Z02";
        if (mpar[INC])   fileout << setw(m) << "E_INC1"  << setw(m) << "E_INC2";
        if (mpar[PA])    fileout << setw(m) << "E_PA1"   << setw(m) << "E_PA2";
        if (mpar[XPOS])  fileout << setw(m) << "E_XPOS1" << setw(m) << "E_XPOS2";
        if (mpar[YPOS])  fileout << setw(m) << "E_YPOS1" << setw(m) << "E_YPOS2";
        if (mpar[VSYS])  fileout << setw(m) << "E_VSYS1" << setw(m) << "E_VSYS2";
        if (mpar[VRAD])  fileout << setw(m) << "E_VRAD1" << setw(m) << "E_VRAD2";
    }
        
    fileout << endl; 
        
    if (verb) { 
        in->pars().setVerbosity(false);
        cout << showpoint << fixed << setprecision(2) << endl ;
        cout << setfill('=') << setw(40) << right << " 3DFIT " << setw(34) << " ";
        cout << setfill(' ') << endl;
    }

    *outr = *inr;
//    global=false;
//    if (global) {
//        Rings<T> *dring = new Rings<T>;
//        *dring = *inr;
//        T minimum=0;
//        T pmin[nfree*dring->nr];
//        if (!minimize(dring, minimum, pmin)) cout << "DIOCANEEEEEE" << endl;

//        int k=0;
//        if (mpar[VROT])  for (int ir=0; ir<inr->nr; ir++) outr->vrot[ir]=pmin[k++];
//        if (mpar[VDISP]) for (int ir=0; ir<inr->nr; ir++) outr->vdisp[ir]=pmin[k++];
//        if (mpar[DENS])  for (int ir=0; ir<inr->nr; ir++) outr->dens[ir]=pmin[k++];
//        if (mpar[Z0])    for (int ir=0; ir<inr->nr; ir++) outr->z0[ir]=pmin[k++];
//        if (mpar[INC])   for (int ir=0; ir<inr->nr; ir++) outr->inc[ir]=pmin[k++];
//        if (mpar[PA])    for (int ir=0; ir<inr->nr; ir++) outr->phi[ir]=pmin[k++];
//        if (mpar[XPOS])  for (int ir=0; ir<inr->nr; ir++) outr->xpos[ir]=pmin[k++];
//        if (mpar[YPOS])  for (int ir=0; ir<inr->nr; ir++) outr->ypos[ir]=pmin[k++];
//        if (mpar[VSYS])  for (int ir=0; ir<inr->nr; ir++) outr->vsys[ir]=pmin[k++];
//        if (mpar[VRAD])  for (int ir=0; ir<inr->nr; ir++) outr->vrad[ir]=pmin[k++];

//        for (int ir=0; ir<inr->nr; ir++) {
//            float radius = ir==0 ?  inr->radsep/4. : inr->radii[ir];
//            double toKpc = KpcPerArc(distance);
//            fileout << setprecision(3) << fixed << left;
//            fileout << setw(m) << radius*toKpc
//                    << setw(m) << radius
//                    << setw(m+1) << outr->vrot[ir]
//                    << setw(m+1) << outr->vdisp[ir]
//                    << setw(m) << outr->inc[ir]
//                    << setw(m) << outr->phi[ir]
//                    << setw(m) << outr->z0[ir]*toKpc*1000
//                    << setw(m) << outr->z0[ir]
//                    << setw(m) << outr->dens[ir]/1E20
//                    << setw(m) << outr->xpos[ir]
//                    << setw(m) << outr->ypos[ir]
//                    << setw(m+1) << outr->vsys[ir]
//                    << setw(m+1) << outr->vrad[ir]
//                    << endl;

//        }
//    }
//    else {
    
    T ***errors = allocate_3D<T>(inr->nr,2,nfree);
    bool fitok[inr->nr];
    double toKpc = KpcPerArc(distance);
    int start_rad = par.STARTRAD<inr->nr ? par.STARTRAD : 0;
    
    int nthreads = in->pars().getThreads();
    
#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    for (int ir=start_rad; ir<inr->nr; ir++) {
        
        if (verb && nthreads==1) {
            time_t t = time(NULL);
            char Time[11] = "          ";
            strftime (Time,11,"[%T]",localtime(&t));
            cout << fixed << setprecision(2)<<"\n Working on ring #"
                 << ir+1 << " at radius " << inr->radii[ir] << " arcsec ("
                 << inr->radii[ir]*toKpc << " Kpc)... " << Time << std::endl;
        }

        Rings<T> *dring = new Rings<T>;
        dring->nr = 2;
        dring->id = ir;
        
        float width1=0, width2=0;
        // Handling the case of a single ring
        if (inr->nr==1) width1 = width2 = inr->radii[0];
        else {
            if (ir==0) width1 = width2 = (inr->radii[1]-inr->radii[0])/2.;
            else if (ir==inr->nr-1) width1 = width2 = (inr->radii[ir]-inr->radii[ir-1])/2.;
            else {
                width1 = (inr->radii[ir]-inr->radii[ir-1])/2.;
                width2 = (inr->radii[ir+1]-inr->radii[ir])/2.;
            }
        }
        
        dring->radii.push_back(max(double(inr->radii[ir]-width1),0.));
        dring->radii.push_back(max(double(inr->radii[ir]+width2),0.));

        for (int i=0; i<dring->nr; i++) {
            dring->vrot.push_back(inr->vrot[ir]);
            dring->vdisp.push_back(inr->vdisp[ir]);
            dring->dens.push_back(inr->dens[ir]);
            dring->z0.push_back(inr->z0[ir]);
            dring->inc.push_back(inr->inc[ir]);
            dring->phi.push_back(inr->phi[ir]);
            dring->xpos.push_back(inr->xpos[ir]);
            dring->ypos.push_back(inr->ypos[ir]);
            dring->vsys.push_back(inr->vsys[ir]);
            dring->vrad.push_back(inr->vrad[ir]);
            dring->vvert.push_back(inr->vvert[ir]);
            dring->dvdz.push_back(inr->dvdz[ir]);
            dring->zcyl.push_back(inr->zcyl[ir]);
        }

        T minimum=0;
        T pmin[nfree];

        fitok[ir] = minimize(dring, minimum, pmin);
        if (!fitok[ir]) continue;

        int k=0;
        if (mpar[VROT])  outr->vrot[ir]=pmin[k++];
        if (mpar[VDISP]) outr->vdisp[ir]=pmin[k++];
        if (mpar[DENS])  outr->dens[ir]=pmin[k++];
        if (mpar[Z0])    outr->z0[ir]=pmin[k++];
        if (mpar[INC])   outr->inc[ir]=pmin[k++];
        if (mpar[PA])    outr->phi[ir]=pmin[k++];
        if (mpar[XPOS])  outr->xpos[ir]=pmin[k++];
        if (mpar[YPOS])  outr->ypos[ir]=pmin[k++];
        if (mpar[VSYS])  outr->vsys[ir]=pmin[k++];
        if (mpar[VRAD])  outr->vrad[ir]=pmin[k++];

        if (verb) {
#pragma omp critical (galfit_outmsg)
{
            if (nthreads>1) {
                cout << "\n Ring #" << ir+1 << " at radius " << inr->radii[ir] << " arcsec ("
                     << inr->radii[ir]*toKpc << " Kpc)... \n";
            }
            
            int m=8, n=11;
            cout << "  Best parameters for ring #" << ir+1
                 << " (fmin = " << setprecision(6) << minimum << "):\n";

            cout << setprecision(2);

            string s;
            s = "    Vrot";
            if (!mpar[VROT]) s += "(f)";
            cout << setw(n) << left << s << setw(3) << right << "= "
                 << setw(m) << outr->vrot[ir] << left << setw(m) << "  km/s";

            s = "        Disp";
            if (!mpar[VDISP]) s += "(f)";
            cout << setw(n+4) << left << s << setw(3) << right << "= "
                 << setw(m-1) << outr->vdisp[ir]
                 << left << setw(m) << "  km/s" << endl;

            s = "    Vrad";
            if (!mpar[VRAD]) s += "(f)";
            cout << setw(n) << left << s << setw(3) << right << "= "
                 << setw(m) << outr->vrad[ir] << left << setw(m) << "  km/s";
            
            s = "        Vsys";
            if (!mpar[VSYS]) s += "(f)";
            cout << setw(n+4) << left << s << setw(3) << right << "= "
                 << setw(m-1) << outr->vsys[ir] << left << setw(m) << "  km/s" << endl;
                
            s = "    Inc";
            if (!mpar[INC]) s += "(f)";
            cout << setw(n) << left << s << setw(3) << right << "= "
                 << setw(m) << outr->inc[ir] << left << setw(m) << "  deg";

            s = "        PA";
            if (!mpar[PA]) s += "(f)";
            cout << setw(n+4) << left << s << setw(3) << right << "= "
                 << setw(m-1) << outr->phi[ir] << left << setw(m) << "  deg" << endl;

            s = "    Xpos";
            if (!mpar[XPOS]) s += "(f)";
            cout << setw(n) << left << s << setw(3) << right << "= "
                 << setw(m) << outr->xpos[ir] << left << setw(m) << "  pix";

            s = "        Ypos";
            if (!mpar[YPOS]) s += "(f)";
            cout << setw(n+4) << left << s << setw(3) << right << "= "
                 << setw(m-1) << outr->ypos[ir] << left << setw(m) << "  pix" << endl;
                
            s = "    Z0";
            if (!mpar[Z0]) s += "(f)";
            cout << setw(n) << left << s << setw(3) << right << "= "
                 << setw(m) << outr->z0[ir]*toKpc << left << setw(m) << "  Kpc";

            //s = "        CD";----
            //if (!mpar[DENS]) s += "(f)";
            //cout << setw(n+4) << left << s << setw(3) << right << "= "
            //     << setw(m-1) << scientific << setprecision(1)
            //     << outr->dens[ir] << left << setw(m) << "  a/cm2";
            
            cout << endl;
        }
}

        if (par.flagERRORS) getErrors(dring,errors[ir],ir,minimum);
        
#pragma omp critical (galfit_write)
{
        // Writing output file. Not ordered if multithread
        fileout << setprecision(3) << fixed << left;
        fileout << setw(m) << outr->radii[ir]*toKpc
                << setw(m+1) << outr->radii[ir]
                << setw(m+1) << outr->vrot[ir]
                << setw(m+1) << outr->vdisp[ir]
                << setw(m) << outr->inc[ir]
                << setw(m) << outr->phi[ir]
                << setw(m) << outr->z0[ir]*toKpc*1000
                << setw(m) << outr->z0[ir]
                << setw(m) << outr->dens[ir]/1E20
                << setw(m) << outr->xpos[ir]
                << setw(m) << outr->ypos[ir]
                << setw(m+1) << outr->vsys[ir]
                << setw(m+1) << outr->vrad[ir];
            
        if (par.flagERRORS) 
            for (int kk=0; kk<nfree; kk++) 
                fileout << setw(m) << errors[ir][0][kk] << setw(m) << errors[ir][1][kk];
        
        fileout << endl;
}
        delete dring;
        
    }
  //  }
         
    fileout.close();

    // If multi-threads rewrite ordered outfile
    if (nthreads>1) {
        fileout.open(fileo.c_str());
    
        fileout << left << setw(m) << "#RAD(Kpc)"
                << setw(m+1)   << "RAD(arcs)"
                << setw(m+1) << "VROT(km/s)"
                << setw(m+1) << "DISP(km/s)"
                << setw(m)   << "INC(deg)" 
                << setw(m)   << "P.A.(deg)" 
                << setw(m)   << "Z0(pc)"
                << setw(m)   << "Z0(arcs)"
                << setw(m)   << "SIG(E20)"
                << setw(m)   << "XPOS(pix)"
                << setw(m)   << "YPOS(pix)"
                << setw(m+1) << "VSYS(km/s)"
                << setw(m+1) << "VRAD(km/s)";
            
        if (par.flagERRORS) {       
            if (mpar[VROT])  fileout << setw(m) << "E_VROT1" << setw(m) << "E_VROT2";
            if (mpar[VDISP]) fileout << setw(m) << "E_DISP1" << setw(m) << "E_DISP2";
            if (mpar[DENS])  fileout << setw(m) << "E_DENS1" << setw(m) << "E_DENS2";
            if (mpar[Z0])    fileout << setw(m) << "E_Z01"   << setw(m) << "E_Z02";
            if (mpar[INC])   fileout << setw(m) << "E_INC1"  << setw(m) << "E_INC2";
            if (mpar[PA])    fileout << setw(m) << "E_PA1"   << setw(m) << "E_PA2";
            if (mpar[XPOS])  fileout << setw(m) << "E_XPOS1" << setw(m) << "E_XPOS2";
            if (mpar[YPOS])  fileout << setw(m) << "E_YPOS1" << setw(m) << "E_YPOS2";
            if (mpar[VSYS])  fileout << setw(m) << "E_VSYS1" << setw(m) << "E_VSYS2";
            if (mpar[VRAD])  fileout << setw(m) << "E_VRAD1" << setw(m) << "E_VRAD2";
        }
        
        fileout << endl; 
        
        for (int ir=start_rad; ir<inr->nr; ir++) {
            if (!fitok[ir]) continue;
            fileout << setprecision(3) << fixed << left;
            fileout << setw(m) << outr->radii[ir]*toKpc
                    << setw(m+1) << outr->radii[ir]
                    << setw(m+1) << outr->vrot[ir]
                    << setw(m+1) << outr->vdisp[ir]
                    << setw(m) << outr->inc[ir]
                    << setw(m) << outr->phi[ir]
                    << setw(m) << outr->z0[ir]*toKpc*1000
                    << setw(m) << outr->z0[ir]
                    << setw(m) << outr->dens[ir]/1E20
                    << setw(m) << outr->xpos[ir]
                    << setw(m) << outr->ypos[ir]
                    << setw(m+1) << outr->vsys[ir]
                    << setw(m+1) << outr->vrad[ir];
            
            if (par.flagERRORS) 
                for (int kk=0; kk<nfree; kk++) 
                    fileout << setw(m) << errors[ir][0][kk] << setw(m) << errors[ir][1][kk];
        
            fileout << endl;
        }
        fileout.close();
    }

    deallocate_3D<T>(errors,inr->nr,2);

    if (verb) {               
        cout << setfill('=') << setw(74) << " " << endl << endl; 
        cout << fixed << setprecision(2) << setfill(' ');
        in->pars().setVerbosity(true);
    }
    
}
template void Galfit<float>::galfit();
template void Galfit<double>::galfit();


template <class T> 
bool Galfit<T>::SecondStage() {
    
    bool isNeeded = mpar[INC] || mpar[PA]   || mpar[Z0] ||
                    mpar[XPOS]|| mpar[YPOS] || mpar[VSYS];
    if (!isNeeded) {par.TWOSTAGE=false; return isNeeded;}
    second = true;
    if (!in->pars().getFlagSlitfit()) func_norm = &Model::Galfit<T>::norm_local;

    int n = outr->nr;
    
    std::vector<T> xx, yy;

    if (anglepar==-1) {
        // Make a bezier interpolation
        std::vector<T> x_bez(n), y_bez_inc(n), y_bez_pa(n), y_bez_z0(n);
        if (mpar[INC]) {
            for (int i=0; i<n; i++) {
                xx.push_back(outr->radii[i]);
                yy.push_back(outr->inc[i]);
            }
            bezier_interp(xx,yy,x_bez,y_bez_inc);
            xx.clear();
            yy.clear();
        }
        if (mpar[PA]) {
            for (int i=0; i<n; i++) {
                xx.push_back(outr->radii[i]);
                yy.push_back(outr->phi[i]);
            }
            bezier_interp(xx,yy,x_bez,y_bez_pa);
            xx.clear();
            yy.clear();
        }
        if (mpar[Z0]) {
            for (int i=0; i<n; i++) {
                xx.push_back(outr->radii[i]);
                yy.push_back(outr->z0[i]);
            }
            bezier_interp(xx,yy,x_bez,y_bez_z0);
            xx.clear();
            yy.clear();
        }

        T xpos_av = findMedian(outr->xpos, n);
        T ypos_av = findMedian(outr->ypos, n);
        T vsys_av = findMedian(outr->vsys, n);

        *inr = *outr;
        for (int i=0; i<n; i++) {
            if (mpar[INC]) inr->inc[i]=y_bez_inc[i];
            if (mpar[PA])  inr->phi[i]=y_bez_pa[i];
            if (mpar[Z0])  inr->z0[i]=y_bez_z0[i];
            inr->xpos[i]=xpos_av;
            inr->ypos[i]=ypos_av;
            inr->vsys[i]=vsys_av;
        }
    }
    else {
        // Make a polynomial fit
        if (n<=anglepar) {
            std::cout << "3DFIT ERROR - Second stage: too few degree of freedom.\n";
            return false;
        }
    
        cout << setprecision(4);
    
        Statistics::Stats<T> stats;
        typename std::vector<T>::iterator where;
        T ww[n], cinc[anglepar], cincout[anglepar];
        T cpa[anglepar], cpaerr[anglepar];
        T cz0[anglepar], cz0err[anglepar];
    
        bool mp[anglepar];
        for (int i=0; i<anglepar;i++) mp[i] = true;
        for (int i=0; i<n; i++) ww[i] = 1;
                
        if (mpar[INC]) {
            for (int i=0; i<n; i++) {
                xx.push_back(outr->radii[i]);
                yy.push_back(outr->inc[i]);
            }
            if (n-2>anglepar) {
                where = std::max_element(yy.begin(), yy.end());
                yy.erase(where);
                xx.erase(xx.begin()+std::distance(yy.begin(), where));
                where = std::min_element(yy.begin(), yy.end());
                yy.erase(where);
                xx.erase(xx.begin()+std::distance(yy.begin(), where));
            }
        
            stats.calculate(&yy[0], yy.size());

            Lsqfit<T> lsq(&xx[0],1,&yy[0],ww,xx.size(),cinc,cincout,mp,anglepar,&polyn,&polynd);
            int nrt = lsq.fit();
            if (nrt<0) {
                std::cout << "3DFIT error: cannot least-square fit the inclination.\n";
                par.TWOSTAGE = false;
                return false;
            }
            if (verb) {
                cout << "  Best parameters for inclination:\n";
                for (int i=0; i<anglepar; i++) {
                    string a = "a"+to_string<int>(i)+" = ";
                    cout << setw(43) << right << a << setw(9) << cinc[i] << setw(4)
                         << " ±" << setw(9) << cincout[i] << endl;
                }
            }
            xx.clear();
            yy.clear();
        }
        if (mpar[PA]) {
            for (int i=0; i<n; i++) {
                xx.push_back(outr->radii[i]);
                yy.push_back(outr->phi[i]);
            }
            if (n-2>anglepar) {
                where = std::max_element(yy.begin(), yy.end());
                yy.erase(where);
                xx.erase(xx.begin()+std::distance(yy.begin(), where));
                where = std::min_element(yy.begin(), yy.end());
                yy.erase(where);
                xx.erase(xx.begin()+std::distance(yy.begin(), where));
            }
            Lsqfit<T> lsq(&xx[0],1,&yy[0],ww,xx.size(),cpa,cpaerr,mp,anglepar,&polyn,&polynd);
            int nrt = lsq.fit();
            if (nrt<0) {
                std::cout << "3DFIT error: cannot least-square fit the position angle.\n";
                par.TWOSTAGE = false;
                return false;
            }
            if (verb) {
                cout << "  Best parameters for position angle:\n";
                for (int i=0; i<anglepar; i++) {
                    std::string a = "a"+to_string<int>(i)+" = ";
                    cout << setw(43) << right << a << setw(9) << cpa[i] << setw(4)
                         << " ±" << setw(9) << cpaerr[i] << endl;
                }
            }
            xx.clear();
            yy.clear();
        }
        if (mpar[Z0]) {
            for (int i=0; i<n; i++) {
                xx.push_back(outr->radii[i]);
                yy.push_back(outr->z0[i]);
            }
            if (n-2>anglepar) {
                where = std::max_element(yy.begin(), yy.end());
                yy.erase(where);
                xx.erase(xx.begin()+std::distance(yy.begin(), where));
                where = std::min_element(yy.begin(), yy.end());
                yy.erase(where);
                xx.erase(xx.begin()+std::distance(yy.begin(), where));
            }
            Lsqfit<T> lsq(&xx[0],1,&yy[0],ww,xx.size(),cz0,cz0err,mp,anglepar,&polyn,&polynd);
            int nrt = lsq.fit();
            if (nrt<0) {
                std::cout << "3DFIT error: cannot least-square fit the thickness.\n";
                par.TWOSTAGE = false;
                return false;
            }
            if (verb) {
                cout << "  Best parameters for thickness:\n";
                for (int i=0; i<anglepar; i++) {
                    std::string a = "a"+to_string<int>(i)+" = ";
                    cout << setw(43) << right << a << setw(9) << cz0[i] << setw(4)
                         << " ±" << setw(9) << cz0err[i] << endl;
                }
            }
            xx.clear();
            yy.clear();
        }

        T xpos_av = findMedian(outr->xpos, n);
        T ypos_av = findMedian(outr->ypos, n);
        T vsys_av = findMedian(outr->vsys, n);
    
        *inr = *outr;
        for (int i=0; i<n; i++) {
            if (mpar[INC]) {
                inr->inc[i]=0;
                for (int j=0; j<anglepar; j++)
                    inr->inc[i] += cinc[j]*std::pow(double(inr->radii[i]),j);
            }
            if (mpar[PA]) {
                inr->phi[i]=0;
                for (int j=0; j<anglepar; j++)
                    inr->phi[i] += cpa[j]*std::pow(double(inr->radii[i]),j);
            }
            if (mpar[Z0]) {
                inr->z0[i]=0;
                for (int j=0; j<anglepar; j++)
                    inr->z0[i] += cz0[j]*std::pow(double(inr->radii[i]),j);
            }
            inr->xpos[i]=xpos_av;
            inr->ypos[i]=ypos_av;
            inr->vsys[i]=vsys_av;
        }
    }

    cout << setprecision(2);
    
    bool oldmpar[MAXPAR];
    for (int i=MAXPAR; i--;) oldmpar[i]=mpar[i];
    int oldnfree =nfree;
    
    mpar[INC] = mpar[PA] = mpar[VSYS] = false;
    mpar[XPOS]= mpar[YPOS] = mpar[Z0] = false;
    nfree = mpar[VROT]+mpar[VDISP]+mpar[VRAD];

    galfit();
    
    for (int i=MAXPAR; i--;) mpar[i]=oldmpar[i];
    nfree = oldnfree;

    return isNeeded;
}
template bool Galfit<float>::SecondStage();
template bool Galfit<double>::SecondStage();


template <class T> 
bool Galfit<T>::setCfield() {
    
    T pixsizeX = fabs(in->Head().Cdelt(0))*arcconv; 
    T pixsizeY = fabs(in->Head().Cdelt(1))*arcconv; 
    
    //Beam Old = {pixsizeX, pixsizeY, 0};
    Beam Old = {0, 0, 0};

    Beam New = {in->Head().Bmaj()*arcconv,
                in->Head().Bmin()*arcconv,
                in->Head().Bpa()};
/*
    if (Old.bmaj<Old.bmin) {
        std::cout << "Old beam major axis < minor axis. Inverting...";
        std::swap(Old.bmaj, Old.bmin);
    }
*/
    bool agreed = ((New.bmaj>=Old.bmaj) && (New.bmin>=Old.bmin)); 
 
    if (!agreed) {
        std::cout << "3DFIT error: new beam smaller than old beam\n";
        return false;
    }
    if (New.bmaj<New.bmin) {
        std::cout << "New beam major axis < minor axis. Inverting...";
        std::swap(New.bmaj, New.bmin);
    }

    // Now calculate the convolution beam;
    double a2  = Old.bmaj/2;
    double b2  = Old.bmin/2;
    double a0  = New.bmaj/2;
    double b0  = New.bmin/2;
    double D0  = a0*a0-b0*b0;
    double D2  = a2*a2-b2*b2;    
    double th2 = Old.bpa*atan(1.)/45.;
    double th0 = New.bpa*atan(1.)/45.;
    double D1  = sqrt(D0*D0+D2*D2-2*D0*D2*cos(2*(th0-th2)));    
    double a1, b1, th1;
    
    double arg = 0.5*(a0*a0+b0*b0-a2*a2-b2*b2+D1); 
    if (arg<0) {
        std::cout << "3DFIT error: unsuitable new beam parameters!\n";
    return false;
    }
    else a1 = sqrt(arg);
      
    arg = 0.5*(a0*a0+b0*b0-a2*a2-b2*b2-D1); 
    if (arg<0) {
        std::cout << "3DFIT error: unsuitable new beam parameters!\n";
    return false;
    }
    else b1 = sqrt(arg);
    
    double nom   = D0*sin(2*th0)-D2*sin(2*th2);
    double denom = D0*cos(2*th0)-D2*cos(2*th2); 
    if (denom==0 && nom==0) th1=0;
    else {
        T twoth1 = atan2(nom,denom);
        th1 = twoth1/2;          
    }
    
    Beam Con = {2*a1, 2*b1, (th1*180./M_PI)-90.};

    // Building the convolution field.
    double phi = Con.bpa*M_PI/180.;
    double cs = cos(phi);
    double sn = sin(phi);  
    double beam[2] = {Con.bmaj, Con.bmin};

    double xr = 0.5*beam[0];
    double yr = 0.5*beam[1];
    double extend = sqrt(-1.0*log(1E-04)/log(2.0));
    xr *= extend;
    yr *= extend;
   
    double x1 = fabs(xr*cs-0.0*sn);
    double y1 = fabs(xr*sn+0.0*cs);
    double x2 = fabs(0.0*cs-yr*sn);
    double y2 = fabs(0.0*sn+yr*cs);
    double x  = (x2>x1 ? x2:x1);
    double y  = (y2>y1 ? y2:y1);
   
    int Xmax = lround(x/pixsizeX); 
    int Ymax = lround(y/pixsizeY); 

    NconX = 2*Xmax+1;    
    NconY = 2*Ymax+1; 
    
    cfield = new double[NconX*NconY];
    cfieldAllocated=true;       
      
    double argfac = -4.0 * log(2.0);
    double totalarea = 0;
    for (int j=-Ymax; j<=Ymax; j++) {
        for (int i=-Xmax; i<=Xmax; i++) {
            int pos = (j+Ymax)*NconX+(i+Xmax);    
            x = i*pixsizeX;
            y = j*pixsizeY;
            xr = x*cs + y*sn;
            yr = -1.0*x*sn + y*cs;
            double argX=0;
            double argY=0;
            if (beam[0]!=0) argX = xr/beam[0];
            if (beam[1]!=0) argY = yr/beam[1];
            double arg = argfac*(argX*argX+argY*argY);
            double c = exp(arg);
            if (c>=1E-04) {
                cfield[pos] = c;
                totalarea += c;
            } 
            else cfield[pos] = 0;
        }
    }
   
    for (int i=0;i<NconX*NconY;i++) cfield[i] = cfield[i]/totalarea;

    return true;
}
template bool Galfit<float>::setCfield();
template bool Galfit<double>::setCfield();


template <class T>
Model::Galmod<T>* Galfit<T>::getModel() {

    Model::Galmod<T> *mod = new Model::Galmod<T>;
    int bhi[2] = {in->DimX(), in->DimY()};
    int blo[2] = {0,0};
    int nv = par.NV;
    if (nv==-1) nv=in->DimZ();
    mod->input(in,bhi,blo,outr,nv,par.LTYPE,1,par.CDENS);
    mod->calculate();
    if (par.SM) mod->smooth();
    return mod;
}
template Model::Galmod<float>* Galfit<float>::getModel();
template Model::Galmod<double>* Galfit<double>::getModel();


template <class T>
void Galfit<T>::setFree() {

    std::string FREE = par.FREE;
    FREE = makelower(FREE);

    int found = FREE.find("vrot");
    if (found<0) mpar[VROT]=false;
    else mpar[VROT]=true;

    found = FREE.find("vdisp");
    if (found<0) mpar[VDISP]=false;
    else mpar[VDISP]=true;

    found = FREE.find("dens");
    if (found<0) mpar[DENS]=false;
    else mpar[DENS]=true;

    found = FREE.find("z0");
    if (found<0) mpar[Z0]=false;
    else mpar[Z0]=true;

    found = FREE.find("inc");
    if (found<0) mpar[INC]=false;
    else mpar[INC]=true;

    found = FREE.find("pa");
    if (found<0) {
        found = FREE.find("phi");
        if (found<0) mpar[PA]=false;
        else mpar[PA]=true;
    }
    else mpar[PA]=true;

    found = FREE.find("xpos");
    if (found<0) mpar[XPOS]=false;
    else mpar[XPOS]=true;

    found = FREE.find("ypos");
    if (found<0) mpar[YPOS]=false;
    else mpar[YPOS]=true;

    found = FREE.find("vsys");
    if (found<0) mpar[VSYS]=false;
    else mpar[VSYS]=true;
    
    found = FREE.find("vrad");
    if (found<0) mpar[VRAD]=false;
    else mpar[VRAD]=true;

    found = FREE.find("all");
    if (found>=0)
        for (int i=0; i<MAXPAR; i++) mpar[i] = true;

    mpar[DENS] = false;
    
    int nfixed=0;
    for (int i=0; i<MAXPAR; i++) nfixed += (1-mpar[i]);
    if (nfixed == MAXPAR) {
        std::cout << " 3DFIT ERROR: NO free parameters!\n";
        std::terminate();
    }
    nfree = MAXPAR-nfixed;
    
}
template void Galfit<float>::setFree();
template void Galfit<double>::setFree();


template <class T>
bool Galfit<T>::AsymmetricDrift(T *rad, T *densprof, T *dispprof, T *inc, int nn) {
    
    // Compute an asymmetric drift correction, following procedure in Iorio+17, sec 4.3
        
    // Fitting dispersion with a third degree polynomial
    int npar1 = 4;
    T cdisp[npar1], cdisperr[npar1], ww[nn];
    bool mp[npar1];
    for (int i=0; i<npar1; i++) mp[i] = true;
    for (int i=0; i<nn; i++) ww[i] = 1;
    Lsqfit<T> lsq1(rad,1,dispprof,ww,nn,cdisp,cdisperr,mp,npar1,&polyn,&polynd);
    int nrt = lsq1.fit();
    if (nrt<0) {
        if (in->pars().isVerbose())  std::cout << "3DFIT ERROR: cannot least-square fit the dispersion for asymmetric drift.\n";
        par.flagADRIFT = false;
        return false;
    }
    
    // Now fitting density*disp2 with a exponential function (line in log space)
    T *fun = new T[nn];
    int npar2 = 2;
    T cfun[npar2], cfunerr[npar2];
    bool mpp[npar2];
    for (int i=0; i<npar2; i++) mpp[i] = true;
    for (int i=0; i<nn; i++) {
        fun[i] = log(dispprof[i]*dispprof[i]*densprof[i]*cos(inc[i]*M_PI/180.));
        ww[i] = 1;
    }
    //Lsqfit<T> lsq2(rad,1,fun,ww,nn,cfun,cfunerr,mpp,npar2,&coreExp,&coreExpd);
    Lsqfit<T> lsq2(rad,1,fun,ww,nn,cfun,cfunerr,mpp,2,&polyn,&polynd);
    nrt = lsq2.fit();
    if (nrt<0) {
        if (in->pars().isVerbose()) std::cout << "3DFIT ERROR: cannot least-square fit the fun for asymmetric drift.\n";
        par.flagADRIFT = false;
        return false;
    }
    
    // Now writing to a text file
    std::ofstream fout(in->pars().getOutfolder()+"asymdrift.txt");
    
    int m=16;
    fout << fixed << setprecision(4);
    fout << "#" << setw(m-1) << "RAD(arcs)"
         << setw(m) << "ASYMDRIFT(km/s)"
         << setw(m) << "DISP_REG(km/s)"
         << setw(m) << "FUN" 
         << setw(m) << "FUN_REG\n";
    
    //T a1 = exp(cfun[0]);
    T a2 = 0;
    T a3 = -1/cfun[1];
    
    for (int i=0; i<nn; i++) {
        T disp_reg = polyn(&rad[i],cdisp,npar1);
        T fun_reg = polyn(&rad[i],cfun,npar2);
        T expn = exp(rad[i]/a3);
        T asdrift = sqrt(rad[i]*disp_reg*disp_reg*expn/(a3*(a2+expn)));
        fout << setw(m) << rad[i] << setw(m) << asdrift
             << setw(m) << disp_reg << setw(m) << fun[i] << setw(m) << fun_reg << std::endl;
    }
    
    fout.close();
    
    delete [] fun;
    
    return true;
    
}
template bool Galfit<float>::AsymmetricDrift(float*,float*,float*,float*,int);
template bool Galfit<double>::AsymmetricDrift(double*,double*,double*,double*,int);



template <class T>
T* Galfit<T>::EstimateInitial(Cube<T> *c, GALFIT_PAR *p){
    
    if (!c->getIsSearched()) c->Search();
    Detection<T> *largest = c->LargestDetection();

    if (largest==NULL) {
        std::cout << " 3DFIT ERROR: No sources detected in the datacube. Cannot fit!!! \n";
        std::terminate();
    }

    ParamGuess<T> *ip = new ParamGuess<T>(c,largest);
    ip->findInitial();

    string pos[2] = {p->XPOS, p->YPOS};
    double *pixs = getCenterCoordinates(pos, c->Head());
    if (p->XPOS!="-1" && p->YPOS!="-1") {
        ip->setXcentre(pixs[0]);
        ip->setYcentre(pixs[1]);
    }
    if (p->PHI!="-1")  ip->setPosang(atof(p->PHI.c_str()));

    ip->fitEllipse();
    if (c->pars().getFlagDebug()) ip->fitIncfromMap();
    if (c->pars().getFlagDebug()) ip->plotGuess();
    
    T *init_par = new T[8];
    init_par[0] = ip->nrings;
    init_par[1] = ip->radsep; 
    init_par[2] = ip->xcentre;
    init_par[3] = ip->ycentre;
    init_par[4] = ip->vsystem;
    init_par[5] = ip->vrot;
    init_par[6] = ip->inclin;
    init_par[7] = ip->posang;
    
    delete ip;
    return init_par;
}
template float* Galfit<float>::EstimateInitial(Cube<float>*,GALFIT_PAR*);
template double* Galfit<double>::EstimateInitial(Cube<double>*,GALFIT_PAR*);


}

