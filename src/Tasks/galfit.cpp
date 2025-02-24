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
#include <Tasks/ringmodel.hh>
#include <Utilities/utils.hh>
#include <Utilities/allocator.hpp>
#include <Utilities/lsqfit.hh>
#include <Utilities/progressbar.hh>
#include <Utilities/paramguess.hh>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Model {


template <class T>
Galfit<T>::~Galfit () {
    if (outDefined) delete outr;
    if (inDefined) delete inr;
    if (line_imDefined) delete line_im;
    if (cfieldAllocated) delete [] cfield; 
    if (chan_noiseAllocated) delete [] chan_noise;
    if (mask2D!=nullptr) delete [] mask2D;
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
    this->mask2D = g.mask2D;

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
    this->reverse   = g.reverse;

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
    
    /// This constructor reads all needed parameters from the Cube object.
    /// Cube object must contain a Param method with all information.

    in   = c;
    par  = c->pars().getParGF();
    verb = c->pars().isVerbose();
    
    c->checkBeam();

    // Read all rings from a BBarolo output file, if given as input
    if (!par.ringfile.empty()) {
        if (par.RADII=="-1" && (par.NRADII==-1 || par.RADSEP==-1))
                par.RADII = "file("+par.ringfile+",2)";
        if (par.VROT=="-1") par.VROT  = "file("+par.ringfile+",3)";
        if (par.VDISP=="-1") par.VDISP = "file("+par.ringfile+",4)";
        if (par.INC=="-1") par.INC   = "file("+par.ringfile+",5)";
        if (par.PHI=="-1") par.PHI   = "file("+par.ringfile+",6)";
        if (par.Z0=="-1") par.Z0    = "file("+par.ringfile+",8)";
        if (par.DENS=="-1") par.DENS  = "file("+par.ringfile+",9)";
        if (par.XPOS=="-1") par.XPOS  = "file("+par.ringfile+",10)";
        if (par.YPOS=="-1") par.YPOS  = "file("+par.ringfile+",11)";
        if (par.VSYS=="-1") par.VSYS  = "file("+par.ringfile+",12)";
        if (par.VRAD=="-1") par.VRAD  = "file("+par.ringfile+",13)";
    }

    // Estimate initial parameters, if needed
    bool toEstimate =  (par.RADII=="-1" && (par.NRADII==-1 || par.RADSEP==-1)) ||
                        par.XPOS=="-1" || par.YPOS=="-1" || par.VSYS=="-1" ||
                        par.VROT=="-1" || par.PHI=="-1"  || par.INC=="-1";

    ParamGuess<T> *ip = nullptr;
    if (toEstimate) {
        ip = EstimateInitial(c,&par);
        if (par.NRADII==-1) par.NRADII = ip->nrings;
        if (par.RADSEP==-1) par.RADSEP = ip->radsep;
        if (par.XPOS=="-1") par.XPOS   = to_string(ip->xcentre);
        if (par.YPOS=="-1") par.YPOS   = to_string(ip->ycentre);
        if (par.VSYS=="-1") par.VSYS   = to_string(ip->vsystem);
        if (par.VROT=="-1") par.VROT   = to_string(ip->vrot);
        if (par.INC=="-1")  par.INC    = to_string(ip->inclin);
        if (par.PHI=="-1")  par.PHI    = to_string(ip->posang);
    }
    
    // Setting defaults for other parameters, if not given
    if (par.VDISP=="-1") par.VDISP = to_string(8.);
    if (par.Z0=="-1")    par.Z0    = to_string(0.);
    if (par.VRAD=="-1")  par.VRAD  = to_string(0.);
    if (par.DENS=="-1")  par.DENS  = to_string(1.);
    if (par.VVERT=="-1") par.VVERT = to_string(0.);
    if (par.DVDZ=="-1")  par.DVDZ  = to_string(0.);
    if (par.ZCYL=="-1")  par.ZCYL  = to_string(0.);

    // Now reading rings with general purpose function
    bool fromfile = false;
    Rings<T> *inR = readRings<T>(par,c->Head(),&fromfile);
        
    if (inR->nr==0) {
        std::cerr << "\n 3DFIT ERROR: The number of radii must be > 0! " << std::endl;
        std::terminate();
    }
    
    if (distance==-1) distance = VeltoDist(fabs(inr->vsys[0]));
    
    // Writing initial rings on screen and in a file
    std::string initfile;
    if (c->pars().getflagGalMod()) initfile = c->pars().getOutfolder()+"rings_model.txt";
    else {
        if (!fromfile && verb) showInitial(inR, std::cout);
        initfile = c->pars().getOutfolder()+"rings_initial.txt";
    }
    printInitial(inR,initfile);
    
    
    // Deciding whether to use reverse fitting based on galaxy inclination
    if (par.REVERSE.find("auto")!=std::string::npos && !c->pars().getflagGalMod()) {
        T incmed = findMedian<T>(&inR->inc[0],inR->nr);
        if (incmed>75) {
            if (verb) {
                //std::cout << "\n 3DFIT WARNING: because the galaxy is highly inclined, I will use a \n"
                //          << " reverse-cumulative fitting scheme. To turn this off, set REVERSE=false.\n";
                std::cout << "\n 3DFIT WARNING: because the galaxy is highly inclined, I recommend to try \n"
                          << " a reverse-cumulative fitting scheme (REVERSE=true) with ring width=beam size.\n";
            }
            reverse = true;
            // INC MUST set to fixed HERE
            std::cout << std::endl;
        }
    }
    else if (par.REVERSE.find("true")!=std::string::npos) reverse=true;
    else reverse = false;
///*
    if (reverse) {
        float beamarcs = c->Head().Bmaj()*arcsconv(c->Head().Cunit(0));
        if (inR->radsep<beamarcs) {
            if (verb) std::cout << " To this end, I will set ring widths >= beam major axis.\n";
            int new_nr = ceil(inR->radii.back()/beamarcs);
            Rings<T> *newr = new Rings<T>;
            newr->radsep = beamarcs;
            for (int i=0; i<new_nr; i++){
                newr->addRing(i*beamarcs+beamarcs/2.,inR->xpos[i],inR->ypos[i],inR->vsys[i],inR->vrot[i],inR->vdisp[i],
                          inR->vrad[i],inR->vvert[i],inR->dvdz[i],inR->zcyl[i],inR->dens[i],inR->z0[i],inR->inc[i],inR->phi[i]);
            }
            inR = newr;
        }
    }
//*/

    // A little trick: fit a 2D model first to get better PA variation
    if (toEstimate && makelower(par.FREE).find("pa")!=std::string::npos && inR->nr>=5 && !c->pars().getflagGalMod()) {

        double topix = c->Head().PixScale()*arcsconv(c->Head().Cunit(0));
        // Initializing rings
        T *radii = new T[inR->nr];
        T *wids  = new T[inR->nr];

        for (int i=0; i<inR->nr; i++) {
            radii[i] = inR->radii[i]/topix;
            wids[i]  = inR->radsep/topix;
        }
        // Initializing a Ringmodel instance
        Ringmodel<T> tr(inR->nr,radii,wids,&inR->vsys[0],&inR->vrot[0],&inR->vrad[0],
                       &inR->phi[0],&inR->inc[0],inR->xpos[0],inR->ypos[0]);

        tr.setfield(ip->Vemap,c->DimX(),c->DimY());
        // Setting free parameters. Order is VSYS, VROT, VEXP, PA, INC, X0, Y0
        // Fitting only VROT and PA
        bool mpar[7] = {false,true,false,true,false,false,false};
        tr.setoption(mpar,3,2,15.);
        // Fitting a tilted-ring model
        tr.ringfit(c->pars().getThreads(),false,false);
        // Writing this model to a file
        std::ofstream fileo(c->pars().getOutfolder()+c->pars().getOutPrefix()+"_2drings.txt");
        tr.printfinal(fileo,c->Head());
        //tr.printfinal(std::cout,c->Head());

        // Updating initial rings
        for (int i=1; i<inR->nr; i++) {
            if (!isNaN(tr.getPosaf(i)) && tr.getVrotf(i)>0) inR->phi[i] = tr.getPosaf(i);
            if (!isNaN(tr.getVrotf(i)) && tr.getVrotf(i)>0) inR->vrot[i] = tr.getVrotf(i);
        }
        inR->phi[0] = inR->phi[1];
        //inR->vrot[0] = inR->vrot[1];

        delete [] radii;
        delete [] wids;
    }

    if (ip!=nullptr) delete ip;

    // Setup all needed parameters
    setup(c,inR,&par);
  
}
template Galfit<float>::Galfit(Cube<float>*);
template Galfit<double>::Galfit(Cube<double>*);


template <class T>
Galfit<T>::Galfit (Cube<T> *c, Rings<T> *inrings, Param *p) {
    // This constructor is used mostly by pyBBarolo and it allows to set
    // separately the Cube, the Rings and all parameters used by the task (including masking etc...)
    c->pars() = *p;
    c->pars().getParGF() = p->getParGF();
    c->pars().checkPars();
    c->Head().setRedshift(p->getRedshift()); 
    c->Head().setWave0(p->getRestwave());

    // Create directory tree if it does not exist
    mkdirp(c->pars().getOutfolder().c_str());

    setup(c,inrings,&c->pars().getParGF());
}
template Galfit<float>::Galfit (Cube<float> *c, Rings<float> *inrings, Param *p);
template Galfit<double>::Galfit (Cube<double> *c, Rings<double> *inrings, Param *p);


template <class T>
void Galfit<T>::setup (Cube<T> *c, Rings<T> *inrings, GALFIT_PAR *p) {
    
    in = c;
    par = *p;
    inr = new Rings<T>;
    *inr = *inrings;
    inDefined = true;
    verb = c->pars().isVerbose();

    // Check that radii are ok.
    for (int ir=0; ir<inr->nr; ir++) {
        if (ir!=0) {
            if (inr->radii[ir]<=inr->radii[ir-1]) {
                std::cerr  << "3DFIT ERROR: Radii not in increasing order.\n";
                std::terminate();
            }
        }
        if (inr->radii[ir]<0) {
            std::cerr << "3DFIT ERROR: Negative radius!!!\n";
            std::terminate();
        }
    }
    
    // Checking that the beam has all information
    in->checkBeam();

    // Setting other 3DFIT variables
    verb = in->pars().isVerbose();
    arcconv = arcsconv(in->Head().Cunit(0));
    distance = par.DISTANCE==-1 ? VeltoDist(fabs(inr->vsys[0])) : par.DISTANCE;
    chan_noise = new float[in->DimZ()];
    chan_noiseAllocated = true;
    for (int z=0; z<in->DimZ(); z++) chan_noise[z]=1;
    if (!in->StatsDef()) in->setCubeStats();
    data_noise = in->stat().getSpread();
    
    wpow = par.WFUNC;
    
    // Read par.FREE and set free parameters
    setFree();
    
    // Creating mask if does not exist and write it in a fitsfile.
    if (!in->MaskAll() || in->pars().getMASK()=="NEGATIVE") in->BlankMask(chan_noise);
    mask = in->Mask();
    mask2D = new bool[in->DimX()*in->DimY()];
    for (int xy=0; xy<in->DimY()*in->DimX(); xy++) {
        mask2D[xy] = 0;
        for (int z=0; z<in->DimZ(); z++)
            mask2D[xy] += mask[xy+z*in->DimX()*in->DimY()];
    }

    // Choosing the normalization functions 
    stringstream ss(par.NORM);
    norms = readVec<string>(ss);
    // Normalizatios for second step + outputs
    if (norms.size()==1) norms.push_back(norms[0]);
    for (auto i=0; i<2; i++) {
        if (norms[i]!="NONE" && norms[i]!="AZIM" && norms[i]!="LOCAL") {
            std::cerr << "3DFIT WARNING: Normalization function not recognized for step " << i+1 << ". Using LOCAL.\n";
            norms[i] = "LOCAL";
        }
    }
    if (norms.size()==2) norms.push_back(norms[1]);

    // Setting limits for fitting parameters
    maxs[VROT]  = *max_element(inr->vrot.begin(),inr->vrot.end())+par.DELTAVROT;
    mins[VROT]  = *min_element(inr->vrot.begin(),inr->vrot.end())-par.DELTAVROT;
    maxs[VDISP] = par.MAXVDISP;
    mins[VDISP] = par.MINVDISP;
    maxs[Z0] = 1000;         // Max scaleheight allowed is 1000 arcs.  
    mins[Z0] = 0.;           // Min scaleheight allowed is 0 arcs.
    maxs[INC] = *max_element(inr->inc.begin(),inr->inc.end())+par.DELTAINC;
    mins[INC] = *min_element(inr->inc.begin(),inr->inc.end())-par.DELTAINC;
    maxs[PA]  = *max_element(inr->phi.begin(),inr->phi.end())+par.DELTAPHI;
    mins[PA]  = *min_element(inr->phi.begin(),inr->phi.end())-par.DELTAPHI; 
    if (mins[VROT]<0)  mins[VROT] = 0;
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
    maxs[VRAD]  = 500;
    mins[VRAD]  = -500;
    
    // Setting the convolution field
    if (par.SM) {
        if (!setCfield()) {
            std::cerr << "3DFIT WARNING: can not set an appropriate convolution "
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

    static int n=0;
    n = n==1 ? 2 : 1;
    std::string fileo = in->pars().getOutfolder()+"rings_final"+to_string(n)+".txt";
    remove(fileo.c_str());
    std::ofstream fout(fileo.c_str());

    if (verb) { 
        in->pars().setVerbosity(false);
        cout << showpoint << fixed << setprecision(2) << endl ;
        cout << setfill('=') << setw(38) << right << " 3DFIT " << setw(32) << " ";
        cout << setfill(' ') << endl;
    }

    *outr = *inr;
    writeHeader(fout,mpar,par.flagERRORS,par.flagBADOUT);

    // Normalizing input data cube such that the maximum value is =10.
    // This helps the convergence of the algorithm and
    // avoid problems with small flux values.
    double scaling = 1;
    if (par.NORMALCUBE) {
        scaling = 10./in->stat().getMax();
        for (auto i=in->NumPix(); i--;) in->Array(i) *= scaling;
        data_noise *= scaling;
    }
    
    T ***errors = allocate_3D<T>(inr->nr,2,nfree);
    bool fitok[inr->nr];
    for (int i=0; i<inr->nr; i++) fitok[i]=false;

    // Selecting the normalization function for current step
    if (norms[n-1]=="NONE") func_norm = &Model::Galfit<T>::norm_none;
    else if (norms[n-1]=="AZIM") func_norm = &Model::Galfit<T>::norm_azim;
    else func_norm = &Model::Galfit<T>::norm_local;

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
    

    bool usereverse = reverse;//(n==2 && reverse) || (!par.TWOSTAGE && reverse);
    if (usereverse) fit_reverse(errors,fitok,fout);
    else fit_straight(errors,fitok,fout);


  //  }

    fout.close();

    // If multi-threads or reverse rewrite ordered outfile
    if (in->pars().getThreads()>1 || usereverse) {
        double toKpc = KpcPerArc(distance);
        fout.open(fileo.c_str());
        writeHeader(fout,mpar,par.flagERRORS, par.flagBADOUT);
        for (int ir=0; ir<inr->nr; ir++) {
            if (fitok[ir] || par.flagBADOUT)
                writeRing(fout,outr,ir,toKpc,nfree,par.flagERRORS,errors,par.flagBADOUT,fitok[ir]);
        }
        fout.close();
    }
    
    // Setting density to 0 if fit is not converged or bad
    bool is2StepNeeded = mpar[INC] || mpar[PA] || mpar[Z0] || mpar[XPOS]|| mpar[YPOS] || mpar[VSYS];
    bool lastround = n==2 || (n==1 && !(par.TWOSTAGE && is2StepNeeded));
    if (lastround) {
        for (int ir=inr->nr; ir--;)
            if (!fitok[ir]) outr->dens[ir]=0;
    }

    deallocate_3D<T>(errors,inr->nr,2);

    if (verb) {
        cout << setfill('=') << setw(69) << "" << endl << endl;
        cout << fixed << setprecision(2) << setfill(' ');
        in->pars().setVerbosity(true);
    }
    
    // Scaling back to original values
    if (par.NORMALCUBE) {
        for (auto i=in->NumPix(); i--;) in->Array(i) /= scaling;
        data_noise /= scaling;
    }

}
template void Galfit<float>::galfit();
template void Galfit<double>::galfit();


template <class T>
void Galfit<T>::fit_straight(T ***errors, bool *fitok, std::ostream &fout) {

    double toKpc = KpcPerArc(distance);
    int start_rad = par.STARTRAD<inr->nr ? par.STARTRAD : 0;
    int nthreads = in->pars().getThreads();

#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    for (int ir=start_rad; ir<inr->nr; ir++) {

        T minimum=0;
        T pmin[nfree];

        if (verb && nthreads==1) {
            time_t t = time(NULL);
            char Time[11] = "          ";
            strftime (Time,11,"[%T]",localtime(&t));
            cout << fixed << setprecision(2)<<"\n Working on ring #"
                 << ir+1 << " at radius " << inr->radii[ir] << " arcsec ("
                 << inr->radii[ir]*toKpc << " Kpc)... " << Time << std::endl;
        }

        Rings<T> *dring = new Rings<T>;
        dring->id = ir;

        float width1=0, width2=0;
        // Handling the case of a single ring
        if (inr->nr==1) width1 = width2 = inr->radsep>0 ? inr->radsep/2. : inr->radii[0]/2.;
        else {
            if (ir==0) width1 = width2 = (inr->radii[1]-inr->radii[0])/2.;
            else if (ir==inr->nr-1) width1 = width2 = (inr->radii[ir]-inr->radii[ir-1])/2.;
            else {
                width1 = (inr->radii[ir]-inr->radii[ir-1])/2.;
                width2 = (inr->radii[ir+1]-inr->radii[ir])/2.;
            }
        }

        T drads[2] = {T(max(double(inr->radii[ir]-width1),0.)), T(max(double(inr->radii[ir]+width2),0.))};

        dring->addRings(2,drads,inr->xpos[ir],inr->ypos[ir],inr->vsys[ir],inr->vrot[ir],inr->vdisp[ir],inr->vrad[ir],
                        inr->vvert[ir],inr->dvdz[ir],inr->zcyl[ir],inr->dens[ir],inr->z0[ir],inr->inc[ir],inr->phi[ir]);

        // Checking that we have enough good pixels to proceed.
        int blo[2], bhi[2];
        getModelSize(dring,blo,bhi);
        double theta=0;
        int nTot=0, nIn=0;
        for (int y=blo[1]; y<=bhi[1]; y++) {
            for (int x=blo[0]; x<=bhi[0]; x++) {
                if (IsIn(x-blo[0],y-blo[1],blo,dring,theta)) {
                    nTot++;
                    if (mask2D[x+y*in->DimX()]>0) nIn++;
                }
            }
        }

        if (nIn<20 && (nIn==0 || nIn<0.7*nTot)) {
            if (verb) {
                std::string msg = " Not enough pixels to fit in ring #"+to_string(ir+1)+". I will skip it.";
                WarningMessage(cout,msg);
            }
            fitok[ir]=false;
            continue;
        }

        // Fitting
        fitok[ir] = minimize(dring, minimum, pmin, nullptr);
        if (!fitok[ir]) {
            if (verb) {
                std::string msg = " Can not achieve convergence in ring #"+to_string(ir+1)+". I'll keep going, but \n"+
                                  " parameters for this ring are wrong! Please, try to change initial \n" +
                                  " conditions and/or the function to minimize.";
                WarningMessage(cout,msg);
            }
            if (par.flagBADOUT) writeRing(fout,outr,ir,toKpc,nfree,par.flagERRORS,errors,par.flagBADOUT,fitok[ir]);
            continue;
        }

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

        // Check that VROT is within limit of cube
        double maxv = fabs(AlltoVel(in->getZphys(in->DimZ()-1),in->Head())-outr->vsys[ir]);
        double minv = fabs(AlltoVel(in->getZphys(0),in->Head())-outr->vsys[ir]);
        double maxvrot = std::max(maxv/sin(outr->inc[ir]*M_PI/180.),minv/sin(outr->inc[ir]*M_PI/180.));
        if (outr->vrot[ir]>maxvrot) {
            if (verb) {
                std::string msg = " Ring #"+to_string(ir+1)+" does not look good. I will ignore it and treat it as a   \n"+
                                  "non-converged fit.";
                WarningMessage(cout,msg);
            }
            fitok[ir] = false;
            if (par.flagBADOUT) writeRing(fout,outr,ir,toKpc,nfree,par.flagERRORS,errors,par.flagBADOUT,fitok[ir]);
            continue;
        }

        if (verb) printRing(cout,outr,ir,minimum,toKpc,mpar,nthreads);

        if (par.flagERRORS) getErrors(dring,errors[ir],ir,minimum);

        writeRing(fout,outr,ir,toKpc,nfree,par.flagERRORS,errors,par.flagBADOUT,fitok[ir]);

        delete dring;
    }
}
template void Galfit<float>::fit_straight(float ***errors, bool *fitok, std::ostream &fout);
template void Galfit<double>::fit_straight(double ***errors, bool *fitok, std::ostream &fout);


template <class T>
void Galfit<T>::fit_reverse(T ***errors, bool *fitok, std::ostream &fout) {

    double toKpc = KpcPerArc(distance);
    int start_rad = par.STARTRAD<inr->nr ? par.STARTRAD : 0;

    for (int ir=inr->nr-1; ir>=start_rad; ir--) {

        T minimum=0;
        T pmin[nfree];

        if (verb) {
            time_t t = time(NULL);
            char Time[11] = "          ";
            strftime (Time,11,"[%T]",localtime(&t));
            cout << fixed << setprecision(2)<<"\n Working on ring #"
                 << ir+1 << " at radius " << inr->radii[ir] << " arcsec ("
                 << inr->radii[ir]*toKpc << " Kpc)... " << Time << std::endl;
        }

        Rings<T> *dring = new Rings<T>;
        dring->id = ir;

        float width1=0, width2=0;
        // Handling the case of a single ring
        if (inr->nr==1) width1 = width2 = inr->radsep>0 ? inr->radsep/2. : inr->radii[0]/2.;
        else {
            if (ir==0) width1 = width2 = (inr->radii[1]-inr->radii[0])/2.;
            else if (ir==inr->nr-1) width1 = width2 = (inr->radii[ir]-inr->radii[ir-1])/2.;
            else {
                width1 = (inr->radii[ir]-inr->radii[ir-1])/2.;
                width2 = (inr->radii[ir+1]-inr->radii[ir])/2.;
            }
        }

        T drads[2] = {T(max(double(inr->radii[ir]-width1),0.)), T(max(double(inr->radii[ir]+width2),0.))};

        dring->addRings(2,drads,inr->xpos[ir],inr->ypos[ir],inr->vsys[ir],inr->vrot[ir],inr->vdisp[ir],inr->vrad[ir],
                        inr->vvert[ir],inr->dvdz[ir],inr->zcyl[ir],inr->dens[ir],inr->z0[ir],inr->inc[ir],inr->phi[ir]);

        if (ir!=inr->nr-1) {
            // Creating a new set of rings fitted so far
            Rings<T> *dring2 = new Rings<T>;
            dring2->id = dring->id+1;
            // Adding previously fitted rings
            int k = ir+1;

            dring2->addRings(inr->nr-ir-1,&outr->radii[k],&outr->xpos[k],&outr->ypos[k],&outr->vsys[k],&outr->vrot[k],&outr->vdisp[k],&outr->vrad[k],
                             &outr->vvert[k],&outr->dvdz[k],&outr->zcyl[k],&outr->dens[k],&outr->z0[k],&outr->inc[k],&outr->phi[k]);

            // Calculating the model so far
            int blo[2], bhi[2], bsize[2];
            getModelSize(outr,blo,bhi);
            int nv = par.NV<0 ? in->DimZ() : par.NV;
            Model::Galmod<T> *modsoFar = new Model::Galmod<T>;
            modsoFar->input(in,bhi,blo,dring2,nv,par.LTYPE,1,par.CDENS);
            modsoFar->calculate();
            //modsoFar->Out()->fitswrite_3d((to_string(ir)+".fits").c_str());
            fitok[ir] = minimize(dring, minimum, pmin, modsoFar);

            delete dring2;
            delete modsoFar;
        }
        else
            fitok[ir] = minimize(dring, minimum, pmin, nullptr);


        if (!fitok[ir]) {
            if (verb) {
                std::string msg = " Can not achieve convergence in ring #"+to_string(ir+1)+". I'll keep going, but \n"+
                                  " parameters for this ring are wrong! Please, try to change initial \n" +
                                  " conditions and/or the function to minimize.";
                WarningMessage(cout,msg);
            }
            if (par.flagBADOUT) writeRing(fout,outr,ir,toKpc,nfree,par.flagERRORS,errors,par.flagBADOUT,fitok[ir]);
            continue;
        }

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

        // Check that VROT is within limit of cube
        double maxv = fabs(AlltoVel(in->getZphys(in->DimZ()-1),in->Head())-outr->vsys[ir]);
        double minv = fabs(AlltoVel(in->getZphys(0),in->Head())-outr->vsys[ir]);
        double maxvrot = std::max(maxv/sin(outr->inc[ir]*M_PI/180.),minv/sin(outr->inc[ir]*M_PI/180.));
        if (outr->vrot[ir]>maxvrot) {
            if (verb) {
                std::string msg = " Ring #"+to_string(ir+1)+" does not look good. I will ignore it and treat it as a   \n"+
                                  "non-converged fit.";
                WarningMessage(cout,msg);
            }
            fitok[ir] = false;
            if (par.flagBADOUT) writeRing(fout,outr,ir,toKpc,nfree,par.flagERRORS,errors,par.flagBADOUT,fitok[ir]);
            continue;
        }

        if (verb) printRing(cout,outr,ir,minimum,toKpc,mpar,1);

        if (par.flagERRORS) getErrors(dring,errors[ir],ir,minimum);

        writeRing(fout,outr,ir,toKpc,nfree,par.flagERRORS,errors,par.flagBADOUT,fitok[ir]);

        delete dring;
    }
}
template void Galfit<float>::fit_reverse(float ***errors, bool *fitok, std::ostream &fout);
template void Galfit<double>::fit_reverse(double ***errors, bool *fitok, std::ostream &fout);


template <class T> 
bool Galfit<T>::SecondStage() {
    
    bool isNeeded = mpar[INC] || mpar[PA]   || mpar[Z0] ||
                    mpar[XPOS]|| mpar[YPOS] || mpar[VSYS];
    if (!isNeeded) {par.TWOSTAGE=false; return isNeeded;}
    second = true;
    // Following line is questionable, although I see a slight improvement in using local
    // in the second step in many cases.
    //if (!in->pars().getFlagSlitfit()) func_norm = &Model::Galfit<T>::norm_local;
    
    int nr = outr->nr;
    
    // Deciding how to regularize the parameters
    std::vector<int> reg = {-3,-3,-3,-3,-3,-3}; // INC PA VSYS XCEN YCEN Z0
    std::vector<std::string> regstr = {"inc","pa","vsys","xpos","ypos","z0"};
    stringstream ss(makelower(par.REGTYPE));
    std::vector<std::string> Polyn = readVec<string>(ss);
    auto getval = [](std::string s){
        int nval = -3;                                     // Auto
        if (isNumber(s))      nval = 1+atoi(s.c_str());    // Polynomial fitting
        else if (s=="bezier") nval = -1;                   // Bezier
        else if (s=="median") nval = -2;                   // Median
        return nval;
    };
    auto keys = splitStrings(Polyn,"=");
    if (Polyn.size()==1 && keys.first.size()==0) {          // If one value, apply to pa and inc
        reg[0] = reg[1] = getval(Polyn[0]);
    }
    else {
        for (int i=0; i<keys.first.size(); i++) {
            for (int j=0; j<reg.size(); j++)
                if (keys.first[i]==regstr[j]) reg[j] = getval(keys.second[i]);
        }
    }

    for (int i=0; i<reg.size(); i++) {
        if (reg[i]==-3) {
            // If 'auto' mode is selected, choosing an appropriate regularization
            if (i<2) {                          // For INC and PA:
                if (nr<=4) reg[i] = 1;          // Constant
                else {
                    // Calculating scatter of inc/pa
                    std::vector<T> myvec = i==0 ? outr->inc : outr->phi;
                    T median = findMedian<T>(&myvec[0],myvec.size());
                    T madfm  = findMADFM(&myvec[0],myvec.size(),median,false);
                    // If scatter is > 3 degrees, it is worth tracing the change
                    if (madfmToSigma(madfm)>3) {
                        if (nr<10) reg[i] = 2;     // Line
                        else reg[i] = -1;          // Bezier
                    }
                    else reg[i] = -2;              // Median
                }
            }
            else reg[i] = -2;                      // For any other: take the Median
        }
    }


    *inr = *outr;
    bool proceed = true;
    if (mpar[INC])   proceed *= regularizeParams(inr->radii,inr->inc,inr->inc,reg[0]);
    if (mpar[PA])    proceed *= regularizeParams(inr->radii,inr->phi,inr->phi,reg[1]);
    if (mpar[VSYS])  proceed *= regularizeParams(inr->radii,inr->vsys,inr->vsys,reg[2]);
    if (mpar[XPOS])  proceed *= regularizeParams(inr->radii,inr->xpos,inr->xpos,reg[3]);
    if (mpar[YPOS])  proceed *= regularizeParams(inr->radii,inr->ypos,inr->ypos,reg[4]);
    if (mpar[Z0])    proceed *= regularizeParams(inr->radii,inr->z0,inr->z0,reg[5]);

    bool oldmpar[MAXPAR];
    for (int i=MAXPAR; i--;) oldmpar[i]=mpar[i];
    int oldnfree =nfree;
    
    mpar[INC] = mpar[PA] = mpar[VSYS] = false;
    mpar[XPOS]= mpar[YPOS] = mpar[Z0] = false;
    nfree = mpar[VROT]+mpar[VDISP]+mpar[VRAD];

    if (proceed) galfit();
    else {
        if (verb) std::cerr << "3DFIT ERROR: Regularization failed."
                               " I will not proceed to second fitting stage. \n";
    }
    
    for (int i=MAXPAR; i--;) mpar[i]=oldmpar[i];
    nfree = oldnfree;

    return isNeeded;
}
template bool Galfit<float>::SecondStage();
template bool Galfit<double>::SecondStage();


template <class T>
bool Galfit<T>::regularizeParams(std::vector<T> x, std::vector<T> y, std::vector<T> &yout, int rtype) {

    int n = x.size();

    // Determinig what parameter we are regularizing
    std::string whichpar;
    if      (&yout==&inr->inc)  whichpar="inclination";
    else if (&yout==&inr->phi)  whichpar="position angle";
    else if (&yout==&inr->vsys) whichpar="systemic velocity";
    else if (&yout==&inr->xpos) whichpar="X center";
    else if (&yout==&inr->ypos) whichpar="Y center";
    else if (&yout==&inr->z0)   whichpar="disk thickness";

    if (rtype==-2) {                            // Take the median
        T val = findMedian<T>(&y[0],n);
        for (int i=n; i--;) yout[i] = val;
    }
    else if (rtype==-1) {                        // Bezier interpolation
        std::vector<T> x_bez(n), y_bez(n);
        if (!bezier_interp(x,y,x,yout)) {
            if (verb) std::cerr << "3DFIT ERROR: cannot find a bezier interpolation for "<< whichpar << ".\n";
            return false;
        }
    }
    else  {                                     // Polynomial interpolation

        if (n<=rtype) {
            if (verb) std::cerr << "3DFIT WARNING - Second stage: too few degree of freedom. Using a constant value.\n";
            rtype = 1;
        }

        std::vector<T> xx = x, yy = y;
        if (n-2>rtype) {
            // Deleting maximum and minimum values for the fit
            auto where = std::max_element(yy.begin(), yy.end());
            yy.erase(where);
            xx.erase(xx.begin()+std::distance(yy.begin(), where));
            where = std::min_element(yy.begin(), yy.end());
            yy.erase(where);
            xx.erase(xx.begin()+std::distance(yy.begin(), where));
        }

        bool mp[rtype];
        for (int i=0; i<rtype;i++) mp[i] = true;
        std::vector<T> ww(n,1);
        T coeff[rtype], coefferr[rtype];

        Lsqfit<T> lsq(&xx[0],1,&yy[0],&ww[0],xx.size(),coeff,coefferr,mp,rtype,&polyn,&polynd);
        if (lsq.fit()<0) {
            if (verb) std::cerr << "3DFIT ERROR: cannot least-square fit " << whichpar << ".\n";
            return false;
        }

        // Filling output vector with polynomial values.
        for (int i=0; i<n; i++) {
            yout[i]=0;
            for (int j=0; j<rtype; j++)
                yout[i] += coeff[j]*std::pow(double(x[i]),j);
        }

        if (verb) {
            // Printing best-fit coefficients
            std::cout << "  Best parameters for " << whichpar << ":\n";
            for (int i=0; i<rtype; i++) {
                std::cout << setprecision(4);
                string a = "a"+to_string<int>(i)+" = ";
                std::cout << setw(43) << right << a << setw(9) << coeff[i]
                          << setw(4) << " ±" << setw(9) << coefferr[i] << std::endl;
            }
        }
    }

    return true;
}
template bool Galfit<float>::regularizeParams(std::vector<float>,std::vector<float>,std::vector<float>&,int);
template bool Galfit<double>::regularizeParams(std::vector<double>,std::vector<double>,std::vector<double>&,int);


template <class T> 
bool Galfit<T>::setCfield() {
    
    T pixsizeX = fabs(in->Head().Cdelt(0))*arcconv; 
    T pixsizeY = fabs(in->Head().Cdelt(1))*arcconv; 
    
    //Beam Old = {pixsizeX, pixsizeY, 0};
    Beam Old = {0, 0, 0};

    Beam New = {in->Head().Bmaj()*3600.,        // Beam always in degrees
                in->Head().Bmin()*3600.,
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

/*
template <class T>
Model::Galmod<T>* Galfit<T>::getModel(Rings<T> *dr) {

    // If no rings are provided, we will use the output rings and built the final model
    if (dr==nullptr) {
        dr = outr;
        // Creating output rings for final Galmod. Moving innermost and outermost ring boundaries.    
        if (dr->nr==1) 
            dr->addRing(dr->radii[0]+dr->radsep/2., dr->xpos[0], dr->ypos[0], dr->vsys[0], dr->vrot[0], dr->vdisp[0], 
                        dr->vrad[0], dr->vvert[0], dr->dvdz[0], dr->zcyl[0], dr->dens[0], dr->z0[0], dr->inc[0], dr->phi[0]);
        else dr->radii[dr->nr-1] += dr->radsep/2.;
        dr->radii[0] = max(double(dr->radii[0]-dr->radsep/2.),0.);
    }

    int bhi[2] = {in->DimX(), in->DimY()};
    int blo[2] = {0,0};
    int nv = par.NV==-1 ? in->DimZ() : par.NV;
    
    Model::Galmod<T> *mod = new Model::Galmod<T>;
    mod->input(in,bhi,blo,dr,nv,par.LTYPE,1,par.CDENS);
    mod->calculate();
    if (par.SM) mod->smooth();
    return mod;
}
template Model::Galmod<float>* Galfit<float>::getModel(Rings<float> *);
template Model::Galmod<double>* Galfit<double>::getModel(Rings<double> *);
*/


template <class T>
Model::Galmod<T>* Galfit<T>::getModel(Rings<T> *dr, int* bhi, int* blo, Model::Galmod<T> *modsoFar, bool finalModel) {

    if (finalModel) {
        // Creating output rings for final Galmod. Moving innermost and outermost ring boundaries.    
        if (dr->nr==1) 
            dr->addRing(dr->radii[0]+dr->radsep/2., dr->xpos[0], dr->ypos[0], dr->vsys[0], dr->vrot[0], dr->vdisp[0], 
                        dr->vrad[0], dr->vvert[0], dr->dvdz[0], dr->zcyl[0], dr->dens[0], dr->z0[0], dr->inc[0], dr->phi[0]);
        else dr->radii[dr->nr-1] += dr->radsep/2.;
        dr->radii[0] = max(double(dr->radii[0]-dr->radsep/2.),0.);
    }
    
    int nv = par.NV==-1 ? in->DimZ() : par.NV;
    int bsize[2] = {bhi[0]-blo[0], bhi[1]-blo[1]};

    Model::Galmod<T> *mod = new Model::Galmod<T>;
    mod->input(in,bhi,blo,dr,nv,par.LTYPE,1,par.CDENS);
    mod->calculate();

    // Adding up the "sofar" model, if requested
    T *modp = mod->Out()->Array();
    if (modsoFar!=nullptr && !finalModel) {
        for (auto i=mod->Out()->NumPix(); i--;) modp[i] += modsoFar->Out()->Array()[i];
    }

    //<<<<< Convolution....
    if (par.SM) {
        if (in->pars().getflagFFT()) Convolve_fft(modp, bsize);
        else Convolve(modp, bsize);
    }

    return mod;

}
template Model::Galmod<float>* Galfit<float>::getModel(Rings<float>*, int*, int*, Model::Galmod<float>*, bool);
template Model::Galmod<double>* Galfit<double>::getModel(Rings<double>*, int*, int*, Model::Galmod<double>*, bool);


template <class T>
void Galfit<T>::setFree() {

    std::string f = makelower(par.FREE);
    
    // Set everything to fixed
    for (int i=0; i<MAXPAR; i++) mpar[i] = false;
    
    // Detect requested free parameters
    if (f.find("vrot")!=std::string::npos) mpar[VROT] = true;
    if (f.find("disp")!=std::string::npos) mpar[VDISP]= true;
    if (f.find("z0")!=std::string::npos)   mpar[Z0]   = true;
    if (f.find("inc")!=std::string::npos)  mpar[INC]  = true;
    if (f.find("pa")!=std::string::npos)   mpar[PA]   = true;
    if (f.find("phi")!=std::string::npos)  mpar[PA]   = true;
    if (f.find("xpos")!=std::string::npos) mpar[XPOS] = true;
    if (f.find("ypos")!=std::string::npos) mpar[YPOS] = true;
    if (f.find("vsys")!=std::string::npos) mpar[VSYS] = true;
    if (f.find("vrad")!=std::string::npos) mpar[VRAD] = true;
    //if (f.find("dens")!=std::string::npos) mpar[DENS] = true;

    if (f.find("all")!=std::string::npos)
        for (int i=0; i<MAXPAR; i++) mpar[i] = true;
    
    int nfixed=0;
    for (int i=0; i<MAXPAR; i++) nfixed += (1-mpar[i]);
    if (nfixed == MAXPAR) {
        std::cerr << " 3DFIT ERROR: NO free parameters!\n";
        std::terminate();
    }
    nfree = MAXPAR-nfixed;
    
}
template void Galfit<float>::setFree();
template void Galfit<double>::setFree();


template <class T>
bool Galfit<T>::AsymmetricDrift(T *rad, T *densprof, T *dispprof, T *rotcur, T *inc, int nn) {
    
    // Compute an asymmetric drift correction, following procedure in Iorio+17, sec 4.3
//    in->pars().setVerbosity(true);
//    for (auto i=0; i<nn; i++) {
 //       std::cout << rad[i] << " " << densprof[i] << " " << dispprof[i] << " " << rotcur[i] << " " << inc[i] << std::endl;
  //  }
    // Fitting dispersion with a polynomial function
    int npar1 = in->pars().getParGF().ADRIFTPOL1 + 1;
    T cdisp[npar1], cdisperr[npar1], ww[nn];
    bool mp[npar1]; 
    if (npar1>0) { // If ADRIFTPOL1<0, just use the density profile
        if (npar1>=nn) npar1 = nn-1;
        for (int i=0; i<npar1; i++) mp[i] = true;
        for (int i=0; i<nn; i++) ww[i] = 1;
        Lsqfit<T> lsq1(rad,1,dispprof,ww,nn,cdisp,cdisperr,mp,npar1,&polyn,&polynd);
        if (lsq1.fit()<0) {
            if (in->pars().isVerbose()) std::cerr << "3DFIT ERROR: cannot least-square fit the dispersion for asymmetric drift.\n";
            par.flagADRIFT = false;
            return false;
        }
    }
    
    // Now fitting log(density*disp2) with a polynomial function 
    T *fun = new T[nn];
    int npar2 = in->pars().getParGF().ADRIFTPOL2 + 1;
    if (npar2>=nn) npar2 = nn-1;
    T cfun[npar2], cfunerr[npar2];
    bool mpp[npar2];
    for (int i=0; i<npar2; i++) mpp[i] = true;
    for (int i=0; i<nn; i++) {
        fun[i] = log(dispprof[i]*dispprof[i]*densprof[i]*cos(inc[i]*M_PI/180.));
        ww[i] = 1;
    }
    //Lsqfit<T> lsq2(rad,1,fun,ww,nn,cfun,cfunerr,mpp,npar2,&coreExp,&coreExpd);
    Lsqfit<T> lsq2(rad,1,fun,ww,nn,cfun,cfunerr,mpp,npar2,&polyn,&polynd);
    if (lsq2.fit()<0) {
        if (in->pars().isVerbose()) std::cerr << "3DFIT ERROR: cannot least-square fit the fun for asymmetric drift.\n";
        par.flagADRIFT = false;
        return false;
    }
    
    // Now writing to a text file
    std::ofstream fout(in->pars().getOutfolder()+"asymdrift.txt");
    
    int m=20;
    fout << fixed << setprecision(4);
    fout << "#" << setw(m-1) << "RAD(arcs)"
         << setw(m) << "VCIRC(km/s)"
         << setw(m) << "ASYMDRIFT_SQ(km/s)2"
         << setw(m) << "DISP_REG(km/s)"
         << setw(m) << "FUN" 
         << setw(m) << "FUN_REG" << std::endl;
        
    for (int i=0; i<nn; i++) {
        // Getting regularized functions
        T disp_reg = npar1>0 ? polyn(&rad[i],cdisp,npar1) : dispprof[i];
        T fun_reg = polyn(&rad[i],cfun,npar2);
        // Derivative of regularized fun
        T fun_der = dpolyn_dx(&rad[i],cfun,npar2);
        // Asymmetric drift correction
        T asdrift_squared = -rad[i]*disp_reg*disp_reg*fun_der;
        T vcirc = sqrt(rotcur[i]*rotcur[i]+asdrift_squared);
        
        fout << setw(m) << rad[i] << setw(m) << vcirc << setw(m) << asdrift_squared 
             << setw(m) << disp_reg << setw(m) << fun[i] << setw(m) << fun_reg << std::endl;
    }
    
    fout.close();
    
    delete [] fun;
    
    return true;
    
}
template bool Galfit<float>::AsymmetricDrift(float*,float*,float*,float*,float*,int);
template bool Galfit<double>::AsymmetricDrift(double*,double*,double*,double*,double*,int);


template <class T>
void Galfit<T>::writeRingFile(std::string filename, Rings<T> *r, T ***errors) {

    std::string fileo = in->pars().getOutfolder()+filename;
    remove(fileo.c_str());
    std::ofstream fout(fileo.c_str());
    double toKpc = KpcPerArc(distance);

    writeHeader(fout,mpar,par.flagERRORS,par.flagBADOUT);
    
    for (int i=0; i<r->nr; i++) {
        writeRing(fout,r,i,toKpc,nfree,par.flagERRORS,errors,false,false);
    }
    
    fout.close();

}
template void Galfit<float>::writeRingFile(std::string, Rings<float>*, float***);
template void Galfit<double>::writeRingFile(std::string, Rings<double>*, double***);


/////////////////////////////////////////////////////////////////////
/// Functions to write GALFIT rings
////////////////////////////////////////////////////////////////////
void writeHeader(std::ostream &fout, bool *mpar, bool writeErrors, bool writeBadRings) {

    int m = 10;
    fout << left << setw(m) << "#RAD(Kpc)"
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

    if (writeErrors) {
        if (mpar[VROT])  fout << setw(m) << "E_VROT1" << setw(m) << "E_VROT2";
        if (mpar[VDISP]) fout << setw(m) << "E_DISP1" << setw(m) << "E_DISP2";
        if (mpar[DENS])  fout << setw(m) << "E_DENS1" << setw(m) << "E_DENS2";
        if (mpar[Z0])    fout << setw(m) << "E_Z01"   << setw(m) << "E_Z02";
        if (mpar[INC])   fout << setw(m) << "E_INC1"  << setw(m) << "E_INC2";
        if (mpar[PA])    fout << setw(m) << "E_PA1"   << setw(m) << "E_PA2";
        if (mpar[XPOS])  fout << setw(m) << "E_XPOS1" << setw(m) << "E_XPOS2";
        if (mpar[YPOS])  fout << setw(m) << "E_YPOS1" << setw(m) << "E_YPOS2";
        if (mpar[VSYS])  fout << setw(m) << "E_VSYS1" << setw(m) << "E_VSYS2";
        if (mpar[VRAD])  fout << setw(m) << "E_VRAD1" << setw(m) << "E_VRAD2";
    }

    if (writeBadRings)
      fout << setw(m) << "FITOK";

    fout << endl;
}


template <class T>
void writeRing(std::ostream &fout, Rings<T> *r, int i, double toKpc, int nfree, bool writeErrors, T ***errors, bool writeBadRings, bool fitOK) {

    int m=10;
#pragma omp critical (galfit_write)
{
    // Writing output file. Not ordered if multithread
    fout << setprecision(3) << fixed << left;

    fout << setw(m) << r->radii[i]*toKpc
         << setw(m+1) << r->radii[i]
         << setw(m+1) << r->vrot[i]
         << setw(m+1) << r->vdisp[i]
         << setw(m) << r->inc[i]
         << setw(m) << r->phi[i]
         << setw(m) << r->z0[i]*toKpc*1000
         << setw(m) << r->z0[i]
         << setw(m) << r->dens[i]/1E20
         << setw(m) << r->xpos[i]
         << setw(m) << r->ypos[i]
         << setw(m+1) << r->vsys[i]
         << setw(m+1) << r->vrad[i];

    if (writeErrors)
        for (int kk=0; kk<nfree; kk++)
            fout << setw(m) << errors[i][0][kk] << setw(m) << errors[i][1][kk];

    if (writeBadRings)
      fout << setw(m) << fitOK;

    fout << endl;
}
}


template <class T>
void printRing(std::ostream &fout, Rings<T> *r, int i, double minimum, double toKpc, bool *mpar, int nthreads){
#pragma omp critical (galfit_outmsg)
{
        if (nthreads>1) {
            fout << "\n Ring #" << i+1 << " at radius " << r->radii[i] << " arcsec ("
                 << r->radii[i]*toKpc << " Kpc)... \n";
        }

        int k=8, n=11;
        fout << "  Best parameters for ring #" << i+1
             << " (fmin = " << scientific << setprecision(3) << minimum << "):\n";

        fout << fixed << setprecision(2);

        string s;
        s = "    Vrot";
        if (!mpar[VROT]) s += "(f)";
        fout << setw(n) << left << s << setw(3) << right << "= "
             << setw(k) << r->vrot[i] << left << setw(k) << "  km/s";

        s = "        Disp";
        if (!mpar[VDISP]) s += "(f)";
        fout << setw(n+4) << left << s << setw(3) << right << "= "
             << setw(k-1) << r->vdisp[i]
             << left << setw(k) << "  km/s" << endl;

        s = "    Vrad";
        if (!mpar[VRAD]) s += "(f)";
        fout << setw(n) << left << s << setw(3) << right << "= "
             << setw(k) << r->vrad[i] << left << setw(k) << "  km/s";

        s = "        Vsys";
        if (!mpar[VSYS]) s += "(f)";
        fout << setw(n+4) << left << s << setw(3) << right << "= "
             << setw(k-1) << r->vsys[i] << left << setw(k) << "  km/s" << endl;

        s = "    Inc";
        if (!mpar[INC]) s += "(f)";
        fout << setw(n) << left << s << setw(3) << right << "= "
             << setw(k) << r->inc[i] << left << setw(k) << "  deg";

        s = "        PA";
        if (!mpar[PA]) s += "(f)";
        fout << setw(n+4) << left << s << setw(3) << right << "= "
             << setw(k-1) << r->phi[i] << left << setw(k) << "  deg" << endl;

        s = "    Xpos";
        if (!mpar[XPOS]) s += "(f)";
        fout << setw(n) << left << s << setw(3) << right << "= "
             << setw(k) << r->xpos[i] << left << setw(k) << "  pix";

        s = "        Ypos";
        if (!mpar[YPOS]) s += "(f)";
        fout << setw(n+4) << left << s << setw(3) << right << "= "
             << setw(k-1) << r->ypos[i] << left << setw(k) << "  pix" << endl;

        s = "    Z0";
        if (!mpar[Z0]) s += "(f)";
        fout << setw(n) << left << s << setw(3) << right << "= "
             << setw(k) << r->z0[i]*toKpc << left << setw(k) << "  Kpc";

        //s = "        CD";----
        //if (!mpar[DENS]) s += "(f)";
        //fout << setw(n+4) << left << s << setw(3) << right << "= "
        //     << setw(k-1) << scientific << setprecision(1)
        //     << outr->dens[ir] << left << setw(k) << "  a/cm2";

        fout << endl;

}
}


void WarningMessage(std::ostream &fout, std::string msg){
#pragma omp critical
{
    fout <<"\n ========================== 3DFIT WARNING ==========================\n"
         << msg << std::endl
         << " ===================================================================\n\n";
}
}

}

