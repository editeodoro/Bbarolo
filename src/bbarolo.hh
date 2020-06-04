/////////////////////////////////////////////
// bbarolo.hh: Contains BB core function
/////////////////////////////////////////////

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


#ifndef BBAROLO_HH_
#define BBAROLO_HH_

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>
#include <Arrays/param.hh>
#include <Arrays/cube.hh>
#include <Arrays/stats.hh>
#include <Tasks/moment.hh>
#include <Tasks/ringmodel.hh>
#include <Tasks/smooth3D.hh>
#include <Tasks/galfit.hh>
#include <Tasks/galwind.hh>
#include <Tasks/ellprof.hh>
#include <Tasks/spacepar.hh>
#include <Utilities/utils.hh>

#ifdef DOUBLE_PRECISION
using BBreal = double;
#else
using BBreal = float;
#endif

bool BBcore (Param *par) {

    // Interface for all BBarolo's tasks

    Cube<BBreal> *c = new Cube<BBreal>;
    c->saveParam(*par);
    
    if (!c->readCube(par->getImageFile())) {
        std::cout << par->getImageFile() << " is not a readable FITS file!\n";
        delete c;
        return false;
    }

    std::string outfolder = c->pars().getOutfolder();
    if (outfolder=="") {
       std::string path = get_currentpath();
            outfolder = path+"/output/"+c->Head().Obname()+"/";
            c->pars().setOutfolder(outfolder);
    }
    mkdirp(outfolder.c_str());
    
    if (par->getCheckCh()) c->CheckChannels();
    
    // Mask making utility ----------------------------------------
    if (par->getMakeMask()) {
        if (par->getMASK().find("LARGEST")!=std::string::npos) c->BlankMask(NULL,true);
        else c->BlankMask(NULL,false);
    }
    // --------------------------------------------------------------
    
    
    // Spatial smoothing utility ------------------------------------------
    if (par->getflagSmooth()) {
        Smooth3D<BBreal> *sm = new Smooth3D<BBreal>;
        sm->cubesmooth(c);
        sm->fitswrite();
        delete sm;
    }
    // --------------------------------------------------------------


    // Spectral smoothing utility -----------------------------------
    if (par->getflagSmoothSpectral()) {
        SpectralSmooth3D<BBreal> *sm = new SpectralSmooth3D<BBreal>(par->getWindowType(),par->getWindowSize());
        sm->smooth(c);
        sm->fitswrite(c);
        delete sm;
    }
    // --------------------------------------------------------------


    // Source finding utility --------------------------------------
    if (par->getflagSearch()) {
        c->search();
        std::ofstream detout((outfolder+"detections.txt").c_str());
        c->printDetections(detout);
        c->writeDetections();
        if (par->getParSE().cubelets) {
            c->writeCubelets();
            c->plotDetections();
        }
    }
    // --------------------------------------------------------------
        
        
    // 3D Cube Fitting task -----------------------------------------
    if (par->getflagGalFit()) {
        Model::Galfit<BBreal> *fit = new Model::Galfit<BBreal>(c);
        fit->galfit();
        if (par->getParGF().TWOSTAGE) fit->SecondStage();
        if (par->getFlagDebug()) fit->writeModel("BOTH",par->getFlagPlots());
        else fit->writeModel(par->getParGF().NORM,par->getFlagPlots());
        delete fit;
    }
    // --------------------------------------------------------------
    
    
    // Cube Model task -----------------------------------------------
    if (par->getflagGalMod()) {
        Model::Galfit<BBreal> *fit = new Model::Galfit<BBreal>(c);
        if (par->getFlagDebug()) fit->writeModel("BOTH",false);
        else fit->writeModel(par->getParGF().NORM,false);
        delete fit;
    }
    //----------------------------------------------------------------

    
    // Full parameter space task --------------------------------------
    if (par->getflagSpace()) {
        Spacepar<BBreal> *sp = new Spacepar<BBreal>(c);
        sp->calculate();
        sp->plotAll_Python();
        delete sp;
    }
    //----------------------------------------------------------------


    // GalWind task --------------------------------------------------
    if (par->getParGW().flagGALWIND) {
        GalWind<BBreal> *w = new GalWind<BBreal>(c);   
        w->compute();
        if (par->getParGF().SM) w->smooth();
        w->writeFITS();
        w->writeMomentMaps();
        w->writePV();
        if (par->getFlagPlots()) w->makePlots();
        delete w;
    }
    //----------------------------------------------------------------


    // 2D tilted-ring fitting task -----------------------------------
    if (par->getFlagRing()) {
        Ringmodel<BBreal> *trmod = new Ringmodel<BBreal>(c);
        trmod->ringfit(c->pars().getThreads(),c->pars().isVerbose(),c->pars().getShowbar());
        std::string fout = c->pars().getOutfolder()+c->Head().Name()+"_2dtrm.txt";
        std::ofstream fileo(fout.c_str());
        trmod->printfinal(fileo,c->Head());
        fileo.close();
        trmod->printfinal(std::cout,c->Head());
        trmod->writeModel(c->pars().getOutfolder()+c->Head().Name()+"_2d_mod.fits",c->Head());
        delete trmod;
    }
    //-----------------------------------------------------------------


    // Moment maps task -----------------------------------------------
    if (par->getMaps()) {
        MomentMap<BBreal> map;
        map.input(c);
        std::string s = outfolder+c->Head().Name();
        bool masking = par->getMASK()=="NONE" ? false : true;
        if (par->getMassDensMap()) {
            map.HIMassDensityMap(masking);
            map.fitswrite_2d((s+"map_massdens.fits").c_str());
        }
        if (par->getTotalMap()) {
            map.ZeroMoment(masking,par->getMapType());
            map.fitswrite_2d((s+"map_0th.fits").c_str());
        }
        if (par->getVelMap()) {
            map.FirstMoment(masking,par->getMapType());
            map.fitswrite_2d((s+"map_1st.fits").c_str());
        }
        if (par->getDispMap()) {
            map.SecondMoment(masking,par->getMapType());
            map.fitswrite_2d((s+"map_2nd.fits").c_str());
        }
        if (par->getRMSMap()) {
            map.RMSMap();
            map.fitswrite_2d((s+"map_RMS.fits").c_str());
        }
        if (par->getGlobProf()) {
            Image2D<BBreal> spectrum;
            spectrum.extractGlobalSpectrum(c);
            std::string specfname = c->pars().getOutfolder()+"spectrum.txt";
            std::ofstream specfile; specfile.open(specfname.c_str());
            for (int i=0; i<c->DimZ(); i++) 
                specfile << c->getZphys(i) << " " << spectrum.Array(i) << std::endl;
        }
    }
    //-----------------------------------------------------------------


    // PVs extraction task --------------------------------------------
    if (par->getFlagPV()) {
        BBreal xpos = par->getXPOS_PV();
        BBreal ypos = par->getYPOS_PV();
        BBreal ang = par->getPA_PV();
        std::string s = outfolder+c->Head().Name();
        Image2D<BBreal> *pv = PositionVelocity(c,xpos,ypos,ang);
        pv->fitswrite_2d((s+"pv_"+to_string(ang,0)+".fits").c_str());
        delete pv;
    }
    //-----------------------------------------------------------------


    // Repixeling task ------------------------------------------------
    if (par->getflagReduce() && !(par->getflagSmooth() || par->getflagSmoothSpectral())) {
        std::string name = c->pars().getOutfolder()+c->Head().Name()+"_red.fits";
        Cube<BBreal> *red = c->Reduce(floor(c->pars().getFactor()));
        red->fitswrite_3d(name.c_str(),true);
        delete red;
    }
    //-----------------------------------------------------------------


    // Kinematics fitting to slit data --------------------------------
    if (par->getFlagSlitfit()) {
        Model::Galfit<BBreal> *sfit = new Model::Galfit<BBreal>;
        sfit->slit_init(c);
        sfit->galfit();
        if (par->getParGF().TWOSTAGE) sfit->SecondStage();
        sfit->writeModel_slit();
    }
    //-----------------------------------------------------------------


    // Radial profile ------------------------------------------------
    if (par->getFlagEllProf()) {
        Tasks::Ellprof<BBreal> *ell = new Tasks::Ellprof<BBreal>(c);
        ell->RadialProfile();
        std::string fout = c->pars().getOutfolder()+c->Head().Name()+"_densprof.txt";
        std::ofstream fileo(fout.c_str());
        ell->printProfile(fileo,ell->getNseg()-1);
        ell->printProfile(std::cout,ell->getNseg()-1);
        ell->writeMap(c->pars().getOutfolder()+c->Head().Name()+"_densmap.fits");
        fileo.close();
        delete ell;
    }
    //-----------------------------------------------------------------
    
    delete c;

    return true;
}


#endif
