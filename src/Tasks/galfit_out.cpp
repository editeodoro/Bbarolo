//------------------------------------------------------------------------
// galfit_out.cpp: Members functions for outputs of the Galfit class.
//------------------------------------------------------------------------

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
#include <cmath>
#include <cfloat>
#include <Arrays/cube.hh>
#include <Arrays/image.hh>
#include <Tasks/galfit.hh>
#include <Tasks/galmod.hh>
#include <Tasks/smooth3D.hh>
#include <Utilities/utils.hh>
#include <Utilities/gnuplot.hh>
#include <Tasks/moment.hh>
#include <Tasks/ellprof.hh>

//#ifdef HAVE_PYTHON
//    #include <Python.h>
//#endif


namespace Model {


template <class T>
void Galfit<T>::writeModel (std::string normtype, bool makeplots) {

    bool verb = in->pars().isVerbose();
    in->pars().setVerbosity(false);
    in->pars().setFluxConvert(false);

    if (verb) std::cout << " Preparing a bunch of cool outputs..." << std::endl;

    std::string outfold = in->pars().getOutfolder();
    std::string prefix = in->pars().getOutPrefix();
    
    // Calculate radial profile along the output rings
    if (verb) std::cout << "    Deriving " << randomAdjective(1) << " radial profile..." << std::flush;

    Tasks::Ellprof<T> ell(in,outr,true);
    ell.RadialProfile();
    
    // If a azimuthal normalization is requested, using the azimuthal profile as input for output rings.
    // This is a little trick to ensure a smoother azimuthal model -> scaling factor avoid very small/large fluxes. 
    float scaling_factor = 1;
    if (normtype=="AZIM" || normtype=="BOTH") {
        double profmin=FLT_MAX;
        for (auto i=0; i<outr->nr; i++) {
            double mean = ell.getSurfDensFaceOn(i);
            if (!isNaN(mean) && profmin>mean && mean>0) profmin = mean;
        }
        while(profmin<0.1) {
            profmin*=10;
            scaling_factor *=10;
        }
        while (profmin>10) {
            profmin /= 10;
            scaling_factor /= 10;
        }
        for (auto i=0; i<outr->nr; i++) {
            if (outr->dens[i]>0) outr->dens[i] = scaling_factor*fabs(ell.getSurfDensFaceOn(i));
            //if (outr->dens[i]==0) outr->dens[i]=profmin;
        }
    }

    std::string dens_out = outfold+"densprof.txt";
    std::ofstream fileo(dens_out.c_str());
    ell.printProfile(fileo,ell.getNseg()-1);
    fileo.close();
    //ell.printProfile(std::cout);
    if (verb) std::cout << " Done." << std::endl;

    // Deleting bad rings, for which I set dens=0 in galfit()
    //for (auto i=0; i<outr->nr; i++) 
    //    if (outr->dens[i]==0) outr->deleteRing(i);

    // Calculate the region within the last ring
    T *ringreg = RingRegion(outr,in->Head());

    if (verb) std::cout << "    Calculating the very last model..." << std::flush;
    int bhi[2] = {in->DimX(), in->DimY()};
    int blo[2] = {0,0};
    Rings<T> *last = new Rings<T>;
    *last = *outr;
    Model::Galmod<T> *mod = getModel(last,bhi,blo,nullptr,true);
    delete last;
    mod->Out()->Head().setMinMax(0.,0.);
    mod->Out()->Head().setName(in->Head().Name()+"mod");
    
    T *outarray = mod->Out()->Array();

    if (normtype=="AZIM" || normtype=="BOTH") {

        // The final model has been build from the azimuthal profile calculated before,
        // thus just need to rescale the model by the same factor used for the profile.
        
        for (auto i=in->NumPix(); i--;) outarray[i] /= scaling_factor;

        if (verb) std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing " << randomAdjective(1) << " azimuthally-normalized model..." << std::flush;
        std::string mfile = outfold+prefix+"mod_azim.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());

        writePVs(mod->Out(),"_azim");
        if (verb) std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing " << randomAdjective(2) << " kinematic maps..." << std::flush;
        writeKinematicMaps(mod->Out(),"_azim");
        if (verb) std::cout << " Done." << std::endl;
    }

    if (normtype=="LOCAL" || normtype=="BOTH") {

        // The final model is normalized pixel-by-pixel, i.e. in each spaxel, the total
        // flux of observation and model is eventually equal.
        for (int y=0; y<in->DimY(); y++) {
            for (int x=0; x<in->DimX(); x++) {
                T factor = 0;
                // Comment next "if" for 2D ring
                //if (!isNaN(ringreg[x+y*in->DimX()])) {
                    T modSum = 0;
                    T obsSum = 0;
                    for (int z=0; z<in->DimZ(); z++) {
                        long Pix = in->nPix(x,y,z);
                        modSum += outarray[Pix];
                        obsSum += in->Array(Pix)*mask[Pix];
                    }
                    if (modSum!=0) factor = obsSum/modSum;
                //}
                for (int z=0; z<in->DimZ(); z++)
                    outarray[in->nPix(x,y,z)] *= factor;
            }
        }
        if (verb && normtype!="BOTH") std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing " << randomAdjective(1) << " locally-normalized model..." << std::flush;
        std::string mfile = outfold+prefix+"mod_local.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());

        writePVs(mod->Out(),"_local");
        if (verb) std::cout << " Done." << std::endl;

        if (verb) {
            std::cout << "    Writing " << std::flush;
            if (normtype=="BOTH") std::cout << "even more " << std::flush;
            std::cout << randomAdjective(2) << " kinematic maps..." << std::flush;
        }
        writeKinematicMaps(mod->Out(),"_local");
        if (verb) std::cout << " Done." << std::endl;
    }

    if (normtype=="NONE") {
        
        // No renormaliza needed, just write the model as it is.
        
        if (verb) std::cout << " Done." << std::endl;
        
        if (verb) std::cout << "    Writing model..." << std::flush;
        std::string mfile = outfold+prefix+"mod_none.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());
        writePVs(mod->Out(),"_none");
        if (verb) std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing " << randomAdjective(1) << " kinematic maps..." << std::flush;
        writeKinematicMaps(mod->Out(),"_none");
        if (verb) std::cout << " Done." << std::endl;
    }
    
    // Adding noise to a model
    if (par.NOISERMS!=0) { 
        if (verb) std::cout << "    Writing " << randomAdjective(1) << " noisy model..." << std::flush;
        mod->addnoise(par.NOISERMS);
        // Writing to FITS file
        std::string mfile = outfold+prefix+"mod_noise.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());
        if (verb) std::cout << " Done." << std::endl;
    }
 
    // Computing asymmetric drift correction
    if (par.flagADRIFT) {
        if (verb) std::cout << "    Computing asymmetric drift correction..." << std::flush;
        T *dens_m = new T[outr->nr];
        int strad = par.STARTRAD<inr->nr ? par.STARTRAD : 0;
        for (auto i=outr->nr; i--;) dens_m[i] = ell.getMedian(i);
        bool ok = AsymmetricDrift(&outr->radii[strad],&dens_m[strad],&outr->vdisp[strad],
                                  &outr->vrot[strad],&outr->inc[strad],outr->nr-strad);
        if (verb) {
            if (ok) std::cout << " Done." << std::endl;
            else std::cout << " Failed." << std::endl;
        }
    }

    // Now plotting everything
    if (makeplots) {
        if (verb) std::cout << "    Producing " << randomAdjective(1) << " plots..." << std::flush;
        int ret = plotParam();
        if (verb) {
            if (ret==0) std::cout << " Done." << std::endl;
            else std::cout << " Something went wrong! Check plotting scripts in the output folder." << std::endl;
        }
    }
    in->pars().setVerbosity(verb);

    if (verb) std::cout << " All done!" << std::endl;

    delete mod;
    
}
template void Galfit<float>::writeModel(std::string, bool);
template void Galfit<double>::writeModel(std::string, bool);


template <class T>
void Galfit<T>::writeOutputs (Cube<T> *mod, Tasks::Ellprof<T> *e, bool makeplots) {
    // Write all the output files and plots
    // Mostly used by BayesianBBarolo.

    std::string suffix = "_none";
    std::string outfold = in->pars().getOutfolder();
    std::string prefix = in->pars().getOutPrefix();

    // Write the density profile
    std::ofstream fileo((outfold+"densprof.txt").c_str());
    e->printProfile(fileo,e->getNseg()-1);
    fileo.close();

    // Write the model cube
    mod->Head().setMinMax(0.,0.);
    mod->Head().setName(prefix+"mod");
    std::string mfile = outfold+prefix+"mod"+suffix+".fits";
    mod->fitswrite_3d(mfile.c_str());

    // Write P-V and kinematic maps
    writePVs(mod,suffix);
    writeKinematicMaps(mod,suffix);

    if (makeplots) plotParam();
}
template void Galfit<float>::writeOutputs(Cube<float>*, Tasks::Ellprof<float>*, bool);
template void Galfit<double>::writeOutputs(Cube<double>*, Tasks::Ellprof<double>*, bool);


template <class T>
void Galfit<T>::writeKinematicMaps(Cube<T> *mod, std::string suffix) {

    std::string outfold = in->pars().getOutfolder();
    std::string prefix  = in->pars().getOutPrefix();

    // Write the total intensity, velocity and dispersion maps of data
    mkdirp((outfold+"maps/").c_str());

    std::vector< MomentMap<T> > allmaps = getAllMoments<T>(in,true,nullptr,"MOMENT");
    allmaps[0].fitswrite_2d((outfold+"maps/"+prefix+"_0mom.fits").c_str());
    allmaps[1].fitswrite_2d((outfold+"maps/"+prefix+"_1mom.fits").c_str());
    allmaps[2].fitswrite_2d((outfold+"maps/"+prefix+"_2mom.fits").c_str());
    
    // Write the total intensity, velocity and dispersion maps of model
    std::vector< MomentMap<T> > modmaps = getAllMoments<T>(mod,false,nullptr,"MOMENT");
    modmaps[0].fitswrite_2d((outfold+"maps/"+prefix+"mod_0mom"+suffix+".fits").c_str());
    modmaps[1].fitswrite_2d((outfold+"maps/"+prefix+"mod_1mom"+suffix+".fits").c_str());
    modmaps[2].fitswrite_2d((outfold+"maps/"+prefix+"mod_2mom"+suffix+".fits").c_str());

}
template void Galfit<float>::writeKinematicMaps(Cube<float>*, std::string);
template void Galfit<double>::writeKinematicMaps(Cube<double>*, std::string);


template <class T>
void Galfit<T>::writePVs(Cube<T> *mod, std::string suffix) {

    // Extracts PVs along major and minor axis and write them in fits.
    std::string outfold = in->pars().getOutfolder()+"pvs/";
    std::string prefix  = in->pars().getOutPrefix();
    
    mkdirp((outfold).c_str());
    
    float meanPA = findMedian(&outr->phi[0], outr->nr);
    float meanXpos = findMedian(&outr->xpos[0], outr->nr);
    float meanYpos = findMedian(&outr->ypos[0], outr->nr);
    float meanPAp90= meanPA+90<360 ? meanPA+90 : meanPA-90;
    
    // Extract pvs of data
    PvSlice<T> *pv_max = PositionVelocity(in,meanXpos,meanYpos,meanPA);
    std::string mfile = outfold+prefix+"_pv_a.fits";
    pv_max->fitswrite_2d(mfile.c_str());
    PvSlice<T> *pv_min = PositionVelocity(in,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+prefix+"_pv_b.fits";
    pv_min->fitswrite_2d(mfile.c_str());

    // Extract pvs of model
    PvSlice<T> *pv_max_m = PositionVelocity(mod,meanXpos,meanYpos,meanPA);
    mfile = outfold+prefix+"mod_pv_a"+suffix+".fits";
    pv_max_m->fitswrite_2d(mfile.c_str());
    PvSlice<T> *pv_min_m = PositionVelocity(mod,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+prefix+"mod_pv_b"+suffix+".fits";
    pv_min_m->fitswrite_2d(mfile.c_str());

    // Extract pvs of mask
    Cube<short> *m = new Cube<short>(in->AxisDim());
    m->saveHead(in->Head());
    m->saveParam(in->pars());
    m->pars().setANTIALIAS(0);
    m->Head().setMinMax(0.,0);
    for (size_t i=in->NumPix(); i--;) m->Array()[i] = short(mask[i]);
    
    PvSlice<short> *pv_max_ma = PositionVelocity(m,meanXpos,meanYpos,meanPA);
    mfile = outfold+prefix+"mask_pv_a.fits";
    pv_max_ma->fitswrite_2d(mfile.c_str());
    PvSlice<short> *pv_min_ma = PositionVelocity(m,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+prefix+"mask_pv_b.fits";
    pv_min_ma->fitswrite_2d(mfile.c_str());
    
    
#ifndef HAVE_PYTHON
#ifdef HAVE_GNUPLOT
    plotPVs_Gnuplot(pv_max,pv_min,pv_max_m,pv_min_m);
#endif
#endif

    delete pv_max; delete pv_min;
    delete pv_max_m; delete pv_min_m;
    delete m;
    
    delete pv_max_ma; delete pv_min_ma;
    
}
template void Galfit<float>::writePVs(Cube<float>*, std::string);
template void Galfit<double>::writePVs(Cube<double>*, std::string);


template <class T>
void Galfit<T>::plotPVs_Gnuplot(Image2D<T> *pva_d, Image2D<T> *pvb_d, Image2D<T> *pva_m, Image2D<T> *pvb_m) {

    /* Plot position-velocity diagrams by using Gnuplot */

    T meanPA = findMean(&outr->phi[0], outr->nr);
    T meanPAp90= meanPA+90<360 ? meanPA+90 : meanPA-90;
    T meanVsys = findMean(&outr->vsys[0], outr->nr);

    std::string outfold = in->pars().getOutfolder();

    /// Gnuplot does not read fitsfile, so write contours in text files.
    /// Processing the PVs of data first.
    Image2D<T> *pv_max = pva_d;
    Image2D<T> *pv_min = pvb_d;
    std::ofstream outpv((outfold+"pv.txt").c_str());
    std::ofstream outpvm((outfold+"pvm.txt").c_str());
    float xmin=FLT_MAX,xmax=-FLT_MAX,xmmin=FLT_MAX,xmmax=-FLT_MAX;
    for (int y=0; y<pv_max->DimY(); y++) {
        for (int x=0;x<pv_max->DimX(); x++) {
            int i = x+y*pv_max->DimX();
            float xphys = (x+1-pv_max->Head().Crpix(0))*pv_max->Head().Cdelt(0)+pv_max->Head().Crval(0);
            float yphys = (y+1-pv_max->Head().Crpix(1))*pv_max->Head().Cdelt(1)+pv_max->Head().Crval(1);
            xphys *=arcconv;
            yphys = AlltoVel(yphys,pv_max->Head());
            if (fabs(xphys)<outr->radii[outr->nr-1]+5*outr->radsep) {
                outpv << xphys << "   " << yphys << "  " << pv_max->Array(i) << endl;
                if (xphys>xmax) xmax=xphys;
                if (xphys<xmin) xmin=xphys;
            }
        }
        for (int x=0;x<pv_min->DimX();x++) {
            int i = x+y*pv_min->DimX();
            float xphys = (x+1-pv_min->Head().Crpix(0))*pv_min->Head().Cdelt(0)+pv_min->Head().Crval(0);
            float yphys = (y+1-pv_min->Head().Crpix(1))*pv_min->Head().Cdelt(1)+pv_min->Head().Crval(1);
            xphys *=arcconv;
            yphys = AlltoVel(yphys,pv_min->Head());
            if (fabs(xphys)<outr->radii[outr->nr-1]+5*outr->radsep) {
                outpvm << xphys << "   " << yphys << "  " << pv_min->Array(i) << endl;
                if (xphys>xmmax) xmmax=xphys;
                if (xphys<xmmin) xmmin=xphys;
            }
        }
        outpv << endl;
        outpvm << endl;
    }
    outpv.close();
    outpvm.close();

    /// Writing in a text file the projected rotation curve.
    outpv.open((outfold+"rcpv.txt").c_str());
    for (int i=par.STARTRAD; i<outr->nr; i++) {
            float vel1 = (outr->vrot[i]*sin(outr->inc[i]*M_PI/180.))+outr->vsys[i];
            float vel2 = outr->vsys[i]-(outr->vrot[i]*sin(outr->inc[i]*M_PI/180.));
            if (meanPA>180.) std::swap(vel1,vel2);
            float radius = outr->radii[i];
            outpv << -radius << "   " << vel1 << endl;
            outpv <<  radius << "   " << vel2 << endl;
    }
    outpv.close();


    /// Processing the PVs of model.
    pv_max = pva_m;
    pv_min = pvb_m;
    outpv.open((outfold+"pv_mod.txt").c_str());
    outpvm.open((outfold+"pvm_mod.txt").c_str());
    for (int y=0; y<pv_max->DimY() ; y++) {
        for (int x=0;x<pv_max->DimX();x++) {
            int i = x+y*pv_min->DimX();
            float xphys = (x+1-pv_max->Head().Crpix(0))*pv_max->Head().Cdelt(0)+pv_max->Head().Crval(0);
            float yphys = (y+1-pv_max->Head().Crpix(1))*pv_max->Head().Cdelt(1)+pv_max->Head().Crval(1);
            yphys = AlltoVel(yphys,pv_max->Head());
            outpv << xphys*arcconv << "   " << yphys << "  " << pv_max->Array(i) << endl;
        }
        for (int x=0;x<pv_min->DimX();x++) {
            int i = x+y*pv_min->DimX();
            float xphys = (x+1-pv_min->Head().Crpix(0))*pv_min->Head().Cdelt(0)+pv_min->Head().Crval(0);
            float yphys = (y+1-pv_min->Head().Crpix(1))*pv_min->Head().Cdelt(1)+pv_min->Head().Crval(1);
            yphys = AlltoVel(yphys,pv_min->Head());
            outpvm << xphys*arcconv << "   " << yphys << "  " << pv_min->Array(i) << endl;
        }
        outpv << endl;
        outpvm << endl;
    }
    outpv.close();
    outpvm.close();


    /// Writing a gnuplot script
    std::string conlevels;
    float sig;
    if (!in->StatsDef()) in->setCubeStats();
    if (in->pars().getParSE().flagUserGrowthT) sig=in->pars().getParSE().threshold;
    else sig = 2.0*in->stat().getSpread();
    int k=0;
    while (sig<in->stat().getMax()) {
        conlevels += to_string(sig)+",";
        sig *= 3;
        k++;
        if (k>10000) break;
    }
    conlevels.erase(conlevels.end()-1);

    float vmin = AlltoVel(in->getZphys(0), in->Head());
    float vmax = AlltoVel(in->getZphys(in->DimZ()-1), in->Head());
    if (vmin>vmax) std::swap(vmin,vmax);
    std::string mfile=outfold+"pv.gnu";
    outpv.open(mfile.c_str());

    outpv << "set contour base"  << endl
          << "set cntrparam levels discrete "<< conlevels << endl
          << "unset surface"  << endl
          << "set table 'cont.tab'" << endl
          << "splot '"<<outfold+"pv.txt'" << endl
          << "set table 'cont_mod.tab'" << endl
          << "splot '"<<outfold+"pv_mod.txt'" << endl
          << "set table 'contm.tab'" << endl
          << "splot '"<<outfold+"pvm.txt'" << endl
          << "set table 'contm_mod.tab'" << endl
          << "splot '"<<outfold+"pvm_mod.txt'" << endl
          << "unset table" << endl
          << "unset key" << endl;

    mfile = outfold+"pv_a_cfr.eps";
    outpv << "set terminal postscript eps enhanced color font 'Helvetica,14'" << endl
          << "set output '" << mfile << "'" << endl
          << "set style line 1 lc rgb '#7F7F7F' lt 2 pt 0 lw 1" << endl
          << "set style line 2 lc rgb '#B22222' lt -1 pt 0 lw 1" << endl
          << "set style line 3 lc rgb '#00008B' lt -1 pt 5 ps 0.6 lw 1" << endl
          << "set xlabel 'Offset [arcsec]'" << endl
          << "set ylabel 'Velocity [km/s]'" << endl
          << "set yrange ["<< to_string(vmin)+":"<< to_string(vmax) << "]" << endl
          << "set yzeroaxis lt -1" << endl
          << "set xrange ["<< to_string(xmin)+":"<< to_string(xmax) << "]" << endl
          << "set y2label ' {/Symbol f} = " << to_string(meanPA) << " deg'" << endl
          << "plot " << to_string(meanVsys) << " ls -1 lc -1, 'cont.tab' w l ls 1, 'cont_mod.tab' w l ls 2, '"<<outfold<<"rcpv.txt' w p ls 3 "<< endl;

    mfile = outfold+"pv_b_cfr.eps";
    outpv << "set output '" << mfile << "'" << endl
          << "set y2label ' {/Symbol f} = " << to_string(meanPAp90) << " deg'" << endl
          << "plot " << to_string(meanVsys) << " ls -1 lc -1, 'contm.tab'w l ls 1, 'contm_mod.tab' w l ls 2"<< endl;

#ifndef HAVE_PYTHON
#ifdef HAVE_GNUPLOT
    Gnuplot gp;
    gp.begin();
    mfile = "load '"+outfold+"pv.gnu'";
    gp.commandln(mfile.c_str());
    gp.end();
    remove ("cont.tab");
    remove ("contm.tab");
    remove ("cont_mod.tab");
    remove ("contm_mod.tab");
#endif
#endif

    remove ((outfold+"pv.txt").c_str());
    remove ((outfold+"pvm.txt").c_str());
    remove ((outfold+"pv_mod.txt").c_str());
    remove ((outfold+"pvm_mod.txt").c_str());
    remove ((outfold+"rcpv.txt").c_str());
    remove ((outfold+"pv.gnu").c_str());


}
template void Galfit<float>::plotPVs_Gnuplot(Image2D<float>*,Image2D<float>*,Image2D<float>*,Image2D<float>*);
template void Galfit<double>::plotPVs_Gnuplot(Image2D<double>*,Image2D<double>*,Image2D<double>*,Image2D<double>*);


template <class T>
int *Galfit<T>::getErrorColumns() {

    int free[nfree];
    int *nc = new int[MAXPAR];
    int k=0;
    for (int nm=0; nm<MAXPAR; nm++) if(mpar[nm]) free[k++]=nm;

    const int err_col=13;

    for (int j=0; j<MAXPAR; j++) {
        if (par.flagERRORS && mpar[j]) {
            nc[j]=err_col;
            for (int i=0; i<nfree; i++) if (free[i]==j) nc[j]+=2*i;
        }
        else nc[j]=-1;
    }

    return nc;
}
template int* Galfit<float>::getErrorColumns();
template int* Galfit<double>::getErrorColumns();


template <class T>
void Galfit<T>::plotPar_Gnuplot () {

    std::string outfold = in->pars().getOutfolder();
    std::string prefix  = in->pars().getOutPrefix();

    mkdirp((outfold+"plotscripts/").c_str());
    std::string mfile = outfold+"plotscripts/gnuscript.gnu";
    std::ofstream gnu(mfile.c_str());

    float xtics = ceil(outr->nr/5.);
    xtics *= outr->radsep;
    
    while (outr->radii.back()/xtics>5.) xtics*=2.;
    while (outr->radii.back()/xtics<2.) xtics/=2.;
    
    int *nc = getErrorColumns();

    /// Setting global option
    gnu << "set terminal postscript eps enhanced color font 'Helvetica,14'" << endl
        << "set output '" << outfold << prefix << "_rc_inc_pa.eps'" << endl
        << "unset key" << endl
        << "set size 0.60, 1" << endl;

    if (par.TWOSTAGE) {
        gnu << "set style line 1 lc rgb '#A9A9A9' lt 4 pt 7 lw 1" << endl
            << "set style line 2 lc rgb '#B22222' lt 9 pt 9 lw 1" << endl;
    }
    else {
        gnu << "set style line 1 lc rgb '#B22222' lt 9 pt 7 lw 1" << endl;
    }

    gnu << "set macros" << endl
        << "XTICS   = 'set xtics " << to_string(xtics) << "; set mxtics 2; set format x \"%g\" '" << endl
        << "NOXTICS = 'unset xlabel; set xtics  " << to_string(xtics) << "; set mxtics 2; set format x '' '" << endl
        << "LABELF  = 'set xlabel font \"Helvetica,13\"; "
        <<            "set ylabel font \"Helvetica,13\" '" << endl
        << "TICSF   = 'set xtics font \"Helvetica,12\"; "
        <<            "set ytics font \"Helvetica,12\" '" << endl
        << "TMARGIN = 'set tmargin at screen 0.95; set bmargin at screen 0.47; "
        <<            "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "MMARGIN = 'set tmargin at screen 0.47; set bmargin at screen 0.27; "
        <<            "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "BMARGIN = 'set tmargin at screen 0.27; set bmargin at screen 0.10; "
        <<            "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "set multiplot layout 3,1 rowsfirst" << endl;

    gnu << "@LABELF" << endl << "@TICSF" << endl;

    /// Plotting rotational velocity
    float maxvel = *max_element(&outr->vrot[0], &outr->vrot[0]+outr->nr);
    maxvel += 0.1*maxvel;
    gnu << "@TMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [-5:" << maxvel << "]" << endl
        << "set ylabel 'V_c  [km/s]'" << endl
        << "set ytics 50" << endl << "set mytics 5" << endl
        << "plot '" << outfold << "rings_final1.txt' ";

    if (nc[VROT]>=0) {
        gnu << "u 2:3:($3+$"+to_string(nc[VROT])+"):($3+$"+to_string(nc[VROT]+1)+") w errorbars ls 1, '"
            << outfold <<"rings_final1.txt' u 2:3 w lp ls 1";
    }
    else gnu << "u 2:3 w lp ls 1";

    if (par.TWOSTAGE) {
        gnu << ", '" << outfold << "rings_final2.txt' ";
        if (par.flagERRORS && mpar[VROT]) {
        gnu << "u 2:3:($3+$13):($3+$14) w errorbars ls 2, '"
            << outfold <<"rings_final2.txt' u 2:3 w lp ls 2";
        }
        else gnu << "u 2:3 w lp ls 2";
    }

    gnu << endl << "set title ''" << endl;
    /// Plotting inclination
    float maxa = *max_element(&outr->inc[0], &outr->inc[0]+outr->nr);
    maxa += 0.1*maxa;
    float mina = *min_element(&outr->inc[0], &outr->inc[0]+outr->nr);
    mina -= 0.1*mina;
    gnu << "@MMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [" << mina << ":" << maxa << "]" << endl
        << "set ylabel 'i [deg]'" << endl
        << "set ytics 5" << endl << "set mytics 5" << endl
        << "plot '" << outfold << "rings_final1.txt' ";

    if (nc[INC]>=0) {
        gnu << "u 2:5:($5+$"+to_string(nc[INC])+"):($5+$"+to_string(nc[INC]+1)+") w errorbars ls 1, '"
            << outfold << "rings_final1.txt' u 2:5 w lp ls 1";
    }
    else gnu << "u 2:5 w lp ls 1";

    if (par.TWOSTAGE)
        gnu << ", '" << outfold << "rings_final2.txt' u 2:5 w lp ls 2";

    gnu << endl;

    /// Plotting position angle
    maxa = *max_element(&outr->phi[0], &outr->phi[0]+outr->nr);
    maxa += 0.1*maxa;
    mina = *min_element(&outr->phi[0], &outr->phi[0]+outr->nr);
    mina -= 0.1*mina;

    gnu << "@BMARGIN" << endl << "@XTICS" << endl
        << "set xlabel 'Radius [arcsec]'" << endl
        << "set yrange [" << mina << ":" << maxa << "]" << endl
        << "set ylabel 'P.A. [deg]'" << endl
        << "set ytics 5" << endl << "set mytics 5" << endl
        << "plot '" << outfold << "rings_final1.txt' ";

    if (nc[PA]>=0) {
        gnu << "u 2:6:($6+$"+to_string(nc[PA])+"):($6+$"+to_string(nc[PA]+1)+") w errorbars ls 1, '"
            << outfold << "rings_final1.txt' u 2:6 w lp ls 1";;
    }
    else gnu << "u 2:6 w lp ls 1";

    if (par.TWOSTAGE)
        gnu << ", '" << outfold << "rings_final2.txt' u 2:6 w lp ls 2";

    gnu << endl;

    gnu << "unset multiplot" << endl;

    gnu << "set output '" << outfold << prefix << "_disp_vsys_z0.eps'" << endl
        << "unset key" << endl
        << "set xlabel 'Radius [arcsec]'" << endl
        << "set xtics 200" << endl
        << "set mxtics 2" << endl
        << "set macros" << endl
        << "TMARGIN = 'set tmargin at screen 0.94; set bmargin at screen 0.66; "
        <<            "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "MMARGIN = 'set tmargin at screen 0.66; set bmargin at screen 0.38; "
        <<            "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "BMARGIN = 'set tmargin at screen 0.38; set bmargin at screen 0.10; "
        <<            "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "set multiplot layout 3,1 rowsfirst" << endl;

    gnu << "@LABELF" << endl << "@TICSF" << endl;

    // Plotting dispersion velocity
    maxa = *max_element(&outr->vdisp[0], &outr->vdisp[0]+outr->nr);
    maxa += 0.1*maxa;
    gnu << "@TMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [0:"<<maxa<<"]\n"
        << "set ylabel '{/Symbol s} [km/s]'\n"
        << "set ytics 5" << endl << "set mytics 5" << endl
        << "plot '"<< outfold << "rings_final1.txt' ";

    if (nc[VDISP]>=0) {
        gnu << "u 2:4:($4+$"+to_string(nc[VDISP])+"):($4+$"+to_string(nc[VDISP]+1)+") w errorbars ls 1, '"
            << outfold <<"rings_final1.txt' u 2:4 w lp ls 1";
    }
    else gnu << "u 2:4 w lp ls 1";

    if (par.TWOSTAGE) {
        gnu << ", '" << outfold << "rings_final2.txt' ";
        if (par.flagERRORS && mpar[VDISP]) {
            gnu << "u 2:4:($3+$15):($3+$16) w errorbars ls 2, '"
                << outfold <<"rings_final2.txt' u 2:4 w lp ls 2";
        }
        else gnu << "u 2:4 w lp ls 2";
    }
    gnu << endl;


    // Plotting systemic velocity
    maxa = *max_element(&outr->vsys[0], &outr->vsys[0]+outr->nr);
    maxa += (0.1*maxa+10);
    mina = *min_element(&outr->vsys[0], &outr->vsys[0]+outr->nr);
    mina -= (0.1*mina+10);
    gnu << "@MMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [" << mina << ":" << maxa << "]" << endl
        << "set ylabel 'V_{sys} [km/s]'" << endl
        << "plot '" << outfold << "rings_final1.txt' ";

    if (nc[VSYS]>=0) {
        gnu << "u 2:12:($12+$"+to_string(nc[VSYS])+"):($12+$"+to_string(nc[VSYS]+1)+") w errorbars ls 1, '"
            << outfold << "rings_final1.txt' u 2:12 w lp ls 1";;
    }
    else gnu << "u 2:12 w lp ls 1";

    if (par.TWOSTAGE)
        gnu << ", '" << outfold << "rings_final2.txt' u 2:12 w lp ls 2";
    gnu << endl;


    // Plotting scale height
    maxa = *max_element(&outr->z0[0], &outr->z0[0]+outr->nr);
    maxa += 0.1*maxa;
    mina = *min_element(&outr->z0[0], &outr->z0[0]+outr->nr);
    mina -= 0.1*mina;
    gnu << "@BMARGIN" << endl << "@XTICS" << endl
        << "set xlabel 'Radius [arcsec]'" << endl
        << "set yrange [" <<mina<<":"<<maxa<<"]\n"
        << "set ylabel 'Scale height [arcsec]'\n"
        << "plot '" << outfold << "rings_final1.txt'";

    if (nc[Z0]>=0) {
        gnu << "u 2:8:($8+$"+to_string(nc[Z0])+"):($8+$"+to_string(nc[Z0]+1)+") w errorbars ls 1, '"
            << outfold << "rings_final1.txt' u 2:8 w lp ls 1";
    }
    else gnu << "u 2:8 w lp ls 1";

    if (par.TWOSTAGE)
        gnu << ", '" << outfold << "rings_final2.txt' u 2:8 w lp ls 2";
    gnu << endl;

    gnu << "unset multiplot" << endl;

    gnu << "set output '" << outfold << prefix << "_xc_yc_cd.eps'" << endl
        << "set multiplot layout 3,1 rowsfirst" << endl;

    gnu << "@LABELF" << endl << "@TICSF" << endl;

    // Plotting xcenter
    maxa = *max_element(&outr->xpos[0], &outr->xpos[0]+outr->nr);
    maxa += 0.1*maxa;
    mina = *min_element(&outr->xpos[0], &outr->xpos[0]+outr->nr);
    mina -= 0.1*mina;
    gnu << "@TMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [" <<mina<<":"<<maxa<<"]\n"
        << "set ylabel 'X_c [pix]'\n"
        //<< "set ytics 5" << endl << "set mytics 5" << endl
        << "plot '" << outfold << "rings_final1.txt' ";

    if (nc[XPOS]>=0) {
        gnu << "u 2:10:($10+$"+to_string(nc[XPOS])+"):($10+$"+to_string(nc[XPOS]+1)+") w errorbars ls 1, '"
            << outfold << "rings_final1.txt' u 2:10 w lp ls 1";;
    }
    else gnu << "u 2:10 w lp ls 1";

    if (par.TWOSTAGE)
        gnu << ", '" << outfold << "rings_final2.txt' u 2:10 w lp ls 2";
    gnu << endl;

    // Plotting ycenter
    maxa = *max_element(&outr->ypos[0], &outr->ypos[0]+outr->nr);
    maxa += 0.1*maxa;
    mina = *min_element(&outr->ypos[0], &outr->ypos[0]+outr->nr);
    mina -= 0.1*mina;
    gnu << "@MMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [" << mina << ":" << maxa << "]" << endl
        << "set ylabel 'Y_c [pix]'" << endl
        //<< "set ytics 5" << endl << "set mytics 5" << endl
        << "plot '" << outfold << "rings_final1.txt' ";

    if (nc[YPOS]>=0) {
        gnu << "u 2:11:($11+$"+to_string(nc[YPOS])+"):($11+$"+to_string(nc[YPOS]+1)+") w errorbars ls 1, '"
            << outfold << "rings_final1.txt' u 2:11 w lp ls 1";;
    }
    else gnu << "u 2:11 w lp ls 1";

    if (par.TWOSTAGE)
        gnu << ", '" << outfold << "rings_final2.txt' u 2:11 w lp ls 2";
    gnu << endl;


    if (mpar[DENS]) {
        maxa = *max_element(&outr->dens[0], &outr->dens[0]+outr->nr);
        maxa += 0.1*maxa;
        mina = *min_element(&outr->dens[0], &outr->dens[0]+outr->nr);
        mina -= 0.1*mina;
        gnu << "@BMARGIN" << endl << "@XTICS" << endl
            << "set xlabel 'Radius [arcsec]'" << endl
            << "set yrange [" <<mina<<":"<<maxa<<"]\n"
            << "set ylabel 'Surface density [10^20 atoms/cm^2]'\n"
            << "plot '" << outfold << "rings_final1.txt' u 2:8 w lp ls 1";

        if (par.TWOSTAGE)
            gnu << ", '" << outfold << "rings_final2.txt' u 2:8 w lp ls 2";
        gnu << endl;
    }

    gnu << "unset multiplot; reset" << endl;
    gnu.close();

#ifndef HAVE_PYTHON
#ifdef HAVE_GNUPLOT
    Gnuplot gp;
    gp.begin();
    mfile = "load '"+outfold+"plotscripts/gnuscript.gnu'";
    gp.commandln(mfile.c_str());
#endif
#endif

}
template void Galfit<float>::plotPar_Gnuplot();
template void Galfit<double>::plotPar_Gnuplot();


template <class T>
int Galfit<T>::plotAll_Python() {
    
    /// This function creates python script for plotting all quantities
    /// (see writeScripts_Python()):
    ///
    /// It needs all output fitsfiles to be in the output directory!
    
    std::vector<std::string> scriptnames = writeScripts_Python();
    
    int returnValue = 0;

#ifdef HAVE_PYTHON
    // This will execute the four scripts concurrently
    std::string cmd;
    for (int i=0; i<scriptnames.size(); i++)
        cmd += "python \""+in->pars().getOutfolder()+"plotscripts/"+scriptnames[i]+"\" > /dev/null 2>&1 & ";
    cmd.erase (cmd.end()-2);
    returnValue = system(cmd.c_str());
    //cmd = "python "+in->pars().getOutfolder()+"plotscripts/plot_all.py > /dev/null 2>&1";
    //returnValue = system(cmd.c_str());
#endif

    return returnValue;    
}
template int Galfit<float>::plotAll_Python();
template int Galfit<double>::plotAll_Python();


template <class T>
std::vector<std::string> Galfit<T>::writeScripts_Python() {
    
    /// This function creates a python script for plotting:
    ///  1) Channel maps
    ///  2) Position-Velocity diagrams along major/minor axis
    ///  3) Output paramenters
    ///  4) Moment maps
    ///  5) Asymmetric drift correction


    std::vector<std::string> scriptnames;
    std::ofstream pyf;

    float crpix3_kms = in->Head().Crpix(2);
    float cdelt3_kms = DeltaVel(in->Head());
    float crval3_kms = AlltoVel(in->Head().Crval(2),in->Head());
    float crval3_kms_pix1 = AlltoVel(in->Head().getZphys(0),in->Head());
    float bmaj = in->Head().Bmaj()/in->Head().PixScale();
    float bmin = in->Head().Bmin()/in->Head().PixScale();
    float bpa  = in->Head().Bpa();
    std::string bunit = "JY/BEAM * KM/S";
    if (FluxtoJy(1,in->Head())==1) bunit = in->Head().Bunit() + " * KM/S";

    float cont = 0;
    int xpos = findMedian(&outr->xpos[0],outr->nr);
    int ypos = findMedian(&outr->ypos[0],outr->nr);
    int xmin=0, ymin=0, zmin=0, disp=0;
    int xmax=in->DimX()-1, ymax=in->DimY()-1, zmax=in->DimZ()-1;
    float vsys_av = findMedian(&outr->vsys[0],outr->nr);

    // Determining contour levels for plots
    if (in->pars().getParGF().PLOTMINCON>0) cont = in->pars().getParGF().PLOTMINCON;
    else {
        if (in->pars().getMASK()=="SEARCH") {
            if (in->pars().getParSE().flagUserGrowthT) cont = in->pars().getParSE().growthThreshold;
            else if (in->pars().getParSE().UserThreshold) cont = in->pars().getParSE().threshold;
            else {
                if (in->pars().getParSE().flagGrowth) cont = in->pars().getParSE().growthCut*in->stat().getSpread();
                else cont = in->pars().getParSE().snrCut*in->stat().getSpread();
            }
        }
        else {
            if (!in->StatsDef()) in->setCubeStats();
            cont = 2.5*in->stat().getSpread();
            //if (in->pars().getParSE().UserThreshold) cont=in->pars().getParSE().threshold;
        }
    }

    // Determining x-y-z extension of plots
    if (in->pars().getMASK()=="SEARCH") {
        // Using largest detection to set the spatial displacement from the center and the spectral range
        Detection<T> *larg = in->getSources()->LargestDetection();
        long ext[4] = {abs(xpos-lround(larg->getXmin()-2*in->Head().Bmaj()/in->Head().PixScale())),
                       abs(xpos-lround(larg->getXmax()+2*in->Head().Bmaj()/in->Head().PixScale())),
                       abs(ypos-lround(larg->getYmin()-2*in->Head().Bmaj()/in->Head().PixScale())),
                       abs(ypos-lround(larg->getYmax()+2*in->Head().Bmaj()/in->Head().PixScale()))};
        disp = *max_element(&ext[0],&ext[0]+4);
        zmin = larg->getZmin()-3;
        zmax = larg->getZmax()+3;
    }
    else if (in->pars().getMASK()=="SMOOTH&SEARCH") {
        // Using mask to set the spatial displacement from the center and the spectral range
        bool *m = in->Mask();
        int Xmax=0,Xmin=in->DimX();
        int Ymax=0,Ymin=in->DimY();
        int Zmax=0,Zmin=in->DimZ();
        for (auto x=0; x<in->DimX(); x++) {
            for (auto y=0; y<in->DimY(); y++) {
                for (auto z=0; z<in->DimZ(); z++) {
                    if (m[in->nPix(x,y,z)]) {
                        if (x>Xmax) Xmax=x; if (y>Ymax) Ymax=y; 
                        if (z>Zmax) Zmax=z; if (x<Xmin) Xmin=x; 
                        if (y<Ymin) Ymin=y; if (z<Zmin) Zmin=z; 
                    }
                }
            }
        }
        long ext[4] = {abs(xpos-lround(Xmin-2*in->Head().Bmaj()/in->Head().PixScale())),
                       abs(xpos-lround(Xmax+2*in->Head().Bmaj()/in->Head().PixScale())),
                       abs(ypos-lround(Ymin-2*in->Head().Bmaj()/in->Head().PixScale())),
                       abs(ypos-lround(Ymax+2*in->Head().Bmaj()/in->Head().PixScale()))};
        disp = *max_element(&ext[0],&ext[0]+4);
        zmin = Zmin-3;
        zmax = Zmax+3;
    }
    else {
        // Using the LOS velocity and model extension to decide plotting boundaries.
        std::vector<T> maxv(outr->nr);
        for (int i=0; i<outr->nr; i++) maxv[i]=outr->vrot[i]*sin(outr->inc[i]*M_PI/180.)+outr->vdisp[i];
        float max_v=*max_element(&maxv[0],&maxv[0]+outr->nr);
    
        double vsys_spec = Vel2Spec(vsys_av,in->Head());
        int z_vsys = lround(in->getZgrid(vsys_spec));
        int disp_v = ceil((1.5*max_v)/fabs(cdelt3_kms));
        
        zmin = z_vsys-2*disp_v>0 ? z_vsys-2*disp_v : 0;
        zmax = z_vsys+2*disp_v<in->DimZ() ? z_vsys+2*disp_v : in->DimZ()-1;
        
        // Spatial boundaries to add 
        disp = fabs((outr->radii.back()/arcconv+2*in->Head().Bmaj())/in->Head().PixScale());
    }

    xmin = xpos-disp>=0 ? xpos-disp : 0;
    xmax = xpos+disp<in->DimX() ? xpos+disp : in->DimX()-1;
    ymin = ypos-disp>=0 ? ypos-disp : 0;
    ymax = ypos+disp<in->DimY() ? ypos+disp : in->DimY()-1;
    if (zmin<0) zmin=0;
    if (zmax>=in->DimZ()) zmax=in->DimZ()-1;


    int *nc = getErrorColumns();

    // All plotting scripts are written in a sub-directory
    mkdirp((in->pars().getOutfolder()+"plotscripts/").c_str());

    /////////////////////////////////////////////////////////////////////////
    /// Script with some general utilities
    /////////////////////////////////////////////////////////////////////////
    pyf.open((in->pars().getOutfolder()+"plotscripts/plot_utils.py").c_str());

    pyf << "###############################################################\n"
        << "#### This script contains some utilities for plotting      ####\n"
        << "###############################################################\n"
        << "import numpy as np \n"
        << "import matplotlib.pyplot as plt \n"
        << std::endl
        << "def defineaxis(nrows,ncols,xlen,ylen,xsep=0,ysep=0,fig_width=8.27,fig_heigth=11.69, fig=None):\n\n"
        << "  if not fig: fig = plt.figure(figsize=(fig_width, fig_heigth)) \n"
        << "  fig_width, fig_heigth = fig.get_size_inches() \n"
        << "  fig_ratio = fig_width / fig_heigth \n"
        << std::endl
        << "  ax = [] \n"
        << "  for i in range(nrows):\n"
        << "    row = []\n"
        << "    y = -(i * (ylen + ysep) * fig_ratio)\n"
        << "    for j in range(ncols):\n"
        << "      x = j * (xlen + xsep)\n"
        << "      row.append(fig.add_axes([x, y, xlen, ylen * fig_ratio])) \n"
        << "    ax.append(row)\n"
        << "  return fig, np.squeeze(ax)\n\n";

    pyf << "def calculate_padding(xsize,ysize):\n"
        << "  n = max(ysize, xsize) \n"
        << "  pad_y = ((n - ysize) // 2, n - ysize - (n - ysize) // 2) \n"
        << "  pad_x = ((n - xsize) // 2, n - xsize - (n - xsize) // 2) \n"
        << "  return (pad_x,pad_y)\n\n";
        
    pyf << "def build_path(x0, y0, rad_pix, pa_deg, nx, ny, oversampling=1):\n\n"
        << "  from scipy.interpolate import make_splprep\n"
        << "  # 1. Ring points \n"
        << "  pa_to_xy = lambda r, pa: (x0-r*np.sin(pa), y0+r*np.cos(pa))\n\n"
        << "  xy_pos = [pa_to_xy(r, np.radians(pa)) for r, pa in zip(rad_pix, pa_deg)]\n"
        << "  xy_neg = [pa_to_xy(-r, np.radians(pa)) for r, pa in zip(rad_pix[::-1], pa_deg[::-1])]\n"
        << "  xy_ring = np.concatenate([xy_neg,[[x0,y0]],xy_pos])\n\n"
        << "  # 2. Spline fit\n"
        << "  spl, u = make_splprep(xy_ring.T, s=0)\n"
        << "  u_inner = np.linspace(0, 1, xy_ring.shape[0]*3)\n"
        << "  x_s, y_s = spl(u_inner)\n\n"
        << "  # 3. Tangents at ends\n"
        << "  t0 = spl(0.0, nu=1)\n"
        << "  t1 = spl(1.0, nu=1)\n"
        << "  t0 /= np.hypot(*t0)\n"
        << "  t1 /= np.hypot(*t1)\n\n"
        << "  # 4. Extend linearly\n"
        << "  n_ext = int(np.maximum(nx,ny)/2)\n"
        << "  s_ext = np.linspace(1, n_ext + 1, n_ext*oversampling)\n"
        << "  xf = x_s[-1] + s_ext * t1[0]\n"
        << "  yf = y_s[-1] + s_ext * t1[1]\n"
        << "  xb = x_s[0] - s_ext * t0[0]\n"
        << "  yb = y_s[0] - s_ext * t0[1]\n\n"
        << "  # 5. Combine full slit and select points within field\n"
        << "  x = np.concatenate([xb[::-1], x_s, xf])\n"
        << "  y = np.concatenate([yb[::-1], y_s, yf])\n"
        << "  mask = (x >= 0) & (x <= nx - 1) & (y >= 0) & (y <= ny - 1)\n"
        << "  x, y = x[mask], y[mask]\n\n"
        << "  # 6. Calculating distances along the path\n"
        << "  icenter = np.argmin((x - x0)**2 + (y - y0)**2)\n"
        << "  ds = np.hypot(np.diff(x), np.diff(y))\n"
        << "  s = np.concatenate([[0.0], np.cumsum(ds)])\n\n"
        << "  # 7. Sorting and interpolating across a uniform grid\n"
        << "  idx = np.argsort(s)\n"
        << "  x, y, s = x[idx], y[idx], s[idx]-s[icenter]\n"
        << "  sp = np.arange(s.min(),s.max(),1/oversampling)\n"
        << "  xp, yp = np.interp(sp,s,x), np.interp(sp,s,y)\n\n"
        << "  return sp,xp,yp\n\n";
       
    pyf << "def extract_pv(array, x, y, order=1, width=1, oversampling=1):\n\n"
        << "  from scipy.ndimage import map_coordinates \n"
        << "  zsize, ysize, xsize = array.shape[:3] \n"
        << "  zs = np.outer(np.arange(0, zsize, 1.0/oversampling), np.ones(len(x)))\n"
        << "  z_len = zs.shape[0]\n\n"
        << "  dx = np.gradient(x)\n"
        << "  dy = np.gradient(y)\n"
        << "  norm = np.hypot(dx, dy)\n"
        << "  nx, ny = dy / norm, -dx / norm\n\n"
        << "  nwidth = 2*width + 1\n"
        << "  offsets = np.linspace(-width/2.0, width/2.0, nwidth)\n\n"
        << "  has_nan = np.any(np.isnan(array))\n"
        << "  samples = []\n"
        << "  for off in offsets:\n"
        << "    xs = x + off * nx\n"
        << "    ys = y + off * ny\n"
        << "    xs_grid = np.outer(np.ones(z_len), xs)\n"
        << "    ys_grid = np.outer(np.ones(z_len), ys)\n\n"
        << "    if has_nan:\n"
        << "      vals = map_coordinates(np.nan_to_num(array),[zs,ys_grid,xs_grid],order=int(order),cval=np.nan)\n"
        << "      bad  = map_coordinates(np.isnan(array).astype(float),[zs,ys_grid,xs_grid],order=0,cval=1.0) > 0\n"
        << "      vals[bad] = np.nan\n"
        << "    else:\n"
        << "      vals = map_coordinates(array,[zs,ys_grid,xs_grid],order=int(order),cval=np.nan)\n\n"
        << "    samples.append(vals)\n"
        << "  return np.nanmean(np.stack(samples, axis=0), axis=0)\n";
    
    pyf.close();


    /////////////////////////////////////////////////////////////////////////
    /// Script to plot best-fit parameters and asymmetric drift correction
    /////////////////////////////////////////////////////////////////////////
    scriptnames.push_back("plot_parameters.py");
    pyf.open((in->pars().getOutfolder()+"plotscripts/"+scriptnames[0]).c_str());

    pyf << "###############################################################\n"
        << "#### This script produces plots of the best-fit parameters ####\n"
        << "###############################################################\n"
        << "import numpy as np \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "lsize = 11 \n"
        << "mpl.rc('xtick',direction='in',top=True) \n"
        << "mpl.rc('ytick',direction='in',right=True) \n"
        << "mpl.rcParams['contour.negative_linestyle'] = 'solid' \n"
        << "plt.rc('font',family='sans-serif',serif='Helvetica',size=10) \n"
        << "params = {'text.usetex': False, 'mathtext.fontset': 'cm', 'mathtext.default': 'regular', 'errorbar.capsize': 0} \n"
        << "plt.rcParams.update(params) \n"
        << std::endl
        << "outprefix = '" << in->pars().getOutPrefix() <<"' \n"
        << "outfolder = '" << in->pars().getOutfolder() <<"' \n"
        << "twostage = " << par.TWOSTAGE << " \n\n";

    pyf << "# Reading in best-fit parameters \n"
        << "rad_sd, surfdens, sd_err = np.genfromtxt(outfolder+'densprof.txt', usecols=(0,3,4),unpack=True) \n"
        << "rad,vrot,disp,inc,pa,z0,xpos,ypos,vsys,vrad = np.genfromtxt(outfolder+'rings_final1.txt',usecols=(1,2,3,4,5,7,9,10,11,12),unpack=True) \n";
    if (outr->nr==1) pyf << "rad,vrot,disp,inc,pa,z0,xpos,ypos,vsys,vrad = np.array([rad]),np.array([vrot]),np.array([disp]),np.array([inc]),np.array([pa]),np.array([z0]),np.array([xpos]),np.array([ypos]),np.array([vsys]),np.array([vrad])\n";
    pyf << "try: err1_l, err1_h = np.zeros(shape=(" << MAXPAR << ",len(rad))), np.zeros(shape=(" << MAXPAR << ",len(rad)))\n"
        << "except: err1_l, err1_h = np.zeros(shape=(" << MAXPAR << ",1)), np.zeros(shape=(" << MAXPAR << ",1))\n"
        << "color=color2='#B22222' \n"
        << "max_vrot,max_vdisp = np.nanmax(vrot),np.nanmax(disp) \n"
        << "max_rad = 1.1*np.nanmax(rad) \n";

    for (int i=0; i<MAXPAR; i++) {
        if (nc[i]>0) {
            pyf << "err1_l[" << i << "], err1_h[" << i << "] = np.genfromtxt(outfolder+'rings_final1.txt',usecols=("
                << nc[i] << "," << nc[i]+1 << "),unpack=True) \n";
        }
    }

    pyf << "\nif twostage: \n"
        << "  rad2, vrot2,disp2,inc2,pa2,z02,xpos2,ypos2,vsys2, vrad2 = np.genfromtxt(outfolder+'rings_final2.txt',usecols=(1,2,3,4,5,7,9,10,11,12),unpack=True)\n";
    if (outr->nr==1) pyf << "  rad2,vrot2,disp2,inc2,pa2,z02,xpos2,ypos2,vsys2,vrad2 = np.array([rad2]),np.array([vrot2]),np.array([disp2]),np.array([inc2]),np.array([pa2]),np.array([z02]),np.array([xpos2]),np.array([ypos2]),np.array([vsys2]),np.array([vrad2])\n";
    pyf << "  err2_l, err2_h = np.zeros(shape=(" << MAXPAR << ",len(rad2))), np.zeros(shape=(" << MAXPAR << ",len(rad2)))\n"
        << "  color='#A0A0A0' \n"
        << "  max_rad = 1.1*np.nanmax(rad2) \n"
        << "  max_vrot,max_vdisp = np.maximum(max_vrot,np.nanmax(vrot2)),np.maximum(max_vdisp,np.nanmax(disp2)) \n";

    for (int i=0; i<MAXPAR; i++) {
        if ((i==0 || i==1 || i==9) && (nc[i]>0)) {
            int ff = i==9 ? 6 : 0;
            pyf << "  err2_l[" << i << "], err2_h[" << i << "] = np.genfromtxt(outfolder+'rings_final2.txt',usecols=("
                << nc[i]-ff << "," << nc[i]+1-ff << "),unpack=True) \n";
        }
    }

    pyf << "\n# Defining figure and axes \n"
        << "fig=plt.figure(figsize=(11,11),dpi=150) \n"
        << "nrows, ncols = 3,3 \n"
        << "xlen, ylen = 0.27, 0.13 \n"
        << "x_sep, y_sep = 0.07,0.015 \n"
        << "bottom_corner = [0.1,0.7] \n"
        << "for i in range (nrows): \n"
        << "  bottom_corner[0], yl = 0.1, ylen \n"
        << "  if i==0: yl *= 1.8 \n"
        << "  for j in range (ncols): \n"
        << "    fig.add_axes([bottom_corner[0],bottom_corner[1],xlen,yl]) \n"
        << "    fig.axes[-1].set_xlim(0,max_rad) \n"
        << "    if i==nrows-1: fig.axes[-1].tick_params(labelbottom=True) \n"
        << "    else: fig.axes[-1].tick_params(labelbottom=False) \n"
        << "    bottom_corner[0]+=xlen+x_sep \n"
        << "  bottom_corner[1]-=(ylen+y_sep) \n\n"
        << "ax = fig.axes \n"
        << std::endl
        << "# Plotting rotation velocity \n"
        << "ax[0].set_ylim(0,1.2*max_vrot) \n"
        << "ax[0].set_ylabel(r'V$_\\mathrm{rot}$ (km/s)', fontsize=lsize) \n"
        << "ax[0].errorbar(rad,vrot, yerr=[-err1_l[0],err1_h[0]],fmt='o', color=color) \n"
        << "if twostage: ax[0].errorbar(rad2,vrot2, yerr=[-err2_l[0],err2_h[0]],fmt='o', color=color2) \n"
        << std::endl
        << "# Plotting velocity dispersion \n"
        << "ax[1].set_ylim(0,1.2*max_vdisp) \n"
        << "ax[1].set_ylabel(r'$\\sigma_\\mathrm{gas}$  (km/s)', fontsize=lsize) \n"
        << "ax[1].errorbar(rad,disp, yerr=[-err1_l[1],err1_h[1]],fmt='o', color=color) \n"
        << "if twostage: ax[1].errorbar(rad2,disp2, yerr=[-err2_l[1],err2_h[1]],fmt='o', color=color2) \n"
        << std::endl
        << "# Plotting surface density \n"
        << "ax[2].set_xlim(0,max_rad) \n"
        << "ax[2].set_ylabel(r'$\\Sigma$ ("<< bunit << ")', fontsize=lsize) \n"
        << "ax[2].errorbar(rad_sd,surfdens, yerr=sd_err,fmt='o', color=color2) \n"
        << std::endl
        << "# Plotting inclination angle \n"
        << "ax[3].set_ylabel('i (deg)', fontsize=lsize) \n"
        << "ax[3].errorbar(rad,inc, yerr=[-err1_l[4],err1_h[4]],fmt='o', color=color) \n"
        << "if twostage: ax[3].errorbar(rad2,inc2,yerr=[-err2_l[4],err2_h[4]], fmt='o-', color=color2) \n"
        << std::endl
        << "# Plotting x-center \n"
        << "ax[4].set_ylabel(r'x$_0$ (pix)', fontsize=lsize) \n"
        << "ax[4].errorbar(rad,xpos, yerr=[-err1_l[6],err1_h[6]],fmt='o', color=color) \n"
        << "if twostage: ax[4].errorbar(rad2,xpos2,yerr=[-err2_l[6],err2_h[6]],fmt='o-', color=color2) \n"
        << std::endl
        << "# Plotting radial velocity \n"
        << "ax[5].set_xlim(0,max_rad) \n"
        << "ax[5].set_ylabel(r'V$_\\mathrm{rad}$ (km/s)', fontsize=lsize) \n"
        << "ax[5].errorbar(rad,vrad, yerr=[-err1_l[9],err1_h[9]],fmt='o', color=color) \n"
        << "if twostage: ax[5].errorbar(rad2,vrad2,yerr=[-err2_l[9],err2_h[9]],fmt='o', color=color2) \n"
        << std::endl
        << "# Plotting position angle \n"
        << "ax[6].set_ylabel(r'$\\phi$ (deg)', fontsize=lsize) \n"
        << "ax[6].set_xlabel('Radius (arcsec)', fontsize=lsize, labelpad=10) \n"
        << "ax[6].errorbar(rad,pa, yerr=[-err1_l[5],err1_h[5]],fmt='o', color=color) \n"
        << "if twostage: ax[6].errorbar(rad2,pa2,yerr=[-err2_l[5],err2_h[5]], fmt='o-', color=color2) \n"
        << std::endl
        << "# Plotting y-center \n"
        << "ax[7].set_ylabel(r'y$_0$ (pix)', fontsize=lsize) \n"
        << "ax[7].set_xlabel('Radius (arcsec)', fontsize=lsize, labelpad=10) \n"
        << "ax[7].errorbar(rad,ypos, yerr=[-err1_l[7],err1_h[7]],fmt='o', color=color) \n"
        << "if twostage: ax[7].errorbar(rad2,ypos2, yerr=[-err2_l[7],err2_h[7]],fmt='o-', color=color2) \n"
        << std::endl
        << "# Plotting systemic velocity \n"
        << "ax[8].set_ylabel(r'v$_\\mathrm{sys}$ (km/s)', fontsize=lsize) \n"
        << "ax[8].set_xlabel('Radius (arcsec)', fontsize=lsize, labelpad=10) \n"
        << "ax[8].errorbar(rad,vsys, yerr=[-err1_l[8],err1_h[8]],fmt='o', color=color) \n"
        << "if twostage: ax[8].errorbar(rad2,vsys2,yerr=[-err2_l[8],err2_h[8]],fmt='o', color=color2) \n"
        << std::endl
        << "fig.savefig(outfolder+'%s_parameters.pdf'%outprefix,bbox_inches='tight') \n";

    if (par.flagADRIFT) {
        pyf << "\n#Asymmetric drift correction\n"
            << "rad_a, vcirc, va2, disp2_a, disp2r, fun, funr = np.genfromtxt(outfolder+'asymdrift.txt',usecols=(0,1,2,3,4,5,6),unpack=True)\n"
            << "if twostage: \n"
            << "  rad, vrot, disp = rad2, vrot2, disp2 \n"
            << "  err1_l, err1_h = err2_l, err2_h\n"
            << "fig=plt.figure(figsize=(10,10), dpi=150) \n"
            << "ax1 = fig.add_axes([0.10,0.1,0.3,0.25])\n"
            << "ax2 = fig.add_axes([0.48,0.1,0.3,0.25])\n"
            << "ax3 = fig.add_axes([0.86,0.1,0.3,0.25])\n"
            << "ax1.set_xlim(0,max_rad)\n"
            << "ax1.set_xlabel('Radius (arcsec)', fontsize=11, labelpad=5)\n"
            << "ax1.set_ylabel('V (km/s)', fontsize=11)\n"
            << "ax1.plot(rad,vrot,'-', color=color2,label=r'V$_\\mathrm{rot}$',zorder=1)\n"
            << "ax1.plot(rad_a,np.sqrt(va2),'--', color='#FABC11',label=r'V$_\\mathrm{A}$',zorder=2)\n"
            << "ax1.errorbar(rad_a,vcirc, yerr=[-err1_l[0],err1_h[0]],fmt='o', color='#43A8D4',label=r'V$_\\mathrm{circ}$',zorder=0)\n"
            << "ax1.legend()\n"
            << "ax2.set_xlim(0,max_rad)\n"
            << "ax2.set_xlabel('Radius (arcsec)', fontsize=11, labelpad=5)\n"
            << "ax2.set_ylabel(r'$\\sigma$ (km/s)', fontsize=11)\n"
            << "ax2.errorbar(rad,disp, yerr=[-err1_l[1],err1_h[1]],fmt='o', color=color2,label=r'$\\sigma_\\mathrm{gas}$',zorder=0)\n"
            << "ax2.plot(rad_a,np.sqrt(disp2r),'-', color='#0FA45A',label=r'$\\sigma_\\mathrm{reg}$',zorder=1)\n"
            << "ax2.legend()\n"
            << "ax3.set_xlim(0,max_rad)\n"
            << "ax3.set_xlabel('Radius (arcsec)', fontsize=11, labelpad=5)\n"
            << "ax3.set_ylabel(r'$\\Sigma_\\mathrm{gas}\\cos(i)$', fontsize=11)\n"
            << "ax3.plot(rad_a,fun,'o', color=color2,label=r'$f_\\mathrm{obs}$',zorder=0)\n"
            << "ax3.plot(rad_a,funr,'-', color='#0FA45A',label=r'$f_\\mathrm{reg}$',zorder=1)\n"
            << "ax3.legend()\n"
            << "fig.savefig(outfolder+'%s_asymmetricdrift.pdf'%outprefix,bbox_inches='tight')\n";
    }

    pyf.close();

    /////////////////////////////////////////////////////////////////////////
    /// Script to plot channel maps of model vs data
    /////////////////////////////////////////////////////////////////////////
    scriptnames.push_back("plot_chanmaps.py");
    pyf.open((in->pars().getOutfolder()+"plotscripts/"+scriptnames.back()).c_str());

    std::string cubefile = in->pars().getImageFile();
    if (in->pars().getFlatContsub()) cubefile = in->pars().getOutfolder()+in->pars().getOutPrefix()+"_contsub.fits";

    pyf << "#####################################################################\n"
        << "#### This script writes a plot of channel maps of model and data ####\n"
        << "#####################################################################\n"
        << "import numpy as np \n"
        << "import os \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "import matplotlib.gridspec as gridspec \n"
        << "from plot_utils import calculate_padding \n"
        << "from astropy.io import fits \n"
        << "from astropy.visualization import PowerStretch \n"
        << "from astropy.visualization.mpl_normalize import ImageNormalize \n"
        << "from astropy.visualization import PercentileInterval \n"
        << "mpl.rc('xtick',direction='in',top=True,labelbottom=False) \n"
        << "mpl.rc('ytick',direction='in',right=True,labelleft=False) \n"
        << "mpl.rcParams['contour.negative_linestyle'] = 'solid' \n"
        << "plt.rc('font',family='sans-serif',serif='Helvetica',size=10) \n"
        << "params = {'text.usetex': False, 'mathtext.fontset': 'cm', 'mathtext.default': 'regular'} \n"
        << "plt.rcParams.update(params) \n"
        << std::endl
        << "gname = '" << in->Head().Name() <<"' \n"
        << "outfolder = '" << in->pars().getOutfolder() <<"' \n"
        << "outprefix = '" << in->pars().getOutPrefix() <<"' \n"
        << "twostage = " << par.TWOSTAGE << " \n"
        << "plotmask = " << par.PLOTMASK << " \n"
        << std::endl
        << "if twostage: xpos,ypos,vsys = np.genfromtxt(outfolder+'rings_final2.txt',usecols=(9,10,11),unpack=True) \n"
        << "else: xpos,ypos,vsys = np.genfromtxt(outfolder+'rings_final1.txt',usecols=(9,10,11),unpack=True) \n"
        << "xcen_m,ycen_m,vsys_m = np.nanmean((xpos,ypos,vsys),axis=1) \n"
        << std::endl
        << "# CHANNEL MAPS: Setting all the needed variables \n"
        << "image = fits.open('" << cubefile << "') \n"
        << "image_mas = fits.open(outfolder+'mask.fits') \n"
        << "xmin, xmax = " << xmin << ", " << xmax << std::endl
        << "ymin, ymax = " << ymin << ", " << ymax << std::endl
        << "zmin, zmax = " << zmin << ", " << zmax << std::endl
        << "data = image[0].data[";
    if (in->Head().NumAx()>3)
        for (int i=0; i<in->Head().NumAx()-3; i++) pyf << "0,";
    pyf << "zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
        << "data_mas = image_mas[0].data[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
        << "head = image[0].header \n"
        << "zsize, ysize, xsize = data.shape \n"
        << "pad_x, pad_y = calculate_padding(xsize,ysize) \n"
        << "cdeltsp=" << in->Head().PixScale()*arcconv << std::endl
        << "cont = " << cont << std::endl
        << "v = np.array([1,2,4,8,16,32,64])*cont \n"
        << "v_neg = [-cont] \n"
        << "interval = PercentileInterval(99.5) \n"
        << "vmax = interval.get_limits(data)[1] \n"
        << "norm = ImageNormalize(vmin=cont, vmax=vmax, stretch=PowerStretch(0.5)) \n"
        << "xcen, ycen = xcen_m-xmin+pad_x[0], ycen_m-ymin+pad_y[0] \n"
        << std::endl
        << "files_mod = [f for f in sorted(os.listdir(outfolder)) if outprefix+'mod' in f]  \n"
        << "if len(files_mod)==0: raise FileNotFoundError('ERROR: no model in output directory') \n\n";

    pyf << "# Beginning channel map plot \n"
        << "for k in range (len(files_mod)): \n"
        << "  image_mod = fits.open(outfolder+files_mod[k]) \n"
        << "  data_mod = image_mod[0].data[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
        << "  fig = plt.figure(figsize=(8.27, 11.69), dpi=150) \n"
        << "  grid = [gridspec.GridSpec(2,5),gridspec.GridSpec(2,5),gridspec.GridSpec(2,5)] \n"
        << "  grid[0].update(top=0.90, bottom=0.645, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
        << "  grid[1].update(top=0.60, bottom=0.345, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
        << "  grid[2].update(top=0.30, bottom=0.045, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
        << std::endl
        << "  num = 0 \n"
        << "  for j in range (0,3): \n"
        << "    for i in range (0,5): \n"
        << "      chan = int(num*(zsize)/15) \n"
        << "      z = np.pad(data[chan,:,:],(pad_y, pad_x),mode='constant',constant_values=np.nan) \n"
        << "      z_mod = np.pad(data_mod[chan,:,:],(pad_y, pad_x),mode='constant',constant_values=np.nan) \n"
        << "      #New matplotlib draws wrong contours when no contours are found. This is a workaround.\n"
        << "      if np.all(z_mod<v[0]): z_mod[:,:] =0\n"
        << "      velo_kms = (chan+1-" << 1-zmin << ")*" << cdelt3_kms << "+" << crval3_kms_pix1 << std::endl
        << "      velo = ' v = ' + str(int(velo_kms-vsys_m)) + ' km/s' \n"
        << "      ax = plt.subplot(grid[j][0,i]) \n"
        << "      ax.set_title(velo, fontsize=10,loc='left') \n"
        << "      ax.imshow(z,origin='lower',cmap = mpl.cm.Greys,norm=norm,aspect='auto',interpolation='none') \n"
        << "      ax.contour(z,v,origin='lower',linewidths=0.7,colors='#00008B') \n"
        << "      ax.contour(z,v_neg,origin='lower',linewidths=0.1,colors='gray') \n"
        << "      ax.plot(xcen,ycen,'x',color='#0FB05A',markersize=7,mew=2) \n"
        << "      if plotmask: \n"
        << "        ax.contour(data_mas[chan],levels=[0],origin='lower',linewidths=2,colors='k') \n"
        << "      if (j==i==0): \n"
        << "        ax.text(0, 1.4, gname, transform=ax.transAxes,fontsize=15,va='center') \n"
        << "        lbar = 0.5*(xmax-xmin)*cdeltsp \n"
        << "        ltex = \"%.0f'' \"%lbar if lbar>10 else \"%.2f'' \"%lbar \n"
        << "        if lbar>600: ltex = \"%.0f' \"%(lbar/60.) \n"
        << "        ax.annotate('', xy=(4.5, 1.4), xycoords='axes fraction', xytext=(5, 1.4),arrowprops=dict(arrowstyle='<->', color='k'))\n"
        << "        ax.text(4.75,1.50,ltex,transform=ax.transAxes,fontsize=11, ha='center')\n"
        << "        bmaj, bmin, bpa = " << bmaj << "/float(xmax-xmin), " << bmin << "/float(ymax-ymin)," << bpa << std::endl
        << "        beam = mpl.patches.Ellipse(xy=(3.5, 1.4), width=bmaj, height=bmin, angle=bpa+90, \\\n"
        << "                                   color='#5605D0', clip_on=False, transform=ax.transAxes, alpha=0.2) \n"
        << "        ax.add_artist(beam) \n"
        << "        ax.text(3.6+bmaj/1.8,1.4,'Beam',transform=ax.transAxes,fontsize=11, ha='left',va='center') \n"
        << "      ax = plt.subplot(grid[j][1,i]) \n"
        << "      ax.tick_params(axis='both',which='both',bottom=True,top=True,labelbottom=False,labelleft=False) \n"
        << "      ax.imshow(z_mod,origin='lower',cmap = mpl.cm.Greys,norm=norm,aspect='auto',interpolation='none') \n"
        << "      ax.contour(z_mod,v,origin='lower',linewidths=0.7,colors='#B22222') \n"
        << "      ax.plot(xcen,ycen,'x',color='#0FB05A',markersize=7,mew=2) \n"
        << "      if (i==0 and j==2): \n"
        << "        clab = r'Contour levels at 2$^n \\, c_{min}$, where $c_{min}$ = %s " << in->Head().Bunit() << " and n = 0,1,..,8 '%cont \n"
        << "        ax.text(0.01,-0.16,clab,transform=ax.transAxes,fontsize=11, ha='left',va='center') \n"
        << "      num = num+1 \n"
        << std::endl
        << "  outfile = '%s_chanmaps'%outprefix \n"
        << "  ntype   = files_mod[k][files_mod[k].rfind('_'):].replace('.fits','') \n"
        << "  fig.savefig(outfolder+outfile+ntype+'.pdf', orientation = 'portrait', format = 'pdf') \n"
        << "  image_mod.close() \n\n"
        << "image.close() \n";

    pyf.close();

    /////////////////////////////////////////////////////////////////////////
    /// Script to plot position-velocity slices along minor/major axis
    /////////////////////////////////////////////////////////////////////////
    float zmin_wcs = AlltoVel(in->getZphys(zmin-0.5),in->Head());
    float zmax_wcs = AlltoVel(in->getZphys(zmax+0.5),in->Head());
    float pa_av = findMedian(&outr->phi[0],outr->nr);
    float pa_min = pa_av+90<360 ? pa_av+90 : pa_av-90;
    //if (zmin_wcs>zmax_wcs) std::swap(zmin_wcs,zmax_wcs);
    bool reverse = (pa_av>=45 && pa_av<225);
    //if (cdelt3_kms<0) reverse = !reverse;

    // OLD VERSION
    //scriptnames.push_back("plot_pvs.py");
    pyf.open((in->pars().getOutfolder()+"plotscripts/plot_pvs_old.py").c_str());

    pyf << "#################################################################################\n"
        << "#### This script writes a plot of position-velocity slices of model and data ####\n"
        << "#################################################################################\n"
        << "import numpy as np \n"
        << "import os \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "from plot_utils import defineaxis \n"
        << "from astropy.io import fits \n"
        << "from astropy.visualization import PowerStretch \n"
        << "from astropy.visualization.mpl_normalize import ImageNormalize \n"
        << "from astropy.visualization import PercentileInterval \n"
        << "mpl.rc('xtick',direction='in') \n"
        << "mpl.rc('ytick',direction='in') \n"
        << "mpl.rcParams['contour.negative_linestyle'] = 'solid' \n"
        << "plt.rc('font',family='sans-serif',serif='Helvetica',size=10) \n"
        << "params = {'text.usetex': False, 'mathtext.fontset': 'cm', 'mathtext.default': 'regular'} \n"
        << "plt.rcParams.update(params) \n"
        << std::endl
        << "gname = '" << in->Head().Name() <<"' \n"
        << "outfolder = '" << in->pars().getOutfolder() <<"' \n"
        << "outprefix = '" << in->pars().getOutPrefix() <<"' \n"
        << "twostage = " << par.TWOSTAGE << " \n"
        << "plotmask = " << par.PLOTMASK << " \n"
        << "zmin, zmax = " << zmin << ", " << zmax << std::endl
        << std::endl
        << "if twostage: rad,vrot,inc,pa,vsys = np.genfromtxt(outfolder+'rings_final2.txt',usecols=(1,2,4,5,11),unpack=True) \n"
        << "else: rad,vrot,inc,pa,vsys = np.genfromtxt(outfolder+'rings_final1.txt',usecols=(1,2,4,5,11),unpack=True) \n"
        << std::endl
        << "files_pva_mod = [f for f in sorted(os.listdir(outfolder+'pvs/')) if outprefix+'mod_pv_a' in f] \n"
        << "files_pvb_mod = [f for f in sorted(os.listdir(outfolder+'pvs/')) if outprefix+'mod_pv_b' in f] \n"
        << "if len(files_pva_mod)==0 or len(files_pvb_mod)==0: raise FileNotFoundError('ERROR: no PV model in output directory') \n\n"
        << std::endl
        << "image_maj     = fits.open(outfolder+'pvs/'+outprefix+'_pv_a.fits') \n"
        << "image_min     = fits.open(outfolder+'pvs/'+outprefix+'_pv_b.fits') \n"
        << "image_mas_maj = fits.open(outfolder+'pvs/'+outprefix+'mask_pv_a.fits') \n"
        << "image_mas_min = fits.open(outfolder+'pvs/'+outprefix+'mask_pv_b.fits') \n"
        << "head = [image_maj[0].header,image_min[0].header] \n"
        << "crpixpv = np.array([head[0]['CRPIX1'],head[1]['CRPIX1']]) \n"
        << "cdeltpv = np.array([head[0]['CDELT1'],head[1]['CDELT1']]) \n"
        << "crvalpv = np.array([head[0]['CRVAL1'],head[1]['CRVAL1']]) \n"
        << "xminpv, xmaxpv = np.floor(crpixpv-1-" << disp << "), np.ceil(crpixpv-1 +"<< disp << ") \n"
        << "if xminpv[0]<0: xminpv[0]=0 \n"
        << "if xminpv[1]<0: xminpv[1]=0 \n"
        << "if xmaxpv[0]>=head[0]['NAXIS1']: xmaxpv[0]=head[0]['NAXIS1']-1 \n"
        << "if xmaxpv[1]>=head[1]['NAXIS1']: xmaxpv[1]=head[1]['NAXIS1']-1 \n"
        << "data_maj = image_maj[0].data[zmin:zmax+1,int(xminpv[0]):int(xmaxpv[0])+1] \n"
        << "data_min = image_min[0].data[zmin:zmax+1,int(xminpv[1]):int(xmaxpv[1])+1] \n"
        << "data_mas_maj = image_mas_maj[0].data[zmin:zmax+1,int(xminpv[0]):int(xmaxpv[0])+1] \n"
        << "data_mas_min = image_mas_min[0].data[zmin:zmax+1,int(xminpv[1]):int(xmaxpv[1])+1] \n"
        << "xmin_wcs = ((xminpv+1-0.5-crpixpv)*cdeltpv+crvalpv)*" << arcconv << std::endl
        << "xmax_wcs = ((xmaxpv+1+0.5-crpixpv)*cdeltpv+crvalpv)*" << arcconv << std::endl
        << "zmin_wcs, zmax_wcs = " << zmin_wcs << ", " << zmax_wcs << std::endl
        << "cont = " << cont << std::endl
        << "v = np.array([1,2,4,8,16,32,64])*cont \n"
        << "v_neg = [-cont] \n"
        << "interval = PercentileInterval(99.5) \n"
        << "vmax = interval.get_limits(data_maj)[1] \n"
        << "norm = ImageNormalize(vmin=cont, vmax=vmax, stretch=PowerStretch(0.5)) \n\n";

    pyf << "radius = np.concatenate((rad,-rad)) \n"
        << "pamaj_av = "<< pa_av << endl
        << "pamin_av = "<< pa_min << endl
        << "costh = np.cos(np.deg2rad(np.abs(pa-pamaj_av))) \n"
        << "vlos1 = vsys+vrot*np.sin(np.deg2rad(inc))*costh \n"
        << "vlos2 = vsys-vrot*np.sin(np.deg2rad(inc))*costh \n";
    if (reverse) pyf << "reverse = True \n";
    else pyf << "reverse = False \n";
    pyf << "if reverse: vlos1, vlos2 = vlos2, vlos1 \n"
        << "vlos = np.concatenate((vlos1,vlos2)) \n"
        << "vsys_m = np.nanmean(vsys) \n"
        << "ext = [[xmin_wcs[0],xmax_wcs[0],zmin_wcs-vsys_m,zmax_wcs-vsys_m],\\" << std::endl
        << "       [xmin_wcs[1],xmax_wcs[1],zmin_wcs-vsys_m,zmax_wcs-vsys_m]] \n"
        << "labsize = 14 \n"
        << "palab = [r'$\\phi = $%i$^\\circ$'%np.round(pamaj_av), r'$\\phi = $%i$^\\circ$'%np.round(pamin_av)] \n\n";

    pyf << "# Beginning PV plot \n"
        << "for k in range (len(files_pva_mod)): \n"
        << "  image_mod_maj = fits.open(outfolder+'pvs/'+files_pva_mod[k]) \n"
        << "  image_mod_min = fits.open(outfolder+'pvs/'+files_pvb_mod[k]) \n"
        << "  data_mod_maj = image_mod_maj[0].data[zmin:zmax+1,int(xminpv[0]):int(xmaxpv[0])+1] \n"
        << "  data_mod_min = image_mod_min[0].data[zmin:zmax+1,int(xminpv[1]):int(xmaxpv[1])+1] \n"
        << "  toplot = [[data_maj,data_min],[data_mod_maj,data_mod_min],[data_mas_maj,data_mas_min]] \n"
        << std::endl
        << "  fig, ax = defineaxis(2,1,0.6,0.42,xsep=0.00,ysep=0.08,fig_width=10,fig_heigth=10) \n\n"
        << "  for i in range (2): \n"
        << "    axis = ax[i] \n"
        << "    axis.tick_params(which='major',length=8, labelsize=labsize) \n"
        << "    axis.set_xlabel('Offset (arcsec)',fontsize=labsize+2) \n"
        << "    axis.set_ylabel(r'$\\mathrm{\\Delta V_{LOS}}$ (km/s)',fontsize=labsize+2) \n"
        << "    axis.text(1, 1.02,palab[i],ha='right',transform=axis.transAxes,fontsize=labsize+4) \n"
        << "    axis2 = axis.twinx() \n"
        << "    axis2.set_xlim([ext[i][0],ext[i][1]]) \n"
        << "    axis2.set_ylim([ext[i][2]+vsys_m,ext[i][3]+vsys_m]) \n"
        << "    axis2.tick_params(which='major',length=8, labelsize=labsize) \n"
        << "    axis2.set_ylabel(r'$\\mathrm{V_{LOS}}$ (km/s)',fontsize=labsize+2) \n"
        << "    axis.imshow(toplot[0][i],origin='lower',cmap = mpl.cm.Greys,norm=norm,extent=ext[i],aspect='auto') \n"
        << "    axis.contour(toplot[0][i],v,origin='lower',linewidths=0.7,colors='#00008B',extent=ext[i]) \n"
        << "    axis.contour(toplot[0][i],v_neg,origin='lower',linewidths=0.1,colors='gray',extent=ext[i]) \n"
        << "    axis.contour(toplot[1][i],v,origin='lower',linewidths=1,colors='#B22222',extent=ext[i]) \n"
        << "    axis.axhline(y=0,color='black') \n"
        << "    axis.axvline(x=0,color='black') \n"
        << "    axis.grid(color='gray', linestyle='--', linewidth=0.3) \n"
        << "    if plotmask: \n"
        << "      axis.contour(toplot[2][i],levels=[0],origin='lower',linewidths=2,colors='k',extent=ext[i]) \n"
        << "    if i==0 : \n"
        << "      axis2.plot(radius,vlos,'yo') \n"
        << "      axis.text(0, 1.1, gname, transform=axis.transAxes,fontsize=22) \n"
        << std::endl
        << "  outfile = '%s_pv'%outprefix \n"
        << "  ntype   = files_pva_mod[k][files_pva_mod[k].rfind('_'):].replace('.fits','') \n"
        << "  fig.savefig(outfolder+outfile+ntype+'.pdf', bbox_inches='tight') \n"
        << "  image_mod_maj.close() \n"
        << "  image_mod_min.close() \n"
        << std::endl
        << "image_maj.close() \n"
        << "image_min.close() \n";

    pyf.close();
    
    // NEW VERSION
    scriptnames.push_back("plot_pvs.py");
    pyf.open((in->pars().getOutfolder()+"plotscripts/"+scriptnames.back()).c_str());

    pyf << "#################################################################################\n"
        << "#### This script writes a plot of position-velocity slices of model and data ####\n"
        << "#################################################################################\n"
        << "import numpy as np \n"
        << "import os \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "from plot_utils import *\n"
        << "from astropy.io import fits \n"
        << "from astropy.visualization import ImageNormalize, PercentileInterval, PowerStretch\n"
        << std::endl
        << "labsize = 18\n"
        << "mpl.rc('xtick',direction='in',top=True) \n"
        << "mpl.rc('ytick',direction='in',right=True) \n"
        << "mpl.rcParams['contour.negative_linestyle'] = 'solid' \n"
        << "plt.rc('font',family='sans-serif',serif='Helvetica',size=labsize) \n"
        << "params = {'text.usetex': False, 'mathtext.fontset': 'cm', 'mathtext.default': 'regular'} \n"
        << "plt.rcParams.update(params) \n"
        << std::endl
        << "gname = '" << in->Head().Name() <<"' \n"
        << "outfolder = '" << in->pars().getOutfolder() <<"' \n"
        << "outprefix = '" << in->pars().getOutPrefix() <<"' \n"
        << "twostage = " << par.TWOSTAGE << " \n"
        << "plotmask = " << par.PLOTMASK << " \n"
        << std::endl
        << "image = fits.open('" << cubefile << "') \n"
        << "image_mas = fits.open(outfolder+'mask.fits') \n"
        << "xmin, xmax = " << xmin << ", " << xmax << std::endl
        << "ymin, ymax = " << ymin << ", " << ymax << std::endl
        << "zmin, zmax = " << zmin << ", " << zmax << std::endl
        << "zmin_wcs, zmax_wcs = " << zmin_wcs << ", " << zmax_wcs << std::endl
        << "data = image[0].data[";
    if (in->Head().NumAx()>3)
        for (int i=0; i<in->Head().NumAx()-3; i++) pyf << "0,";
    pyf << "zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
        << "data_mas = image_mas[0].data[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
        << "zsize, ysize, xsize = data.shape \n"
        << "pad_x, pad_y = calculate_padding(xsize,ysize) \n"
        << "cdeltsp=" << in->Head().PixScale()*arcconv << std::endl
        << "cont = " << cont << std::endl
        << "v = np.array([1,2,4,8,16,32,64])*cont \n"
        << "v_neg = [-cont] \n"
        << "mom0 = fits.open(outfolder+'/maps/'+outprefix+'_0mom.fits')[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "mom1 = fits.open(outfolder+'/maps/'+outprefix+'_1mom.fits')[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "mom0p = np.pad(mom0,(pad_y, pad_x),mode='constant',constant_values=np.nan) \n"
        << "mom1p = np.pad(mom1,(pad_y, pad_x),mode='constant',constant_values=np.nan) \n"
        << "extm = [0,mom0p.shape[1]-1,0,mom0p.shape[0]-1]\n"
        << "norm0 = ImageNormalize(mom0, interval=PercentileInterval(98.0), stretch=PowerStretch(0.5))\n"
        << "norm1 = ImageNormalize(mom1, interval=PercentileInterval(99.5))\n"
        << std::endl
        << "frile = outfolder + ('rings_final2.txt' if twostage else 'rings_final1.txt')\n"
        << "rings = np.genfromtxt(frile,usecols=(1,2,3,4,5,9,10,11),unpack=True)\n"
        << "rad,vrot,vdisp,inc,pa,xpos,ypos,vsys = rings\n"
        << "xcen_m,ycen_m,vsys_m,inc_m,pa_m = np.nanmean((xpos,ypos,vsys,inc,pa),axis=1) \n"
        << "xcen_m, ycen_m = xcen_m-xmin, ycen_m-ymin \n"
        << "proj_vmax = np.nanmax(vrot)*np.sin(np.radians(inc_m))\n"
        << "max_vdisp = np.nanmax(vdisp)\n"
        << "radii = np.concatenate((rad,-rad)) \n"
        << "vlos1 = vrot*np.sin(np.deg2rad(inc))\n"
        << "vlos = np.concatenate((vlos1,-vlos1)) \n\n"
        << "files_mod = [f for f in sorted(os.listdir(outfolder)) if outprefix+'mod' in f] \n"
        << "if len(files_mod)==0: raise FileNotFoundError('ERROR: no model in output directory') \n\n"
        << std::endl
        << "def plot_pv(ax,x0,y0,theta,**kwargs):\n\n"
        << "  s,x,y = build_path(x0, y0, rad/cdeltsp, theta, xsize, ysize)\n"
        << "  pv_data = extract_pv(data,x,y,order=0)\n"
        << "  pv_mod  = extract_pv(data_mod,x,y)\n"
        << "  norm = ImageNormalize(vmin=cont, vmax=np.percentile(data, 99.8), stretch=PowerStretch(0.5)) \n"
        << "  #norm = ImageNormalize(data, interval=PercentileInterval(99.9), stretch=PowerStretch(1.0))\n\n"
        << "  s_arcsec = s * cdeltsp\n"
        << "  ext = [s_arcsec[0], s_arcsec[-1], zmin_wcs-vsys_m,zmax_wcs-vsys_m]\n\n"
        << "  ax.imshow(pv_data,origin='lower',cmap='Greys',norm=norm,aspect='auto',interpolation='nearest',extent=ext)\n"
        << "  ax.contour(pv_data,v,origin='lower',linewidths=1.0,colors='#00008B',extent=ext) \n"
        << "  ax.contour(pv_data,v_neg,origin='lower',linewidths=0.2,colors='gray',extent=ext) \n"
        << "  ax.contour(pv_mod,v,origin='lower',linewidths=1.5,colors='#B22222',extent=ext) \n"
        << "  ax.axhline(y=0,color='black') \n"
        << "  ax.axvline(x=0,color='black') \n"
        << "  ax.grid(color='gray', linestyle='--', linewidth=0.3) \n\n"
        << "  if plotmask:\n"
        << "    pv_mask = extract_pv(data_mas,x,y) \n"
        << "    ax.contour(pv_mask,levels=[0],origin='lower',linewidths=3.5,colors='k',extent=ext)\n\n"
        << "  ax.set_xlim(-1.3*rad[-1],1.3*rad[-1])\n"
        << "  #ax.set_ylim(-1.1*(proj_vmax+2*max_vdisp),1.1*(proj_vmax+2*max_vdisp))\n"
        << "  ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=6))\n"
        << "  ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=6))\n"
        << "  ax.tick_params(labelbottom=True,labelleft=True)\n"
        << "  ax.tick_params(which='major',length=8, labelsize=labsize) \n"
        << "  ax.set_ylabel(r'$\\mathrm{\\Delta V_{LOS}}$ (km/s)',fontsize=labsize+2) \n"
        << "  ax.set_xlabel('Offset (arcsec)',fontsize=labsize+2)\n\n"
        << "  ax2 = ax.twinx()\n"
        << "  ax2.set_ylim([ext[2]+vsys_m,ext[3]+vsys_m]) \n"
        << "  ax2.tick_params(which='major',length=8, labelsize=labsize) \n"
        << "  ax2.set_ylabel(r'$\\mathrm{V_{LOS}}$ (km/s)',fontsize=labsize+2) \n\n"
        << "  return s,x,y\n\n"
        << std::endl
        << "# Beginning position-velocity plot\n"
        << "for k in range (len(files_mod)): \n"
        << "  image_mod = fits.open(outfolder+files_mod[k]) \n"
        << "  data_mod = image_mod[0].data[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1]\n\n"
        << "  fig, axes = defineaxis(2,1,1,0.7,xsep=0,ysep=0.1,fig_width=9,fig_heigth=9)\n"
        << "  axes = np.ravel(axes)\n\n"
        << "  # Major axis\n"
        << "  s,x,y = plot_pv(axes[0],xcen_m,ycen_m,pa)\n"
        << "  axes[0].text(0, 1.05, gname, transform=axes[0].transAxes,fontsize=32) \n\n"
        << "  fig.canvas.draw()\n"
        << "  renderer = fig.canvas.get_renderer() \n"
        << "  bbox = axes[0].get_tightbbox(renderer).transformed(fig.transFigure.inverted())\n"
        << "  pos = axes[0].get_position()\n"
        << "  ax_m = fig.add_axes([bbox.x1+0.14,pos.y1-0.4,0.4,0.4])\n\n"
        << "  ax_m.tick_params(labelbottom=False,labelleft=False,top=False,bottom=False,left=False,right=False) \n"
        << "  ax_m.imshow(mom0p,origin='lower',cmap='Greys',aspect='auto',interpolation='nearest',norm=norm0,extent=extm)\n"
        << "  ax_m.imshow(mom1p,origin='lower',cmap='seismic',aspect='auto',interpolation='nearest',norm=norm1,extent=extm, alpha=0.4)\n"
        << "  ax_m.plot(xcen_m+pad_x[0],ycen_m+pad_y[0],'x',color='y',markersize=7,mew=1.5,zorder=10)\n"
        << "  ax_m.annotate('',xy=(0.8,-0.1),xycoords='axes fraction',xytext=(1.0, -0.1),arrowprops=dict(arrowstyle='<->', color='k', lw=1.5))\n"
        << "  ax_m.text(0.9,-0.16,f\"{0.2*xsize*cdeltsp:.0f}''\",transform=ax_m.transAxes,fontsize=12, ha='center')\n"
        << "  ax_m.axvline(xcen_m+pad_x[0],ls='dotted',c='k',alpha=0.2)\n"
        << "  ax_m.axhline(ycen_m+pad_y[1],ls='dotted',c='k',alpha=0.2)\n"
        << "  _,ax_m_xmax = ax_m.get_xlim()\n"
        << "  _,ax_m_ymax = ax_m.get_ylim()\n"
        << "  m = (x >= 0) & (x < ax_m_xmax) & (y >= 0) & (y < ax_m_ymax)\n"
        << "  ax_m.plot(x[m]+pad_x[0],y[m]+pad_y[0],'-',lw=1.5,c='#B22222') \n"
        << "  ax_m.plot(x[~m]+pad_x[0],y[~m]+pad_y[0],'-',lw=1.5,c='#B22222',alpha=0.3) \n\n"
        << "  # Plotting rotation curve \n"
        << "  r = np.sqrt((x-xcen_m)**2+(y-ycen_m)**2)\n"
        << "  r[s<0] *= -1\n"
        << "  x_rad = np.interp(radii/cdeltsp,r,x)\n"
        << "  y_rad = np.interp(radii/cdeltsp,r,y)\n"
        << "  s_rad = np.interp(radii/cdeltsp,r,s)\n"
        << "  axes[0].plot(s_rad*cdeltsp,vlos,'o',c='gold',ms=10,alpha=0.8,mec='cornsilk') \n"
        << "  ax_m.scatter(x_rad+pad_x[0],y_rad+pad_y[0],c='gold',s=10,alpha=0.8) \n\n"
        << "  # Minor axis\n"
        << "  s,x,y = plot_pv(axes[1],xcen_m,ycen_m,pa+90)\n"
        << "  m = (x >= 0) & (x < ax_m_xmax) & (y >= 0) & (y < ax_m_ymax)\n"
        << "  ax_m.plot(x[m]+pad_x[0],y[m]+pad_y[0],'--',lw=1.5,c='navy') \n"
        << "  ax_m.plot(x[~m]+pad_x[0],y[~m]+pad_y[0],'--',lw=1.5,c='navy',alpha=0.3) \n\n"
        << "  xran_min = rad[-1]*np.sin(np.radians(inc_m))\n"
        << "  axes[1].set_xlim(-1.3*xran_min,1.3*xran_min)\n\n"
        << "  outfile = '%s_pv'%outprefix\n"
        << "  ntype   = files_mod[k][files_mod[k].rfind('_'):].replace('.fits','') \n"
        << "  fig.savefig(outfolder+outfile+ntype+'.pdf', bbox_inches='tight') \n"
        << "  image_mod.close() \n\n"
        << "image.close() \n"
        << "image_mas.close() \n";

    pyf.close();

    /////////////////////////////////////////////////////////////////////////
    /// Script to plot position-velocity slices along several cuts
    /////////////////////////////////////////////////////////////////////////
    scriptnames.push_back("plot_pvslices.py");
    pyf.open((in->pars().getOutfolder()+"plotscripts/"+scriptnames.back()).c_str());

    pyf << "#######################################################################\n"
        << "#### This script writes a plot of P-V slices across several cuts   ####\n"
        << "#######################################################################\n"
        << "import numpy as np \n"
        << "import os \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "from plot_utils import *\n"
        << "from astropy.io import fits \n"
        << "from astropy.visualization import ImageNormalize, PercentileInterval, PowerStretch\n"
        << std::endl
        << "labsize = 13\n"
        << "mpl.rc('xtick',direction='in',top=True) \n"
        << "mpl.rc('ytick',direction='in',right=True) \n"
        << "mpl.rcParams['contour.negative_linestyle'] = 'solid' \n"
        << "plt.rc('font',family='sans-serif',serif='Helvetica',size=labsize) \n"
        << "params = {'text.usetex': False, 'mathtext.fontset': 'cm', 'mathtext.default': 'regular'} \n"
        << "plt.rcParams.update(params) \n"
        << std::endl
        << "gname = '" << in->Head().Name() <<"' \n"
        << "outfolder = '" << in->pars().getOutfolder() <<"' \n"
        << "outprefix = '" << in->pars().getOutPrefix() <<"' \n"
        << "twostage = " << par.TWOSTAGE << " \n"
        << "plotmask = " << par.PLOTMASK << " \n"
        << std::endl
        << "image = fits.open('" << cubefile << "') \n"
        << "image_mas = fits.open(outfolder+'mask.fits') \n"
        << "xmin, xmax = " << xmin << ", " << xmax << std::endl
        << "ymin, ymax = " << ymin << ", " << ymax << std::endl
        << "zmin, zmax = " << zmin << ", " << zmax << std::endl
        << "zmin_wcs, zmax_wcs = " << zmin_wcs << ", " << zmax_wcs << std::endl
        << "data = image[0].data[";
    if (in->Head().NumAx()>3)
        for (int i=0; i<in->Head().NumAx()-3; i++) pyf << "0,";
    pyf << "zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
        << "data_mas = image_mas[0].data[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
        << "zsize, ysize, xsize = data.shape \n"
        << "pad_x, pad_y = calculate_padding(xsize,ysize) \n"
        << "cdeltsp=" << in->Head().PixScale()*arcconv << std::endl
        << "cont = " << cont << std::endl
        << "v = np.array([1,2,4,8,16,32,64])*cont \n"
        << "v_neg = [-cont] \n"
        << "mom0 = fits.open(outfolder+'/maps/'+outprefix+'_0mom.fits')[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "mom0p = np.pad(mom0,(pad_y, pad_x),mode='constant',constant_values=np.nan) \n"
        << "extm = [0,mom0p.shape[1]-1,0,mom0p.shape[0]-1]\n"
        << "norm = ImageNormalize(data, interval=PercentileInterval(99.9), stretch=PowerStretch(1.0))\n"
        << "norm0 = ImageNormalize(mom0, interval=PercentileInterval(99.0), stretch=PowerStretch(0.5))\n"
        << std::endl
        << "frile = outfolder + ('rings_final2.txt' if twostage else 'rings_final1.txt')\n"
        << "rings = np.genfromtxt(frile,usecols=(1,2,3,4,5,9,10,11),unpack=True)\n"
        << "rad,vrot,vdisp,inc,pa,xpos,ypos,vsys = rings\n"
        << "xcen_m,ycen_m,vsys_m,inc_m,pa_m = np.nanmean((xpos,ypos,vsys,inc,pa),axis=1) \n"
        << "xcen_m,ycen_m = xcen_m-xmin, ycen_m-ymin \n"
        << "proj_vmax = np.nanmax(vrot)*np.sin(np.radians(inc_m))\n"
        << "max_vdisp = np.nanmax(vdisp)\n"
        << std::endl
        << "s_maj,x_maj,y_maj = build_path(xcen_m, ycen_m, rad/cdeltsp, np.full(len(rad),pa_m), xsize, ysize)\n"
        << "s_min,x_min,y_min = build_path(xcen_m, ycen_m, rad/cdeltsp, np.full(len(rad),pa_m+90), xsize, ysize)\n"
        << "s_maj = s_maj*cdeltsp\n"
        << "s_min = s_min*cdeltsp\n"
        << "maj_offsets = np.linspace(0,rad[-1],3,endpoint=True)\n"
        << "min_offsets = np.linspace(0,rad[-1]*np.cos(np.radians(inc_m)),3,endpoint=True)\n"
        << "maj_offsets = np.concatenate([-maj_offsets[1:][::-1],maj_offsets,])\n"
        << "min_offsets = np.concatenate([-min_offsets[1:][::-1],min_offsets,])\n"
        << "x_maj_s = np.interp(min_offsets, s_min, x_min)\n"
        << "y_maj_s = np.interp(min_offsets, s_min, y_min)\n"
        << "x_min_s = np.interp(maj_offsets, s_maj, x_maj)\n"
        << "y_min_s = np.interp(maj_offsets, s_maj, y_maj)\n\n"
        << "files_mod = [f for f in sorted(os.listdir(outfolder)) if outprefix+'mod' in f] \n"
        << "if len(files_mod)==0: raise FileNotFoundError('ERROR: no model in output directory') \n\n"
        << "cmap = plt.cm.PuBu_r \n"
        << "colors = cmap(np.linspace(0, 1, 13)) \n"
        <<  std::endl << std::endl
        << "def plot_pv(i,ax,ax_m,x0,y0,theta,**kwargs):\n\n"
        << "  s,x,y = build_path(x0, y0, rad/cdeltsp, np.full(len(rad),theta), xsize, ysize)\n"
        << "  pv_data = extract_pv(data,x,y)\n"
        << "  pv_mod  = extract_pv(data_mod,x,y)\n"
        << "  s_arcsec = s * cdeltsp\n"
        << "  ext = [s_arcsec[0], s_arcsec[-1], zmin_wcs-vsys_m,zmax_wcs-vsys_m]\n\n"
        << "  ax.imshow(pv_data,origin='lower',cmap='Blues',norm=norm,aspect='auto',interpolation='nearest',extent=ext)\n"
        << "  ax.contour(pv_data,v,origin='lower',linewidths=0.7,colors='#00008B',extent=ext) \n"
        << "  ax.contour(pv_data,v_neg,origin='lower',linewidths=0.1,colors='gray',extent=ext) \n"
        << "  ax.contour(pv_mod,v,origin='lower',linewidths=1.2,colors='#B22222',extent=ext) \n"
        << "  ax.axhline(y=0,color='black') \n"
        << "  ax.axvline(x=0,color='black') \n"
        << "  ax.grid(color='gray', linestyle='--', linewidth=0.3) \n\n"
        << "  if plotmask:\n"
        << "    pv_mask = extract_pv(data_mas,x,y) \n"
        << "    ax.contour(pv_mask,levels=[0],origin='lower',linewidths=1.5,colors='k',extent=ext)\n\n"
        << "  ax.set_xlim(-1.3*rad[-1],1.3*rad[-1])\n"
        << "  #ax.set_ylim(-1.1*(proj_vmax+2*max_vdisp),1.1*(proj_vmax+2*max_vdisp))\n"
        << "  ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5))\n"
        << "  ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5))\n"
        << "  if i==0 or i==5: \n"
        << "    ax.tick_params(labelbottom=True,labelleft=True)\n"
        << "    ax.set_ylabel(r'$\\mathrm{\\Delta V_{LOS}}$ (km/s)',fontsize=labsize+2) \n"
        << "  else: ax.tick_params(labelbottom=True,labelleft=False)\n"
        << "  ax.set_xlabel('Offset (arcsec)',fontsize=labsize)\n"
        << "  ax_m.plot(x+pad_x[0],y+pad_y[0],**kwargs) \n\n";
    
    pyf << "# Beginning position-velocity plot \n"
        << "for k in range (len(files_mod)): \n"
        << "  image_mod = fits.open(outfolder+files_mod[k]) \n"
        << "  data_mod = image_mod[0].data[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1]\n\n"
        << "  fig, axes = defineaxis(4,5,0.3,0.3,xsep=0.01,ysep=0.1,fig_width=10,fig_heigth=10)\n"
        << "  axes = np.ravel(axes)\n\n"
        << "  pos = axes[4].get_position()\n"
        << "  ax_m1 = fig.add_axes([pos.x1+0.03,pos.y0-0.2,0.3,0.3])\n"
        << "  pos = axes[14].get_position()\n"
        << "  ax_m2 = fig.add_axes([pos.x1+0.03,pos.y0-0.2,0.3,0.3])\n"
        << std::endl
        << "  # Plotting moment map\n"
        << "  pm = np.pad(mom0,(pad_y, pad_x),mode='constant',constant_values=np.nan) \n"
        << "  for ax in [ax_m1,ax_m2]:\n"
        << "    ax.set_prop_cycle(color=colors)\n"
        << "    ax.tick_params(labelbottom=False,labelleft=False,top=False,bottom=False,left=False,right=False)\n"
        << "    ax.imshow(pm,origin='lower',cmap='Greys',aspect='auto',interpolation='nearest',norm=norm0,extent=extm)\n"
        << "    ax.set_title('Slices')\n"
        << "    ax.plot(xcen_m+pad_x[0],ycen_m+pad_y[0],'x',color='yellow',markersize=7,mew=1.5,zorder=12) \n"
        << "    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5))\n"
        << "    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5))\n\n"
        << "  pos = ax_m2.get_position()\n"
        << "  ax_m2.set_position([pos.x0,pos.y0 - 0.07,pos.width,pos.height])\n"
        << std::endl
        << "  dtheta = 18\n"
        << "  for i,ax in enumerate(axes[0:10]):\n"
        << "    theta = pa_m + i*dtheta  \n"
        << "    ls = 'solid' if i<5 else 'dashed'\n"
        << "    plot_pv(i,ax,ax_m1,xcen_m,ycen_m,theta,ls=ls,label=rf'{i*dtheta:.0f}$^\\circ$')\n"
        << "    ax.text(0.03,1.03,rf\"$\\theta$ = {i*dtheta:.0f}$^\\circ$\",transform=ax.transAxes,fontsize=labsize)\n"
        << "    theta += 15.\n"
        << std::endl
        << "  for i,ax in enumerate(axes[10:]):\n"
        << "    pos = ax.get_position()\n"
        << "    ax.set_position([pos.x0,pos.y0 - 0.05,pos.width,pos.height])\n"
        << std::endl
        << "    if i<5:\n"
        << "      xi,yi,off,theta = x_maj_s[i],y_maj_s[i],min_offsets[i],pa_m\n"
        << "      ls = 'solid'\n"
        << "    else:\n"
        << "      xi,yi,off,theta = x_min_s[i-5],y_min_s[i-5],maj_offsets[i-5],pa_m + 90\n"
        << "      ls = 'dashed'\n"
        << std::endl
        << "    plot_pv(i,ax,ax_m2,xi,yi,theta,ls=ls,label=rf\"{off:.1f}''\")\n"
        << "    ax.text(0.03,1.03,rf\"offset = {off:.1f}''\",transform=ax.transAxes,fontsize=labsize)\n"
        << "    ax_m2.plot(xi+pad_x[0],yi+pad_y[0],'o',c='royalblue',zorder=11)\n"
        << "  leg_opts = dict(ncol=2,loc='upper left',bbox_to_anchor=(0.0, -0.05),frameon=False)\n"
        << "  ax_m1.legend(**leg_opts,handlelength=3.3)\n"
        << "  ax_m2.legend(**leg_opts,handlelength=2)\n"
        << std::endl
        << "  outfile = '%s_pvslices'%outprefix \n"
        << "  ntype   = files_mod[k][files_mod[k].rfind('_'):].replace('.fits','') \n"
        << "  fig.savefig(outfolder+outfile+ntype+'.pdf', bbox_inches='tight') \n"
        << "  image_mod.close() \n"
        << std::endl
        << "image.close() \n"
        << "image_mas.close() \n";
    
    pyf.close();

    /////////////////////////////////////////////////////////////////////////
    /// Script to plot kinematic maps of model vs data
    /////////////////////////////////////////////////////////////////////////
    scriptnames.push_back("plot_kinmaps.py");
    pyf.open((in->pars().getOutfolder()+"plotscripts/"+scriptnames.back()).c_str());

    pyf << "#######################################################################\n"
        << "#### This script writes a plot of kinematic maps of model and data ####\n"
        << "#######################################################################\n"
        << "import numpy as np \n"
        << "import os \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "from matplotlib.colorbar import ColorbarBase \n"
        << "from plot_utils import calculate_padding, defineaxis \n"
        << "from astropy.io import fits \n"
        << "from astropy.visualization import PercentileInterval \n"
        << "from copy import copy \n"
        << "mpl.rc('xtick',direction='in') \n"
        << "mpl.rc('ytick',direction='in') \n"
        << "plt.rc('font',family='sans-serif',serif='Helvetica',size=10) \n"
        << "params = {'text.usetex': False, 'mathtext.fontset': 'cm', 'mathtext.default': 'regular'} \n"
        << "plt.rcParams.update(params) \n"
        << std::endl
        << "gname = '" << in->Head().Name() <<"' \n"
        << "outfolder = '" << in->pars().getOutfolder() <<"' \n"
        << "outprefix = '" << in->pars().getOutPrefix() <<"' \n"
        << "twostage = " << par.TWOSTAGE << " \n"
        << "xmin, xmax = " << xmin << ", " << xmax << std::endl
        << "ymin, ymax = " << ymin << ", " << ymax << std::endl
        << std::endl
        << "# Opening maps and retrieving intensity map units\n"
        << "f0 = fits.open(outfolder+'/maps/'+outprefix+'_0mom.fits') \n"
        << "f1 = fits.open(outfolder+'/maps/'+outprefix+'_1mom.fits') \n"
        << "f2 = fits.open(outfolder+'/maps/'+outprefix+'_2mom.fits') \n"
        << "bunit = f0[0].header['BUNIT'] \n"
        << "bunit = bunit.replace(' ', '').lower() \n"
        << "# Now plotting moment maps \n"
        << "mom0 = f0[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "mom1 = f1[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "mom2 = f2[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "maskmap = np.copy(mom1) \n"
        << "maskmap[mom1==mom1] = 1 \n"
        << "ysize, xsize = mom0.shape \n"
        << "pad_x, pad_y = calculate_padding(xsize,ysize) \n"
        << std::endl
        << "if twostage: rad,inc,pa,xpos,ypos,vsys = np.genfromtxt(outfolder+'rings_final2.txt',usecols=(1,4,5,9,10,11),unpack=True) \n"
        << "else: rad,inc,pa,xpos,ypos,vsys = np.genfromtxt(outfolder+'rings_final1.txt',usecols=(1,4,5,9,10,11),unpack=True) \n"
        << "xcen_m,ycen_m,inc_m,pa_m,vsys_m=np.nanmean((xpos,ypos,inc,pa,vsys),axis=1) \n"
        << "xcen, ycen = xcen_m-xmin+pad_x[0], ycen_m-ymin+pad_y[0] \n"
        << std::endl
        << "files_mod0 = [f for f in sorted(os.listdir(outfolder+'maps/')) if outprefix+'mod_0mom' in f] \n"
        << "files_mod1 = [f for f in sorted(os.listdir(outfolder+'maps/')) if outprefix+'mod_1mom' in f] \n"
        << "files_mod2 = [f for f in sorted(os.listdir(outfolder+'maps/')) if outprefix+'mod_2mom' in f] \n"
        << "if len(files_mod0)==0 or len(files_mod1)==0 or len(files_mod2)==0 : raise FileNotFoundError('ERROR: no model maps in output directory') \n\n"
        << std::endl
        << "cmaps = [plt.get_cmap('Spectral_r'),plt.get_cmap('coolwarm'),plt.get_cmap('PuOr_r')] \n"
        << "barlab = ['Intensity ('+bunit+')', r'V$_\\mathrm{LOS}$ (km/s)', r'$\\sigma_\\mathrm{obs}$ (km/s)'] \n"
        << "barlab2 = [r'I$_\\mathrm{res}$ ('+bunit+')', r'V$_\\mathrm{res}$ (km/s)', r'$\\sigma_\\mathrm{res}$ (km/s)'] \n"
        << "titles = ['DATA', 'MODEL','RESIDUALS'] \n"
        << "mapname = ['MOMENT 0TH', 'MOMENT 1ST', 'MOMENT 2ND'] \n"
        << "x = np.arange(0,xmax-xmin,0.1) \n"
        << "y = np.tan(np.radians(pa_m-90))*(x-xcen)+ycen \n"
        << "ext = [0,xmax-xmin+pad_x[0]+pad_x[1],0, ymax-ymin+pad_y[0]+pad_y[1]] \n"
        << "rad_pix = rad/"<< in->Head().PixScale()*arcconv << std::endl
        << "try: nr = len(rad_pix) \n"
        << "except: nr = 1 \n"
        << "interval = PercentileInterval(99.5) \n"
        << std::endl
        << "for k in range (len(files_mod0)): \n"
        << "  mom0_mod = fits.open(outfolder+'/maps/'+files_mod0[k])[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "  mom1_mod = fits.open(outfolder+'/maps/'+files_mod1[k])[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "  mom2_mod = fits.open(outfolder+'/maps/'+files_mod2[k])[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "  to_plot = [[mom0,mom1-vsys_m,mom2],[mom0_mod,mom1_mod-vsys_m,mom2_mod],[mom0-mom0_mod,mom1-mom1_mod,mom2-mom2_mod]] \n"
        << std::endl
        << "  nrows, ncols, x_len, y_len = 3, 3, 0.2, 0.2 \n"
        << "  fig, ax = defineaxis(nrows,ncols,x_len,y_len,xsep=0.00,ysep=0.08,fig_width=11,fig_heigth=11) \n\n"
        << "  for i in range (ax.shape[0]): \n"
        << "    cmap = copy(cmaps[i]) \n"
        << "    cmap.set_bad('w',1.) \n"
        << "    vmin, vmax = interval.get_limits(to_plot[1][i]) \n"
        << "    vmin, vmax = (-1.1*np.nanmax(vmax),1.1*np.nanmax(vmax)) if i==1 else (vmin,vmax) \n"
        << "    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax) \n"
        << "    cbax = fig.add_axes([ax[i][0].get_position().x0,ax[i][0].get_position().y0-0.025,2*x_len-0.003,0.02]) \n"
        << "    cb1 = ColorbarBase(cbax, orientation='horizontal', cmap=cmap, norm=norm) \n"
        << "    cb1.set_label(barlab[i],fontsize=13) \n"
        << "    for j in range (ax.shape[1]): \n"
        << "      axis = ax[i][j] \n"
        << "      axis.tick_params(labelbottom=False,labelleft=False,right=True,top=True) \n"
        << "      axis.set_xlim(ext[0],ext[1]) \n"
        << "      axis.set_ylim(ext[2],ext[3]) \n"
        << "      if j==2: \n"
        << "        vmax = np.nanmax(interval.get_limits(to_plot[j][i])) \n"
        << "        norm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax) \n"
		<< "      pp = np.pad(to_plot[j][i]*maskmap,(pad_y, pad_x),mode='constant',constant_values=np.nan) \n"
		<< "      axis.imshow(pp,origin='lower',cmap=cmap,norm=norm,aspect='auto',extent=ext,interpolation='nearest') \n"
        << "      axis.plot(xcen,ycen,'x',color='#000000',markersize=7,mew=1.5,zorder=10) \n"
        << std::endl
        << "      if i==0: \n"
        << "        axis.text(0.5,1.05,titles[j],ha='center',transform=axis.transAxes,fontsize=15) \n"
        << "        axis.plot(x,y,'--',color='k',linewidth=1) \n"
        << "        if j!=2 and nr>3:  \n"
        << "          axmaj = rad_pix[-1] \n"
        << "          axmin = axmaj*np.cos(np.radians(inc_m))  \n"
        << "          posa = np.radians(pa_m-90)  \n"
        << "          t = np.linspace(0,2*np.pi,100)  \n"
        << "          xt = xcen+axmaj*np.cos(posa)*np.cos(t)-axmin*np.sin(posa)*np.sin(t)  \n"
        << "          yt = ycen+axmaj*np.sin(posa)*np.cos(t)+axmin*np.cos(posa)*np.sin(t)  \n"
        << "          axis.plot(xt,yt,'-',c='k',lw=0.8)  \n"
        << "      elif i==1: \n"
        << "        axis.plot(x,y,'--',color='k',linewidth=1) \n"
        << "        if nr<10: \n"
        << "          x_pix = rad_pix*np.cos(np.radians(pa_m-90)) \n"
        << "          y_pix = rad_pix*np.sin(np.radians(pa_m-90)) \n"
        << "          axis.scatter(x_pix+xcen,y_pix+ycen,c='grey',s=12) \n"
        << "          axis.scatter(xcen-x_pix,ycen-y_pix,c='grey',s=12) \n"
        << "        if nr>5 and not all(np.diff(pa)==0): \n"
        << "          x_pix = rad_pix*np.cos(np.radians(pa-90)) \n"
        << "          y_pix = rad_pix*np.sin(np.radians(pa-90)) \n"
        << "          axis.plot(xcen-x_pix,ycen-y_pix,'-',color='grey',lw=1) \n"
        << "          axis.plot(x_pix+xcen,y_pix+ycen,'-',color='grey',lw=1) \n"
        << "        if j!=2: \n"
        << "          cmax = np.nanmax(np.abs(to_plot[0][i])) \n"
        << "          levels = np.linspace(0.166 * cmax, cmax, 6) \n"
        << "          axis.contour(pp,levels=[0],colors='forestgreen',origin='lower',extent=ext) \n"
        << "          axis.contour(pp,levels=-levels[::-1],colors='navy',linewidths=0.7,origin='lower',extent=ext,linestyles='solid') \n"
        << "          axis.contour(pp,levels=levels,colors='darkred',linewidths=0.7,origin='lower',extent=ext) \n"
        << "      if j==0: axis.text(-0.12,0.5,mapname[i],va='center',rotation=90,transform=axis.transAxes,fontsize=15) \n\n"
        << "    cbax = fig.add_axes([ax[i][2].get_position().x0+0.003,ax[i][2].get_position().y0-0.025,x_len-0.003,0.02]) \n"
        << "    cb2 = ColorbarBase(cbax, orientation='horizontal', cmap=cmap, norm=norm) \n"
        << "    cb2.set_label(barlab2[i],fontsize=13) \n"
        << "    cb2.ax.locator_params(nbins=3) \n"
        << "    for c in [cb1,cb2]: \n"
        << "      c.solids.set_edgecolor('face') \n"
        << "      c.outline.set_linewidth(0) \n"
        << std::endl
        << "  outfile = '%s_maps'%outprefix \n"
        << "  ntype   = files_mod0[k][files_mod0[k].rfind('_'):].replace('.fits','') \n"
        << "  fig.savefig(outfolder+outfile+ntype+'.pdf', bbox_inches = 'tight') \n\n";

    pyf.close();

    /////////////////////////////////////////////////////////////////////////
    /// One script to rule them all
    /////////////////////////////////////////////////////////////////////////
    pyf.open((in->pars().getOutfolder()+"plotscripts/plot_all.py").c_str());
    pyf << "########################################################################\n"
        << "#### This script simply calls all other python scripts for plotting ####\n"
        << "########################################################################\n"
        << "import os \n"
        << std::endl
        << "scriptdir = '"<< in->pars().getOutfolder() << "plotscripts/' \n"
        << "cmd = '' \n"
        << std::endl
        << "for f in os.listdir(scriptdir): \n"
        << "  if '.py' in f and f!='plot_all.py' and f!='plot_utils.py': \n"
        << "    cmd += 'python \"%s/%s\" & '%(scriptdir,f) \n"
        << std::endl
        << "os.system(cmd[:-2]) \n";

    pyf.close();
    
    return scriptnames;

}
template std::vector<std::string> Galfit<float>::writeScripts_Python();
template std::vector<std::string> Galfit<double>::writeScripts_Python();
//*/

template <class T>
void Galfit<T>::printDetails (Rings<T> *dr, T fmin, long pix, std::ostream& str) {

    int m=7, n=9;

    bool details = true;
    if (details) {
        str << endl << setfill('-') << setw(80) << " " << setfill(' ') << endl;
        str << setw(n) << right << "Fmin";
        str << setw(n) << right << "Pixs";
        if (mpar[VSYS])  str << setw(m+1) << right << "VSYS";
        if (mpar[XPOS])  str << setw(m) << right << "XPOS";
        if (mpar[YPOS])  str << setw(m) << right << "YPOS";
        if (mpar[VROT])  str << setw(m) << right << "VROT";
        if (mpar[VDISP]) str << setw(m-1) << right << "DISP";
        if (mpar[INC])   str << setw(m-1) << right << "INC";
        if (mpar[PA])    str << setw(m) << right << "PA";
        if (mpar[Z0])    str << setw(m-1) << right << "Z0";
        str << endl << setfill('-') << setw(80) << " " << setfill(' ') << endl;
    }

    str << setw(n) << scientific << setprecision(2) << right << fmin;
    str << setw(n) << fixed << right << pix;

    str << setprecision(1);
    if (mpar[VSYS])  str << setw(m+1) << right << dr->vsys.back();
    if (mpar[XPOS])  str << setw(m) << right << dr->xpos.back();
    if (mpar[YPOS])  str << setw(m) << right << dr->ypos.back();
    if (mpar[VROT])  str << setw(m) << right << dr->vrot.back();
    if (mpar[VDISP]) str << setw(m-1) << right << dr->vdisp.back();
    if (mpar[INC])   str << setw(m-1) << right << dr->inc.back();
    if (mpar[PA])    str << setw(m) << right << dr->phi.back();
    if (mpar[Z0])    str << setw(m-1) << right << dr->z0.back();

    str << endl;

}
template void Galfit<float>::printDetails (Rings<float> *, float, long, std::ostream&);
template void Galfit<double>::printDetails (Rings<double> *, double, long, std::ostream&);


template <class T>
void Galfit<T>::showInitial (Rings<T> *inr, std::ostream& Stream) {

    int m=9;
    int n=12;
    GALFIT_PAR &pp = in->pars().getParGF();
    Stream << showpoint << fixed << setprecision(2) << endl;

    Stream << setfill('=') << setw(44) << right << " Initial parameters " << setw(25) << " ";
    Stream << setfill(' ') << endl;

    Stream << "    (i) input by the user    (d) default    (e) estimated by me" << endl << endl;

    string s;
    s = "   Fitting #" + to_string<int>(inr->nr);
    if (pp.NRADII==-1) s += "(e) ";
    else s += "(i) ";
    s += "rings of width " + to_string<T>(inr->radsep);
    if (pp.RADSEP==-1) s += "(e) ";
    else s += "(i) ";
    s += "arcsec";
    Stream << s << endl << endl;


    s = "    Xpos";
    if (pp.XPOS=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->xpos[0] << left << setw(m) << "  pix";


    s = "        Ypos";
    if (pp.YPOS=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n+4) << left << s << setw(3) << right << "= "
         << setw(m-1) << inr->ypos[0]
         << left << setw(m) << "  pix" << endl;

    s = "    Vsys";
    if (pp.VSYS=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->vsys[0] << left << setw(m) << "  km/s";

    s = "        Vrot";
    if (pp.VROT=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n+4) << left << s << setw(3) << right << "= "
         << setw(m-1) << inr->vrot[0]
         << left << setw(m) << "  km/s" << endl;

    s = "    Inc";
    if (pp.INC=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->inc[0] << left << setw(m) << "  deg";

    s = "        PA";
    if (pp.PHI=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n+4) << left << s << setw(3) << right << "= "
           << setw(m-1) << inr->phi[0] << left << setw(m) << "  deg" << endl;

    s = "    Z0";
    if (pp.Z0=="-1") s += " (d)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->z0[0] << left << setw(m) << "  arcs";

    s = "        Disp";
    if (pp.VDISP=="-1") s += " (d)";
    else s += " (i)";
    Stream << setw(n+4) << left << s << setw(3) << right << "= "
         << setw(m-1) << inr->vdisp[0]
         << left << setw(m) << "  km/s" << endl;
    
    s = "    Vrad";
    if (pp.VRAD=="-1") s += " (d)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->vrad[0] << left << setw(m) << "  km/s";

    Stream   << endl << endl;

    Stream   << endl;

}
template void Galfit<float>::showInitial(Rings<float>*,std::ostream&);
template void Galfit<double>::showInitial(Rings<double>*,std::ostream&);


template <class T>
void Galfit<T>::printInitial (Rings<T> *inr, std::string outfile) {

    int m=11;
    std::ofstream initout(outfile.c_str());    
    initout << "#" << setfill('=');
    initout << setw(66) << right << " Initial parameters " << setw(46) << " " << endl;
    initout << setfill(' ') << setprecision(3) << fixed;
    initout << left << setw(m) << "#RAD(arcs)"
            << setw(m) << "VROT(km/s)"
            << setw(m) << "DISP(km/s)"
            << setw(m) << "INC(deg)"
            << setw(m) << "P.A.(deg)"
            << setw(m) << "Z0(arcs)"
            << setw(m) << "SIG(E20)"
            << setw(m) << "XPOS(pix)"
            << setw(m) << "YPOS(pix)"
            << setw(m) << "VSYS(km/s)" << endl;
    
    for (int i=0; i<inr->nr; i++) {
        initout << setw(m) << inr->radii[i]
            << setw(m) << inr->vrot[i]
            << setw(m) << inr->vdisp[i]
            << setw(m) << inr->inc[i]
            << setw(m) << inr->phi[i]
            << setw(m) << inr->z0[i]
            << setw(m) << inr->dens[i]
            << setw(m) << inr->xpos[i]
            << setw(m) << inr->ypos[i]
            << setw(m) << inr->vsys[i] << endl;
    }
    initout.close();

}
template void Galfit<float>::printInitial(Rings<float>*,std::string);
template void Galfit<double>::printInitial(Rings<double>*,std::string);


}

