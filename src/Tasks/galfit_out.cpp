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

    if (verb) std::cout << " Preparing a bunch of cool outputs..." << std::endl;

    std::string outfold = in->pars().getOutfolder();
    std::string object = in->Head().Name();

    // Get the total intensity, velocity and dispersion maps of data
    mkdirp((outfold+"maps/").c_str());
    if (verb) std::cout << "    Extracting "<< randomAdjective(2) << " 2D maps..." << std::flush;
    std::vector< MomentMap<T> > allmaps = getAllMoments<T>(in,true,nullptr,"MOMENT");
    allmaps[0].fitswrite_2d((outfold+"maps/"+object+"_0mom.fits").c_str());
    allmaps[1].fitswrite_2d((outfold+"maps/"+object+"_1mom.fits").c_str());
    allmaps[2].fitswrite_2d((outfold+"maps/"+object+"_2mom.fits").c_str());
    if (verb) std::cout << " Done." << std::endl;

    // Calculate the total flux inside last ring in data
    T *ringreg = RingRegion(outr,in->Head());
    float totflux_data=0, totflux_model=0;
    for (auto i=0; i<in->DimX()*in->DimY(); i++) {
        if (!isNaN(ringreg[i])) {
            for (auto z=0; z<in->DimZ(); z++) {
                size_t npix = i+z*in->DimY()*in->DimX();
                totflux_data += in->Array(npix)*in->Mask()[npix];
            }
        }
    }

    // Calculate radial profile along the output rings
    if (verb) std::cout << "    Deriving " << randomAdjective(1) << " radial profile..." << std::flush;
    T meanPA = findMean(&outr->phi[0], outr->nr);
    int nseg = 1;
    float segments[4] = {0, 360., 0., 0};
    if (par.SIDE=="A") {
        nseg = 2;
        segments[2]=-90;
        segments[3]=90;
    }
    else if(par.SIDE=="R") {
        nseg = 2;
        segments[2]=90;
        segments[3]=-90;
    }
    if (meanPA>180) std::swap(segments[2], segments[3]);

    Tasks::Ellprof<T> ell(&allmaps[0],outr,nseg,segments);
    ell.RadialProfile();
    
    if (normtype=="AZIM" || normtype=="BOTH") {
        double profmin=FLT_MAX;
        for (auto i=0; i<outr->nr; i++) {
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
        for (auto i=0; i<outr->nr; i++) {
            if (outr->dens[i]>0) outr->dens[i]=factor*fabs(ell.getMean(i))*1E20;
            //if (outr->dens[i]==0) outr->dens[i]=profmin*1E20;
        }
    }
    
    std::string dens_out = outfold+"densprof.txt";
    std::ofstream fileo(dens_out.c_str());
    ell.printProfile(fileo,nseg-1);
    fileo.close();
    //ell.printProfile(std::cout);
    if (verb) std::cout << " Done." << std::endl;

    // Deleting bad rings, for which I set dens=0 in galfit()
    //for (auto i=0; i<outr->nr; i++) 
    //    if (outr->dens[i]==0) outr->deleteRing(i);

    if (verb) std::cout << "    Calculating the very last model..." << std::flush;
    Model::Galmod<T> *mod = getModel();
    mod->Out()->Head().setMinMax(0.,0.);
    mod->Out()->Head().setName(object+"mod");
    
    T *outarray = mod->Out()->Array();

    if (normtype=="AZIM" || normtype=="BOTH") {

        // The final model has been build from the azimuthal profile calculated before,
        // thus just need to rescale the model to the total flux of data inside last ring.

        // Calculate total flux of model within last ring
        totflux_data = totflux_model = 0;
        
        /*
        for (auto i=0; i<in->DimX()*in->DimY(); i++) {
            if (!isNaN(ringreg[i])) {
                for (auto z=0; z<in->DimZ(); z++) {
                    long npix = i+z*in->DimY()*in->DimX();
                    totflux_model += outarray[npix]*in->Mask(npix);
                    totflux_data  += in->Array(npix)*in->Mask(npix);
                }
            }
            //else {
            //      for (size_t z=0; z<in->DimZ(); z++)
            //             outarray[i+z*in->DimY()*in->DimX()]=0;
            //}
        }
        */
        
        //////////////////////////////////////////////////////////////////
        // New normalization method, taking only pixels in a given range 
        // (see start, stop below)
        std::vector<T> d, m;
        for (auto i=0; i<in->DimX()*in->DimY(); i++) {
            if (!isNaN(ringreg[i])) {
                for (auto z=0; z<in->DimZ(); z++) {
                    long npix = i+z*in->DimY()*in->DimX();
                    if (in->Mask(npix) && in->Array(npix)>0) {
                        d.push_back(in->Array(npix));
                        m.push_back(outarray[npix]);
                    }
                }
            }
        }
        
        std::sort (d.begin(), d.end()); 
        std::sort (m.begin(), m.end()); 
        int start = 0.6*d.size();
        int stop  = 0.99*d.size();
        
        for (int i=start; i<stop; i++) {
            totflux_data  += d[i];
            totflux_model += m[i];
        }
        //*/
        ////////////////////////////////////////////////////////////////////////
        

        double factor = totflux_data/totflux_model;
        for (auto i=in->NumPix(); i--;) outarray[i] *= factor;
        if (verb) std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing " << randomAdjective(1) << " azimuthally-normalized model..." << std::flush;
        std::string mfile = outfold+object+"mod_azim.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());
        writePVs(mod->Out(),"_azim");
        if (verb) std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing " << randomAdjective(2) << " kinematic maps..." << std::flush;
        std::vector< MomentMap<T> > modmaps = getAllMoments<T>(mod->Out(),false,nullptr,"MOMENT");
        modmaps[0].fitswrite_2d((outfold+"maps/"+object+"_azim_0mom.fits").c_str());
        modmaps[1].fitswrite_2d((outfold+"maps/"+object+"_azim_1mom.fits").c_str());
        modmaps[2].fitswrite_2d((outfold+"maps/"+object+"_azim_2mom.fits").c_str());

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
        std::string mfile = outfold+object+"mod_local.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());

        writePVs(mod->Out(),"_local");
        if (verb) std::cout << " Done." << std::endl;

        if (verb) {
            std::cout << "    Writing " << std::flush;
            if (normtype=="BOTH") std::cout << "even more " << std::flush;
            std::cout << randomAdjective(2) << " kinematic maps..." << std::flush;
        }
        std::vector< MomentMap<T> > modmaps = getAllMoments<T>(mod->Out(),false,nullptr,"MOMENT");
        modmaps[0].fitswrite_2d((outfold+"maps/"+object+"_local_0mom.fits").c_str());
        modmaps[1].fitswrite_2d((outfold+"maps/"+object+"_local_1mom.fits").c_str());
        modmaps[2].fitswrite_2d((outfold+"maps/"+object+"_local_2mom.fits").c_str());

        if (verb) std::cout << " Done." << std::endl;
    }

    if (normtype=="NONE") {
        
        // The final model has been build from an input density profile in cm^-2
        // Here I renormalize to have to the integral of the input density profile.
        // Output cube will have units of Jy/beam for HI data.
        
        // Current flux model is assumed to be in JY/beam
        mod->Out()->Head().setBunit("JY/BEAM");
        // Getting current profile
        MomentMap<T> *dmap = new MomentMap<T>;
        dmap->input(mod->Out());
        dmap->ZeroMoment(false);
        Tasks::Ellprof<T> ell(dmap,outr,nseg,segments);
        ell.RadialProfile();
        delete dmap;
        
        // Calculating integral of input density profile and current
        double totmass_req=0, totmass_curr=0;
        for (int i=0; i<outr->nr; i++) {
            totmass_req  += outr->dens[i];                                  // In cm^-2
            totmass_curr += ell.getSurfDensFaceOn(i)/in->Head().BeamArea(); // In JY*KM/S/pc2
        }
        
        // Converting everything to Msun/pc2
        const double arctorad = 1/3600.*M_PI/180.;
        totmass_req  = 3.0856*3.0856*8.41185687e-22*totmass_req;
        totmass_curr = 2.36E-07*totmass_curr/(arctorad*arctorad);

        // Re-normalization
        for (auto i=in->NumPix(); i--;) outarray[i] *= totmass_req/totmass_curr;
        if (verb) std::cout << " Done." << std::endl;
        
        if (verb) std::cout << "    Writing model..." << std::flush;
        std::string mfile = outfold+object+"mod_nonorm.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());
        writePVs(mod->Out(),"_nonorm");
        if (verb) std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing " << randomAdjective(1) << " kinematic maps..." << std::flush;
        std::vector< MomentMap<T> > modmaps = getAllMoments<T>(mod->Out(),false,nullptr,"MOMENT");
        modmaps[0].fitswrite_2d((outfold+"maps/"+object+"_nonorm_0mom.fits").c_str());
        modmaps[1].fitswrite_2d((outfold+"maps/"+object+"_nonorm_1mom.fits").c_str());
        modmaps[2].fitswrite_2d((outfold+"maps/"+object+"_nonorm_2mom.fits").c_str());
        if (verb) std::cout << " Done." << std::endl;
        
    }
    
    // Adding noise to a model
    if (par.NOISERMS!=0) { 
        if (verb) std::cout << "    Writing " << randomAdjective(1) << " noisy model..." << std::flush;
        size_t nPix = mod->Out()->NumPix();
        // Getting gaussian noise
        T *noise = SimulateNoise<T>(par.NOISERMS,nPix);
        double fac = 1.;
        if (par.SM) {
            // Smoothing the noise
            Beam obeam = {0., 0., 0};    
            Beam nbeam = {in->Head().Bmaj()*arcconv,in->Head().Bmin()*arcconv,in->Head().Bpa()};
            Smooth3D<T> *sm = new Smooth3D<T>;
            sm->setUseBlanks(false);
            sm->smooth(in, obeam, nbeam, noise, noise);
            // Rescaling the noise to match the required rms
            T stdd = findStddev(noise,nPix);
            fac = par.NOISERMS/stdd;
            delete sm;
        }
        // Add noise to the model
        for (size_t i=nPix; i--;) mod->Out()->Array(i) += (noise[i]*fac);    
        // Writing to FITS file
        std::string mfile = outfold+object+"mod_noise.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());
        delete [] noise;
        if (verb) std::cout << " Done." << std::endl;
    }

    
    // Computing asymmetric drift correction
    if (par.flagADRIFT) {
        if (verb) std::cout << "    Computing asymmetric drift correction..." << std::flush;
        T *dens_m = new T[outr->nr];
        int strad = par.STARTRAD<inr->nr ? par.STARTRAD : 0;
        for (int i=0; i<outr->nr; i++) dens_m[i] = ell.getMedian(i);
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
void Galfit<T>::writePVs(Cube<T> *mod, std::string suffix) {

    // Extracts PVs along major and minor axis and write them in fits.
    std::string outfold = in->pars().getOutfolder()+"pvs/";
    std::string object = in->Head().Name();
    
    mkdirp((outfold).c_str());
    
    float meanPA = findMedian(&outr->phi[0], outr->nr);
    float meanXpos = findMedian(&outr->xpos[0], outr->nr);
    float meanYpos = findMedian(&outr->ypos[0], outr->nr);
    float meanPAp90= meanPA+90<360 ? meanPA+90 : meanPA-90;
    
    
    // Extract pvs of data
    PvSlice<T> *pv_max = PositionVelocity(in,meanXpos,meanYpos,meanPA);
    std::string mfile = outfold+object+"_pv_a.fits";
    pv_max->fitswrite_2d(mfile.c_str());
    PvSlice<T> *pv_min = PositionVelocity(in,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+object+"_pv_b.fits";
    pv_min->fitswrite_2d(mfile.c_str());


    // Extract pvs of data
    PvSlice<T> *pv_max_m = PositionVelocity(mod,meanXpos,meanYpos,meanPA);
    mfile = outfold+object+"mod_pv_a"+suffix+".fits";
    pv_max_m->fitswrite_2d(mfile.c_str());
    PvSlice<T> *pv_min_m = PositionVelocity(mod,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+object+"mod_pv_b"+suffix+".fits";
    pv_min_m->fitswrite_2d(mfile.c_str());


    // Extract pvs of mask
    Cube<short> *m = new Cube<short>(in->AxisDim());
    m->saveHead(in->Head());
    m->saveParam(in->pars());
    m->pars().setANTIALIAS(0);
    m->Head().setMinMax(0.,0);
    for (size_t i=in->NumPix(); i--;) m->Array()[i] = short(mask[i]);
    
    PvSlice<short> *pv_max_ma = PositionVelocity(m,meanXpos,meanYpos,meanPA);
    mfile = outfold+object+"mask_pv_a.fits";
    pv_max_ma->fitswrite_2d(mfile.c_str());
    PvSlice<short> *pv_min_ma = PositionVelocity(m,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+object+"mask_pv_b.fits";
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
    std::string object = in->Head().Name();

    mkdirp((in->pars().getOutfolder()+"plotscripts/").c_str());
    std::string mfile = outfold+"plotscripts/gnuscript.gnu";
    std::ofstream gnu(mfile.c_str());

    float xtics = ceil(outr->nr/5.);
    xtics *= outr->radsep;
    
    while (outr->radii.back()/xtics>5.) xtics*=2.;
    while (outr->radii.back()/xtics<2.) xtics/=2.;
    
    int *nc = getErrorColumns();

    /// Setting global option
    gnu << "set terminal postscript eps enhanced color font 'Helvetica,14'" << endl
        << "set output '" << outfold << object << "_rc_inc_pa.eps'" << endl
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

    gnu << "set output '" << outfold << object << "_disp_vsys_z0.eps'" << endl
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
        << "plot '"<<in->pars().getOutfolder()<<"rings_final1.txt' ";

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

    gnu << "set output '" << outfold << object << "_xc_yc_cd.eps'" << endl
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
        maxa += 0.1*maxa/1.E20;
        mina = *min_element(&outr->dens[0], &outr->dens[0]+outr->nr);
        mina -= 0.1*mina/1.E20;
        gnu << "@BMARGIN" << endl << "@XTICS" << endl
            << "set xlabel 'Radius [arcsec]'" << endl
            << "set yrange [" <<mina<<":"<<maxa<<"]\n"
            << "set ylabel 'Surface density [10^20 atoms/cm^2]'\n"
            << "plot '"<<in->pars().getOutfolder()<<"rings_final1.txt' u 2:8 w lp ls 1";

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

    /// This function creates and runs a python script for plotting:
    ///  1) Channel maps
    ///  2) Position-Velocity diagrams along major/minor axis
    ///  3) Output paramenters
    ///  4) Moment maps
    ///  5) Asymmetric drift correction
    ///
    /// It needs all output fitsfiles to be in the output directory!

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
        << "gname = '" << in->Head().Name() <<"' \n"
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
        << "\trad2, vrot2,disp2,inc2,pa2,z02,xpos2,ypos2,vsys2, vrad2 = np.genfromtxt(outfolder+'rings_final2.txt',usecols=(1,2,3,4,5,7,9,10,11,12),unpack=True)\n";
    if (outr->nr==1) pyf << "\trad2,vrot2,disp2,inc2,pa2,z02,xpos2,ypos2,vsys2,vrad2 = np.array([rad2]),np.array([vrot2]),np.array([disp2]),np.array([inc2]),np.array([pa2]),np.array([z02]),np.array([xpos2]),np.array([ypos2]),np.array([vsys2]),np.array([vrad2])\n";
    pyf << "\terr2_l, err2_h = np.zeros(shape=(" << MAXPAR << ",len(rad2))), np.zeros(shape=(" << MAXPAR << ",len(rad2)))\n"
        << "\tcolor='#A0A0A0' \n"
        << "\tmax_rad = 1.1*np.nanmax(rad2) \n"
        << "\tmax_vrot,max_vdisp = np.maximum(max_vrot,np.nanmax(vrot2)),np.maximum(max_vdisp,np.nanmax(disp2)) \n";

    for (int i=0; i<MAXPAR; i++) {
        if ((i==0 || i==1 || i==9) && (nc[i]>0)) {
            int ff = i==9 ? 6 : 0;
            pyf << "\terr2_l[" << i << "], err2_h[" << i << "] = np.genfromtxt(outfolder+'rings_final2.txt',usecols=("
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
        << "\tbottom_corner[0], yl = 0.1, ylen \n"
        << "\tif i==0: yl *= 1.8 \n"
        << "\tfor j in range (ncols): \n"
        << "\t\tfig.add_axes([bottom_corner[0],bottom_corner[1],xlen,yl]) \n"
        << "\t\tfig.axes[-1].set_xlim(0,max_rad) \n"
        << "\t\tif i==nrows-1: fig.axes[-1].tick_params(labelbottom=True) \n"
        << "\t\telse: fig.axes[-1].tick_params(labelbottom=False) \n"
        << "\t\tbottom_corner[0]+=xlen+x_sep \n"
        << "\tbottom_corner[1]-=(ylen+y_sep) \n\n"
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
        << "fig.savefig(outfolder+'%s_parameters.pdf'%gname,bbox_inches='tight') \n";

    if (par.flagADRIFT) {
        pyf << "\n#Asymmetric drift correction\n"
            << "rad_a, vcirc, va2, dispr, fun, funr = np.genfromtxt(outfolder+'asymdrift.txt',usecols=(0,1,2,3,4,5),unpack=True)\n"
            << "if twostage: \n"
            << "\trad, vrot, disp = rad2, vrot2, disp2 \n"
            << "\terr1_l, err1_h = err2_l, err2_h\n"
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
            << "ax2.plot(rad_a,dispr,'-', color='#0FA45A',label=r'$\\sigma_\\mathrm{reg}$',zorder=1)\n"
            << "ax2.legend()\n"
            << "ax3.set_xlim(0,max_rad)\n"
            << "ax3.set_xlabel('Radius (arcsec)', fontsize=11, labelpad=5)\n"
            << "ax3.set_ylabel(r'$f = \\log(\\sigma_\\mathrm{gas}^2\\Sigma_\\mathrm{gas}\\cos(i)$)', fontsize=11)\n"
            << "ax3.plot(rad_a,fun,'o', color=color2,label=r'$f_\\mathrm{obs}$',zorder=0)\n"
            << "ax3.plot(rad_a,funr,'-', color='#0FA45A',label=r'$f_\\mathrm{reg}$',zorder=1)\n"
            << "ax3.legend()\n"
            << "fig.savefig(outfolder+'%s_asymmetricdrift.pdf'%gname,bbox_inches='tight')\n";
    }

    pyf.close();

    /////////////////////////////////////////////////////////////////////////
    /// Script to plot channel maps of model vs data
    /////////////////////////////////////////////////////////////////////////
    scriptnames.push_back("plot_chanmaps.py");
    pyf.open((in->pars().getOutfolder()+"plotscripts/"+scriptnames[1]).c_str());

    std::string cubefile = in->pars().getImageFile();
    if (in->pars().getFlatContsub()) cubefile = in->pars().getOutfolder()+in->Head().Name()+"_contsub.fits";

    pyf << "#####################################################################\n"
        << "#### This script writes a plot of channel maps of model and data ####\n"
        << "#####################################################################\n"
        << "import numpy as np \n"
        << "import os \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "import matplotlib.gridspec as gridspec \n"
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
        << "zsize= data.shape[0] \n"
        << "cdeltsp=" << in->Head().PixScale()*arcconv << std::endl
        << "cont = " << cont << std::endl
        << "v = np.array([1,2,4,8,16,32,64])*cont \n"
        << "v_neg = [-cont] \n"
        << "interval = PercentileInterval(99.5) \n"
        << "vmax = interval.get_limits(data)[1] \n"
        << "norm = ImageNormalize(vmin=cont, vmax=vmax, stretch=PowerStretch(0.5)) \n"
        << "xcen, ycen = xcen_m-xmin, ycen_m-ymin  \n"
        << std::endl
        << "files_mod = [] \n"
        << "for thisFile in sorted(os.listdir(outfolder)): \n"
        << "\tif 'mod_azim.fits'  in thisFile: files_mod.append(thisFile) \n"
        << "\tif 'mod_local.fits' in thisFile: files_mod.append(thisFile) \n"
        << "if len(files_mod)==0: exit('ERROR: no model in output directory') \n\n";

    pyf << "# Beginning channel map plot \n"
        << "for k in range (len(files_mod)): \n"
        << "\timage_mod = fits.open(outfolder+files_mod[k]) \n"
        << "\tdata_mod = image_mod[0].data[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
        << "\tfig = plt.figure(figsize=(8.27, 11.69), dpi=150) \n"
        << "\tgrid = [gridspec.GridSpec(2,5),gridspec.GridSpec(2,5),gridspec.GridSpec(2,5)] \n"
        << "\tgrid[0].update(top=0.90, bottom=0.645, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
        << "\tgrid[1].update(top=0.60, bottom=0.345, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
        << "\tgrid[2].update(top=0.30, bottom=0.045, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
        << std::endl
        << "\tnum = 0 \n"
        << "\tfor j in range (0,3): \n"
        << "\t\tfor i in range (0,5): \n"
        << "\t\t\tchan = int(num*(zsize)/15) \n"
        << "\t\t\tz = data[chan,:,:] \n"
        << "\t\t\tz_mod = data_mod[chan,:,:] \n"
        << "\t\t\t#New matplotlib draws wrong contours when no contours are found. This is a workaround.\n"
        << "\t\t\tif np.all(z_mod<v[0]): z_mod[:,:] =0\n"
        << "\t\t\tvelo_kms = (chan+1-" << 1-zmin << ")*" << cdelt3_kms << "+" << crval3_kms_pix1 << std::endl
        << "\t\t\tvelo = ' v = ' + str(int(velo_kms-vsys_m)) + ' km/s' \n"
        << "\t\t\tax = plt.subplot(grid[j][0,i]) \n"
        << "\t\t\tax.set_title(velo, fontsize=10,loc='left') \n"
        << "\t\t\tax.imshow(z,origin='lower',cmap = mpl.cm.Greys,norm=norm,aspect='auto',interpolation='none') \n"
        << "\t\t\tax.contour(z,v,origin='lower',linewidths=0.7,colors='#00008B') \n"
        << "\t\t\tax.contour(z,v_neg,origin='lower',linewidths=0.1,colors='gray') \n"
        << "\t\t\tax.plot(xcen,ycen,'x',color='#0FB05A',markersize=7,mew=2) \n"
        << "\t\t\tif plotmask: \n"
        << "\t\t\t\tax.contour(data_mas[chan],levels=[0],origin='lower',linewidths=2,colors='k') \n"
        << "\t\t\tif (j==i==0): \n"
        << "\t\t\t\tax.text(0, 1.4, gname, transform=ax.transAxes,fontsize=15,va='center') \n"
        << "\t\t\t\tlbar = 0.5*(xmax-xmin)*cdeltsp \n"
        << "\t\t\t\tltex = \"%.0f'' \"%lbar if lbar>10 else \"%.2f'' \"%lbar \n"
        << "\t\t\t\tif lbar>600: ltex = \"%.0f' \"%(lbar/60.) \n"
        << "\t\t\t\tax.annotate('', xy=(4.5, 1.4), xycoords='axes fraction', xytext=(5, 1.4),arrowprops=dict(arrowstyle='<->', color='k'))\n"
        << "\t\t\t\tax.text(4.75,1.50,ltex,transform=ax.transAxes,fontsize=11, ha='center')\n"
        << "\t\t\t\tbmaj, bmin, bpa = " << bmaj << "/float(xmax-xmin), " << bmin << "/float(ymax-ymin)," << bpa << std::endl
        << "\t\t\t\tbeam = mpl.patches.Ellipse(xy=(3.5, 1.4), width=bmaj, height=bmin, angle=bpa+90, \\\n"
        << "\t\t\t\t                           color='#5605D0', clip_on=False, transform=ax.transAxes, alpha=0.2) \n"
        << "\t\t\t\tax.add_artist(beam) \n"
        << "\t\t\t\tax.text(3.6+bmaj/1.8,1.4,'Beam',transform=ax.transAxes,fontsize=11, ha='left',va='center') \n"
        << "\t\t\tax = plt.subplot(grid[j][1,i]) \n"
        << "\t\t\tax.tick_params(axis='both',which='both',bottom=True,top=True,labelbottom=False,labelleft=False) \n"
        << "\t\t\tax.imshow(z_mod,origin='lower',cmap = mpl.cm.Greys,norm=norm,aspect='auto',interpolation='none') \n"
        << "\t\t\tax.contour(z_mod,v,origin='lower',linewidths=0.7,colors='#B22222') \n"
        << "\t\t\tax.plot(xcen,ycen,'x',color='#0FB05A',markersize=7,mew=2) \n"
        << "\t\t\tif (i==0 and j==2): \n"
        << "\t\t\t\tclab = r'Contour levels at 2$^n \\, c_{min}$, where $c_{min}$ = %s " << in->Head().Bunit() << " and n = 0,1,..,8 '%cont \n"
        << "\t\t\t\tax.text(0.01,-0.16,clab,transform=ax.transAxes,fontsize=11, ha='left',va='center') \n"
        << "\t\t\tnum = num+1 \n"
        << std::endl
        << "\toutfile = '%s_chanmaps'%gname \n"
        << "\tif ('azim' in files_mod[k]): outfile += '_azim' \n"
        << "\telif ('local' in files_mod[k]): outfile += '_local' \n"
        << "\tfig.savefig(outfolder+outfile+'.pdf', orientation = 'portrait', format = 'pdf') \n"
        << "\timage_mod.close() \n\n"
        << "image.close() \n";

    pyf.close();

    /////////////////////////////////////////////////////////////////////////
    /// Script to plot position-velocity slices of model vs data
    /////////////////////////////////////////////////////////////////////////
    float zmin_wcs = AlltoVel(in->getZphys(zmin-0.5),in->Head());
    float zmax_wcs = AlltoVel(in->getZphys(zmax+0.5),in->Head());
    float pa_av = findMedian(&outr->phi[0],outr->nr);
    float pa_min = pa_av+90<360 ? pa_av+90 : pa_av-90;
    //if (zmin_wcs>zmax_wcs) std::swap(zmin_wcs,zmax_wcs);
    bool reverse = (pa_av>=45 && pa_av<225);
    //if (cdelt3_kms<0) reverse = !reverse;


    scriptnames.push_back("plot_pvs.py");
    pyf.open((in->pars().getOutfolder()+"plotscripts/"+scriptnames[2]).c_str());

    pyf << "#################################################################################\n"
        << "#### This script writes a plot of position-velocity slices of model and data ####\n"
        << "#################################################################################\n"
        << "import numpy as np \n"
        << "import os \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
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
        << "twostage = " << par.TWOSTAGE << " \n"
        << "plotmask = " << par.PLOTMASK << " \n"
        << "zmin, zmax = " << zmin << ", " << zmax << std::endl
        << std::endl
        << "if twostage: rad,vrot,inc,pa,vsys = np.genfromtxt(outfolder+'rings_final2.txt',usecols=(1,2,4,5,11),unpack=True) \n"
        << "else: rad,vrot,inc,pa,vsys = np.genfromtxt(outfolder+'rings_final1.txt',usecols=(1,2,4,5,11),unpack=True) \n"
        << std::endl
        << "files_pva_mod, files_pvb_mod = [], [] \n"
        << "for thisFile in sorted(os.listdir(outfolder+'pvs/')): \n"
        << "\tif 'pv_a_azim.fits' in thisFile: files_pva_mod.append(thisFile) \n"
        << "\tif 'pv_b_azim.fits' in thisFile: files_pvb_mod.append(thisFile) \n"
        << "\tif 'pv_a_local.fits' in thisFile: files_pva_mod.append(thisFile) \n"
        << "\tif 'pv_b_local.fits' in thisFile: files_pvb_mod.append(thisFile) \n"
        << std::endl
        << "image_maj     = fits.open(outfolder+'pvs/'+gname+'_pv_a.fits') \n"
        << "image_min     = fits.open(outfolder+'pvs/'+gname+'_pv_b.fits') \n"
        << "image_mas_maj = fits.open(outfolder+'pvs/'+gname+'mask_pv_a.fits') \n"
        << "image_mas_min = fits.open(outfolder+'pvs/'+gname+'mask_pv_b.fits') \n"
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
        << "\timage_mod_maj = fits.open(outfolder+'pvs/'+files_pva_mod[k]) \n"
        << "\timage_mod_min = fits.open(outfolder+'pvs/'+files_pvb_mod[k]) \n"
        << "\tdata_mod_maj = image_mod_maj[0].data[zmin:zmax+1,int(xminpv[0]):int(xmaxpv[0])+1] \n"
        << "\tdata_mod_min = image_mod_min[0].data[zmin:zmax+1,int(xminpv[1]):int(xmaxpv[1])+1] \n"
        << "\ttoplot = [[data_maj,data_min],[data_mod_maj,data_mod_min],[data_mas_maj,data_mas_min]] \n"
        << std::endl
        << "\tfig = plt.figure(figsize=(10,10), dpi=150) \n"
        << "\tx_len, y_len, y_sep = 0.6, 0.42, 0.08 \n"
        << "\tax, bottom_corner = [], [0.1,0.7] \n"
        << "\tfor i in range (2): \n"
        << "\t\tbottom_corner[0], axcol = 0.1, [] \n"
        << "\t\tax.append(fig.add_axes([bottom_corner[0],bottom_corner[1],x_len,y_len])) \n"
        << "\t\tbottom_corner[1]-=(y_len+y_sep) \n"
        << std::endl
        << "\tfor i in range (2): \n"
        << "\t\taxis = ax[i] \n"
        << "\t\taxis.tick_params(which='major',length=8, labelsize=labsize) \n"
        << "\t\taxis.set_xlabel('Offset (arcsec)',fontsize=labsize+2) \n"
        << "\t\taxis.set_ylabel(r'$\\mathrm{\\Delta V_{LOS}}$ (km/s)',fontsize=labsize+2) \n"
        << "\t\taxis.text(1, 1.02,palab[i],ha='right',transform=axis.transAxes,fontsize=labsize+4) \n"
        << "\t\taxis2 = axis.twinx() \n"
        << "\t\taxis2.set_xlim([ext[i][0],ext[i][1]]) \n"
        << "\t\taxis2.set_ylim([ext[i][2]+vsys_m,ext[i][3]+vsys_m]) \n"
        << "\t\taxis2.tick_params(which='major',length=8, labelsize=labsize) \n"
        << "\t\taxis2.set_ylabel(r'$\\mathrm{V_{LOS}}$ (km/s)',fontsize=labsize+2) \n"
        << "\t\taxis.imshow(toplot[0][i],origin='lower',cmap = mpl.cm.Greys,norm=norm,extent=ext[i],aspect='auto') \n"
        << "\t\taxis.contour(toplot[0][i],v,origin='lower',linewidths=0.7,colors='#00008B',extent=ext[i]) \n"
        << "\t\taxis.contour(toplot[0][i],v_neg,origin='lower',linewidths=0.1,colors='gray',extent=ext[i]) \n"
        << "\t\taxis.contour(toplot[1][i],v,origin='lower',linewidths=1,colors='#B22222',extent=ext[i]) \n"
        << "\t\taxis.axhline(y=0,color='black') \n"
        << "\t\taxis.axvline(x=0,color='black') \n"
        << "\t\taxis.grid(color='gray', linestyle='--', linewidth=0.3) \n"
        << "\t\tif plotmask: \n"
        << "\t\t\taxis.contour(toplot[2][i],levels=[0],origin='lower',linewidths=2,colors='k',extent=ext[i]) \n"
        << "\t\tif i==0 : \n"
        << "\t\t\taxis2.plot(radius,vlos,'yo') \n"
        << "\t\t\taxis.text(0, 1.1, gname, transform=axis.transAxes,fontsize=22) \n"
        << std::endl
        << "\toutfile = '%s_pv'%gname \n"
        << "\tif ('azim' in files_pva_mod[k]): outfile += '_azim' \n"
        << "\telif ('local' in files_pva_mod[k]): outfile += '_local' \n"
        << "\tfig.savefig(outfolder+outfile+'.pdf', bbox_inches='tight') \n"
        << "\timage_mod_maj.close() \n"
        << "\timage_mod_min.close() \n"
        << std::endl
        << "image_maj.close() \n"
        << "image_min.close() \n";

    pyf.close();


    /////////////////////////////////////////////////////////////////////////
    /// Script to plot kinematic maps of model vs data
    /////////////////////////////////////////////////////////////////////////
    scriptnames.push_back("plot_kinmaps.py");
    pyf.open((in->pars().getOutfolder()+"plotscripts/"+scriptnames[3]).c_str());

    pyf << "#######################################################################\n"
        << "#### This script writes a plot of kinematic maps of model and data ####\n"
        << "#######################################################################\n"
        << "import numpy as np \n"
        << "import os \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "from matplotlib.colorbar import ColorbarBase \n"
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
        << "twostage = " << par.TWOSTAGE << " \n"
        << "xmin, xmax = " << xmin << ", " << xmax << std::endl
        << "ymin, ymax = " << ymin << ", " << ymax << std::endl
        << std::endl
        << "if twostage: rad,inc,pa,xpos,ypos,vsys = np.genfromtxt(outfolder+'rings_final2.txt',usecols=(1,4,5,9,10,11),unpack=True) \n"
        << "else: rad,inc,pa,xpos,ypos,vsys = np.genfromtxt(outfolder+'rings_final1.txt',usecols=(1,4,5,9,10,11),unpack=True) \n"
        << "xcen_m,ycen_m,inc_m,pa_m,vsys_m=np.nanmean((xpos,ypos,inc,pa,vsys),axis=1) \n"
        << "xcen, ycen = xcen_m-xmin, ycen_m-ymin \n"
        << std::endl
        << "# Opening maps and retrieving intensity map units\n"
        << "f0 = fits.open(outfolder+'/maps/'+gname+'_0mom.fits') \n"
        << "f1 = fits.open(outfolder+'/maps/'+gname+'_1mom.fits') \n"
        << "f2 = fits.open(outfolder+'/maps/'+gname+'_2mom.fits') \n"
        << "bunit = f0[0].header['BUNIT'] \n"
        << "bunit = bunit.replace(' ', '').lower() \n"
        << "# Now plotting moment maps \n"
        << "mom0 = f0[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "mom1 = f1[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "mom2 = f2[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "maskmap = np.copy(mom1) \n"
        << "maskmap[mom1==mom1] = 1 \n"
        << std::endl
        << "files_mod0, files_mod1, files_mod2 = [], [], [] \n"
        << "for thisFile in sorted(os.listdir(outfolder+'/maps/')): \n"
        << "\tif 'azim_0mom.fits' in thisFile: files_mod0.append(thisFile) \n"
        << "\tif 'azim_1mom.fits' in thisFile: files_mod1.append(thisFile) \n"
        << "\tif 'azim_2mom.fits' in thisFile: files_mod2.append(thisFile) \n"
        << "\tif 'local_0mom.fits' in thisFile: files_mod0.append(thisFile) \n"
        << "\tif 'local_1mom.fits' in thisFile: files_mod1.append(thisFile) \n"
        << "\tif 'local_2mom.fits' in thisFile: files_mod2.append(thisFile) \n"
        << std::endl
        << "cmaps = [plt.get_cmap('Spectral_r'),plt.get_cmap('RdBu_r',25),plt.get_cmap('PuOr_r')] \n"
        << "barlab = ['Intensity ('+bunit+')', r'V$_\\mathrm{LOS}$ (km/s)', r'$\\sigma$ (km/s)'] \n"
        << "barlab2 = [r'I$_\\mathrm{res}$ ('+bunit+')', r'V$_\\mathrm{res}$ (km/s)', r'$\\sigma_\\mathrm{res}$ (km/s)'] \n"
        << "titles = ['DATA', 'MODEL','RESIDUAL'] \n"
        << "mapname = ['INTENSITY', 'VELOCITY', 'DISPERSION'] \n"
        << "x = np.arange(0,xmax-xmin,0.1) \n"
        << "y = np.tan(np.radians(pa_m-90))*(x-xcen)+ycen \n"
        << "ext = [0,xmax-xmin,0, ymax-ymin] \n"
        << "rad_pix = rad/"<< in->Head().PixScale()*arcconv << std::endl
        << "try: nr = len(rad_pix) \n"
        << "except: nr = 1 \n"
        << "interval = PercentileInterval(99.5) \n"
        << std::endl
        << "for k in range (len(files_mod0)): \n"
        << "\tmom0_mod = fits.open(outfolder+'/maps/'+files_mod0[k])[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "\tmom1_mod = fits.open(outfolder+'/maps/'+files_mod1[k])[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "\tmom2_mod = fits.open(outfolder+'/maps/'+files_mod2[k])[0].data[ymin:ymax+1,xmin:xmax+1] \n"
        << "\tto_plot = [[mom0,mom1-vsys_m,mom2],[mom0_mod,mom1_mod-vsys_m,mom2_mod],[mom0-mom0_mod,mom1-mom1_mod,mom2-mom2_mod]] \n"
        << std::endl
        << "\tfig=plt.figure(figsize=(11,11), dpi=150) \n"
        << "\tnrows, ncols = 3, 3 \n"
        << "\tx_len, y_len = 0.2, 0.2 \n"
        << "\tx_sep, y_sep = 0.00,0.08 \n"
        << "\tax, bottom_corner = [], [0.1,0.7] \n"
        << "\tfor i in range (nrows): \n"
        << "\t\tbottom_corner[0], axcol = 0.1, [] \n"
        << "\t\tfor j in range (ncols): \n"
        << "\t\t\taxcol.append(fig.add_axes([bottom_corner[0],bottom_corner[1],x_len,y_len])) \n"
        << "\t\t\tbottom_corner[0]+=x_len+x_sep \n"
        << "\t\tax.append(axcol) \n"
        << "\t\tbottom_corner[1]-=(y_len+y_sep) \n"
        << std::endl
        << "\tfor i in range (nrows): \n"
        << "\t\tcmap = copy(cmaps[i]) \n"
        << "\t\tcmap.set_bad('w',1.) \n"
        << "\t\tvmin, vmax = interval.get_limits(to_plot[1][i]) \n"
        << "\t\tvmin, vmax = (-1.1*np.nanmax(vmax),1.1*np.nanmax(vmax)) if i==1 else (vmin,vmax) \n"
        << "\t\tnorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax) \n"
        << "\t\tcbax = fig.add_axes([ax[i][0].get_position().x0,ax[i][0].get_position().y0-0.025,2*x_len-0.003,0.02]) \n"
        << "\t\tcb1 = ColorbarBase(cbax, orientation='horizontal', cmap=cmap, norm=norm) \n"
        << "\t\tcb1.set_label(barlab[i],fontsize=13) \n"
        << "\t\tfor j in range (ncols): \n"
        << "\t\t\taxis = ax[i][j] \n"
        << "\t\t\taxis.tick_params(labelbottom=False,labelleft=False,right=True,top=True) \n"
        << "\t\t\taxis.set_xlim(ext[0],ext[1]) \n"
        << "\t\t\taxis.set_ylim(ext[2],ext[3]) \n"
        << "\t\t\tif j==2: \n"
        << "\t\t\t\tvmax = np.nanmax(interval.get_limits(to_plot[j][i])) \n"
        << "\t\t\t\tnorm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax) \n"
        << "\t\t\taxis.imshow(to_plot[j][i]*maskmap,origin='lower',cmap=cmap,norm=norm,aspect='auto',extent=ext,interpolation='nearest') \n"
        << "\t\t\taxis.plot(xcen,ycen,'x',color='#000000',markersize=7,mew=1.5) \n"
        << std::endl
        << "\t\t\tif i==0: \n"
        << "\t\t\t\taxis.text(0.5,1.05,titles[j],ha='center',transform=axis.transAxes,fontsize=15) \n"
        << "\t\t\t\taxis.plot(x,y,'--',color='k',linewidth=1) \n"
        << "\t\t\t\tif j!=2 and nr>3:  \n"
        << "\t\t\t\t\taxmaj = rad_pix[-1] \n"
        << "\t\t\t\t\taxmin = axmaj*np.cos(np.radians(inc_m))  \n"
        << "\t\t\t\t\tposa = np.radians(pa_m-90)  \n"
        << "\t\t\t\t\tt = np.linspace(0,2*np.pi,100)  \n"
        << "\t\t\t\t\txt = xcen+axmaj*np.cos(posa)*np.cos(t)-axmin*np.sin(posa)*np.sin(t)  \n"
        << "\t\t\t\t\tyt = ycen+axmaj*np.sin(posa)*np.cos(t)+axmin*np.cos(posa)*np.sin(t)  \n"
        << "\t\t\t\t\taxis.plot(xt,yt,'-',c='k',lw=0.8)  \n"
        << "\t\t\telif i==1: \n"
        << "\t\t\t\taxis.plot(x,y,'--',color='k',linewidth=1) \n"
        << "\t\t\t\tif nr<10: \n"
        << "\t\t\t\t\tx_pix = rad_pix*np.cos(np.radians(pa_m-90)) \n"
        << "\t\t\t\t\ty_pix = rad_pix*np.sin(np.radians(pa_m-90)) \n"
        << "\t\t\t\t\taxis.scatter(x_pix+xcen,y_pix+ycen,c='grey',s=12) \n"
        << "\t\t\t\t\taxis.scatter(xcen-x_pix,ycen-y_pix,c='grey',s=12) \n"
        << "\t\t\t\tif nr>5 and not all(np.diff(pa)==0): \n"
        << "\t\t\t\t\tx_pix = rad_pix*np.cos(np.radians(pa-90)) \n"
        << "\t\t\t\t\ty_pix = rad_pix*np.sin(np.radians(pa-90)) \n"
        << "\t\t\t\t\taxis.plot(xcen-x_pix,ycen-y_pix,'-',color='grey',lw=1) \n"
        << "\t\t\t\t\taxis.plot(x_pix+xcen,y_pix+ycen,'-',color='grey',lw=1) \n"
        << "\t\t\t\tif j!=2: axis.contour(to_plot[j][i]*maskmap,levels=[0],colors='green',origin='lower',extent=ext) \n"
        << "\t\t\tif j==0: axis.text(-0.12,0.5,mapname[i],va='center',rotation=90,transform=axis.transAxes,fontsize=15) \n\n"
        << "\t\tcbax = fig.add_axes([ax[i][2].get_position().x0+0.003,ax[i][2].get_position().y0-0.025,x_len-0.003,0.02]) \n"
        << "\t\tcb2 = ColorbarBase(cbax, orientation='horizontal', cmap=cmap, norm=norm) \n"
        << "\t\tcb2.set_label(barlab2[i],fontsize=13) \n"
        << "\t\tcb2.ax.locator_params(nbins=3) \n"
        << "\t\tfor c in [cb1,cb2]: \n"
        << "\t\t\tc.solids.set_edgecolor('face') \n"
        << "\t\t\tc.outline.set_linewidth(0) \n"
        << std::endl
        << "\toutfile = '%s_maps'%gname \n"
        << "\tif ('azim' in files_mod0[k]): outfile += '_azim' \n"
        << "\telif ('local' in files_mod0[k]): outfile += '_local' \n"
        << "\tfig.savefig(outfolder+outfile+'.pdf', bbox_inches = 'tight') \n\n";

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
        << "\tif '.py' in f and f!='plot_all.py': \n"
        << "\t\tcmd += 'python \"%s/%s\" & '%(scriptdir,f) \n"
        << std::endl
        << "os.system(cmd[:-2]) \n";

    pyf.close();


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
            << setw(m) << inr->dens[i]/1.E20
            << setw(m) << inr->xpos[i]
            << setw(m) << inr->ypos[i]
            << setw(m) << inr->vsys[i] << endl;
    }
    initout.close();

}
template void Galfit<float>::printInitial(Rings<float>*,std::string);
template void Galfit<double>::printInitial(Rings<double>*,std::string);


template <class T>
void Galfit<T>::DensityProfile (T *surf_dens, int *count) {

    // -STRONZATEEE-----------------------------------------------------------

    if (!in->getIsSearched()) in->search();
    Detection<T> *obj = in->pObject(0);
    obj->calcFluxes(obj->getPixelSet(in->Array(), in->AxisDim()));
    obj->calcWCSparams(in->Head());
    obj->calcIntegFlux(in->DimZ(), obj->getPixelSet(in->Array(), in->AxisDim()), in->Head(), false, in->pars().getFluxConvert());
    obj->setMass(2.365E5*obj->getIntegFlux()*distance*distance);
    T *surf_bright_faceon = new T[outr->nr];
    T *mass_surf_dens = new T[outr->nr];
    float surf_bright_tot = 0;
    float Area = fabs(in->Head().Cdelt(0))*arcsconv(in->Head().Cunit(0))*fabs(in->Head().Cdelt(1))*arcsconv(in->Head().Cunit(1));
    for (int i=0;i<outr->nr;i++) {
        surf_bright_faceon[i] = surf_dens[i]/Area;
        surf_bright_faceon[i] *= cos(outr->inc[i]*M_PI/180.)/count[i];
        surf_bright_tot += surf_dens[i];
        //surf_dens[i]/=count[i];
        //cout << surf_bright_tot << endl;
    }
    //cout << scientific<< obj->getMass() << std::endl << surf_bright_tot;

    std::ofstream fileout((in->pars().getOutfolder()+"surface_dens.txt").c_str());

    int m=18;
    fileout << "#" << setw(m-1) << "RADIUS" << setw(m) << "SURFBRIGHT" << setw(m) << "SURFDENS" << std::endl;
    fileout << "#" << setw(m-1) << "(arcsec)" << setw(m) << "("+in->Head().Bunit()+"/arc2)" << setw(m) << "(Msun/pc2)" << std::endl;

    for (int i=0;i<outr->nr;i++) {
        mass_surf_dens[i]=obj->getMass()*surf_bright_faceon[i]/surf_bright_tot/(4.848*4.848*distance*distance);
        fileout << setw(m) << outr->radii[i] << setw(m) << surf_bright_faceon[i] << setw(m) << mass_surf_dens[i] << std::endl;
    }

    delete [] surf_bright_faceon;
    delete [] mass_surf_dens;
}
template void Galfit<float>::DensityProfile (float*, int *);
template void Galfit<double>::DensityProfile (double*, int *);


}

