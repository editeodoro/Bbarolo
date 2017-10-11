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
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <Arrays/cube.hh>
#include <Arrays/image.hh>
#include <Utilities/galfit.hh>
#include <Utilities/galmod.hh>
#include <Utilities/utils.hh>
#include <Utilities/gnuplot.hh>
#include <Utilities/moment.hh>
#include <Utilities/ellprof.hh>

//#ifdef HAVE_PYTHON
//    #include <Python.h>
//#endif

#define VROT  0
#define VDISP 1
#define DENS  2
#define Z0    3
#define INC   4
#define PA    5
#define XPOS  6
#define YPOS  7
#define VSYS  8
#define VRAD  9
#define MAXPAR 10

namespace Model {


template <class T>
void Galfit<T>::writeModel (std::string normtype) {

    bool verb = in->pars().isVerbose();
    in->pars().setVerbosity(false);

    if (verb) std::cout << " Preparing a bunch of cool outputs..." << std::endl;

    std::string outfold = in->pars().getOutfolder();
    std::string object = in->Head().Name();

    // Get the total intensity, velocity and dispersion maps of data
    mkdirp((outfold+"maps/").c_str());
    if (verb) std::cout << "    Extracting modest intensity map..." << std::flush;
    MomentMap<T> *totalmap = new MomentMap<T>;
    totalmap->input(in);
    totalmap->SumMap(true);
    totalmap->fitswrite_2d((outfold+"maps/"+object+"_0mom.fits").c_str());
    MomentMap<T> *field = new MomentMap<T>;
    field->input(in);
    field->FirstMoment(true);
    field->fitswrite_2d((outfold+"maps/"+object+"_1mom.fits").c_str());
    field->SecondMoment(true);
    field->fitswrite_2d((outfold+"maps/"+object+"_2mom.fits").c_str());
    if (verb) std::cout << " Done." << std::endl;

    // Calculate the total flux inside last ring in data
    T *ringreg = getFinalRingsRegion();
    float totflux_data=0, totflux_model=0;
    for (size_t i=0; i<in->DimX()*in->DimY(); i++) {
        if (!isNaN(ringreg[i]) && !isNaN(totalmap->Array(i))) {
            totflux_data += totalmap->Array(i);
        }
    }

    float mass = 2.365E5*FluxtoJy(totflux_data,in->Head())*fabs(DeltaVel<float>(in->Head()))*distance*distance;

    // Calculate radial profile along the output rings
    if (verb) std::cout << "    Deriving trustful radial profile..." << std::flush;
    T meanPA = findMean(&outr->phi[0], outr->nr);
    int nseg = 1;
    float segments[4] = {0, 360., 0., 0};
    if (in->pars().getSIDE()!="A") {
        nseg = 2;
        segments[2]=-90;
        segments[3]=90;
    }
    else if(in->pars().getSIDE()=="R") {
        nseg = 2;
        segments[2]=90;
        segments[3]=-90;
    }
    if (meanPA>180) std::swap(segments[2], segments[3]);

    Tasks::Ellprof<T> ell(totalmap,outr,nseg,segments);
    ell.setOptions(mass,distance);  //To set the mass and the distance
    ell.RadialProfile();
    std::string dens_out = outfold+"densprof.txt";
    std::ofstream fileo;
    fileo.open(dens_out.c_str());
    ell.printProfile(fileo,nseg-1);
    fileo.close();
    //ell.printProfile(std::cout);
    delete totalmap;

    if (normtype=="AZIM" || normtype=="BOTH") {
        double profmin=FLT_MAX;
        for (size_t i=0; i<outr->nr; i++) {
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
        for (size_t i=0; i<outr->nr; i++) {
            outr->dens[i]=factor*fabs(ell.getMean(i))*1E20;
            if (outr->dens[i]==0) outr->dens[i]=profmin*1E20;
        }
    }
    if (verb) std::cout << " Done." << std::endl;

    if (verb) std::cout << "    Calculating the very last model..." << std::flush;
    Model::Galmod<T> *mod = getModel();
    mod->Out()->Head().setMinMax(0.,0.);
    mod->Out()->Head().setName(object+"mod");

    T *outarray = mod->Out()->Array();

    if (normtype=="AZIM" || normtype=="BOTH") {

        // The final model has been build from the azimuthal profile calculated before,
        // thus just need to rescale the model to the total flux of data inside last ring.

        // Calculate total flux of model within last ring
        for (size_t i=0; i<in->DimX()*in->DimY(); i++) {
            if (!isNaN(ringreg[i])) {
                for (size_t z=0; z<in->DimZ(); z++)
                    totflux_model += outarray[i+z*in->DimY()*in->DimX()];
            }
            else {
                 // for (size_t z=0; z<in->DimZ(); z++)
                   //      outarray[i+z*in->DimY()*in->DimX()]=0;
              }
        }

        double factor = totflux_data/totflux_model;
        for (int i=0; i<in->NumPix(); i++) outarray[i] *= factor;
        if (verb) std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing useful azimuthally-normalized model..." << std::flush;
        std::string mfile = outfold+object+"mod_azim.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());
        writePVs(mod->Out(),"_azim");
        if (verb) std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing shameful kinematic maps..." << std::flush;
        MomentMap<T> *map = new MomentMap<T>;
        map->input(mod->Out());
        map->SumMap(false);
        map->fitswrite_2d((outfold+"maps/"+object+"_azim_0mom.fits").c_str());
        map->FirstMoment(false);
        map->fitswrite_2d((outfold+"maps/"+object+"_azim_1mom.fits").c_str());
        map->SecondMoment(false);
        map->fitswrite_2d((outfold+"maps/"+object+"_azim_2mom.fits").c_str());
        delete map;
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

        if (verb) std::cout << "    Writing awesome locally-normalized model..." << std::flush;
        std::string mfile = outfold+object+"mod_local.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());

        writePVs(mod->Out(),"_local");
        if (verb) std::cout << " Done." << std::endl;

        if (verb) {
            std::cout << "    Writing " << std::flush;
            if (normtype=="BOTH") std::cout << "even more " << std::flush;
            std::cout << "shameful kinematic maps..." << std::flush;
        }
        MomentMap<T> *map = new MomentMap<T>;
        map->input(mod->Out());
        map->SumMap(false);
        map->fitswrite_2d((outfold+"maps/"+object+"_local_0mom.fits").c_str());
        map->FirstMoment(false);
        map->fitswrite_2d((outfold+"maps/"+object+"_local_1mom.fits").c_str());
        map->SecondMoment(false);
        map->fitswrite_2d((outfold+"maps/"+object+"_local_2mom.fits").c_str());
        delete map;
        if (verb) std::cout << " Done." << std::endl;
    }

    if (normtype=="NONE") {
        
        // Calculate total flux of model within last ring
        for (size_t i=0; i<in->DimX()*in->DimY(); i++) {
            if (!isNaN(ringreg[i])) {
                for (size_t z=0; z<in->DimZ(); z++)
                    totflux_model += outarray[i+z*in->DimY()*in->DimX()];
            }
        }

        double factor = totflux_data/totflux_model;
        if (totflux_data==0) factor=1;
        for (int i=0; i<in->NumPix(); i++) outarray[i] *= factor;
        if (verb) std::cout << " Done." << std::endl;
        
        if (verb) std::cout << "    Writing model..." << std::flush;
        std::string mfile = outfold+object+"mod_nonorm.fits";
        mod->Out()->fitswrite_3d(mfile.c_str());
        writePVs(mod->Out(),"_nonorm");
        if (verb) std::cout << " Done." << std::endl;

        if (verb) std::cout << "    Writing shameful kinematic maps..." << std::flush;
        MomentMap<T> *map = new MomentMap<T>;
        map->input(mod->Out());
        map->SumMap(false);
        map->fitswrite_2d((outfold+"maps/"+object+"_nonorm_0mom.fits").c_str());
        map->FirstMoment(false);
        map->fitswrite_2d((outfold+"maps/"+object+"_nonorm_1mom.fits").c_str());
        map->SecondMoment(false);
        map->fitswrite_2d((outfold+"maps/"+object+"_nonorm_2mom.fits").c_str());
        delete map;
        if (verb) std::cout << " Done." << std::endl;
        
/*////////        TO BE REMOVED ////////////////////////////////////////////////////////
        map = new MomentMap<T>;
        map->input(mod->Out());
        map->SumMap(false);
        Tasks::Ellprof<T> ell(map,outr,nseg,segments);
        ell.setOptions(mass,distance);  //To set the mass and the distance
        ell.RadialProfile();
        std::string dens_out = outfold+"densprofmod.txt";
        std::ofstream fileo;
        fileo.open(dens_out.c_str());
        ell.printProfile(fileo,nseg-1);
        fileo.close();
//////////////////////////////////////////////////////////////////////////////////////////*/
        
        
    }

    // Now plotting everything
    if (verb) std::cout << "    Writing creative plots..." << std::flush;
    int ret = plotParam();
    if (verb) {
        if (ret==0) std::cout << " Done." << std::endl;
        else std::cout << " Something went wrong! Check pyscript.py in the output folder." << std::endl;
    }

    in->pars().setVerbosity(verb);

    if (verb) std::cout << " All done!" << std::endl;

    delete mod;

}
template void Galfit<float>::writeModel(std::string);
template void Galfit<double>::writeModel(std::string);


template <class T>
void Galfit<T>::writePVs(Cube<T> *mod, std::string suffix) {

    // Extracts PVs along major and minor axis and write them in fits.
    T meanPA = findMean(&outr->phi[0], outr->nr);
    T meanXpos = findMean(&outr->xpos[0], outr->nr);
    T meanYpos = findMean(&outr->ypos[0], outr->nr);
    T meanPAp90= meanPA+90<360 ? meanPA+90 : meanPA-90;
    std::string outfold = in->pars().getOutfolder();
    std::string object = in->Head().Name();

    Image2D<T> *pv_max = PositionVelocity(in,meanXpos,meanYpos,meanPA);
    std::string mfile = outfold+object+"_pv_a.fits";
    pv_max->fitswrite_2d(mfile.c_str());
    Image2D<T> *pv_min = PositionVelocity(in,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+object+"_pv_b.fits";
    pv_min->fitswrite_2d(mfile.c_str());

    Image2D<T> *pv_max_m = PositionVelocity(mod,meanXpos,meanYpos,meanPA);
    mfile = outfold+object+"mod_pv_a"+suffix+".fits";
    pv_max_m->fitswrite_2d(mfile.c_str());
    Image2D<T> *pv_min_m = PositionVelocity(mod,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+object+"mod_pv_b"+suffix+".fits";
    pv_min_m->fitswrite_2d(mfile.c_str());

#ifndef HAVE_PYTHON
#ifdef HAVE_GNUPLOT
    plotPVs_Gnuplot(pv_max,pv_min,pv_max_m,pv_min_m);
#endif
#endif

    delete pv_max; delete pv_min;
    delete pv_max_m; delete pv_min_m;
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
    for (int i=in->pars().getStartRad(); i<outr->nr; i++) {
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
    if (in->pars().getFlagUserThreshold()) sig=in->pars().getThreshold();
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
    mfile = "load '"+outfold+"gnuscript.gnu'";
    if (!in->pars().getflagGalMod()) gp.commandln(mfile.c_str());
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
        if (flagErrors && mpar[j]) {
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

    std::ofstream gnu;
    std::string mfile = outfold+"gnuscript.gnu";
    gnu.open(mfile.c_str());

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

    if (in->pars().getTwoStage()) {
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
        << "plot '" << outfold << "ringlog1.txt' ";

    if (nc[VROT]>=0) {
        gnu << "u 2:3:($3+$"+to_string(nc[VROT])+"):($3+$"+to_string(nc[VROT]+1)+") w errorbars ls 1, '"
            << outfold <<"ringlog1.txt' u 2:3 w lp ls 1";
    }
    else gnu << "u 2:3 w lp ls 1";

    if (in->pars().getTwoStage()) {
        gnu << ", '" << outfold << "ringlog2.txt' ";
        if (flagErrors && mpar[VROT]) {
        gnu << "u 2:3:($3+$13):($3+$14) w errorbars ls 2, '"
            << outfold <<"ringlog2.txt' u 2:3 w lp ls 2";
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
        << "plot '" << outfold << "ringlog1.txt' ";

    if (nc[INC]>=0) {
        gnu << "u 2:5:($5+$"+to_string(nc[INC])+"):($5+$"+to_string(nc[INC]+1)+") w errorbars ls 1, '"
            << outfold << "ringlog1.txt' u 2:5 w lp ls 1";
    }
    else gnu << "u 2:5 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:5 w lp ls 2";


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
        << "plot '" << outfold << "ringlog1.txt' ";

    if (nc[PA]>=0) {
        gnu << "u 2:6:($6+$"+to_string(nc[PA])+"):($6+$"+to_string(nc[PA]+1)+") w errorbars ls 1, '"
            << outfold << "ringlog1.txt' u 2:6 w lp ls 1";;
    }
    else gnu << "u 2:6 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:6 w lp ls 2";

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
        << "plot '"<<in->pars().getOutfolder()<<"ringlog1.txt' ";

    if (nc[VDISP]>=0) {
        gnu << "u 2:4:($4+$"+to_string(nc[VDISP])+"):($4+$"+to_string(nc[VDISP]+1)+") w errorbars ls 1, '"
            << outfold <<"ringlog1.txt' u 2:4 w lp ls 1";
    }
    else gnu << "u 2:4 w lp ls 1";

    if (in->pars().getTwoStage()) {
        gnu << ", '" << outfold << "ringlog2.txt' ";
        if (flagErrors && mpar[VDISP]) {
            gnu << "u 2:4:($3+$15):($3+$16) w errorbars ls 2, '"
                << outfold <<"ringlog2.txt' u 2:4 w lp ls 2";
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
        << "plot '" << outfold << "ringlog1.txt' ";

    if (nc[VSYS]>=0) {
        gnu << "u 2:12:($12+$"+to_string(nc[VSYS])+"):($12+$"+to_string(nc[VSYS]+1)+") w errorbars ls 1, '"
            << outfold << "ringlog1.txt' u 2:12 w lp ls 1";;
    }
    else gnu << "u 2:12 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:12 w lp ls 2";
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
        << "plot '" << outfold << "ringlog1.txt'";

    if (nc[Z0]>=0) {
        gnu << "u 2:8:($8+$"+to_string(nc[Z0])+"):($8+$"+to_string(nc[Z0]+1)+") w errorbars ls 1, '"
            << outfold << "ringlog1.txt' u 2:8 w lp ls 1";
    }
    else gnu << "u 2:8 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:8 w lp ls 2";
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
        << "plot '" << outfold << "ringlog1.txt' ";

    if (nc[XPOS]>=0) {
        gnu << "u 2:10:($10+$"+to_string(nc[XPOS])+"):($10+$"+to_string(nc[XPOS]+1)+") w errorbars ls 1, '"
            << outfold << "ringlog1.txt' u 2:10 w lp ls 1";;
    }
    else gnu << "u 2:10 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:10 w lp ls 2";
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
        << "plot '" << outfold << "ringlog1.txt' ";

    if (nc[YPOS]>=0) {
        gnu << "u 2:11:($11+$"+to_string(nc[YPOS])+"):($11+$"+to_string(nc[YPOS]+1)+") w errorbars ls 1, '"
            << outfold << "ringlog1.txt' u 2:11 w lp ls 1";;
    }
    else gnu << "u 2:11 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:11 w lp ls 2";
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
            << "plot '"<<in->pars().getOutfolder()<<"ringlog1.txt' u 2:8 w lp ls 1";

        if (in->pars().getTwoStage())
            gnu << ", '" << outfold << "ringlog2.txt' u 2:8 w lp ls 2";
        gnu << endl;
    }

    gnu << "unset multiplot; reset" << endl;
    gnu.close();

#ifndef HAVE_PYTHON
#ifdef HAVE_GNUPLOT
    Gnuplot gp;
    gp.begin();
    mfile = "load '"+outfold+"gnuscript.gnu'";
    gp.commandln(mfile.c_str());
#endif
#endif

}
template void Galfit<float>::plotPar_Gnuplot();
template void Galfit<double>::plotPar_Gnuplot();


template <class T>
int Galfit<T>::plotAll_Python() {

    /* This function creates and runs a python script for plotting:
     *   1) Channel maps
     *   2) Position-Velocity diagrams along major/minor axis
     *   3) Output paramenters
     *
     * It needs all output fitsfiles to be in the output directory!
     */

    std::string scriptname = "pyscript.py";
    std::ofstream py_file((in->pars().getOutfolder()+scriptname).c_str());

    float crpix3_kms = in->Head().Crpix(2);
    float cdelt3_kms = DeltaVel<float>(in->Head());
    float crval3_kms = AlltoVel(in->Head().Crval(2),in->Head());

    float cont = 0;
    int xpos = findMedian(&outr->xpos[0],outr->nr);
    int ypos = findMedian(&outr->ypos[0],outr->nr);
    int xmin=0, ymin=0, zmin=0, disp=0;
    int xmax=in->DimX()-1, ymax=in->DimY()-1, zmax=in->DimZ()-1;
    float vsys_av = findMedian(&outr->vsys[0],outr->nr);

    int free[nfree];
    for (int nm=0, k=0; nm<MAXPAR; nm++) {
        if (mpar[nm]) free[k++]=nm;
    }

    if (in->pars().getMASK()=="SEARCH") {
        if (in->pars().getFlagUserGrowthThreshold()) cont = in->pars().getGrowthThreshold();
        else if (in->pars().getFlagUserThreshold()) cont = in->pars().getThreshold();
        else {
            if (in->pars().getFlagGrowth()) cont = in->pars().getGrowthCut()*in->stat().getSpread();
            else cont = in->pars().getCut()*in->stat().getSpread();
        }
        Detection<T> *larg = in->LargestDetection();
        long ext[4] = {abs(xpos-lround(larg->getXmin()-2*in->Head().Bmaj()/in->Head().PixScale())),
                       abs(xpos-lround(larg->getXmax()+2*in->Head().Bmaj()/in->Head().PixScale())),
                       abs(ypos-lround(larg->getYmin()-2*in->Head().Bmaj()/in->Head().PixScale())),
                       abs(ypos-lround(larg->getYmax()+2*in->Head().Bmaj()/in->Head().PixScale()))};
        disp = *max_element(&ext[0],&ext[0]+4);
        zmin = larg->getZmin()-3;
        zmax = larg->getZmax()+3;
        if (zmin<0) zmin=0;
        if (zmax>=in->DimZ()) zmax=in->DimZ()-1;
    }
    else {
        if (!in->StatsDef()) in->setCubeStats();
        cont = 2.5*in->stat().getSpread();
        if (in->pars().getFlagUserThreshold()) cont=in->pars().getThreshold();
        std::vector<T> maxv(outr->nr);
        for (int i=0; i<outr->nr; i++) maxv[i]=outr->vrot[i]*sin(outr->inc[i]*M_PI/180.)+outr->vdisp[i];
       //float max_vrot = *max_element(&outr->vrot[0],&outr->vrot[0]+outr->nr);
        float max_v=*max_element(&maxv[0],&maxv[0]+outr->nr);
        disp = fabs((outr->radii.back()/arcconv+2*in->Head().Bmaj())/in->Head().PixScale());
        int z_vsys = (vsys_av-crval3_kms)/cdelt3_kms+crpix3_kms-1;
        //int disp_v = ceil((1.5*max_vrot)*sin(inc_av*M_PI/180.)/fabs(DeltaVel<float>(in->Head())));
        int disp_v = ceil((1.5*max_v)/fabs(DeltaVel<float>(in->Head())));
        zmin = z_vsys-2*disp_v>0 ? z_vsys-2*disp_v : 0;
        zmax = z_vsys+2*disp_v<in->DimZ() ? z_vsys+2*disp_v : in->DimZ()-1;
    }

    xmin = xpos-disp>=0 ? xpos-disp : 0;
    xmax = xpos+disp<in->DimX() ? xpos+disp : in->DimX()-1;
    ymin = ypos-disp>=0 ? ypos-disp : 0;
    ymax = ypos+disp<in->DimY() ? ypos+disp : in->DimY()-1;

    int *nc = getErrorColumns();

    py_file << "import numpy as np \n"
            << "import os \n"
            << "import matplotlib \n"
            << "import matplotlib.pyplot as plt \n"
            << "import matplotlib.gridspec as gridspec \n"
            << "from matplotlib.colorbar import ColorbarBase \n"
            << "from astropy.io import fits \n"
            << "from astropy.visualization import LinearStretch, PowerStretch \n"
            << "from astropy.visualization.mpl_normalize import ImageNormalize \n"
            << "from astropy.visualization import PercentileInterval \n\n";


    py_file << "# PARAMETERS: plotting the fit parameters \n"
            << "outfolder = '" << in->pars().getOutfolder() <<"' \n"
            << "file1 = outfolder+'ringlog1.txt' \n"
            << "file2 = outfolder+'ringlog2.txt' \n"
            << "filesb= outfolder+'densprof.txt' \n"
            << "twostage=False \n";

    if (in->pars().getTwoStage()) py_file << "twostage=True \n";

    py_file << "rad,vrot,disp,inc,pa,z0,xpos,ypos,vsys= np.genfromtxt(file1,skip_header=1,usecols=(1,2,3,4,5,7,9,10,11),unpack=True) \n"
            << "err1_l, err1_h = np.zeros(shape=(" << MAXPAR << ",len(rad))), np.zeros(shape=(" << MAXPAR << ",len(rad)))\n"
            << "color=color2='#B22222' \n"
            << "max_vrot,max_vdisp,max_inc,max_pa=np.max(vrot),np.max(disp),np.max(inc),np.max(pa) \n"
            << "max_z0,max_xpos,max_ypos,max_vsys=np.max(z0),np.max(xpos),np.max(ypos),np.max(vsys) \n"
            << "max_rad = 1.1*np.max(rad) \n";

    for (int i=0; i<MAXPAR; i++) {
        if (nc[i]>0) {
            py_file << "err1_l[" << i << "], err1_h[" << i << "] = np.genfromtxt(file1,skip_header=1,usecols=("
                    << nc[i] << "," << nc[i]+1 << "),unpack=True) \n";
        }
    }

    py_file << "\nif twostage==True: \n"
            << "\trad2, vrot2,disp2,inc2,pa2,z02,xpos2,ypos2,vsys2=np.genfromtxt(file2,skip_header=1,usecols=(1,2,3,4,5,7,9,10,11),unpack=True)\n"
            << "\terr2_l, err2_h = np.zeros(shape=(" << MAXPAR << ",len(rad2))), np.zeros(shape=(" << MAXPAR << ",len(rad2)))\n"
            << "\tcolor='#A0A0A0' \n"
            << "\tmax_vrot,max_vdisp,max_inc,max_pa=np.maximum(max_vrot,np.max(vrot2)),np.maximum(max_vdisp,np.max(disp2)),np.maximum(max_inc,np.max(inc2)),np.maximum(max_pa,np.max(pa2)) \n"
            << "\tmax_z0,max_xpos,max_ypos,max_vsys=np.maximum(max_z0,np.max(z02)),np.maximum(max_xpos,np.max(xpos2)),np.maximum(max_ypos,np.max(ypos2)),np.maximum(max_vsys,np.max(vsys2)) \n";

    for (int i=0; i<2; i++) {
        if (nc[i]>0) {
            py_file << "\terr2_l[" << i << "], err2_h[" << i << "] = np.genfromtxt(file2,skip_header=1,usecols=("
                    << nc[i] << "," << nc[i]+1 << "),unpack=True) \n";
        }
    }

    py_file << std::endl
            << "rad_sd, surfdens, sd_err = np.genfromtxt(filesb, usecols=(0,2,3),unpack=True) \n";

    py_file << "\nfig1=plt.figure(figsize=(11.69,8.27), dpi=150)  \n"
            << "plt.rc('font',family='sans-serif',serif='Helvetica',size=10)  \n"
            << "params = {'text.usetex': False, 'mathtext.fontset': 'cm', 'mathtext.default': 'regular', 'errorbar.capsize': 0} \n"
            << "plt.rcParams.update(params) \n"
            << "fig_ratio = 11.69/8.27 \n"
            << "nrows, ncols = 3,3 \n"
            << "x_axis_length, y_axis_length = 0.27, 0.13 \n"
            << "x_sep, y_sep = 0.07,0.015 \n"
            << "ax, bottom_corner = [], [0.1,0.7]\n"
            << "for i in range (nrows): \n"
            << "\tbottom_corner[0], axcol, ylen = 0.1, [], y_axis_length \n"
            << "\tif i==0: ylen *= 1.8 \n"
            << "\tfor j in range (ncols): \n"
            << "\t\taxcol.append(fig1.add_axes([bottom_corner[0],bottom_corner[1],x_axis_length,ylen*fig_ratio])) \n"
            << "\t\tbottom_corner[0]+=x_axis_length+x_sep \n"
            << "\tax.append(axcol) \n"
            << "\tbottom_corner[1]-=(y_axis_length+y_sep)*fig_ratio \n";

    py_file << std::endl
            << "axis=ax[0][0]  \n"
            << "axis.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on')  \n"
            << "axis.set_xlim(0,max_rad)  \n"
            << "axis.set_ylim(0,1.2*max_vrot)  \n"
            << "axis.set_ylabel('v$_\\mathrm{rot}$ (km/s)', fontsize=14)  \n"
            << "axis.errorbar(rad,vrot, yerr=[err1_l[0],-err1_h[0]],fmt='o', color=color)  \n"
            << "if twostage==True: axis.errorbar(rad2,vrot2, yerr=[err2_l[0],-err2_h[0]],fmt='o', color=color2) \n";

    py_file << std::endl
            << "axis=ax[1][0]  \n"
            << "axis.set_xlim(0,max_rad)  \n"
            << "axis.set_ylabel('i (deg)', fontsize=14)  \n"
            << "axis.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') \n"
            << "axis.errorbar(rad,inc, yerr=[err1_l[4],-err1_h[4]],fmt='o', color=color)  \n"
            << "if twostage==True: axis.errorbar(rad2,inc2,yerr=[err2_l[4],-err2_h[4]], fmt='o-', color=color2) \n";

    py_file << std::endl
            << "axis=ax[2][0]  \n"
            << "axis.set_xlim(0,max_rad)  \n"
            << "axis.set_ylabel('$\\phi$ (deg)', fontsize=14)  \n"
            << "axis.set_xlabel('Radius (arcsec)', fontsize=14, labelpad=10) \n"
            << "axis.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='on',labelleft='on')  \n"
            << "axis.errorbar(rad,pa, yerr=[err1_l[5],-err1_h[5]],fmt='o', color=color)  \n"
            << "if twostage==True: axis.errorbar(rad2,pa2,yerr=[err2_l[5],-err2_h[5]], fmt='o-', color=color2)  \n";

    py_file << std::endl
            << "axis=ax[0][1]  \n"
            << "axis.set_xlim(0,max_rad)  \n"
            << "axis.set_ylim(0,1.2*max_vdisp)  \n"
            << "axis.set_ylabel('$\\sigma_\\mathrm{gas}$  (km/s)', fontsize=14)  \n"
            << "axis.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') \n"
            << "axis.errorbar(rad,disp, yerr=[err1_l[1],-err1_h[1]],fmt='o', color=color)  \n"
            << "if twostage==True: axis.errorbar(rad2,disp2, yerr=[err2_l[1],-err2_h[1]],fmt='o', color=color2)  \n"   ;

    py_file << std::endl
            << "axis=ax[1][1]  \n"
            << "axis.set_xlim(0,max_rad)  \n"
            << "axis.set_ylabel('x$_0$ (pix)', fontsize=14)  \n"
            << "axis.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on')   \n"
            << "axis.errorbar(rad,xpos, yerr=[err1_l[6],-err1_h[6]],fmt='o', color=color)  \n"
            << "if twostage==True: axis.errorbar(rad2,xpos2,yerr=[err2_l[6],-err2_h[6]],fmt='o-', color=color2)  \n";

    py_file << std::endl
            << "axis=ax[2][1]  \n"
            << "axis.set_xlim(0,max_rad)  \n"
            << "axis.set_ylabel('y$_0$ (pix)', fontsize=14)  \n"
            << "axis.set_xlabel('Radius (arcsec)', fontsize=14, labelpad=10) \n"
            << "axis.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='on',labelleft='on')  \n"
            << "axis.errorbar(rad,ypos, yerr=[err1_l[7],-err1_h[7]],fmt='o', color=color)  \n"
            << "if twostage==True: axis.errorbar(rad2,ypos2, yerr=[err2_l[7],-err2_h[7]],fmt='o-', color=color2) \n";

    py_file << std::endl
            << "axis=ax[0][2]  \n"
            << "axis.set_xlim(0,max_rad)  \n"
            << "axis.set_ylabel('$\\Sigma}$ ("+in->Head().Bunit()+"/arcs$^2$)', fontsize=14)  \n"
            << "axis.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on')  \n"
            << "axis.errorbar(rad_sd,surfdens, yerr=sd_err,fmt='o', color=color2)  \n";

    py_file << std::endl
            << "axis=ax[1][2]  \n"
            << "axis.set_xlim(0,max_rad)  \n"
            << "axis.set_ylabel('z$_0$ (arcs)', fontsize=14)  \n"
            << "axis.tick_params(axis='both',which='both',bottom='off',top='on',labelbottom='off',labelleft='on')  \n"
            << "axis.errorbar(rad,z0, yerr=[err1_l[3],-err1_h[3]],fmt='o', color=color)  \n"
            << "if twostage==True: axis.errorbar(rad2,z02,yerr=[err2_l[3],-err2_h[3]],fmt='o-', color=color2)  \n";

    py_file << std::endl
            << "axis=ax[2][2]  \n"
            << "axis.set_xlim(0,max_rad) \n"
            << "axis.set_ylabel('v$_\\mathrm{sys}$ (km/s)', fontsize=14) \n"
            << "axis.set_xlabel('Radius (arcsec)', fontsize=14, labelpad=10) \n"
            << "axis.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='on',labelleft='on')  \n"
            << "axis.errorbar(rad,vsys, yerr=[err1_l[8],-err1_h[8]],fmt='o', color=color)  \n"
            << "if twostage==True: axis.errorbar(rad2,vsys2,yerr=[err2_l[8],-err2_h[8]],fmt='o-', color=color2) \n";

    py_file << std::endl
            << "plt.savefig(outfolder+'plot_parameters.pdf', orientation = 'landscape', format = 'pdf',bbox_inches='tight') \n"
            << std::endl << std::endl;


    py_file << "# CHANNEL MAPS: Setting all the needed variables \n"
            << "image = fits.open('" << in->pars().getImageFile() << "') \n"
            << "xmin, xmax = " << xmin << ", " << xmax << std::endl
            << "ymin, ymax = " << ymin << ", " << ymax << std::endl
            << "zmin, zmax = " << zmin << ", " << zmax << std::endl
            << "imagedata = image[0].data[";
    if (in->Head().NumAx()>3)
        for (int i=0; i<in->Head().NumAx()-3; i++) py_file << "0,";
    py_file << "zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
            << "head = image[0].header \n"
            << "zsize=imagedata[:,0,0].size \n"
            << "cont = " << cont << std::endl
            << "v=[cont,2*cont,4*cont,8*cont,16*cont,32*cont,64*cont] \n"
            << "v_neg = [-cont] \n"
            << "interval = PercentileInterval(99.5) \n"
            << "vmax = interval.get_limits(imagedata)[1] \n"
            << "norm = ImageNormalize(vmin=cont, vmax=vmax, stretch=PowerStretch(0.5)) \n\n";

    py_file << "files_mod, files_pva_mod, files_pvb_mod, typ = [], [], [], [] \n"
            << "for thisFile in os.listdir(outfolder): \n"
            << "\tif 'mod_azim.fits'  in thisFile: files_mod.append(thisFile) \n"
            << "\tif 'pv_a_azim.fits' in thisFile: files_pva_mod.append(thisFile) \n"
            << "\tif 'pv_b_azim.fits' in thisFile: files_pvb_mod.append(thisFile) \n"
            << std::endl
            << "if len(files_mod)==1: typ.append('AZIM') \n"
            << std::endl
            << "for thisFile in os.listdir(outfolder): \n"
            << "\tif 'mod_local.fits'  in thisFile: files_mod.append(thisFile) \n"
            << "\tif 'pv_a_local.fits' in thisFile: files_pva_mod.append(thisFile) \n"
            << "\tif 'pv_b_local.fits' in thisFile: files_pvb_mod.append(thisFile) \n"
            << std::endl
            << "if len(files_mod)==2: typ.append('LOCAL') \n"
            << "elif (len(files_mod)==1 and len(typ)==0): typ.append('LOCAL') \n"
            << "else: exit() \n\n";

    py_file << "# Beginning channel map plot \n"
            << "for k in range (len(files_mod)): \n"
            << "\timage_mod = fits.open(outfolder+files_mod[k]) \n"
            << "\timagedata_mod = image_mod[0].data[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1] \n"
            << "\tplt.figure(figsize=(8.27, 11.69), dpi=100) \n"
            << "\tgrid = [gridspec.GridSpec(2,5),gridspec.GridSpec(2,5),gridspec.GridSpec(2,5)] \n"
            << "\tgrid[0].update(top=0.95, bottom=0.695, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
            << "\tgrid[1].update(top=0.65, bottom=0.395, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
            << "\tgrid[2].update(top=0.35, bottom=0.095, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
            << "\tmatplotlib.rcParams['contour.negative_linestyle'] = 'solid' \n"
            << std::endl
            << "\tnum = 0 \n"
            << "\tfor j in range (0,3): \n"
            << "\t\tfor i in range (0,5): \n"
            << "\t\t\tchan = int(num*(zsize)/15) \n"
            << "\t\t\tz = imagedata[chan,:,:] \n"
            << "\t\t\tz_mod = imagedata_mod[chan,:,:] \n"
            << "\t\t\tvelo_kms = (chan+1-" << crpix3_kms-zmin << ")*" << cdelt3_kms << "+" << crval3_kms << std::endl
            << "\t\t\tvelo = ' v = ' + str(int(velo_kms)) + ' km/s' \n"
            << "\t\t\tplt.subplot(grid[j][0,i]) \n"
            << "\t\t\tplt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='off') \n"
            << "\t\t\tplt.title(velo, fontsize=8,loc='left') \n"
            << "\t\t\tplt.imshow(z,origin='lower',cmap = matplotlib.cm.Greys,norm=norm,aspect='auto') \n"
            << "\t\t\tplt.contour(z,v,origin='lower',linewidths=0.7,colors='#00008B') \n"
            << "\t\t\tplt.contour(z,v_neg,origin='lower',linewidths=0.1,colors='gray') \n"
            << "\t\t\tplt.subplot(grid[j][1,i]) \n"
            << "\t\t\tplt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='off') \n"
            << "\t\t\tplt.imshow(z_mod,origin='lower',cmap = matplotlib.cm.Greys,norm=norm,aspect='auto') \n"
            << "\t\t\tplt.contour(z_mod,v,origin='lower',linewidths=0.7,colors='#B22222') \n"
            << "\t\t\tnum = num+1 \n"
            << std::endl
            << "\toutfile = 'plot_chanmaps.pdf' \n"
            << "\tif (typ[k]=='AZIM'): outfile = 'plot_chanmaps_azim.pdf' \n"
            << "\tif (typ[k]=='LOCAL'): outfile = 'plot_chanmaps_local.pdf' \n"
            << "\tplt.savefig(outfolder+outfile, orientation = 'portrait', format = 'pdf') \n"
            << "\timage_mod.close() \n\n"
            << "image.close() \n";

    py_file << std::endl << std::endl;

    float zmin_wcs = AlltoVel(in->getZphys(zmin),in->Head());
    float zmax_wcs = AlltoVel(in->getZphys(zmax),in->Head());
    int pa_av = lround(findMedian(&outr->phi[0],outr->nr));
    int pa_min = pa_av+90<360 ? pa_av+90 : pa_av-90;
    if (zmin_wcs>zmax_wcs) std::swap(zmin_wcs,zmax_wcs);
    bool reverse = (pa_av>=45 && pa_av<225);
    if (cdelt3_kms<0) reverse = !reverse;

    py_file << "# Now plotting the position-velocity diagrams \n"
            << "image_maj     = fits.open('"<< in->pars().getOutfolder() << in->Head().Name() << "_pv_a.fits') \n"
            << "image_min     = fits.open('"<< in->pars().getOutfolder() << in->Head().Name() << "_pv_b.fits') \n"
            << "head = [image_maj[0].header,image_min[0].header] \n"
            << "crpixpv = np.array([head[0]['CRPIX1'],head[1]['CRPIX1']]) \n"
            << "cdeltpv = np.array([head[0]['CDELT1'],head[1]['CDELT1']]) \n"
            << "crvalpv = np.array([head[0]['CRVAL1'],head[1]['CRVAL1']]) \n"
            << "xminpv, xmaxpv = np.floor(crpixpv-1-" << disp << "), np.ceil(crpixpv-1 +"<< disp << ") \n"
            << "if xminpv[0]<0: xminpv[0]=0 \n"
            << "if xminpv[1]<0: xminpv[1]=0 \n"
            << "if xmaxpv[0]>=head[0]['NAXIS1']: xmaxpv[0]=head[0]['NAXIS1']-1 \n"
            << "if xmaxpv[1]>=head[1]['NAXIS1']: xmaxpv[1]=head[1]['NAXIS1']-1 \n"
            << "imagedata_maj = image_maj[0].data[zmin:zmax+1,int(xminpv[0]):int(xmaxpv[0])+1] \n"
            << "imagedata_min = image_min[0].data[zmin:zmax+1,int(xminpv[1]):int(xmaxpv[1])+1] \n"
            << "xmin_wcs = ((xminpv+1-crpixpv)*cdeltpv+crvalpv)*" << arcconv << std::endl
            << "xmax_wcs = ((xmaxpv+1-crpixpv)*cdeltpv+crvalpv)*" << arcconv << std::endl
            << "zmin_wcs, zmax_wcs = " << zmin_wcs << ", " << zmax_wcs << std::endl;

    py_file << std::endl
            << "radius = np.concatenate((rad,-rad)) \n"
            << "vrotation, inclin, vsystem, posang = vrot, inc, vsys, pa  \n"
            << "if twostage==True:\n "
            << "\tradius, vrotation, inclin, vsystem, posang = np.concatenate((rad2,-rad2)), vrot2, inc2, vsys2, pa2 \n"
            << "vlos1 = vrotation*np.sin(np.deg2rad(inclin))+vsystem \n"
            << "vlos2 = vsystem-vrotation*np.sin(np.deg2rad(inclin)) \n";
    if (reverse) py_file << "reverse = True \n";
    else py_file << "reverse = False \n";
    py_file << "if reverse==True: vlos1, vlos2 = vlos2, vlos1 \n"
            << "vlos = np.concatenate((vlos1,vlos2)) \n"
            << "vsysmean, pamean = np.nanmean(vsystem), np.nanmean(posang) \n"
            << "ext = [[xmin_wcs[0],xmax_wcs[0],zmin_wcs-vsysmean,zmax_wcs-vsysmean],\\" << std::endl
            << "       [xmin_wcs[1],xmax_wcs[1],zmin_wcs-vsysmean,zmax_wcs-vsysmean]] \n"
            << "labsize = 15 \n"
            << "palab = ['$\\phi = $" << pa_av <<"$^\\circ$', '$\\phi = $" << pa_min << "$^\\circ$'] \n\n";
    
    py_file << "# Beginning PV plot \n"
            << "for k in range (len(files_mod)): \n"
            << "\timage_mod_maj = fits.open(outfolder+files_pva_mod[k]) \n"
            << "\timage_mod_min = fits.open(outfolder+files_pvb_mod[k]) \n"
            << "\timagedata_mod_maj = image_mod_maj[0].data[zmin:zmax+1,int(xminpv[0]):int(xmaxpv[0])+1] \n"    
            << "\timagedata_mod_min = image_mod_min[0].data[zmin:zmax+1,int(xminpv[1]):int(xmaxpv[1])+1] \n"
            << "\ttoplot = [[imagedata_maj,imagedata_min],[imagedata_mod_maj,imagedata_mod_min]] \n"
            << std::endl
            << "\tfig=plt.figure(figsize=(11.69,8.27), dpi=150) \n"
            << "\tfig_ratio = 11.69/8.27 \n"
            << "\tx_len, y_len, y_sep = 0.6, 0.42, 0.08 \n"
            << "\tax, bottom_corner = [], [0.1,0.7] \n"
            << "\tfor i in range (2): \n"
            << "\t\tbottom_corner[0], axcol = 0.1, [] \n"
            << "\t\tax.append(fig.add_axes([bottom_corner[0],bottom_corner[1],x_len,y_len*fig_ratio])) \n"
            << "\t\tbottom_corner[1]-=(y_len+y_sep)*fig_ratio \n"
            << std::endl
            << "\tfor i in range (2): \n"
            << "\t\taxis = ax[i] \n"
            << "\t\taxis.tick_params(which='major',length=8, labelsize=labsize) \n"
            << "\t\taxis.set_xlabel('Offset (arcsec)',fontsize=labsize+2) \n"
            << "\t\taxis.set_ylabel('$\\mathrm{\\Delta V_{LOS}}$ (km/s)',fontsize=labsize+2) \n"
            << "\t\taxis.text(1, 1.02,palab[i],ha='right',transform=axis.transAxes,fontsize=labsize+4) \n"
            << "\t\taxis2 = axis.twinx() \n"
            << "\t\taxis2.set_xlim([ext[i][0],ext[i][1]]) \n"
            << "\t\taxis2.set_ylim([ext[i][2]+vsysmean,ext[i][3]+vsysmean]) \n"
            << "\t\taxis2.tick_params(which='major',length=8, labelsize=labsize) \n"
            << "\t\taxis2.set_ylabel('$\\mathrm{V_{LOS}}$ (km/s)',fontsize=labsize+2) \n"
            << "\t\taxis.imshow(toplot[0][i],origin='lower',cmap = matplotlib.cm.Greys,norm=norm,extent=ext[i],aspect='auto') \n"
            << "\t\taxis.contour(toplot[0][i],v,origin='lower',linewidths=0.7,colors='#00008B',extent=ext[i],aspect='auto') \n"
            << "\t\taxis.contour(toplot[0][i],v_neg,origin='lower',linewidths=0.1,colors='gray',extent=ext[i],aspect='auto') \n"
            << "\t\taxis.contour(toplot[1][i],v,origin='lower',linewidths=0.7,colors='#B22222',extent=ext[i],aspect='auto') \n"
            << "\t\taxis.axhline(y=0,color='black') \n"
            << "\t\taxis.axvline(x=0,color='black') \n"
            << "\t\tif i==0 :axis2.plot(radius,vlos,'yo') \n"
            << std::endl
            << "\toutfile = 'plot_pv.pdf' \n"
            << "\tif (typ[k]=='AZIM'): outfile = 'plot_pv_azim.pdf' \n"
            << "\tif (typ[k]=='LOCAL'): outfile = 'plot_pv_local.pdf' \n"
            << "\tplt.savefig(outfolder+outfile, bbox_inches='tight') \n"
            << "\timage_mod_maj.close() \n"
            << "\timage_mod_min.close() \n"
            << std::endl
            << "image_maj.close() \n"
            << "image_min.close() \n";
    


    py_file << std::endl << "# Now plotting moment maps \n"
            << "mom0 = fits.open(outfolder+'/maps/"<< in->Head().Name() << "_0mom.fits')[0].data[ymin:ymax+1,xmin:xmax+1] \n"
            << "mom1 = fits.open(outfolder+'/maps/"<< in->Head().Name() << "_1mom.fits')[0].data[ymin:ymax+1,xmin:xmax+1] \n"
            << "mom2 = fits.open(outfolder+'/maps/"<< in->Head().Name() << "_2mom.fits')[0].data[ymin:ymax+1,xmin:xmax+1] \n"
            << std::endl
            << "files_mod0, files_mod1, files_mod2 = [], [], [] \n"
            << "for thisFile in os.listdir(outfolder+'/maps/'): \n"
            << "\tif 'azim_0mom.fits' in thisFile: files_mod0.append(thisFile) \n"
            << "\tif 'azim_1mom.fits' in thisFile: files_mod1.append(thisFile) \n"
            << "\tif 'azim_2mom.fits' in thisFile: files_mod2.append(thisFile) \n"
            << std::endl
            << "for thisFile in os.listdir(outfolder+'/maps/'): \n"
            << "\tif 'local_0mom.fits' in thisFile: files_mod0.append(thisFile) \n"
            << "\tif 'local_1mom.fits' in thisFile: files_mod1.append(thisFile) \n"
            << "\tif 'local_2mom.fits' in thisFile: files_mod2.append(thisFile) \n"
            << std::endl
            << "norm0 = ImageNormalize(vmin=np.nanmin(mom0), vmax=np.nanmax(mom0), stretch=LinearStretch()) \n"
            << "norm1 = ImageNormalize(vmin=np.nanmin(mom1), vmax=np.nanmax(mom1), stretch=LinearStretch()) \n"
            << "norm2 = ImageNormalize(vmin=np.nanmin(mom2), vmax=np.nanmax(mom2), stretch=LinearStretch()) \n"
            << "norm = [norm0, norm1, norm2] \n"
            << "cmaps = [matplotlib.cm.jet,matplotlib.cm.jet,matplotlib.cm.jet] \n"
            << "barlab = ['Intensity (" << in->Head().Bunit() <<")', 'V$_\\mathrm{LOS}$ (km/s)', '$\\sigma$ (km/s)'] \n"
            << "titles = ['DATA', 'MODEL'] \n"
            << "mapname = ['INTENSITY', 'VELOCITY', 'DISPERSION'] \n"
            << std::endl
            << "xcen, ycen, phi = [np.nanmean(xpos)-xmin,np.nanmean(ypos)-ymin,np.nanmean(pa)] \n"
            << "if twostage==True: xcen, ycen, phi = [np.nanmean(xpos2),np.nanmean(ypos2),np.nanmean(pa2)] \n"
            << "x = np.arange(0,xmax-xmin,0.1) \n"
            << "y = np.tan(np.radians(phi-90))*(x-xcen)+ycen \n"
            << "ext = [0,xmax-xmin,0, ymax-ymin] \n"
            << std::endl
            << "for k in range (len(files_mod0)): \n"
            << "\tmom0_mod = fits.open(outfolder+'/maps/'+files_mod0[k])[0].data[ymin:ymax+1,xmin:xmax+1] \n"
            << "\tmom1_mod = fits.open(outfolder+'/maps/'+files_mod1[k])[0].data[ymin:ymax+1,xmin:xmax+1] \n"
            << "\tmom2_mod = fits.open(outfolder+'/maps/'+files_mod2[k])[0].data[ymin:ymax+1,xmin:xmax+1] \n"
            << "\tto_plot = [[mom0,mom1,mom2],[mom0_mod,mom1_mod,mom2_mod]] \n"
            << std::endl
            << "\tfig=plt.figure(figsize=(11.69,8.27), dpi=150) \n"
            << "\tfig_ratio = 11.69/8.27 \n"
            << "\tnrows, ncols = 3, 2 \n"
            << "\tx_len, y_len = 0.2, 0.2 \n"
            << "\tx_sep, y_sep = 0.00,0.02 \n"
            << "\tax, ax_cb, bottom_corner = [], [], [0.1,0.7] \n"
            << "\tfor i in range (nrows): \n"
            << "\t\tbottom_corner[0], axcol = 0.1, [] \n"
            << "\t\tfor j in range (ncols): \n"
            << "\t\t\taxcol.append(fig.add_axes([bottom_corner[0],bottom_corner[1],x_len,y_len*fig_ratio])) \n"
            << "\t\t\tbottom_corner[0]+=x_len+x_sep \n"
            << "\t\tax.append(axcol) \n"
            << "\t\tax_cb.append(fig.add_axes([bottom_corner[0]+0.01,bottom_corner[1],0.02,y_len*fig_ratio])) \n"
            << "\t\tbottom_corner[1]-=(y_len+y_sep)*fig_ratio \n"
            << std::endl
            << "\tfor i in range (nrows): \n"
            << "\t\tfor j in range (ncols): \n"
            << "\t\t\taxis = ax[i][j] \n"
            << "\t\t\taxis.tick_params(labelbottom='off',labelleft='off') \n"
            << "\t\t\taxis.set_xlim(ext[0],ext[1]) \n"
            << "\t\t\taxis.set_ylim(ext[2],ext[3]) \n"
            << "\t\t\taxis.imshow(to_plot[j][i],origin='lower',cmap=cmaps[i],norm=norm[i],aspect='auto',extent=ext) \n"
            << "\t\t\taxis.plot(xcen,ycen,'x',color='#000000',markersize=7,mew=1.5) \n"
            << std::endl
            << "\t\t\tcb = ColorbarBase(ax_cb[i], orientation='vertical', cmap=cmaps[i], norm=norm[i]) \n"
            << "\t\t\tcb.solids.set_edgecolor('face') \n"
            << "\t\t\tcb.set_label(barlab[i],fontsize=13) \n"
            << std::endl
            << "\t\t\tif i==0: axis.text(0.5,1.05,titles[j],ha='center',transform=axis.transAxes,fontsize=15) \n"
            << "\t\t\telif i==1: axis.plot(x,y,color='#808080',linewidth=2) \n"
            << "\t\t\tif j==0: axis.text(-0.1,0.5,mapname[i],va='center',rotation=90,transform=axis.transAxes,fontsize=15) \n"
            << std::endl
            << "\tif (typ[k]=='AZIM'): outfile = 'plot_maps_azim.pdf' \n"
            << "\tif (typ[k]=='LOCAL'): outfile = 'plot_maps_local.pdf' \n"
            << "\tplt.savefig(outfolder+outfile, bbox_inches = 'tight') \n";

    py_file.close();

#ifdef HAVE_PYTHON
    std::string cmd = "python "+in->pars().getOutfolder()+scriptname+" > /dev/null 2>&1";
    return system(cmd.c_str());
    //Py_Initialize();
    //std::string cmd = "import subprocess, os\n";
    //cmd += "FNULL = open(os.devnull, 'w')\n";
    //cmd += "subprocess.call(['python','"+in->pars().getOutfolder()+scriptname+"'],stdout=FNULL, stderr=FNULL)\n";
    //int ret = PyRun_SimpleString(cmd.c_str());
    //Py_Finalize();
    //return ret;    
#endif

    return -1;
}
template int Galfit<float>::plotAll_Python();
template int Galfit<double>::plotAll_Python();


template <class T>
void Galfit<T>::printDetails (Rings<T> *dr, T fmin, long pix, std::ostream& str) {


    int m=7, n=9;

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
    Stream << showpoint << fixed << setprecision(2) << endl;

    Stream << setfill('=') << setw(44) << right << " Initial parameters " << setw(25) << " ";
    Stream << setfill(' ') << endl;

    Stream << "    (i) input by the user    (d) default    (e) estimated by me" << endl << endl;

    string s;
    s = "   Fitting #" + to_string<int>(inr->nr);
    if (in->pars().getNRADII()==-1) s += "(e) ";
    else s += "(i) ";
    s += "rings of width " + to_string<T>(inr->radsep);
    if (in->pars().getRADSEP()==-1) s += "(e) ";
    else s += "(i) ";
    s += "arcsec";
    Stream << s << endl << endl;


    s = "    Xpos";
    if (in->pars().getXPOS()=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->xpos[0] << left << setw(m) << "  pix";


    s = "        Ypos";
    if (in->pars().getYPOS()=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n+4) << left << s << setw(3) << right << "= "
         << setw(m-1) << inr->ypos[0]
         << left << setw(m) << "  pix" << endl;

    s = "    Vsys";
    if (in->pars().getVSYS()=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->vsys[0] << left << setw(m) << "  km/s";

    s = "        Vrot";
    if (in->pars().getVROT()=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n+4) << left << s << setw(3) << right << "= "
         << setw(m-1) << inr->vrot[0]
         << left << setw(m) << "  km/s" << endl;

    s = "    Inc";
    if (in->pars().getINC()=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->inc[0] << left << setw(m) << "  deg";

    s = "        PA";
    if (in->pars().getPHI()=="-1") s += " (e)";
    else s += " (i)";
    Stream << setw(n+4) << left << s << setw(3) << right << "= "
           << setw(m-1) << inr->phi[0] << left << setw(m) << "  deg" << endl;

    double toKpc = KpcPerArc(distance);
    s = "    Z0";
    if (in->pars().getZ0()=="-1") s += " (d)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->z0[0]*toKpc*1000 << left << setw(m) << "  pc";

    s = "        Disp";
    if (in->pars().getVDISP()=="-1") s += " (d)";
    else s += " (i)";
    Stream << setw(n+4) << left << s << setw(3) << right << "= "
         << setw(m-1) << inr->vdisp[0]
         << left << setw(m) << "  km/s" << endl;
    
    s = "    Vrad";
    if (in->pars().getVRAD()=="-1") s += " (d)";
    else s += " (i)";
    Stream << setw(n) << left << s << setw(3) << right << "= "
           << setw(m) << inr->vrad[0] << left << setw(m) << "  km/s";

    Stream   << endl << endl;

    Stream   << endl;

}
template void Galfit<float>::showInitial(Rings<float>*,std::ostream&);
template void Galfit<double>::showInitial(Rings<double>*,std::ostream&);


template <class T>
void Galfit<T>::printInitial (Rings<T> *inr) {

    int m=11;
    std::ofstream initout;
    initout.open((in->pars().getOutfolder()+"initial_rings.txt").c_str());
    initout << setfill('=');
    initout << setw(66) << right << " Initial parameters " << setw(46) << " " << endl;
    initout << setfill(' ');

    initout << left << setw(m) << "RAD(arcs)"
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
template void Galfit<float>::printInitial(Rings<float>*);
template void Galfit<double>::printInitial(Rings<double>*);


template <class T>
void Galfit<T>::DensityProfile (T *surf_dens, int *count) {

    // -STRONZATEEE-----------------------------------------------------------

    if (!in->getIsSearched()) in->Search();
    Detection<T> *obj = in->pObject(0);
    obj->calcFluxes(obj->getPixelSet(in->Array(), in->AxisDim()));
    obj->calcWCSparams(in->Head());
    obj->calcIntegFlux(in->DimZ(), obj->getPixelSet(in->Array(), in->AxisDim()), in->Head());
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

#undef VROT
#undef VDISP
#undef DENS
#undef Z0
#undef INC
#undef PA
#undef XPOS
#undef YPOS
#undef VSYS
#undef VRAD
#undef MAXPAR
