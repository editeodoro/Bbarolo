//------------------------------------------------------------------------
// galfit_out.cpp: Members functions for outputs of the Galfit class.
//------------------------------------------------------------------------

/*-----------------------------------------------------------------------
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 Bbarp;p is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with Bbarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning Bbarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "../Arrays/cube.hh"
#include "../Arrays/image.hh"
#include "galfit.hh"
#include "galmod.hh"
#include "utils.hh"
#include "gnuplot.hh"
#include "moment.hh"

#define VROT  0
#define VDISP 1
#define DENS  2
#define Z0    3
#define INC   4
#define PA    5
#define XPOS  6
#define YPOS  7
#define VSYS  8
#define MAXPAR 9

namespace Model {

template <class T>
void Galfit<T>::writeModel_norm() {

    bool verbosity = in->pars().isVerbose();
    if (verbosity && !in->pars().getflagGalMod()) {
        in->pars().setVerbosity(false);
        std::cout << " Preparing a bunch of cool outputs..." << std::flush;
    }

    std::string outfold = in->pars().getOutfolder();
    std::string object = in->Head().Name();

    Model::Galmod<T> *mod = getModel();

    Cube<T> *res = new Cube<T>(in->AxisDim());
    res->saveHead(in->Head());
    res->saveParam(in->pars());

///*	Anello 2D --------------------------------------------------------------------
    T *surf_dens = new T[outr->nr];
    int *count = new int[outr->nr];
    for (int i=0;i<outr->nr;i++) {surf_dens[i]=0; count[i]=0;}

    for (int y=0; y<in->DimY(); y++) {
        for (int x=0; x<in->DimX(); x++) {
            T modSum = 0;
            T obsSum = 0;
            T factor = 0;
            for (int z=0; z<in->DimZ(); z++) {
                long Pix = in->nPix(x,y,z);
                modSum += mod->Out()->Array(Pix);
                obsSum += in->Array(Pix)*mask[Pix];
            }

            for (int ir=0;ir<outr->nr;ir++) {
                double pixScale = ((fabs(in->Head().Cdelt(0))*arcconv)+
                                    (fabs(in->Head().Cdelt(1))*arcconv))/2.;
                double inc = outr->inc[ir]*M_PI/180.;
                double phi = outr->phi[ir]*M_PI/180.;
                double radsep1=0, radsep2=0;
                if (ir==0) {
                    radsep1 = 0;
                    radsep2 = (outr->radii[1]-outr->radii[0])/2;
                }
                else if (ir==outr->nr-1) {
                    radsep1 = (outr->radii[ir]-outr->radii[ir-1])/2;
                    radsep2 = 0;
                }
                else {
                    radsep1 = (outr->radii[ir]-outr->radii[ir-1])/2;
                    radsep2 = (outr->radii[ir+1]-outr->radii[ir])/2;
                }
                double r1 = (outr->radii[ir]-radsep1)/pixScale;
                double r2 = (outr->radii[ir]+radsep2)/pixScale;
                double x0 = outr->xpos[ir];
                double y0 = outr->ypos[ir];
                double xr =  -(x-x0)*sin(phi)+(y-y0)*cos(phi);
                double yr = (-(x-x0)*cos(phi)-(y-y0)*sin(phi))/cos(inc);
                double r = sqrt(xr*xr+yr*yr);
                bool isin = ir==outr->nr-1 ? r>=r1 && r<=r2+sqrt(in->Head().BeamArea()/M_PI) : r>=r1 && r<=r2;
                if (!isin) continue;
                surf_dens[ir] += obsSum;
                count[ir]++;
            }

            if (modSum!=0) factor = obsSum/modSum;
            else continue;
            for (int z=0; z<in->DimZ(); z++) {
                long Pix = in->nPix(x,y,z);
                mod->Out()->Array()[Pix] *= factor;

            }
        }
    }

    for (int i=0; i<in->NumPix(); i++) {
        T model = mod->Out()->Array(i);
        if (mask[i] || model!=0) res->Array()[i] = (in->Array(i)-model);
        else res->Array()[i] = 0;
    }
    // --------------------------------------------------------------------------------
//*/


    /*// ANELLO 3D ----------------------------------------------------------------------
    T *ringreg = getFinalRingsRegion();
    for (int y=0; y<in->DimY(); y++) {
        for (int x=0; x<in->DimX(); x++) {
            T factor = 0;
            if (!isNaN(ringreg[x+y*in->DimX()])) {
                T modSum = 0;
                T obsSum = 0;
                for (int z=0; z<in->DimZ(); z++) {
                    long Pix = in->nPix(x,y,z);
                    modSum += mod->Out()->Array(Pix);
                    obsSum += in->Array(Pix)*mask[Pix];
                }
                if (modSum!=0) factor = obsSum/modSum;
            }
            for (int z=0; z<in->DimZ(); z++) {
                long Pix = in->nPix(x,y,z);
                mod->Out()->Array()[Pix] *= factor;
            }
        }
    }
    //------------------------------------------------------------------------------------
 */


    std::string mfile = outfold+object+"mod.fits";
    mod->Out()->Head().setDataMax(0.);
    mod->Out()->Head().setDataMin(0.);

    mod->Out()->fitswrite_3d(mfile.c_str());
    mfile = outfold+object+"res.fits";
    //res->fitswrite_3d(mfile.c_str());

    plotParam(mod->Out(), res);
    plotChanMaps();

    if (verbosity && !in->pars().getflagGalMod()) {
        in->pars().setVerbosity(true);
        std::cout << " OK!" << std::endl;
    }

    delete mod;
    delete res;


}
template void Galfit<float>::writeModel_norm();
template void Galfit<double>::writeModel_norm();


template <class T>
void Galfit<T>::writeModel_azim() {

    bool verbosity = in->pars().isVerbose();
    if (verbosity) {
        in->pars().setVerbosity(false);
        std::cout << " Preparing a bunch of cool outputs..." << std::flush;
    }

    std::string outfold = in->pars().getOutfolder();
    std::string object = in->Head().Name();

    Model::Galmod<T> *mod = getModel();
    Cube<T> *out = mod->Out();

    Cube<T> *res = new Cube<T>(in->AxisDim());
    res->saveHead(in->Head());
    res->saveParam(in->pars());

    T *surf_dens = new T[outr->nr];
    T *map_mod = new T[out->DimX()*out->DimY()];

    int *count = new int[outr->nr];
    float *rmap  = new float[out->DimX()*out->DimY()];
    for (int i=0;i<out->DimX()*out->DimY();i++) rmap[i]=-1;
    for (int i=0;i<outr->nr;i++) {
        surf_dens[i]=0;
        count[i]=0;
    }

  /*
    /// ANELLO 2D -------------------------------------------------------------
    for (int y=0; y<in->DimY(); y++) {
        for (int x=0; x<in->DimX(); x++) {
            map_mod[x+y*out->DimX()]=0;
            T obsSum = 0;
            for (int z=0; z<in->DimZ(); z++) {
                long Pix = in->nPix(x,y,z);
                map_mod[x+y*out->DimX()] += mod->Out()->Array(Pix);
                obsSum += in->Array(Pix)*mask[Pix];
            }

            for (int ir=0;ir<outr->nr;ir++) {
                double pixScale = ((fabs(in->Head().Cdelt(0))*arcconv)+
                                    (fabs(in->Head().Cdelt(1))*arcconv))/2.;
                double inc = outr->inc[ir]*M_PI/180.;
                double phi = outr->phi[ir]*M_PI/180.;
                double radsep1=0, radsep2=0;
                if (ir==0) {
                    radsep1 = 0;
                    radsep2 = (inr->radii[1]-inr->radii[0])/2;
                }
                else if (ir==inr->nr-1) {
                    radsep1 = (inr->radii[ir]-inr->radii[ir-1])/2;
                    radsep2 = 0;
                }
                else {
                    radsep1 = (outr->radii[ir]-inr->radii[ir-1])/2;
                    radsep2 = (outr->radii[ir+1]-inr->radii[ir])/2;
                }
                double r1 = (outr->radii[ir]-radsep1)/pixScale;
                double r2 = (outr->radii[ir]+radsep2)/pixScale;
                double x0 = outr->xpos[ir];
                double y0 = outr->ypos[ir];
                double xr =  -(x-x0)*sin(phi)+(y-y0)*cos(phi);
                double yr = (-(x-x0)*cos(phi)-(y-y0)*sin(phi))/cos(inc);
                double r = sqrt(xr*xr+yr*yr);
                //bool isin = r>=r1 && r<=r2;
                bool isin = ir==outr->nr-1 ? r>=r1 && r<=r2+sqrt(in->Head().BeamArea()/M_PI) : r>=r1 && r<=r2;
                if (!isin) continue;

                rmap[x+y*out->DimX()] = ir;

                if (obsSum!=0) {surf_dens[ir]+=obsSum; count[ir]++;}

            }

        }
    }
    //--------------------------------------------------------------------------
    */

///*
    // ANELLO 3D --------------------------------------------------------------
    T *ringreg = getFinalRingsRegion();
    for (int ir=0;ir<outr->nr;ir++) {
        Rings<T> *dr = new Rings<T>;
        dr->nr = 2;
        double radsep1, radsep2;
        if (ir==0) {
            radsep1 = (outr->radii[1]-outr->radii[0])/2;
            dr->radii.push_back(outr->radii[0]);
            dr->radii.push_back(outr->radii[0]+radsep1);
        }
        else if (ir==outr->nr-1){
            radsep2 = (outr->radii[ir]-outr->radii[ir-1])/2;
            dr->radii.push_back(outr->radii[ir]-radsep2);
            dr->radii.push_back(outr->radii[ir]+radsep2);
        }
        else {
            radsep1 = (outr->radii[ir]-outr->radii[ir-1])/2;
            radsep2 = (outr->radii[ir+1]-outr->radii[ir])/2;
            dr->radii.push_back(outr->radii[ir]-radsep1);
            dr->radii.push_back(outr->radii[ir]+radsep2);
        }

        for (int i=0; i<dr->nr; i++) {
            dr->vrot.push_back(outr->vrot[ir]);
            dr->vdisp.push_back(outr->vdisp[ir]);
            dr->dens.push_back(outr->dens[ir]);
            dr->z0.push_back(outr->z0[ir]);
            dr->inc.push_back(outr->inc[ir]);
            dr->phi.push_back(outr->phi[ir]);
            dr->xpos.push_back(outr->xpos[ir]);
            dr->ypos.push_back(outr->ypos[ir]);
            dr->vsys.push_back(outr->vsys[ir]);
        }

        int bhi[2] = {in->DimX(), in->DimY()};
        int blo[2] = {0,0};
        std::vector<PixelInfo::Pixel<T> > *ppix = getRingRegion(dr,bhi,blo);
        typename std::vector<PixelInfo::Pixel<T> >::iterator pix;

        for(pix=ppix->begin();pix<ppix->end();pix++) {
            long x = pix->getX();
            long y = pix->getY();
            map_mod[x+y*out->DimX()]=0;
            T obsSum = 0;
            for (int z=0; z<in->DimZ(); z++) {
                    long Pix = in->nPix(x,y,z);
                    map_mod[x+y*out->DimX()] += mod->Out()->Array(Pix);
                    obsSum += in->Array(Pix)*mask[Pix];
            }
            rmap[x+y*out->DimX()]=ir;
            if (obsSum!=0) {surf_dens[ir]+=obsSum; count[ir]++;}
        }

        delete ppix;
        delete dr;
    }
    delete [] ringreg;

    // -----------------------------------------------------------------------------
//*/

    for (int i=0;i<outr->nr;i++) surf_dens[i]/=count[i];

    //FitsWrite_2D("maremma.fits", rmap, in->DimX(),in->DimY());

    for (int y=0; y<out->DimY(); y++) {
        for (int x=0; x<out->DimX(); x++) {
            int radius = rmap[x+y*out->DimX()];
            T factor;
            if (radius==-1 || map_mod[x+y*out->DimX()]==0) factor=0;
            else factor = surf_dens[radius]/map_mod[x+y*out->DimX()];
            for (int z=0; z<in->DimZ(); z++) {
                long Pix = out->nPix(x,y,z);
                mod->Out()->Array()[Pix] *= factor;
            }
        }
    }

    for (int y=0; y<out->DimY(); y++) {
        for (int x=0; x<out->DimX(); x++) {
            rmap[x+y*out->DimX()] = 0;
            for (int z=0; z<in->DimZ(); z++) {
                long Pix = out->nPix(x,y,z);				\
                rmap[x+y*out->DimX()] +=  mod->Out()->Array()[Pix];
            }
        }
    }

    Image2D<float> *m = new Image2D<float>(in->AxisDim());
    m->saveHead(in->Head());
    m->saveParam(in->pars());
    m->setImage(rmap, in->AxisDim());
    std::string outname = outfold+object+"mod_azim_map.fits";
    //m->fitswrite_2d(outname.c_str());
    delete m;

    std::string mfile = outfold+object+"mod_azim.fits";
    out->Head().setDataMax(0.);
    out->Head().setDataMin(0.);
    mod->Out()->fitswrite_3d(mfile.c_str());
    //res->fitswrite_3d("DIOBONO0.fits");

    plotParam(mod->Out(), res);
    plotChanMaps();

    if (verbosity) {
        in->pars().setVerbosity(true);
        std::cout << " OK!" << std::endl;
    }


    delete [] rmap;
    delete [] surf_dens;
    delete [] map_mod;
    delete [] count;
    delete mod;
    delete res;

}
template void Galfit<float>::writeModel_azim();
template void Galfit<double>::writeModel_azim();


template <class T>
void Galfit<T>::plotParam (Cube<T> *mod, Cube<T> *res) {

    std::string outfold = in->pars().getOutfolder();
    std::string object = in->Head().Name();

    T meanPA = findMean(&outr->phi[0], outr->nr);
    T meanXpos = findMean(&outr->xpos[0], outr->nr);
    T meanYpos = findMean(&outr->ypos[0], outr->nr);
    T meanPAp90= meanPA+90<360 ? meanPA+90 : meanPA-90;
    T meanVsys = findMean(&outr->vsys[0], outr->nr);

    Image2D<T> *pv_max = PositionVelocity(in,meanXpos,meanYpos,meanPA);
    std::string mfile = outfold+object+"_pv_a.fits";
    pv_max->fitswrite_2d(mfile.c_str());
    Image2D<T> *pv_min = PositionVelocity(in,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+object+"_pv_b.fits";
    pv_min->fitswrite_2d(mfile.c_str());


    std::ofstream outpv((outfold+"pv.txt").c_str());
    std::ofstream outpvm((outfold+"pvm.txt").c_str());
    float xmin=1.E10,xmax=0,xmmin=1.E10,xmmax=0;
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

    outpv.open((outfold+"rcpv.txt").c_str());
    for (int i=in->pars().getStartRad(); i<outr->nr; i++) {
            float vel1 = (outr->vrot[i]*sin(outr->inc[i]*M_PI/180.))+outr->vsys[i];
            float vel2 = outr->vsys[i]-(outr->vrot[i]*sin(outr->inc[i]*M_PI/180.));
            if (meanPA<90 || meanPA>270) std::swap(vel1,vel2);
            float radius = outr->radii[i];
            if (i==0) radius += (outr->radsep/4.);
            outpv << -radius << "   " << vel1 << endl;
            outpv <<  radius << "   " << vel2 << endl;
    }
    outpv.close();


    pv_max = PositionVelocity(mod,meanXpos,meanYpos,meanPA);
    mfile = outfold+object+"mod_pv_a.fits";
    pv_max->fitswrite_2d(mfile.c_str());
    pv_min = PositionVelocity(mod,meanXpos,meanYpos,meanPAp90);
    mfile = outfold+object+"mod_pv_b.fits";
    pv_min->fitswrite_2d(mfile.c_str());

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

//	pv_max = PositionVelocity(res,meanXpos,meanYpos,meanPA);
//	mfile = outfold+object+"res_pv_a.fits";
//	pv_max->fitswrite_2d(mfile.c_str());
//	pv_min = PositionVelocity(res,meanXpos,meanYpos,meanPAp90);
//	mfile = outfold+object+"res_pv_b.fits";
//	pv_min->fitswrite_2d(mfile.c_str());

    delete pv_max;
    delete pv_min;

    int free[nfree];
    int k;
    for (int nm=0, k=0; nm<MAXPAR; nm++) {
        if (mpar[nm]) free[k++]=nm;
    }


    const int err_col=13;
    std::ofstream gnu;
    mfile = outfold+"gnuscript.gnu";
    gnu.open(mfile.c_str());


    float xtics = lround(outr->nr/5.);
    xtics *= outr->radsep;
    while (outr->radii.back()/xtics>5) xtics*=2;
    while (outr->radii.back()/xtics<2) xtics/=2;

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
        <<			  "set ylabel font \"Helvetica,13\" '" << endl
        << "TICSF	= 'set xtics font \"Helvetica,12\"; "
        <<			  "set ytics font \"Helvetica,12\" '" << endl
        << "TMARGIN = 'set tmargin at screen 0.95; set bmargin at screen 0.47; "
        << 			  "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "MMARGIN = 'set tmargin at screen 0.47; set bmargin at screen 0.27; "
        << 			  "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
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

    if (flagErrors && mpar[VROT]) {
        int nc=err_col;
        for (int i=0; i<nfree; i++) if (free[i]==VROT) nc+=2*i;
        gnu << "u 2:3:($3+$"+to_string(nc)+"):($3+$"+to_string(nc+1)+") w errorbars ls 1, '"
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

    if (flagErrors && mpar[INC]) {
        int nc=err_col;
        for (int i=0; i<nfree; i++) if (free[i]==INC) nc+=2*i;
        gnu << "u 2:5:($5+$"+to_string(nc)+"):($5+$"+to_string(nc+1)+") w errorbars ls 1, '"
            << outfold << "ringlog1.txt' u 2:5 w lp ls 1";
    }
    else gnu << "u 2:5 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:5 w lp ls 2";


    gnu	<< endl;

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

    if (flagErrors && mpar[PA]) {
        int nc=err_col;
        for (int i=0; i<nfree; i++) if (free[i]==PA) nc+=2*i;
        gnu << "u 2:6:($6+$"+to_string(nc)+"):($6+$"+to_string(nc+1)+") w errorbars ls 1, '"
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
        << 			  "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "MMARGIN = 'set tmargin at screen 0.66; set bmargin at screen 0.38; "
        << 			  "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
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

    if (flagErrors && mpar[VDISP]) {
        int nc=err_col;
        for (int i=0; i<nfree; i++) if (free[i]==VDISP) nc+=2*i;
        gnu << "u 2:4:($4+$"+to_string(nc)+"):($4+$"+to_string(nc+1)+") w errorbars ls 1, '"
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
    gnu	<< endl;



    // Plotting systemic velocity
    maxa = *max_element(&outr->vsys[0], &outr->vsys[0]+outr->nr);
    maxa += (0.1*maxa+10);
    mina = *min_element(&outr->vsys[0], &outr->vsys[0]+outr->nr);
    mina -= (0.1*mina+10);
    gnu << "@MMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [" << mina << ":" << maxa << "]" << endl
        << "set ylabel 'V_{sys} [km/s]'" << endl
        //<< "set ytics 5" << endl << "set mytics 5" << endl
        << "plot '" << outfold << "ringlog1.txt' u 2:12 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:12 w lp ls 2";
    gnu	<< endl;


    // Plotting scale height
    maxa = *max_element(&outr->z0[0], &outr->z0[0]+outr->nr);
    maxa += 0.1*maxa;
    mina = *min_element(&outr->z0[0], &outr->z0[0]+outr->nr);
    mina -= 0.1*mina;
    gnu << "@BMARGIN" << endl << "@XTICS" << endl
        << "set xlabel 'Radius [arcsec]'" << endl
        << "set yrange [" <<mina<<":"<<maxa<<"]\n"
        << "set ylabel 'Scale height [arcsec]'\n"
        << "plot '"<<in->pars().getOutfolder()<<"ringlog1.txt' u 2:8 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:8 w lp ls 2";
    gnu	<< endl;

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
        << "plot '"<<in->pars().getOutfolder()<<"ringlog1.txt' u 2:10 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:10 w lp ls 2";
    gnu	<< endl;

    // Plotting ycenter
    maxa = *max_element(&outr->ypos[0], &outr->ypos[0]+outr->nr);
    maxa += 0.1*maxa;
    mina = *min_element(&outr->ypos[0], &outr->ypos[0]+outr->nr);
    mina -= 0.1*mina;
    gnu << "@MMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [" << mina << ":" << maxa << "]" << endl
        << "set ylabel 'Y_c [pix]'" << endl
        //<< "set ytics 5" << endl << "set mytics 5" << endl
        << "plot '" << outfold << "ringlog1.txt' u 2:11 w lp ls 1";

    if (in->pars().getTwoStage())
        gnu << ", '" << outfold << "ringlog2.txt' u 2:11 w lp ls 2";
    gnu	<< endl;


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
        gnu	<< endl;
    }

    gnu << "unset multiplot; reset" << endl;
    gnu.close();

    /// Plotting pv contours
    std::string conlevels;
    float sig;
    if (!in->StatsDef()) in->setCubeStats();
    if (in->pars().getFlagUserThreshold()) sig=in->pars().getThreshold();
    else sig = 2.0*in->stat().getSpread();
    k=0;
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
    mfile=outfold+"pv.gnu";
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

#ifdef HAVE_GNUPLOT
    Gnuplot gp;
    gp.begin();
    mfile = "load '"+outfold+"gnuscript.gnu'";
    if (!in->pars().getflagGalMod()) gp.commandln(mfile.c_str());
    mfile = "load '"+outfold+"pv.gnu'";
//    gp.commandln(mfile.c_str());
    gp.end();
    remove ("cont.tab");
    remove ("contm.tab");
    remove ("cont_mod.tab");
    remove ("contm_mod.tab");
#endif

    remove ((outfold+"pv.txt").c_str());
    remove ((outfold+"pvm.txt").c_str());
    remove ((outfold+"pv_mod.txt").c_str());
    remove ((outfold+"pvm_mod.txt").c_str());
    remove ((outfold+"rcpv.txt").c_str());
    remove ((outfold+"pv.gnu").c_str());

}
template void Galfit<float>::plotParam(Cube<float>*, Cube<float>*);
template void Galfit<double>::plotParam(Cube<double>*, Cube<double>*);


template <class T>
void Galfit<T>::plotChanMaps() {

    std::ofstream py_file((in->pars().getOutfolder()+"maps.py").c_str());

    float crpix3_kms = in->Head().Crpix(2);
    float cdelt3_kms = DeltaVel<float>(in->Head());
    float crval3_kms = AlltoVel(in->Head().Crval(2),in->Head());

    float cont = 0;
    int xpos = findMedian(&outr->xpos[0],outr->nr);
    int ypos = findMedian(&outr->ypos[0],outr->nr);
    int xmin=0, ymin=0, zmin=0, disp=0;
    int xmax=in->DimX()-1, ymax=in->DimY()-1, zmax=in->DimZ()-1;
    float vsys_av = findMedian(&outr->vsys[0],outr->nr);

    if (in->pars().getMASK()=="SEARCH") {
        if (in->pars().getFlagUserGrowthThreshold()) cont = in->pars().getGrowthThreshold();
        else if (in->pars().getFlagUserThreshold()) cont = in->pars().getThreshold();
        else {
            if (in->pars().getFlagGrowth()) cont = in->pars().getGrowthCut()*in->stat().getSpread();
            else cont = in->pars().getCut()*in->stat().getSpread();
        }
        Detection<T> *larg = in->LargestDetection();
        int ext[4] = {fabs(xpos-lround(larg->getXmin()-2*in->Head().Bmaj()/in->Head().PixScale())),
                      fabs(xpos-lround(larg->getXmax()+2*in->Head().Bmaj()/in->Head().PixScale())),
                      fabs(ypos-lround(larg->getYmin()-2*in->Head().Bmaj()/in->Head().PixScale())),
                      fabs(ypos-lround(larg->getYmax()+2*in->Head().Bmaj()/in->Head().PixScale()))};
        disp = *max_element(&ext[0],&ext[0]+4);
        zmin = larg->getZmin();
        zmax = larg->getZmax();
    }
    else {
        if (!in->StatsDef()) in->setCubeStats();
        cont = 2.5*in->stat().getSpread();
        float inc_av = findMedian(&outr->inc[0],outr->nr);
        std::vector<T> maxv(outr->nr);
        for (int i=0; i<outr->nr; i++) maxv[i]=outr->vrot[i]*sin(outr->inc[i]*M_PI/180.)+outr->vdisp[i];
       //float max_vrot = *max_element(&outr->vrot[0],&outr->vrot[0]+outr->nr);
        float max_v=*max_element(&maxv[0],&maxv[0]+outr->nr);
        disp = fabs((outr->radii.back()/arcconv+2*in->Head().Bmaj())/in->Head().PixScale());
        int z_vsys = (vsys_av-crval3_kms)/cdelt3_kms+crpix3_kms-1;
        //int disp_v = ceil((1.5*max_vrot)*sin(inc_av*M_PI/180.)/fabs(DeltaVel<float>(in->Head())));
        int disp_v = ceil((1.5*max_v)/fabs(DeltaVel<float>(in->Head())));
        zmin = z_vsys-disp_v>0 ? z_vsys-disp_v : 0;
        zmax = z_vsys+disp_v<in->DimZ() ? z_vsys+disp_v : in->DimZ()-1;
    }

    xmin = xpos-disp>0 ? xpos-disp : 0;
    xmax = xpos+disp<in->DimX() ? xpos+disp : in->DimX()-1;
    ymin = ypos-disp>0 ? ypos-disp : 0;
    ymax = ypos+disp<in->DimY() ? ypos+disp : in->DimY()-1;
//    std::string range = "["+to_string(zmin)+":"+to_string(zmax)+","+to_string(ymin)+":"+to_string(ymax)+","+to_string(xmin)+":"+to_string(xmax)+"]";

    py_file << "import matplotlib \n"
            << "import matplotlib.pyplot as plt \n"
            << "import matplotlib.gridspec as gridspec \n"
            << "from astropy.io import fits \n"
            << "from astropy.visualization import PowerStretch \n"
            << "from astropy.visualization.mpl_normalize import ImageNormalize \n"
            << "from astropy.visualization import PercentileInterval \n\n";

    py_file << "# Setting all the needed variables \n"
            << "image = fits.open('" << in->pars().getImageFile() << "') \n"
            << "image_mod = fits.open('" << in->pars().getOutfolder() << in->Head().Name() << "mod.fits') \n"
            << "xmin = " << xmin << std::endl << "xmax = " << xmax << std::endl
            << "ymin = " << ymin << std::endl << "ymax = " << ymax << std::endl
            << "zmin = " << zmin << std::endl << "zmax = " << zmax << std::endl
            << "imagedata = image[0].data[";
    if (in->Head().NumAx()>3)
        for (int i=0; i<in->Head().NumAx()-3; i++) py_file << "0,";
    py_file << "zmin:zmax,ymin:ymax,xmin:xmax] \n"
            << "imagedata_mod = image_mod[0].data[zmin:zmax,ymin:ymax,xmin:xmax] \n"
            << "head = image[0].header \n"
            << "zsize=imagedata[:,0,0].size \n"
            << "cont = " << cont << std::endl
            << "v=[cont,2*cont,4*cont,8*cont,16*cont,32*cont,64*cont] \n"
            << "v_neg = [-cont] \n"
            << "interval = PercentileInterval(99.5) \n"
            << "vmax = interval.get_limits(imagedata)[1] \n"
            << "norm = ImageNormalize(vmin=cont, vmax=vmax, stretch=PowerStretch(0.5)) \n\n";

    py_file << "# Beginning plot \n"
            << "plt.figure(figsize=(8.27, 11.69), dpi=100) \n"
            << "grid = [gridspec.GridSpec(2,5),gridspec.GridSpec(2,5),gridspec.GridSpec(2,5)] \n"
            << "grid[0].update(top=0.95, bottom=0.695, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
            << "grid[1].update(top=0.65, bottom=0.395, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
            << "grid[2].update(top=0.35, bottom=0.095, left=0.05, right=0.95, wspace=0.0, hspace=0.0) \n"
            << "matplotlib.rcParams['contour.negative_linestyle'] = 'solid' \n"
            << "plt.rc('font',family='sans-serif',serif='Helvetica') \n\n";

    py_file << "num = 0 \n"
            << "for j in range (0,3): \n"
            << "\tfor i in range (0,5): \n"
            << "\t\tchan = num*(zsize)/15 \n"
            << "\t\tz = imagedata[chan,:,:] \n"
            << "\t\tz_mod = imagedata_mod[chan,:,:] \n"
            << "\t\tvelo_kms = (chan+1-" << crpix3_kms-zmin << ")*" << cdelt3_kms << "+" << crval3_kms << std::endl
            << "\t\tvelo = ' v = ' + str(int(velo_kms)) + ' km/s' \n"
            << "\t\tplt.subplot(grid[j][0,i]) \n"
            << "\t\tplt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='off') \n"
            << "\t\tplt.title(velo, fontsize=8,loc='left') \n"
            << "\t\tplt.imshow(z,origin='lower',cmap = matplotlib.cm.Greys,norm=norm) \n"
            << "\t\tplt.contour(z,v,origin='lower',linewidths=0.7,colors='#00008B') \n"
            << "\t\tplt.contour(z,v_neg,origin='lower',linewidths=0.1,colors='gray') \n"
            << "\t\tplt.subplot(grid[j][1,i]) \n"
            << "\t\tplt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='off') \n"
            << "\t\tplt.imshow(z_mod,origin='lower',cmap = matplotlib.cm.Greys,norm=norm ) \n"
            << "\t\tplt.contour(z_mod,v,origin='lower',linewidths=0.7,colors='#B22222') \n"
            << "\t\tnum = num+1 \n\n"
            << "plt.savefig('" << in->pars().getOutfolder() << "chanmaps.pdf', orientation = 'portrait', format = 'pdf') \n"
            << "image.close() \n"
            << "image_mod.close() \n";

    py_file << std::endl << std::endl;

    float zmin_wcs = AlltoVel(in->getZphys(zmin),in->Head());
    float zmax_wcs = AlltoVel(in->getZphys(zmax),in->Head());
    int pa_av = lround(findMedian(&outr->phi[0],outr->nr));
    int pa_min = pa_av+90<360 ? pa_av+90 : pa_av-90;
    if (zmin_wcs>zmax_wcs) std::swap(zmin_wcs,zmax_wcs);

    py_file << "# Now plotting the position-velocity diagrams \n"
            << "image_maj     = fits.open('"<< in->pars().getOutfolder() << in->Head().Name() << "_pv_a.fits') \n"
            << "image_mod_maj = fits.open('"<< in->pars().getOutfolder() << in->Head().Name() << "mod_pv_a.fits') \n"
            << "image_min     = fits.open('"<< in->pars().getOutfolder() << in->Head().Name() << "_pv_b.fits') \n"
            << "image_mod_min = fits.open('"<< in->pars().getOutfolder() << in->Head().Name() << "mod_pv_b.fits') \n"
            << "head = image_maj[0].header \n"
            << "xmin = head['CRPIX1']-1 - " << disp << std::endl
            << "xmax = head['CRPIX1']-1 + " << disp << std::endl
            << "if xmin<0: xmin=0 \n"
            << "if xmax>=head['NAXIS1']: xmax=head['NAXIS1']-1 \n"
            << "imagedata_min = image_min[0].data[zmin:zmax,xmin:xmax] \n"
            << "imagedata_mod_min = image_mod_min[0].data[zmin:zmax,xmin:xmax] \n"
            << "imagedata_maj = image_maj[0].data[zmin:zmax,xmin:xmax] \n"
            << "imagedata_mod_maj = image_mod_maj[0].data[zmin:zmax,xmin:xmax] \n"
            << "xmin_wcs = ((xmin+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1'])*" << arcconv << std::endl
            << "xmax_wcs = ((xmax+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1'])*" << arcconv << std::endl
            << "zmin_wcs = " << zmin_wcs << std::endl
            << "zmax_wcs = " << zmax_wcs << std::endl;

    std::string radius = "radius = [";
    std::string vrot = "vrot = [";
    for (int i=in->pars().getStartRad(); i<outr->nr; i++) {
        float vel1 = (outr->vrot[i]*sin(outr->inc[i]*M_PI/180.))+outr->vsys[i];
        float vel2 = outr->vsys[i]-(outr->vrot[i]*sin(outr->inc[i]*M_PI/180.));
        float rad = outr->radii[i];
        bool reverse = (pa_av>0 && pa_av<90) || (pa_av>180 && pa_av<270);
        if (reverse) std::swap(vel1,vel2);
        if (i==0) rad += (outr->radsep/4.);
        radius += (to_string(rad,1)+","+to_string(-rad,1)+",");
        vrot += (to_string(vel1,1)+","+to_string(vel2,1)+",");
    }
    radius.erase(radius.size()-1);
    vrot.erase(vrot.size()-1);
    radius += "]";
    vrot += "]";

    py_file << radius << std::endl << vrot << std::endl
            << "ext = [xmin_wcs,xmax_wcs,zmin_wcs,zmax_wcs] \n"
            << "plt.figure(figsize=(8.27, 11.69), dpi=100) \n"
            << "grid = gridspec.GridSpec(2,1) \n"
            << "grid.update(top=0.95, bottom=0.05, left=0.10, right=0.90) \n"
            << "plt.subplot(grid[0,0]) \n"
            << "plt.xlabel('Offset (arcsec)',fontsize=14) \n"
            << "plt.ylabel('V$_\\mathrm{LOS}$ (km/s)',fontsize=14) \n"
            << "plt.title('$\\phi = " << pa_av <<"^\\circ$',loc='right',fontsize=20) \n"
            << "plt.imshow(imagedata_maj,origin='lower',cmap = matplotlib.cm.Greys,norm=norm,extent=ext,aspect='auto') \n"
            << "plt.contour(imagedata_maj,v,origin='lower',linewidths=0.7,colors='#00008B',extent=ext,aspect='auto') \n"
            << "plt.contour(imagedata_maj,v_neg,origin='lower',linewidths=0.1,colors='gray',extent=ext,aspect='auto') \n"
            << "plt.contour(imagedata_mod_maj,v,origin='lower',linewidths=0.7,colors='#B22222',extent=ext,aspect='auto') \n"
            << "plt.axhline(y=" << vsys_av << ",color='black') \n"
            << "plt.axvline(x=0,color='black') \n"
            << "plt.plot(radius,vrot,'yo') \n"
            << "plt.subplot(grid[1,0]) \n"
            << "plt.xlabel('Offset (arcsec)',fontsize=14) \n"
            << "plt.ylabel('V$_\\mathrm{LOS}$ (km/s)',fontsize=14) \n"
            << "plt.title('$\\phi = " << pa_min << "^\\circ$',loc='right',fontsize=20) \n"
            << "plt.imshow(imagedata_min,origin='lower',cmap = matplotlib.cm.Greys,norm=norm,extent=ext,aspect='auto') \n"
            << "plt.contour(imagedata_min,v,origin='lower',linewidths=0.7,colors='#00008B',extent=ext,aspect='auto') \n"
            << "plt.contour(imagedata_min,v_neg,origin='lower',linewidths=0.1,colors='gray',extent=ext,aspect='auto') \n"
            << "plt.contour(imagedata_mod_min,v,origin='lower',linewidths=0.7,colors='#B22222',extent=ext,aspect='auto') \n"
            << "plt.axhline(y=" << vsys_av << ",color='black') \n"
            << "plt.axvline(x=0,color='black') \n"
            << "plt.savefig('" << in->pars().getOutfolder() << "pv_diagrams.pdf', orientation = 'portrait', format = 'pdf') \n"
            << "image_maj.close() \n"
            << "image_min.close() \n"
            << "image_mod_maj.close() \n"
            << "image_mod_min.close() \n";

    py_file.close();

#ifdef HAVE_PYTHON
    std::string cmd = "python -W ignore "+in->pars().getOutfolder()+"maps.py";
    system(cmd.c_str());
#endif

}
template void Galfit<float>::plotChanMaps();
template void Galfit<double>::plotChanMaps();


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
#undef MAXPAR
