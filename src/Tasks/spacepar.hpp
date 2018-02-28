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

#ifndef SPACEPAR_HH_
#define SPACEPAR_HH_

#include <iostream> 
#include <cfloat> 
#include <Arrays/cube.hh>
#include <Tasks/galfit.hh>
#include <Utilities/utils.hh>
#include <Utilities/gnuplot.hh>
#include <Utilities/progressbar.hh>

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


using namespace Model;

template <class T>
class Spacepar : public Galfit<T>
{   
public:
    Spacepar(Cube<T> *c);
    ~Spacepar(){if (parspAll) delete parspace;}
    
    void calculate();
    void plotAll_Python();
    
private:    

    Cube<T> *parspace;  // A cube object with all parameter spaces
    bool    parspAll;   // Wheter parspace has been allocated.
    std::string p1;     // Type of parameter 1
    std::string p2;     // Type of parameter 2
    T   minp1;          // Minimum of parameter 1
    T   minp2;          // Minimum of parameter 2
    T   maxp1;          // Maximum of parameter 1
    T   maxp2;          // Maximum of parameter 2
    T   delp1;          // Step size for parameter 1
    T   delp2;          // Step size for parameter 2
    
    void plotmin_Gnuplot(int n);
    void setOutSpacepar(int p1_nsteps, int p2_nsteps);
};


template <class T>
Spacepar<T>::Spacepar(Cube<T> *c) : Galfit<T>::Galfit(c) {
    
    p1 = Galfit<T>::in->pars().getP1();
    p2 = Galfit<T>::in->pars().getP2();
    p1 = makeupper(p1);
    p2 = makeupper(p2);

    if (p1!="VROT" && p1!="VDISP" && p1!="Z0"   && p1!="INC" &&
        p1!="PA"   && p1!="XPOS"  && p1!="YPOS" && p1!="VSYS"&& p1!="VRAD") {
        std::cout << "SPACEPAR ERROR: Unknown parameter 1 "
                  << p1 << std::endl;
        std::terminate();
    }
    if (p2!="VROT" && p2!="VDISP" && p2!="Z0"   && p2!="INC" &&
        p2!="PA"   && p2!="XPOS"  && p2!="YPOS" && p2!="VSYS" && p2!="VRAD") {
        std::cout << "SPACEPAR ERROR: Unknown parameter 2 "
                  << p2 << std::endl;
        std::terminate();
    }
    if (p1==p2) {
        std::cout << "SPACEPAR ERROR: p1=p2!!" << std::endl;
        std::terminate();
    }
    
    for (int i=0; i<MAXPAR; i++) Galfit<T>::mpar[i]=false;

    minp1 = Galfit<T>::in->pars().getP1p(0);
    minp2 = Galfit<T>::in->pars().getP2p(0);
    maxp1 = Galfit<T>::in->pars().getP1p(1);
    maxp2 = Galfit<T>::in->pars().getP2p(1);
    delp1 = Galfit<T>::in->pars().getP1p(2);
    delp2 = Galfit<T>::in->pars().getP2p(2);
    
    // Changing default GALFIT limits for parameters
    Galfit<T>::maxs[VROT] = 10000;
    Galfit<T>::mins[VROT] = 0;
    Galfit<T>::maxs[VDISP] = 10000;
    Galfit<T>::mins[VDISP] = 0;
    Galfit<T>::maxs[Z0] = 10000;
    Galfit<T>::mins[Z0] = 0;
    Galfit<T>::maxs[INC] = 90;
    Galfit<T>::mins[INC] = 0;
    Galfit<T>::maxs[PA] = 360;
    Galfit<T>::mins[PA] = 0;
    Galfit<T>::maxs[XPOS] = Galfit<T>::in->DimX()-1;
    Galfit<T>::mins[XPOS] = 0;
    Galfit<T>::maxs[YPOS] = Galfit<T>::in->DimY()-1;
    Galfit<T>::mins[YPOS] = 0;
    Galfit<T>::maxs[VRAD] = 10000;
    Galfit<T>::mins[VRAD] = -10000;

    std::string cunit2 = Galfit<T>::in->Head().Cunit(2);
    Galfit<T>::maxs[VSYS] = AlltoVel(Galfit<T>::in->getZphys(Galfit<T>::in->DimZ()-1), Galfit<T>::in->Head());
    Galfit<T>::mins[VSYS] = AlltoVel(Galfit<T>::in->getZphys(0), Galfit<T>::in->Head());
    
    parspAll = false;
}


template <class T>
void Spacepar<T>::calculate() {
    
    bool verb = Galfit<T>::in->pars().isVerbose();
    bool debug = Galfit<T>::in->pars().getFlagDebug();
    
    cout << showpoint << fixed << setprecision(2) << endl ;
    if (verb) { 
        Galfit<T>::in->pars().setVerbosity(false);
        cout << setfill('=') << setw(40) << right << " SPACEPAR " << setw(34) << " ";
        cout << setfill(' ') << endl << endl;
        cout << " Exploring parameter space: \n"
             << "  " << left << setw(5) << p1 << " in [" << minp1 << ", " << maxp1 << "] with delta = " << delp1 << std::endl
             << "  " << left << setw(5) << p2 << " in [" << minp2 << ", " << maxp2 << "] with delta = " << delp2;
        cout << endl << endl;
    }
     
    std::string filename = "./output/ring.dat";
    std::ofstream file;
        
    // Filling vectors with steps in the two parameters
    int p1_nsteps = abs(maxp1-minp1)/delp1 + 1 ;
    int p2_nsteps = abs(maxp2-minp2)/delp2 + 1 ;
    T *p1steps = new T[p1_nsteps];
    T *p2steps = new T[p2_nsteps];
    for (int i=0; i<p1_nsteps;i++) p1steps[i] = minp1+i*delp1;
    for (int i=0; i<p2_nsteps;i++) p2steps[i] = minp2+i*delp2;
    
    // Setting output FITS File
    setOutSpacepar(p1_nsteps,p2_nsteps);

    int start_rad = Galfit<T>::par.STARTRAD < Galfit<T>::inr->nr ? Galfit<T>::par.STARTRAD : 0;    
    // Starting loop over rings
    for (int ir=start_rad; ir<Galfit<T>::inr->nr; ir++) {
        
        file.open(filename.c_str());
        if (verb) 
            cout << "\n Working on ring #" << ir+1 << " at radius " 
                 << Galfit<T>::inr->radii[ir] << " arcsec: \n";  
        
        T p1min=0, p2min=0;
        T funmin=FLT_MAX;
        
        ProgressBar bar("   Calculating models...",true);
        bar.setShowbar(Galfit<T>::in->pars().getShowbar());

        int nthreads = Galfit<T>::in->pars().getThreads();
        // Starting loops over parameters
#pragma omp parallel num_threads(nthreads) shared(funmin,p1min,p2min)
{
        // NB: THIS LOOP IN PARALLEL CAUSES SOMETIMES SEGFAULT (in fftw_execute)
        if(verb) bar.init(p1_nsteps);
#pragma omp for 
        for (int i=0; i<p1_nsteps; i++) {
            if (verb) bar.update(i+1);
            for (int j=0; j<p2_nsteps; j++) {

                Rings<T> *dring = new Rings<T>;
                Rings<T> *inri = Galfit<T>::inr;
                
                dring->nr = 2;
                float width1=0, width2=0;
                if (ir==0) width1 = width2 = (inri->radii[1]-inri->radii[0])/2.;
                else if (ir==inri->nr-1) width1 = width2 = (inri->radii[ir]-inri->radii[ir-1])/2.;
                else {
                    width1 = (inri->radii[ir]-inri->radii[ir-1])/2.;
                    width2 = (inri->radii[ir+1]-inri->radii[ir])/2.;
                }

                dring->radii.push_back(max(double(inri->radii[ir]-width1),0.));
                dring->radii.push_back(max(double(inri->radii[ir]+width2),0.));

                for (int k=0; k<dring->nr; k++) {
                    if (p1=="VROT") dring->vrot.push_back(p1steps[i]);
                    else if (p2=="VROT") dring->vrot.push_back(p2steps[j]);
                    else dring->vrot.push_back(inri->vrot[ir]);
                    
                    if (p1=="VDISP") dring->vdisp.push_back(p1steps[i]);
                    else if (p2=="VDISP") dring->vdisp.push_back(p2steps[j]);
                    else dring->vdisp.push_back(inri->vdisp[ir]);
                    
                    if (p1=="Z0") dring->z0.push_back(p1steps[i]);
                    else if (p2=="Z0") dring->z0.push_back(p2steps[j]);
                    else dring->z0.push_back(inri->z0[ir]);
                    
                    if (p1=="INC") dring->inc.push_back(p1steps[i]);
                    else if (p2=="INC") dring->inc.push_back(p2steps[j]);
                    else dring->inc.push_back(inri->inc[ir]);
                    
                    if (p1=="PA") dring->phi.push_back(p1steps[i]);
                    else if (p2=="PA") dring->phi.push_back(p2steps[j]);
                    else dring->phi.push_back(inri->phi[ir]);
                    
                    if (p1=="XPOS") dring->xpos.push_back(p1steps[i]);
                    else if (p2=="XPOS") dring->xpos.push_back(p2steps[j]);
                    else dring->xpos.push_back(inri->xpos[ir]);
                    
                    if (p1=="YPOS") dring->ypos.push_back(p1steps[i]);
                    else if (p2=="YPOS") dring->ypos.push_back(p2steps[j]);
                    else dring->ypos.push_back(inri->ypos[ir]);
                    
                    if (p1=="VSYS") dring->vsys.push_back(p1steps[i]);
                    else if (p2=="VSYS") dring->vsys.push_back(p2steps[j]);
                    else dring->vsys.push_back(inri->vsys[ir]);

                    if (p1=="VRAD") dring->vrad.push_back(p1steps[i]);
                    else if (p2=="VRAD") dring->vrad.push_back(p2steps[j]);
                    else dring->vrad.push_back(inri->vrad[ir]);
                    
                    dring->dens.push_back(inri->dens[ir]);
                    dring->vvert.push_back(inri->vvert[ir]);
                    dring->dvdz.push_back(inri->dvdz[ir]);
                    dring->zcyl.push_back(inri->zcyl[ir]);
                }
                
                T par[MAXPAR];
                par[VROT] = dring->vrot.back();
                par[VDISP] = dring->vdisp.back();
                par[DENS] = dring->dens.back();
                par[Z0] = dring->z0.back();
                par[INC] = dring->inc.back();
                par[PA] = dring->phi.back();
                par[XPOS] = dring->xpos.back();
                par[YPOS] = dring->ypos.back();
                par[VSYS] = dring->vsys.back();
                par[VRAD] = dring->vrad.back();
                     
                T minfunc = Galfit<T>::func3D(dring, par);
                parspace->Array(i,j,ir) = minfunc;
                delete dring;
                
#pragma omp critical (spacepar_min) 
{
                if (minfunc<funmin) {
                    funmin = minfunc;
                    p1min = p1steps[i];
                    p2min = p2steps[j];
                }
                file << p1steps[i] << "  " << p2steps[j] << "  " << minfunc << endl;
}
            }
#pragma omp critical (spacepar_endl) 
            file << endl;
        }
}
        file.close();
        
        if (debug) plotmin_Gnuplot(ir);

        if (verb) {
            bar.fillSpace("OK.\n");
            std::cout << "   Minimum (fmin=" << setprecision(6) << funmin << ") at: "
                      << p1 << " = " << setprecision(1) << p1min << ", " 
                      << p2 << " = " << p2min << std::endl;
        }
    
    }

    remove(filename.c_str());
    Galfit<T>::in->pars().setVerbosity(verb);
    
    delete [] p1steps;
    delete [] p2steps;
}


template <class T>
void Spacepar<T>::plotmin_Gnuplot (int nr) {

#ifdef HAVE_GNUPLOT
    Gnuplot gp;
    gp.begin(); 
    gp.commandln("set terminal png enhanced");
    //gp.commandln("set terminal postscript eps color enhanced");
    gp.commandln("set pm3d map");
    gp.commandln("set size square");
    gp.commandln("unset key");
    std::string cmd = "set xrange ["+to_string(minp1)+":"+to_string(maxp1)+"]";
    gp.commandln(cmd.c_str());
    cmd = "set yrange ["+to_string(minp2)+":"+to_string(maxp2)+"]"; 
    gp.commandln(cmd.c_str());
    cmd = Galfit<T>::in->pars().getOutfolder()+"r"+to_string<int>(nr+1)+".png";
    gp.commandln(("set output '"+cmd+"'").c_str());
    if       (p1=="VROT")   cmd="set xlabel 'V_{rot}  [Km/s]'";
    else if (p1=="VDISP") cmd="set xlabel 'Dispersion [km/s]'";
    else if (p1=="Z0")  cmd="set xlabel 'Scale height [arcs]'";
    else if (p1=="INC")     cmd="set xlabel 'Inclination [degree]'";
    else if (p1=="PA")  cmd="set xlabel 'Position angle [degree]'";
    else if (p1=="XPOS")    cmd="set xlabel 'Xcenter [pixel]'";
    else if (p1=="YPOS")    cmd="set xlabel 'Ycenter [pixel]'";
    else if (p1=="VSYSP")   cmd="set xlabel 'Systemic velocity [km/s]'";
    else if (p1=="VRAD")   cmd="set xlabel 'Radial velocity [km/s]'";
    gp.commandln(cmd.c_str());
    if       (p2=="VROT")  cmd="set ylabel 'V_rot  [Km/s]'";
    else if (p2=="VDISP") cmd="set ylabel 'Dispersion [km/s]'";
    else if (p2=="Z0")  cmd="set ylabel 'Scale height [arcs]'";
    else if (p2=="INC")     cmd="set ylabel 'Inclination [degree]'";
    else if (p2=="PA")  cmd="set ylabel 'Position angle [degree]'";
    else if (p2=="XPOS")    cmd="set xlabel 'Xcenter [pixel]'";
    else if (p2=="YPOS")    cmd="set xlabel 'Ycenter [pixel]'";
    else if (p2=="VSYS")   cmd="set xlabel 'Systemic velocity [km/s]'";
    else if (p2=="VRAD")   cmd="set xlabel 'Radial velocity [km/s]'";
    gp.commandln(cmd.c_str());
    
    cmd = "set title 'Ring #"+to_string(nr+1)+" ("+
          to_string(Galfit<T>::inr->radii[nr])+" arcs)'";
    gp.commandln(cmd.c_str());
    
    cmd = "splot './output/ring.dat'";
    gp.commandln(cmd.c_str());

    gp.end();
#endif  
}

template <class T>
void Spacepar<T>::plotAll_Python() {
    
    if (!parspAll) return;
    
    bool verb = Galfit<T>::in->pars().isVerbose();
    
    // Writing fitsfile with parameter space
    std::string outfold = Galfit<T>::in->pars().getOutfolder();
    std::string object  = Galfit<T>::in->Head().Name();
    std::string fname   = outfold+"spacepar_"+object+".fits";
    parspace->fitswrite_3d(fname.c_str());
    
    // Writing python script for output plots
    std::string scriptname = "spacepar.py";
    std::string fout = outfold+"spacepar_"+object+".pdf";
    std::ofstream py_file((outfold+scriptname).c_str());
    
    py_file << "# A script to produce output plots for SPACEPAR task.\n"
            << "import numpy as np\n" 
            << "import matplotlib\n"
            << "import matplotlib.pyplot as plt\n"
            << "from astropy.io import fits\n" 
            << "matplotlib.rc('xtick',direction='in')\n" 
            << "matplotlib.rc('ytick',direction='in')\n" 
            << "matplotlib.rc('font',family='sans-serif',serif='Helvetica',size=10)\n\n"

            << "f = fits.open('"<< fname << "')[0]\n"
            << "d, h = f.data, f.header\n"
            << "p1, p2 = h['CTYPE1'], h['CTYPE2']\n"
            << "p1u, p2u, ru = h['CUNIT1'], h['CUNIT2'], h['CUNIT3']\n"
            << "p1ran = (np.arange(0,d.shape[2])+1-h['CRPIX1'])*h['CDELT1']+h['CRVAL1']\n"
            << "p2ran = (np.arange(0,d.shape[1])+1-h['CRPIX2'])*h['CDELT2']+h['CRVAL2']\n"
            << "rings = (np.arange(0,d.shape[0])+1-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']\n\n"

            << "nrad = d.shape[0]\n"
            << "ncols = 4\n"
            << "nrows = int(np.ceil(nrad/float(ncols)))\n"
            << "ext = [p1ran[0],p1ran[-1],p2ran[0],p2ran[-1]]\n"
            << "cmap = plt.get_cmap('nipy_spectral') #plt.get_cmap('gnuplot')\n\n"

            << "fig = plt.figure(figsize=(8,8))\n"
            << "x_axis_len, y_axis_len = 0.27, 0.27 \n"
            << "x_sep, y_sep = 0.07, 0.08 \n"
            << "count = 0\n"
            << "axis, bottom_corner = [], [0.1,0.7]\n"
            << "for i in range (nrows):\n"
            << "\tbottom_corner[0] = 0.1\n"
            << "\tfor j in range (ncols):\n"
            << "\t\tif (count>=nrad): break\n"
            << "\t\taxis.append(fig.add_axes([bottom_corner[0],bottom_corner[1],x_axis_len,y_axis_len]))\n"
            << "\t\tbottom_corner[0]+=x_axis_len+x_sep\n"
            << "\t\tcount += 1\n" 
            << "\tbottom_corner[1]-=(y_axis_len+y_sep)\n\n"

            << "for i in range (nrad):\n"
            << "\tnr = int(i/ncols) + 1\n"
            << "\tnc = i - (nr-1)*ncols + 1\n"   
            << "\ttoplot = d[i]\n"
            << "\ta = np.unravel_index(np.argmin(toplot),toplot.shape)\n"
            << "\tp1min, p2min = p1ran[a[1]], p2ran[a[0]]\n"
            << "\tax = axis[i]\n"
            << "\tax.set_xlim(ext[0],ext[1])\n"
            << "\tax.set_ylim(ext[2],ext[3])\n"
            << "\tax.imshow(toplot,origin='lower',extent=ext,aspect='auto',cmap=cmap)\n"
            << "\tax.plot(p1min,p2min,'x',mew=2,ms=8,c='w')\n"
            << "\tradstr = 'R = %s %s'%(rings[i],ru)\n"
            << "\tminstr = 'min = (%.1f %s, %.1f %s)'%(p1min,p1u,p2min,p2u)\n"
            << "\tax.text(0.01,1.1,radstr,transform=ax.transAxes)\n"
            << "\tax.text(0.01,1.03,minstr,transform=ax.transAxes)\n"
            << "\tif nc==1: ax.set_ylabel(p2+' ('+p2u+')')\n"
            << "\tif (nr==nrows) or (nr==nrows-1 and nrad%ncols!=0 and nc>nrad%ncols):\n" 
            << "\t\tax.set_xlabel(p1+' ('+p1u+')')\n\n"
                
            << "fig.savefig('"+fout+"',bbox_inches='tight')\n";
    
#ifdef HAVE_PYTHON
    std::string cmd = "python "+outfold+scriptname;//+" > /dev/null 2>&1";
    if (Galfit<T>::in->pars().getFlagPlots()) {
        if (verb) std::cout << "\n\nProducing stunning parameter space plots..." << std::flush;
        int ret = system(cmd.c_str());
        if (verb) {
            if (ret==0) std::cout << " Done.\n";
            else std::cout << " Something went wrong! Check spacepar.py in the output folder.\n";
        }
    }
#endif
    
}

template <class T>
void Spacepar<T>::setOutSpacepar(int p1_nsteps, int p2_nsteps) {
    
    // Allocate output cube containing all parameter spaces
    if (parspAll) delete parspace;
    int start_rad = Galfit<T>::par.STARTRAD < Galfit<T>::inr->nr ? Galfit<T>::par.STARTRAD : 0;
    int dims[3] = {p1_nsteps,p2_nsteps,Galfit<T>::inr->nr};
    parspace = new Cube<T>(dims);
    Header h;
    h.setNumAx(3);
    h.setDimAx(0,dims[0]);
    h.setDimAx(1,dims[1]);
    h.setDimAx(2,dims[2]);
    h.setCrpix(0,1);
    h.setCrpix(1,1);
    h.setCrpix(2,1);
    h.setCrval(0,minp1);
    h.setCrval(1,minp2);
    h.setCrval(2,Galfit<T>::inr->radii[start_rad]);
    h.setCdelt(0,delp1);
    h.setCdelt(1,delp2);
    h.setCdelt(2,Galfit<T>::inr->radsep);
    h.setCtype(0,p1);
    h.setCtype(1,p2);
    h.setCtype(2,"RING_RAD");
    h.setCunit(0,"km/s");
    if (p1=="INC" || p1=="PA") h.setCunit(0,"deg");
    else if (p1=="X0" || p1=="Y0") h.setCunit(0,"pix");
    else if (p1=="Z0") h.setCunit(0,"arcsec");   
    h.setCunit(1,"km/s");
    if (p2=="INC" || p2=="PA") h.setCunit(1,"deg");
    else if (p2=="X0" || p2=="Y0") h.setCunit(1,"pix");
    else if (p2=="Z0") h.setCunit(1,"arcsec");
    h.setCunit(2,"arcsec");
    h.setBtype("FMIN_"+to_string<int>(Galfit<T>::par.FTYPE));
    h.setBunit(Galfit<T>::in->Head().Bunit());
    parspace->saveHead(h);
    parspAll = true;
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

#endif
