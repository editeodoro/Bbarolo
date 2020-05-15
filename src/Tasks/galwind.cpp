// -----------------------------------------------------------------------
// galwind.cpp: Functions of the GalWind class.
// -----------------------------------------------------------------------

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
#include <string>
#include <cmath>
#include <vector>
#include <Tasks/galwind.hh>
#include <Tasks/moment.hh>
#include <Tasks/smooth3D.hh>
#include <Tasks/galmod.hh>
#include <Utilities/utils.hh>
#include <Utilities/progressbar.hh>


template <class T>
GalWind<T>::GalWind(Cube<T> *c, T x0, T y0, T pa, T inc, T disp, T dens, T vsys, T vw, 
                    T openang, T htot, int denstype, int ntot, int cdens, int nv, int nthreads) {
            
    in          = c;
    par.XPOS    = to_string(x0);
    par.YPOS    = to_string(y0);
    par.PHI     = to_string(pa);
    par.INC     = to_string(inc);
    par.VDISP   = to_string(disp);
    par.DENS    = to_string(dens);
    par.VSYS    = to_string(vsys);
    par.VWIND   = to_string(vw);
    par.OPENANG = to_string(openang);
    par.HTOT    = htot;
    par.DENSTYPE= denstype;
    par.NTOT    = ntot;
    par.CDENS   = cdens;
    par.NV      = nv;
    in->pars().setThreads(nthreads);

    in->checkBeam();
}


template <class T>
GalWind<T>& GalWind<T>::operator= (const GalWind<T>& gw) {
    
    if(this==&gw) return *this;
    
    this->in  = gw.in;
    this->par = gw.par;
    if(this->outDefined) delete this->out;
    this->outDefined = gw.outDefined;
    if(this->outDefined) *this->out = *gw.out;
    
    return *this;
}


template <class T>
bool GalWind<T>::compute() {
    
    // Allocate output datacube
    if (outDefined) delete out;
    out = new Cube<T>(in->AxisDim());
    out->saveHead(in->Head());
    outDefined = true;
    for (size_t i=0; i<out->NumPix(); i++) out->Array(i) = 0;
    
    if (par.WTYPE==0) return compute_cylindrical();
    else return compute_spherical();
    
}


template <class T>
bool GalWind<T>::compute_spherical() {
    
    // This function computes the model a spherical geometry, i.e.
    // decomposing the outflow in spherical shells. Quantities are defined 
    // for each shell at increasing radius from the galactic center.
    
    int nv = par.NV>1 ? par.NV : 20;
    
    // Fixed parameters
    float PA   = atof(par.PHI.c_str());  
    float inc  = atof(par.INC.c_str());   
    float Vsys = atof(par.VSYS.c_str());
    // Getting the center position 
    string pos[2] = {par.XPOS, par.YPOS};
    double *pixs = getCenterCoordinates(pos, in->Head());
    float x0  = pixs[0], y0=pixs[1];
    
    // Parameters that can vary shell by shell
    std::vector<T> Radii, Vwind, Vrot, Dens, Vdisp, Oang;
    bool radii_b = getDataColumn(Radii,par.RADII);
    bool vwind_b = getDataColumn(Vwind,par.VWIND);
    bool vdisp_b = getDataColumn(Vdisp,par.VDISP);
    bool dens_b  = getDataColumn(Dens,par.DENS);
    bool opa_b   = getDataColumn(Oang,par.OPENANG);
    bool vrot_b  = getDataColumn(Vrot,par.OPENANG);
    
    // Setting number of shells and radsep
    size_t nr = par.NTOT;
    double radsep = par.HTOT/nr;    
    if (radii_b) {
        for (unsigned i=1; i<Radii.size()-1; i++)
            radsep += Radii[i+1]-Radii[i];
        radsep/=(Radii.size()-2);
        nr = Radii.size();
    }

    // Initializing shells
    Shells<T> *inS = new Shells<T>;
    inS->ns = nr;
    inS->sep = radsep;
    for (int i=0; i<inS->ns; i++) {
        if (radii_b) inS->radii.push_back(Radii[i]);
        else inS->radii.push_back(i*radsep+radsep/2.);
        if (vrot_b) inS->vrot.push_back(Vrot[i]);
        else inS->vrot.push_back(atof(par.VROT.c_str()));
        if (vdisp_b) inS->vdisp.push_back(Vdisp[i]);
        else inS->vdisp.push_back(atof(par.VDISP.c_str()));
        if (vwind_b) inS->vsph.push_back(Vwind[i]);
        else inS->vsph.push_back(atof(par.VWIND.c_str()));
        if (dens_b) inS->dens.push_back(Dens[i]*1.E20);
        else inS->dens.push_back(atof(par.VWIND.c_str())*1.E20);
        if (opa_b) inS->openang.push_back(Oang[i]);
        else inS->openang.push_back(atof(par.OPENANG.c_str()));
        inS->xpos.push_back(x0);
        inS->ypos.push_back(y0);
        inS->inc.push_back(inc);
        inS->pa.push_back(PA);
        inS->vsys.push_back(Vsys);
    }
    
    int Blo[2] = {0,0};
    int Bup[2] = {in->DimX(),in->DimY()};
    
    Model::Galmod_wind<T> *w = new Model::Galmod_wind<T>;
    w->input(in, Bup, Blo, inS, nv, par.LTYPE, 1, par.CDENS);
    w->calculate();
    
    for (size_t i=0; i<in->NumPix(); i++) 
        out->Array(i) += w->Out()->Array(i);
    
    delete inS;
    delete w;
     
    return true;
}


template <class T>
bool GalWind<T>::compute_cylindrical() {
    
    // This function computes the model using Lelli's approach, i.e.
    // decomposing the outflow in cylinders. Quantities are defined for 
    // each cylinder at increasing height from the galactic plane.
    
    // Collecting needed parameters
    float arcconv = arcsconv(in->Head().Cunit(0));
    float PixSize = in->Head().PixScale()*arcconv;
    int nv = par.NV>1 ? par.NV : 20;
    float OpA  = atof(par.OPENANG.c_str());
    float z0 = par.HTOT/(2*par.NTOT);  //half-height of each cylinder in arcsec
    float zpix = z0/PixSize;           //half-height of each cylinder in M_PIxels
    float width = (2*z0)*tan(OpA/2.*M_PI/180.); //width of each ring (Rout/k)
    
    // Fixed parameter
    float PA   = atof(par.PHI.c_str());  
    float inc  = atof(par.INC.c_str());   
    float Vsys = atof(par.VSYS.c_str());
    // Getting the center position 
    string pos[2] = {par.XPOS, par.YPOS};
    double *pixs = getCenterCoordinates(pos, in->Head());
    float x0  = pixs[0], y0=pixs[1];
    
    // Parameters that can vary cylinder by cylinder
    std::vector<T> Vwind, Dens, Vdisp, Vrot;
    bool vwind_b = getDataColumn(Vwind,par.VWIND);
    bool vdisp_b = getDataColumn(Vdisp,par.VDISP);
    bool dens_b  = getDataColumn(Dens,par.DENS);
    bool vrot_b  = getDataColumn(Dens,par.DENS);
    for (size_t i=0; i<par.NTOT; i++) {
        // If only 1 value is given, fill vectors with this value
        if (!vwind_b) Vwind.push_back(atof(par.VWIND.c_str()));
        if (!vdisp_b) Vdisp.push_back(atof(par.VDISP.c_str()));
        if (!dens_b)  Dens.push_back(atof(par.DENS.c_str()));
        if (!vrot_b)  Vrot.push_back(atof(par.VROT.c_str()));
        
    }
    // Checking array dimensions. They must be equal to NTOT
    if (Vwind.size()!=par.NTOT || Vdisp.size()!=par.NTOT || 
        Dens.size()!=par.NTOT || Vrot.size()!=par.NTOT) {
        std::cerr << " GALWIND ERROR: VWIND, VDISP and DENS must have size=NTOT!\n";
        return false;
    }
     

    // Starting progress bar
    bool verb = in->pars().isVerbose();
    if (verb) in->pars().setVerbosity(false);
    ProgressBar bar(" Generating outflow model... ",false);
    bar.setShowbar(in->pars().getShowbar());
    
    int nthreads = in->pars().getThreads();
#pragma omp parallel num_threads(nthreads)
{
    // Allocating Galmod objects for the two cylinders
    Model::Galmod<T> *con1 = new Model::Galmod<T>;
    Model::Galmod<T> *con2 = new Model::Galmod<T>;
    if (verb) bar.init(par.NTOT);
    
#pragma omp for schedule(dynamic)
    // Start main loop
    for (size_t k=1; k<(par.NTOT+1); k++){
        if (verb) bar.update(k);

        // Assign distribution and kinematics to cylinder k
        float Rin  = par.DENSTYPE==2 ? (k-1)*width : 0;         // inner ring
        float Rout = float(k)*width;                            // outer ring
        float dens = par.DENSTYPE==0 ? Dens[k-1]/float(k*k) : Dens[k-1];  // density
        
        // Assign center to cylinder k
        float xtmp = x0 + (2.*k-1)*zpix*cos(PA*M_PI/180.)*sin(inc*M_PI/180.);
        float ytmp = y0 + (2.*k-1)*zpix*sin(PA*M_PI/180.)*sin(inc*M_PI/180.); 
        
        // Initializing rings
        Rings<T> *r = new Rings<T>;
        r->nr = par.DENSTYPE==2 ? 2 : k+1;
        r->radsep = (Rout-Rin)/(r->nr-1);
        for (int i=0; i<r->nr; i++) {
            r->radii.push_back(Rin+i*r->radsep);
            float F = par.DENSTYPE==2 ? 1 : float(i)/float(r->nr-1);
            r->vrad.push_back(-Vwind[k-1]*sin(F*OpA/2.*M_PI/180.));  
            r->vvert.push_back(Vwind[k-1]*cos(F*OpA/2.*M_PI/180.)); 
            r->xpos.push_back(xtmp); 
            r->ypos.push_back(ytmp);  
            r->vsys.push_back(Vsys);  
            r->vdisp.push_back(Vdisp[k-1]);
            r->inc.push_back(inc); 
            r->phi.push_back(PA);  
            r->z0.push_back(z0);   
            r->dens.push_back(dens*1E20); 
            r->vrot.push_back(Vrot[k-1]);  
            r->dvdz.push_back(0.);  
            r->zcyl.push_back(0.);  
        }

        // Build first cone
        con1->input(in, r, nv, par.LTYPE, 1, par.CDENS, -k);
        con1->calculate();
    
        // Build second cone
        xtmp = x0 - (2.*k-1)*zpix*cos(PA*M_PI/180.)*sin(inc*M_PI/180.);
        ytmp = y0 - (2.*k-1)*zpix*sin(PA*M_PI/180.)*sin(inc*M_PI/180.); 
        
        for (int i=0; i<r->nr; i++) {
            r->xpos[i]  = xtmp; 
            r->ypos[i]  = ytmp;
            r->vvert[i] = -r->vvert[i];
        }
        
        con2->input(in, r, nv, par.LTYPE, 1, par.CDENS, -k);
        con2->calculate();

        for (size_t i=0; i<in->NumPix(); i++) 
            out->Array(i) += (con1->Out()->Array(i)+con2->Out()->Array(i));

        delete r;
    }
    
    delete con1;
    delete con2;
}

    if (verb) bar.fillSpace(" Done.\n");
    
    in->pars().setVerbosity(verb);
    return true;
}


template <class T>
bool GalWind<T>::smooth(bool scalefac) {
        
    if (outDefined && par.SM) {
        // Smoothing
        float arcconv = arcsconv(in->Head().Cunit(0)); 
        Beam oldbeam = {0., 0., 0};
        Beam newbeam = {in->Head().Bmaj()*arcconv, in->Head().Bmin()*arcconv, in->Head().Bpa()};

        Smooth3D<T> *smoothed = new Smooth3D<T>;    
        smoothed->setUseScalefac(scalefac);
        smoothed->setUseBlanks(false);
        smoothed->smooth(out, oldbeam, newbeam, out->Array(), out->Array());

        for (size_t i=0; i<out->NumPix(); i++)
            if (out->Array(i)<1.E-12) out->Array()[i] = 0.;

        delete smoothed;
        return true;
    }
    else return false;
    
    
}


template <class T>
bool GalWind<T>::writeFITS(std::string fname, bool fullHead) {
    
    if (outDefined) {
        std::string fn = fname;
        if (fn=="") fn = in->pars().getOutfolder()+in->Head().Name()+"_wind.fits";
    
        return out->fitswrite_3d(fn.c_str(),fullHead);
    }
    else return false;

}


template <class T>
bool GalWind<T>::writeMomentMaps() {
    
    if (outDefined) {
        std::string fn = in->pars().getOutfolder()+in->Head().Name()+"_wind";
        MomentMap<T> *map = new MomentMap<T>;
        map->input(out);
        map->ZeroMoment(false);
        map->fitswrite_2d((fn+"_0mom.fits").c_str());
        map->FirstMoment(false);
        map->fitswrite_2d((fn+"_1mom.fits").c_str());
        map->SecondMoment(false);
        map->fitswrite_2d((fn+"_2mom.fits").c_str());
        delete map;
        return true;
    }
    else return false;
}


template <class T>
bool GalWind<T>::writePV(std::string fname){
    
    if (outDefined) {
        std::string fn = fname;
        if (fn=="") fn = in->pars().getOutfolder()+in->Head().Name()+"_wind_pv.fits";
        float ang = atof(par.PHI.c_str())-90;
        string pos[2] = {par.XPOS, par.YPOS};
        double *pixs = getCenterCoordinates(pos, in->Head());
        float x0  = pixs[0], y0=pixs[1];
        // Extract pv
        Image2D<T> *pv_max = PositionVelocity(out,x0,y0,ang);
        return pv_max->fitswrite_2d(fn.c_str());
    }
    else return false;
    
}


template <class T>
bool GalWind<T>::makePlots(){
    
    if (!outDefined) return false;
    
    std::string scriptname = "pyscript_w.py";
    std::ofstream pyf((in->pars().getOutfolder()+scriptname).c_str());

    float PA   = atof(par.PHI.c_str());  
    float Vsys = atof(par.VSYS.c_str());
    // Getting the center position 
    string pos[2] = {par.XPOS, par.YPOS};
    double *pixs = getCenterCoordinates(pos, in->Head());
    float x0  = pixs[0], y0=pixs[1];
    
    pyf << "import numpy as np \n"
        << "import os \n"
        << "import matplotlib as mpl \n"
        << "import matplotlib.pyplot as plt \n"
        << "from matplotlib.colorbar import ColorbarBase \n"
        << "from astropy.io import fits \n"
        << "from astropy.visualization import LinearStretch, PowerStretch \n"
        << "from astropy.visualization.mpl_normalize import ImageNormalize \n"
        << "from astropy.visualization import PercentileInterval \n"
        << "mpl.style.use('classic')\n"
        << "mpl.rc('font',family='Times New Roman') \n"
        << "params = {'text.usetex': False, 'mathtext.fontset': 'stix'}\n"
        << "plt.rcParams.update(params)\n\n";
    
    pyf << "# PARAMETERS: plotting the fit parameters \n"
        << "gname = '" << in->Head().Name() <<"' \n"
        << "outfolder = '" << in->pars().getOutfolder() <<"' \n\n";
            
    pyf << "nc = 7 \n"
        << "xlen, ylen = 0.15, 0.15 \n"
        << "xsep, ysep = 0.01, 0.01 \n"
        << "cmap = plt.get_cmap('Greys') \n"
        << "interval = PercentileInterval(99) \n";

    pyf << "# Plotting channel maps\n"
        << "fig = plt.figure(figsize=(10, 10))\n" 
        << "ax = []\n"
        << "for jj in range (nc):\n"
        << "\tax.append(fig.add_axes([jj*(xlen+xsep),1,xlen,ylen]))\n\n"
        << "fitc = fits.open(outfolder+gname+'_wind.fits')[0]\n"
        << "d, h = fitc.data, fitc.header\n"
        << "lims = interval.get_limits(d) \n"
        << "norm = ImageNormalize(vmin=lims[0], vmax=lims[1], stretch=PowerStretch(1))\n"
        << "cont = np.array([1,2,4,16,64,256,1024])*0.04 \n"
        << "crpix3_kms = " << in->Head().Crpix(2) << std::endl
        << "cdelt3_kms = " << DeltaVel<float>(in->Head()) << std::endl
        << "crval3_kms = " << AlltoVel(in->Head().Crval(2),in->Head()) << "\n"    
        << "velos = (np.arange(0,d.shape[0])+1-crpix3_kms)*cdelt3_kms+crval3_kms\n\n";
        
    pyf << "for k in range(nc): \n"
        << "\ta = ax[k] \n"
        << "\tchan = int(k*(d.shape[0])/nc) \n"           
        << "\ttopl = d[chan]\n"     
        << "\ta.tick_params(top=False,right=False,bottom=False,left=False,labelbottom=False,labelleft=False) \n"
        << "\tif np.all(topl<cont[0]): topl[:,:] = 0 \n"
        << "\ta.imshow(topl,origin='lower',cmap=cmap,norm=norm) \n"
        << "\ta.contour(topl,levels=cont,colors='#009999') \n"
        << "\ttext = r'$V_\\mathrm{LOS}$ = %.1f km/s'%(velos[chan]) \n"
        << "\ta.text(0.5,1.05,text,transform=a.transAxes,ha='center')\n"
        << "fig.savefig(outfolder+gname+'_wind.pdf',bbox_inches='tight')\n\n";
            
    pyf << "# Plotting kinematic maps\n"
        << "fig = plt.figure(figsize=(10, 10))\n" 
        << "ax, axcb = [], [] \n"
        << "for jj in range (4):\n"
        << "\tif jj==3: ax.append(fig.add_axes([jj*(xlen+xsep),1,1.5*xlen,ylen]))\n"
        << "\telse: ax.append(fig.add_axes([jj*(xlen+xsep),1,xlen,ylen]))\n"
        << "\tif jj<3:\n"
        << "\t\taxcb.append(fig.add_axes([jj*(xlen+xsep),1-0.019,xlen,0.015]))\n\n"
        << "pv   = fits.open(outfolder+gname+'_wind_pv.fits')[0].data\n"
        << "mom0 = fits.open(outfolder+gname+'_wind_0mom.fits')[0].data\n"
        << "mom1 = fits.open(outfolder+gname+'_wind_1mom.fits')[0].data\n"
        << "mom2 = fits.open(outfolder+gname+'_wind_2mom.fits')[0].data\n"
        << "mom0 /= np.nanmax(mom0)\n"
        << "topl = [mom0,mom1,mom2,pv]\n"
        << "cmaps = [plt.get_cmap('afmhot'),plt.get_cmap('Spectral_r'),plt.get_cmap('viridis'),plt.get_cmap('Greys')]\n"
        << "norms = [ImageNormalize(vmin=0,vmax=1),ImageNormalize(vmin=np.nanmin(mom1),vmax=np.nanmax(mom1)),ImageNormalize(vmin=np.nanmin(mom2),vmax=np.nanmax(mom2))]\n"
        << "titls = ['Intensity','Velocity','Dispersion','Position-Velocity']\n"
        << "labs  = ['Normalized $I$', r'$V_\\mathrm{LOS}$ (km/s)', r'$\\sigma$ (km/s)', 'Offset (arbitrary units)']\n"
        << "for k in range(len(topl)):\n"
        << "\ta = ax[k]\n"
        << "\ta.tick_params(top=False,right=False,bottom=False,left=False,labelbottom=False,labelleft=False)\n"
        << "\ta.text(0.5,1.1,titls[k],transform=a.transAxes,ha='center',va='center')\n"
        << "\tif k==1:\n"            
        << "\t\txcen, ycen, phi = "<< x0 << "," << y0 << "," << PA << std::endl
        << "\t\tif phi==90: a.axvline(xcen,c='k',ls='--',alpha=0.7) \n"
        << "\t\telse: \n"
        << "\t\t\txx = np.arange(d.shape[-1])\n"
        << "\t\t\tyy = np.tan(np.radians(phi))*(xx-xcen)+ycen \n"
        << "\t\t\ta.plot(xx,yy,c='k',ls='--',alpha=0.7)\n\n"
        
        << "\tif k==3: \n"
        << "\t\text = [-0.99,0.99 ,velos[0],velos[-1]]\n"
        << "\t\ta.set_ylim(ext[2]-2*abs(cdelt3_kms),ext[3]+2*abs(cdelt3_kms))\n"
        << "\t\ta.tick_params(top=True,left=True,bottom=True,right=True,labelright='on')\n"
        << "\t\ta.tick_params(axis='x', pad=7)\n"
        << "\t\ta.imshow(topl[k],origin='lower',cmap=cmaps[k],aspect='auto',extent=ext)\n"
        << "\t\ta.contour(topl[k],levels=cont,colors='#769BCB',extent=ext)\n"
        << "\t\ta.axhline(" << Vsys << ",c='k',ls='-',alpha=0.5)\n"
        << "\t\ta.axvline(0,c='k',ls='-',alpha=0.5)\n"
        << "\t\ta.tick_params(labelbottom=True) \n"
        << "\t\ta.text(1.25, 0.5,r'$V_\\mathrm{LOS}$ (km/s)',transform=a.transAxes,ha='center',va='center',rotation=90) \n"
        << "\telse:\n"
        << "\t\tm = a.imshow(topl[k],origin='lower',cmap=cmaps[k],aspect='auto',norm=norms[k]) \n"
        << "\t\tfor axis in ['top','bottom','left','right']: axcb[k].spines[axis].set_linewidth(0) \n"
        << "\t\tcbar = plt.colorbar(mappable=m,cax=axcb[k],orientation='horizontal',norm=norms[k])\n"
        << "\t\tcbar.solids.set_edgecolor('face')\n"
        << "\t\tcbar.outline.set_linewidth(0)\n"
        << "\t\tcbar.locator = mpl.ticker.MaxNLocator(nbins=5)\n"
        << "\t\tcbar.update_ticks()\n"
        << "\ta.text(0.5,-0.35,labs[k],transform=a.transAxes,ha='center')\n"
        << "fig.savefig(outfolder+gname+'_wind_maps.pdf',bbox_inches='tight')\n\n";
    
    pyf.close();

#ifdef HAVE_PYTHON
    // Now plotting everything
    if (in->pars().isVerbose()) std::cout << " Producing " << randomAdjective(1) << " plots..." << std::flush;
    std::string cmd = "python "+in->pars().getOutfolder()+scriptname+" > /dev/null 2>&1";
    int ret = system(cmd.c_str());
    if (in->pars().isVerbose()) {
        if (ret==0) std::cout << " Done." << std::endl;
        else std::cout << " Something went wrong! Check pyscript.py in the output folder." << std::endl;
    }
    
#endif

    return true;
}


// Explicit instantiation of the class
template class GalWind<float>;
template class GalWind<double>;
