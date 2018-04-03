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
    par.OPENANG = openang;
    par.HTOT    = htot;
    par.DENSTYPE= denstype;
    par.NTOT    = ntot;
    par.CDENS   = cdens;
    par.NV      = nv;
    in->pars().setThreads(nthreads);
    
    in->checkBeam();
}
template GalWind<float>::GalWind(Cube<float>*,float,float,float,float,float,float,float,float,float,float,int,int,int,int,int);
template GalWind<double>::GalWind(Cube<double>*,double,double,double,double,double,double,double,double,double,double,int,int,int,int,int);


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
template GalWind<float>& GalWind<float>::operator=(const GalWind<float>&);
template GalWind<double>& GalWind<double>::operator=(const GalWind<double>&);


template <class T>
bool GalWind<T>::compute() {
    
    // Collecting needed parameters
    float arcconv = arcsconv(in->Head().Cunit(0));
    float PixSize = in->Head().PixScale()*arcconv;
    int nv = par.NV>1 ? par.NV : 20;
    float z0 = par.HTOT/(2*par.NTOT);  //half-height of each cylinder in arcsec
    float zpix = z0/PixSize;           //half-height of each cylinder in M_PIxels
    float width = (2*z0)*tan(par.OPENANG/2.*M_PI/180.); //width of each ring (Rout/k)
        
    // Fixed parameter
    float PA   = atof(par.PHI.c_str());  
    float inc  = atof(par.INC.c_str());   
    float Vsys = atof(par.VSYS.c_str());  
    // Getting the center position 
    string pos[2] = {par.XPOS, par.YPOS};
    double *pixs = getCenterCoordinates(pos, in->Head());
    float x0  = pixs[0], y0=pixs[1];
    
    // Parameters that can vary cylinder by cylinder
    std::vector<T> Vwind, Dens, Vdisp;
    bool vwind_b = getDataColumn(Vwind,par.VWIND);
    bool vdisp_b = getDataColumn(Vdisp,par.VDISP);
    bool dens_b  = getDataColumn(Dens,par.DENS);
    for (size_t i=0; i<par.NTOT; i++) {
        // If only 1 value is given, fill vectors with this value
        if (!vwind_b) Vwind.push_back(atof(par.VWIND.c_str()));
        if (!vdisp_b) Vdisp.push_back(atof(par.VDISP.c_str()));
        if (!dens_b)  Dens.push_back(atof(par.DENS.c_str()));
    }
    // Checking array dimensions. They must be equal to NTOT
    if (Vwind.size()!=par.NTOT || Vdisp.size()!=par.NTOT || Dens.size()!=par.NTOT) {
        std::cerr << " GALWIND ERROR: VWIND, VDISP and DENS must have size=NTOT!\n";
        return false;
    }
     
    // Allocate output datacube
    if (outDefined) delete out;
    out = new Cube<T>(in->AxisDim());
    out->saveHead(in->Head());
    outDefined = true;
    for (size_t i=0; i<out->NumPix(); i++) out->Array(i) = 0;

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
    
#pragma omp for
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
            r->vrad.push_back(-Vwind[k-1]*sin(F*par.OPENANG/2.*M_PI/180.));  
            r->vvert.push_back(Vwind[k-1]*cos(F*par.OPENANG/2.*M_PI/180.)); 
            r->xpos.push_back(xtmp); 
            r->ypos.push_back(ytmp);  
            r->vsys.push_back(Vsys);  
            r->vdisp.push_back(Vdisp[k-1]);
            r->inc.push_back(inc); 
            r->phi.push_back(PA);  
            r->z0.push_back(z0);   
            r->dens.push_back(dens*1E20); 
            r->vrot.push_back(0.);  
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
template bool GalWind<float>::compute();
template bool GalWind<double>::compute();


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
template bool GalWind<float>::smooth(bool);
template bool GalWind<double>::smooth(bool);


template <class T>
bool GalWind<T>::writeFITS(std::string fname, bool fullHead) {
    
    if (outDefined) {
        std::string fn = fname;
        if (fn=="") fn = in->pars().getOutfolder()+in->Head().Name()+"_wind.fits";
    
        return out->fitswrite_3d(fn.c_str(),fullHead);
    }
    else return false;

}
template bool GalWind<float>::writeFITS(std::string, bool);
template bool GalWind<double>::writeFITS(std::string, bool);


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
template bool GalWind<float>::writeMomentMaps();
template bool GalWind<double>::writeMomentMaps();

