//----------------------------------------------------------
// ringmodel.cpp: Member functions for Ringmodel class
//----------------------------------------------------------

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
#include <cfloat>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <fstream>
#include <Tasks/ringmodel.hh>
#include <Arrays/cube.hh>
#include <Arrays/rings.hh>
#include <Arrays/param.hh>
#include <Tasks/moment.hh>
#include <Utilities/lsqfit.hh>
#include <Utilities/utils.hh>
#include <Utilities/progressbar.hh>


#define MAXPAR  7
#define nint(x)     (x>0.0?(int)(x+0.5):(int)(x-0.5))


void Ringmodel::defaults() {
    
    tol = 0.001;
    cor[0] = -1;
    cor[1] = -1;    
    
    mask = new bool[MAXPAR];
    for (int i=0; i<MAXPAR; i++) mask[i] = true;
    
    side = 3;
    wpow = 2;
    thetaf = 15;
    
    allAllocated = false;
}


Ringmodel::Ringmodel() {
    
    defaults();
    
}


Ringmodel::Ringmodel(int nrings) {
    
    defaults(); 
    
    nrad  = nrings;
    rads  = new float [nrings]; 
    wids  = new float [nrings];                 
    vsysi = new float [nrings];             
    vsysf = new float [nrings];             
    vsyse = new float [nrings];                 
    vroti = new float [nrings];             
    vrotf = new float [nrings];         
    vrote = new float [nrings];             
    vexpi = new float [nrings];             
    vexpf = new float [nrings];             
    vexpe = new float [nrings];                 
    posai = new float [nrings];         
    posaf = new float [nrings];     
    posae = new float [nrings];     
    incli = new float [nrings];     
    inclf = new float [nrings];     
    incle = new float [nrings];         
    xposf = new float [nrings];         
    xpose = new float [nrings];                 
    yposf = new float [nrings];     
    ypose = new float [nrings];         

    npts  = new int [nrings];           
    chis  = new float [nrings];     
    
    elp = allocate_2D<float>(nrings, 4);    
    
    allAllocated = true;
    
}


Ringmodel::Ringmodel (Cube<float> *c)  {

    in = c;
    GALFIT_PAR p = in->pars().getParGF();
    
    // Reading ring inputs
    Rings<float> *r = readRings<float>(p,c->Head());
    int nrings = r->nr-1;
    float *widths = new float[nrings];
    float *radii = new float[nrings];
    for (int i=0; i<nrings; i++) {
        radii[i]  = r->radii[i]/(arcsconv(c->Head().Cunit(0))*c->Head().PixScale());
        widths[i] = (r->radii[i+1]-r->radii[i])/(arcsconv(c->Head().Cunit(0))*c->Head().PixScale());
    }

    set(r->nr-1,radii,widths,&r->vsys[0],&r->vrot[0],&r->vrad[0],
        &r->phi[0],&r->inc[0],r->xpos[0],r->ypos[0]);
    

    bool mpar[MAXPAR];
    for (int i=0; i<MAXPAR; i++) mpar[i]=false;
    std::string FREE = p.FREE;
    FREE = makelower(FREE);

    int found = FREE.find("vsys");
    if (found<0) mpar[VSYS]=false;
    else mpar[VSYS]=true;

    found = FREE.find("vrot");
    if (found<0) mpar[VROT]=false;
    else mpar[VROT]=true;

    found = FREE.find("pa");
    if (found<0) mpar[PA]=false;
    else mpar[PA]=true;

    found = FREE.find("inc");
    if (found<0) mpar[INC]=false;
    else mpar[INC]=true;

    found = FREE.find("xpos");
    if (found<0) mpar[X0]=false;
    else mpar[X0]=true;

    found = FREE.find("ypos");
    if (found<0) mpar[Y0]=false;
    else mpar[Y0]=true;

    found = FREE.find("all");
    if (found>=0)
        for (int i=0; i<MAXPAR; i++) mpar[i] = true;

    int hside;
    if (p.SIDE=="R") hside=1;
    else if (p.SIDE=="A") hside=2;
    else hside = 3;

    int wfunc = p.WFUNC;

    setoption (mpar,hside,wfunc,15.);

    // Calculating 1st moment map
    MomentMap<float> map;
    map.input(c);
    if (c->Head().NumAx()>2) map.FirstMoment(true);
    else {
        for (int i=0; i<c->NumPix(); i++) map.Array(i) = c->Array(i);
        map.setHead(1);
    }
    map.fitswrite_2d((c->pars().getOutfolder()+c->Head().Name()+"map_1st.fits").c_str());

    int boxlow[2] = {0,0};
    int boxup[2] = {map.DimX()-1, map.DimY()-1};
    setfield(map.Array(),map.DimX(),map.DimY(),boxup,boxlow);
}


Ringmodel::Ringmodel (int nrings, float *radii, float *widths, float *vsys, float *vrot, 
                      float *vexp, float *posang, float *incl, float xcenter, float ycenter)  {
    
    set(nrings, radii, widths, vsys, vrot, vexp, posang, incl, xcenter, ycenter);
    
}


Ringmodel::Ringmodel (int nrings, float *radii, float widths, float vsys, float vrot, 
                      float vexp, float posang, float incl, float xcenter, float ycenter)  {
    
    set(nrings, radii, widths, vsys, vrot, vexp, posang, incl, xcenter, ycenter);
    
}


void Ringmodel::set (int nrings, float *radii, float *widths, float *vsys, float *vrot, 
                      float *vexp, float *posang, float *incl, float xcenter, float ycenter) {
                          
    defaults();
    
    nrad  = nrings;
    rads  = radii;
    wids  = widths;
    vsysi = vsys;
    vroti = vrot; 
    vexpi = vexp;
    posai = posang;
    incli = incl;
    xposi = xcenter;
    yposi = ycenter;
    vsysf = new float [nrings];             
    vsyse = new float [nrings];                             
    vrotf = new float [nrings];         
    vrote = new float [nrings];                         
    vexpf = new float [nrings];             
    vexpe = new float [nrings];                         
    posaf = new float [nrings];     
    posae = new float [nrings];         
    inclf = new float [nrings];     
    incle = new float [nrings];         
    xposf = new float [nrings];         
    xpose = new float [nrings];                 
    yposf = new float [nrings];     
    ypose = new float [nrings];         

    npts  = new int [nrings];           
    chis  = new float [nrings];     
    
    elp = allocate_2D<float>(nrings, 4);        
    
    allAllocated = true;
    
}


void Ringmodel::set (int nrings, float *radii, float widths, float vsys, float vrot, 
                      float vexp, float posang, float incl, float xcenter, float ycenter) {

    defaults();
    
    nrad  = nrings;
    xposi = xcenter;
    yposi = ycenter;
    
    rads  = radii;
    wids  = new float [nrings];                 
    vsysi = new float [nrings];                             
    vroti = new float [nrings];                             
    vexpi = new float [nrings];                             
    posai = new float [nrings];             
    incli = new float [nrings]; 
    
    for (int i=0; i<nrad; i++) {
        
        wids  [i] = widths;
        vsysi [i] = vsys;
        vroti [i] = vrot; 
        vexpi [i] = vexp;
        posai [i] = posang;
        incli [i] = incl;

    }   

    vsysf = new float [nrings];             
    vsyse = new float [nrings];                             
    vrotf = new float [nrings];         
    vrote = new float [nrings];                         
    vexpf = new float [nrings];             
    vexpe = new float [nrings];                         
    posaf = new float [nrings];     
    posae = new float [nrings];         
    inclf = new float [nrings];     
    incle = new float [nrings];         
    xposf = new float [nrings];         
    xpose = new float [nrings];                 
    yposf = new float [nrings];     
    ypose = new float [nrings];         

    npts  = new int [nrings];           
    chis  = new float [nrings];     
    
    elp = allocate_2D<float>(nrings, 4);        
    
    allAllocated = true;

}


void Ringmodel::setoption (bool *maskpar, int hside, int wfunc, float freeangle) {
    
  /// This function sets some options for fitting.
  ///
  /// \param maskpar    Which parameter do you want to fix.
  /// \param hside      Which half of velocity field:
  ///                   1 = Receding half.
  ///                   2 = Approching half.
  ///                   3 = Both halves.
  /// \param wfunc      Which weighting function:
  ///                   1 = Uniform.
  ///                   2 = Cosine.
  ///                   3 = Cosine-squared.
  /// \param freeangle  The angle around mionr axis within 
  ///                   witch radial velocities are discarded.
    
    side = hside;   
    if ((side!=1) && (side!=2) && (side!=3)) {
        std::cout << "Not allowed half of galaxy. Setting to both.\n";
        side = 3;
    }
    
    thetaf = freeangle;
    
    switch(wfunc) {
        case 1: 
            wpow = 0;
            break;       
        case 2: 
            wpow = 1;
            break;
        case 3: 
            wpow = 2;
            break;
        default: 
            std::cout << "Not allowed weighting function. Setting to COSINE.\n";
            wpow = 1;
            break;
         
    }

    // We allow only fitting of systemic velocity and centre position
    // when both halves of the galaxy are used.
    
    for (int i=0; i<MAXPAR; i++) mask[i] = maskpar[i];
    
    if (side!=3) {
        mask[VSYS] = mask[X0] = mask[Y0] = 0;
    }

    int nfixed = 0;
    for (int i=0; i<MAXPAR; i++) nfixed += (1-mask[i]);
    
    if (nfixed == MAXPAR) std::cout << "NO free parameters!\n";
    

}
    

void Ringmodel::setfield (float *Array, int xsize, int ysize, int *boxup, int *boxlow) {
    
    int npoints;
    
    for (int i=0; i<2; i++) {
        bup[i] = boxup[i];
        blo[i] = boxlow[i];
    }
    
    npoints = (bup[1]-blo[1]+1)*(bup[0]-blo[0]+1);
    
    vfield = new float [npoints];
    fieldAllocated = true;
    
    for (int i=blo[0]; i<=bup[0]; i++) {
        for (int j=blo[1]; j<=bup[1]; j++) {
            float v = Array[i+j*xsize];
            int velpix = (i-blo[0]) + (j-blo[1])*(bup[0]-blo[0]+1);
            vfield[velpix]= v;
        }
        
    }

}


void Ringmodel::ringfit() {
    
  /// This function makes a loop over all concentric rings 
  /// and do the fit.
  /// It is the calling function for Ringmodel class.
   

    if (fieldAllocated && allAllocated) {
        
        ProgressBar bar(" Fitting 2D tilted-ring model... ", true);
        bar.setShowbar(in->pars().getShowbar());        
        bool verb = in->pars().isVerbose();
        int nthreads = in->pars().getThreads();

#pragma omp parallel num_threads(nthreads)
{
        if (verb) bar.init(nrad);
#pragma omp for
        for (int ir = 0; ir<nrad; ir++) {        
            if (verb) bar.update(ir+1);
    
            int n;
            float   e[MAXPAR];
            float   p[MAXPAR];
            float   q = 0.0;

            p[VSYS] = vsysi[ir];
            p[VROT] = vroti[ir];
            p[VEXP] = vexpi[ir];
            p[PA]   = posai[ir];
            p[INC]  = incli[ir];
            p[X0]   = xposi;
            p[Y0]   = yposi;
            float ri = rads[ir] - 0.5 * wids[ir];
            float ro = rads[ir] + 0.5 * wids[ir];
            if ( ri < 0.0 ) ri = 0.0;
        
            if (rotfit(ri, ro, p, e, n, q)>0) {
            
                vsysf[ir] = p[0];
                if (e[0] < 999.99) vsyse[ir] = e[0];
                else vsyse[ir] = 999.99;
            
                vrotf[ir] = p[1];
                if (e[1] < 999.99) vrote[ir] = e[1];
                else vrote[ir] = 999.99;
            
                vexpf[ir] = p[2];
                if (e[2] < 999.99) vexpe[ir] = e[2];
                else vexpe[ir] = 999.99;
            
                posaf[ir] = p[3];
                if (e[3] < 999.99) posae[ir] = e[3];
                else posae[ir] = 999.99;
            
                inclf[ir] = p[4];
                if (e[4] < 999.99) incle[ir] = e[4];
                else incle[ir] = 999.99;
            
                xposf[ir] = p[5];
                if (e[5] < 999.99) xpose[ir] = e[5];
                else xpose[ir] = 999.99;
            
                yposf[ir] = p[6];
                if (e[6] < 999.99) ypose[ir] = e[6];
                else ypose[ir] = 999.99;
            
                elp[ir][0] = elp4[0];
                elp[ir][1] = elp4[1];
                elp[ir][2] = elp4[2];
                elp[ir][3] = elp4[3];
                npts[ir] = n;
                chis[ir] = q;
            } 
            else {
                // Fit did not succeed
                vsysf[ir]=vrotf[ir]=vexpf[ir]=posaf[ir]=inclf[ir]=xposf[ir]=yposf[ir]=log(-1);
                vsyse[ir]=vrote[ir]=vexpe[ir]=posae[ir]=incle[ir]=xpose[ir]=ypose[ir]=log(-1);
            }
        }
}    
        if (verb) bar.fillSpace("Done.\n");
    
    }
    else {
        std::cout << "2DFIT ERROR: Arrays are not allocated!\n";
    }
    
    
}


int Ringmodel::rotfit (float ri, float ro, float *p, float *e, int &n, float &q) {
    
  /// This function does a least squares fit to the radial velocity field.
  ///
  /// \param  ri        Inner radius of ring.
  /// \param  ro        Outer radius of ring.
  /// \param  p         Fitted parameters of ring.
  /// \param  e         Errors in parameters.
  /// \param  n         Number of points in the fit.
  /// \param  q         Chi-squared.
  ///
  /// \return           Error or success.    

    int     i, h;
    int     ier = 0;            // Error return. 
    int     nfr;                // Number of free parameters.
    int     nrt;                // Return code from lsqfit.
    int     stop;               // Stop ? 
    int     t = 100;            // Max. number of iterations.
    float   b[MAXPAR];          // Partial derivatives.
    float   chi;                // Old chi-squared.
    float   df[MAXPAR];         // Difference vector.
    float   eps[MAXPAR];        // Stop criterium.
    float   flip;               // Direction.
    float   lab = 0.001;        // Mixing parameter.        
    
    std::vector<float> x,y,w;           // (x,y) position, f(x,y) and w(x,y).
    float *pf  = new float [MAXPAR];    // Intermediate results.

    for (nfr=0, i=0; i<MAXPAR; i++) {
        eps[i] = 0.1;                   // Convergence criterium.
        nfr += mask[i];
    }
   
    h = 0;                              
    n = getdat(x, y, w, p, ri, ro, q, nfr);

    stop = 0;
    do {                    
        int npar = MAXPAR;                              // Number of parameters.
        int xdim = 2;                                   // Function is two-dimensional.
        h += 1;             
        chi = q;                                        // Save chi-squared.
        for (i=0; i<MAXPAR; i++) {                      // Loop to save initial estimates.
            pf[i] = p[i];
        }
        
        Lsqfit<float> lsqfit(&x[0], xdim, &y[0], &w[0], n, pf, e, mask, npar, &func, &derv, tol, t, lab);
        nrt = lsqfit.fit();

        if (nrt<0) break;                               // Stop because of error.
        for (i=0; i<MAXPAR; i++)                        // Calculate difference vector.
            df[i] = pf[i] - p[i];
        
        flip = 1.0;                                     // Factor for inner loop.
        while (1) {                                     // Inner loop. 
            for (i=0; i<MAXPAR; i++)                    // Calculate new parameters.
                pf[i] = flip * df[i] + p[i];
         
            if (pf[INC] > 90.0)  pf[INC] -= 180.0;      // In case inclination > 90.
            
            n = getdat(x, y, w, pf, ri, ro, q, nfr);
         
            if (q<chi) {                                // Better fit.
                for (i=0; i<MAXPAR; i++)                // Save new parameters.
                    p[i] = pf[i];
                break;  
            } 
            else {
                if ((2*h) > t) {
                    for (stop=1, i=0; i<MAXPAR; i++) {
                        stop = (stop && (fabs(flip * df[i]) < eps[i]));
                    }
                } 
                else {
                    if (q == chi && chi == 0.0)
                        stop = 1;
                    else
                        stop = ((fabs(q-chi)/chi) < tol);
                }
                if (stop) {
                    q = chi;
                    break;
                }
            }
            if (flip > 0.0) flip *= -1.0;
            else flip *= -0.5;
         
        }
    } while (!stop && h<t);
   
   
    // Find out why we quit fitting and printing errors
    
    if (stop)                   // Good fit:  ier = number of big loops.
        ier = h;                
    else if (nrt < 0)           // Error from lsqfit: ier = return code of lsqfit.
        ier = nrt;
    else if ( h == t )          // Maximum number of iterations.
        ier = -4;
   
    switch (ier) {
        
        case -1:
            std::cout << "Error fitting ring model: Too many free parameters!\n";
            break;
      
        case -2: 
            std::cout << "Error fitting ring model: No free parameters!\n";
            break;
      
        case -3: 
            std::cout << "Error fitting ring model: Not enough degrees of freedom!\n";
            break;
      
        case -4: 
            std::cout << "Error fitting ring model: Maximum number of iterations too small!\n";
            break;
        
        case -5: 
            std::cout << "Error fitting ring model: Diagonal of matrix contains zeroes!\n";
            break;
        
        case -6: 
            std::cout << "Error fitting ring model: Deter. of the coeff. matrix is zero!\n";
            break;
        
        case -7: 
            std::cout << "Error fitting ring model: Square root of negative number!\n";
            break;
      
        default: 
            break;
      
    }
        
    // Calculate ellipse parameters.
   
    if (ier==1 && cor[0]>-1 && cor[1]>-1 ) {
        int   maxpar = MAXPAR;
        float a11 = 0.0;
        float a12 = 0.0;
        float a22 = 0.0;
        float sigma2 = 0.0;

        for (i=0; i < n; i++) {
            derv( &x[2*i], p, b, maxpar);
            a11 = a11 + w[i] * b[cor[0]] * b[cor[0]];
            a22 = a22 + w[i] * b[cor[1]] * b[cor[1]];
            a12 = a12 + w[i] * b[cor[0]] * b[cor[1]];
            sigma2 = sigma2+w[i]*std::pow(double(y[i]-func(&x[2*i], p, maxpar)), double(2.0));
        }
        sigma2 = sigma2 / (float) (n);
        elp4[0] = a11;
        elp4[1] = a12;
        elp4[2] = a22;
        elp4[3] = sigma2;
    }
    
    delete [] pf;
    
    return ier;
}


int Ringmodel::getdat (std::vector<float> &x, std::vector<float> &y, std::vector<float> &w, 
                       float *p, float ri, float ro, float &q, int nfr) {
    
  /// The function selects the data from velocity field
  /// and calculates differences.
  ///
  /// \param  x         Sky coordinates of pixel inside ring.
  /// \param  y         Radial velocities.
  /// \param  w         Weights of radial velocities.
  /// \param  p         Parameters of ring.
  /// \param  ri        Inner radius of ring.
  /// \param  ro        Outer radius of ring.
  /// \param  q         Chi-squared.
  /// \param  nfr       Degrees of freedom.
  ///
  /// \return           Number of points in ring.
  
    int     n=0;                        // Return value = number of points.
    const double F = M_PI/180.;
    
    // Reset variables
    x.clear();
    y.clear();
    w.clear();
    q    = 0.0;
    // Definition of parameters.
    float phi  = p[PA];       
    float inc  = p[INC];              
    float x0   = p[X0];           
    float y0   = p[Y0];           
    float free = fabs(sin(F*thetaf)); // Free angle.
    float sinp = sin(F*phi);      
    float cosp = cos(F*phi);  
    float sini = sin(F*inc);      
    float cosi = cos(F*inc);     
    float a    = sqrt(1.0-cosp*cosp*sini*sini);
    float b    = sqrt(1.0-sinp*sinp*sini*sini); 
    int llo    = std::max(blo[0], nint(x0-a*ro)-1);
    int lup    = std::min(bup[0], nint(x0+a*ro)+1);
    int mlo    = std::max(blo[1], nint(y0-b*ro)-1);
    int mup    = std::min(bup[1], nint(y0+b*ro)+1);

    if ((llo > lup) || (mlo > mup)) {
        std::cout << "Ring is outside the map!"<<std::endl;
        q = FLT_MAX;
        return 0;
    }
   
    int nlt = bup[0]-blo[0]+1;                     // Number of pixels in X.
   
    for (int rx=llo; rx<lup; rx++) {
        for (int ry=mlo; ry<mup; ry++) {
            int   ip = (ry-blo[1])*nlt+ (rx-blo[0]);
            float v  = vfield[ip];                 // Radial velocity at this position.
            
            if (v==v) {
                float xr = (-(rx-x0)*sinp + (ry-y0)*cosp);
                float yr = (-(rx-x0)*cosp - (ry-y0)*sinp)/cosi;
                float r = sqrt(xr*xr+yr*yr);
                float theta = 0.;
                if (r>=0.1) theta = atan2(yr, xr)/F;   
                float costh = fabs(cos(F*theta));
                    
                if (r>ri && r<ro && costh>free) {      // If we are inside the ring.

                    bool use = true;
                    if (side==1) use = (fabs(theta)<=90.0);      //< Receding half.
                    if (side==2) use = (fabs(theta)>=90.0);      //< Approaching half. 

                    if (use) {                          // Load data point ? 
                        n += 1;
                        float xx[2] = {float(rx),float(ry)};
                        float vz = func (xx, p, MAXPAR);
                        float s = v - vz;          // Corrected difference
                        float wi = std::pow(costh, wpow); // Weight of this point.
                        x.push_back(rx);           // Load X-coordinate.
                        x.push_back(ry);           // Load Y-coordinate.
                        y.push_back(v);            // Load LOS velocity.
                        w.push_back(wi);           // Load weight.
                        q += s*s*wi;               // Calculate chi-squared.
                        
                    }
                }
            }
        }   
    }
    
    if (n > nfr) q = sqrt(q/(float)(n-nfr));     // Enough data points ? 
    else q = FLT_MAX;       
    
    return n;
    
}


void Ringmodel::print (std::ostream& Stream) {
        
    if (allAllocated && fieldAllocated) {
        
        using namespace std;
        
        int m = 11;
        Stream  << showpoint << fixed;
        Stream  << endl << setfill('-');
        
        Stream << showpoint << fixed << setfill('-') << endl;
        Stream << setw(70) << " Initial values for fitting parameters " 
               << setw(34) << "  " << endl << endl;
                    
        Stream << setfill(' ');
    
        Stream << setw(m) << right << "Nrads"    << setw(m) << right << "Rmax " 
               << setw(m) << right << "Vsys "    << setw(m) << right << "Vrot " 
               << setw(m) << right << "Vexp "    << setw(m) << right << "P.A."
               << setw(m) << right << "Incl."    << setw(m) << right << "Xcenter"
               << setw(m) << right << "Ycenter"<< endl;
    
        Stream << setw(m) << right  << "[#] "  << setw(m) << right << "[pix]" 
               << setw(m) << right << "[KM/S]" << setw(m) << right << "[KM/S]"
               << setw(m) << right << "[KM/S]" << setw(m) << right << "[deg]"
               << setw(m) << right << "[deg]"    << setw(m) << right << "[pix] "
               << setw(m) << right << "[pix] " << endl << endl << endl;
        
        Stream << setfill(' ');
            
        Stream << setw(m) << right << setprecision(2) << nrad 
               << setw(m) << right << setprecision(2) << rads[nrad-1] 
               << setw(m) << right << setprecision(2) << vsysi[0]
               << setw(m) << right << setprecision(2) << vroti[0] 
               << setw(m) << right << setprecision(2) << vexpi[0]
               << setw(m) << right << setprecision(2) << posai[0]
               << setw(m) << right << setprecision(2) << incli[0]
               << setw(m) << right << setprecision(2) << xposi
               << setw(m) << right << setprecision(2) << yposi;
        
        Stream << endl << endl << endl;                 
        Stream << endl << endl << endl << endl << endl << setfill('-');
    
        
        
        Stream  << "  Results from fitting the tilted rings model:\n\n\n";
        
        for (int i=0; i<nrad; i++) {
            
            Stream  << setw(64) << " Fitted parameters for ring # " << i+1
                    << " " << setw(37) << " " << endl << endl;
                    
            Stream  << setfill(' ');
    
            Stream  << setw(m) << right << "Radius" 
                    << setw(m) << right << "Width" 
                    << setw(m) << right << "Vsys "
                    << setw(m) << right << "Vrot "  
                    << setw(m) << right << "Vexp " 
                    << setw(m) << right << "P.A."
                    << setw(m) << right << "Incl."
                    << setw(m) << right << "Xcenter"
                    << setw(m) << right << "Ycenter";
    
            Stream  << endl;
    
            Stream  << setw(m) << right  << "[pix] " 
                    << setw(m) << right << "[pix]" 
                    << setw(m) << right << "[KM/S]"
                    << setw(m) << right << "[KM/S]"
                    << setw(m) << right << "[KM/S]" 
                    << setw(m) << right << "[deg]"
                    << setw(m) << right << "[deg]"
                    << setw(m) << right << "[pix] "
                    << setw(m) << right << "[pix] ";
    
            Stream  << endl;
    
            Stream  << endl << endl;
        
            Stream  << setfill(' ');
        
                
            Stream  << setw(m) << right << setprecision(2) << rads [i] 
                    << setw(m) << right << setprecision(2) << wids [i] 
                    << setw(m) << right << setprecision(2) << vsysf[i]
                    << setw(m) << right << setprecision(2) << vrotf[i] 
                    << setw(m) << right << setprecision(2) << vexpf[i]
                    << setw(m) << right << setprecision(2) << posaf[i]
                    << setw(m) << right << setprecision(2) << inclf[i]
                    << setw(m) << right << setprecision(2) << xposf[i]
                    << setw(m) << right << setprecision(2) << yposf[i];
                        
        
            Stream  << endl << endl;
            
            
            Stream  << "\n  Parameters hold fixed: ";
        
            for (int j=0; j<MAXPAR; j++) 
                if (!mask[j]) {
                    
                    if (j==VSYS) Stream << "  Vsys = " << vsysi[i];
                    if (j==VROT) Stream << "  Vrot = " << vroti[i];
                    if (j==VEXP) Stream << "  Vexp = " << vexpi[i];
                    if (j==PA) Stream << "  P.A. = " << posai[i];
                    if (j==INC) Stream << "  Inclination = " << incli[i];
                    if (j==X0) Stream << "  Xcenter = " << xposi;
                    if (j==Y0) Stream << "  Ycenter = " << yposi;
            }
            
            Stream <<"\n  Reduced chi-squared: " << chis[i];
            
            Stream  << endl << endl << endl << endl << setfill('-');
        }
    }
}

void Ringmodel::printfinal (std::ostream& Stream) {

    int m=10;
    Stream  << fixed << setprecision(2);
    Stream  << "#" << setw(m) << right << "Radius"
            << setw(m) << right << "Radius"
            << setw(m) << right << "Vsys "
            << setw(m) << right << "Vrot "
            << setw(m) << right << "Vexp "
            << setw(m) << right << "P.A."
            << setw(m) << right << "Incl."
            << setw(m) << right << "Xcenter"
            << setw(m) << right << "Ycenter" << endl;

    Stream  << "#" << setw(m) << right  << "[pix] "
            << setw(m) << right << "[arcs]"
            << setw(m) << right << "[KM/S]"
            << setw(m) << right << "[KM/S]"
            << setw(m) << right << "[KM/S]"
            << setw(m) << right << "[deg]"
            << setw(m) << right << "[deg]"
            << setw(m) << right << "[pix] "
            << setw(m) << right << "[pix] " << endl;

    for (int i=0; i<nrad; i++) {

        Stream  << setw(m) << right << setprecision(2) << rads[i]
                << setw(m) << right << setprecision(2) << rads[i]*in->Head().PixScale()*arcsconv(in->Head().Cunit(0))
                << setw(m) << right << setprecision(2) << vsysf[i]
                << setw(m) << right << setprecision(2) << vrotf[i]
                << setw(m) << right << setprecision(2) << vexpf[i]
                << setw(m) << right << setprecision(2) << posaf[i]
                << setw(m) << right << setprecision(2) << inclf[i]
                << setw(m) << right << setprecision(2) << xposf[i]
                << setw(m) << right << setprecision(2) << yposf[i] << std::endl;

    }
}


void Ringmodel::writeModel (std::string fname) {
    
    int dim[2] = {in->DimX(),in->DimY()};
    Image2D<float> model(dim);
    model.saveHead(in->Head());
    for (int i=model.NumPix(); i--;) model[i] = log(-1);

    const double F = M_PI/180.;
    float p[MAXPAR], e[MAXPAR];
    float ot = thetaf;
    thetaf = 0;
    ///*
    for (int ir=0; ir<nrad; ir++) {
        p[X0]   = xposf[ir];
        p[Y0]   = yposf[ir];
        p[PA]   = posaf[ir];
        p[INC]  = inclf[ir];
        p[VEXP] = vexpf[ir];
        p[VROT] = vrotf[ir];
        p[VSYS] = vsysf[ir];
                
        float ri = rads[ir] - 0.5 * wids[ir];
        float ro = rads[ir] + 0.5 * wids[ir];
        if (ri<0.0) ri = 0.0;
        float q = 0.0;
        std::vector<float> x,y,w;
        int n = getdat(x,y,w,p,ri,ro,q,0);
        
        for (int j=1; j<=n; j++) {
            float xx[2] = {x[2*j-2],x[2*j-1]};
            if (xx[0]<model.DimX() && xx[1]<model.DimY()) {
                float vv = func(xx,p,MAXPAR);
                size_t pp = xx[0]+xx[1]*model.DimX();
                if (model[pp]!=model[pp]) model[pp] = vv;
                else model[pp] = (model[pp]+vv)/2.;
            }
        }
    }
    //*/
    // Making a last loop for filling holes in the model field;
    p[X0]   = xposf[nrad-1];
    p[Y0]   = yposf[nrad-1];
    p[PA]   = posaf[nrad-1];
    p[INC]  = inclf[nrad-1];
    p[VEXP] = vexpf[nrad-1];
    p[VROT] = vrotf[nrad-1];
    p[VSYS] = vsysf[nrad-1];
    float rl = rads[nrad-1]+0.5*wids[nrad-1];
    float q = 0;
    
    std::vector<float> x,y,w;
    int n = getdat(x,y,w,p,0,rl,q,0);
    
    for (int j=1; j<=n; j++) {
        float xx[2] = {x[2*j-2],x[2*j-1]};
        if (xx[0]<model.DimX() && xx[1]<model.DimY()) {
            size_t pp = xx[0]+xx[1]*model.DimX();
            if (model[pp]!=model[pp]) {
                float xr = (-(xx[0]-p[X0])*sin(F*p[PA]) + (xx[1]-p[Y0])*cos(F*p[PA]));
                float yr = (-(xx[0]-p[X0])*cos(F*p[PA]) - (xx[1]-p[Y0])*sin(F*p[PA]))/cos(F*p[INC]);
                float R  = sqrt(xr*xr+yr*yr);
                // Finding closest radius 
                float bdif = 1E18;
                int ind = 0;
                for (int j=0; j<nrad; j++) {
                    float diff = fabs(R-rads[j]);
                    if (diff<bdif) {
                        bdif = diff;
                        ind  = j;
                    }
                }
                p[X0]   = xposf[ind];
                p[Y0]   = yposf[ind];
                p[PA]   = posaf[ind];
                p[INC]  = inclf[ind];
                p[VEXP] = vexpf[ind];
                p[VROT] = vrotf[ind];
                p[VSYS] = vsysf[ind];
                model[pp] = func(xx, p, MAXPAR);
            }
        }
    }
    
    model.fitswrite_2d(fname.c_str());
    thetaf = ot;
}


std::ostream& operator<< (std::ostream& Stream, Ringmodel& r) {
    
    r.print(Stream);
    return Stream;
}


float func (float *c, float *p, int npar) {
  
  /// The function calculates radial velocity from rotation curve.
  /// 
  /// \param  c     Grid position in plane of galaxy. Dimension = 2.
  /// \param  p     List of parameters of ring.
  /// \param  npar  Number of parameters.
  ///
  /// \return       The radial velocity in requested point.
  

    float   vs, vc, vr;                     // Parameters of velocity field. 
    float   x, y;                           // Sky coordinates.
    float   cost1, sint1; 
    float   x1, y1, r;
   
    static float phi= 0.0, inc = 0.0;       // Saved parameters (static type).
    static float cosp1 = 1.0;
    static float sinp1 = 0.0;
    static float cosi1 = 1.0;
    static float sini1 = 0.0;

    vs = p[VSYS];                           // Systemic velocity.
    vc = p[VROT];                           // Circular velocity.
    vr = p[VEXP];                           // Expansion velocity.
    
    if (p[PA] != phi) {                     // Position angle.
        phi = p[PA];                            
        float factor = M_PI/180.;
        cosp1 = cos(factor*phi);                        
        sinp1 = sin (factor*phi);                           
    }
   
    if (p[INC] != inc) {                    // Inclination.
        inc = p[INC];   
        float factor = M_PI/180.;           
        cosi1 = cos(factor*inc);                            
        sini1 = sin (factor*inc);               
    }
   
    x = c[0] - p[X0];                       // Calculate X. 
    y = c[1] - p[Y0];                       // Calculate Y.
    x1 = (-x*sinp1 + y*cosp1);              // X in plane of galaxy.
    y1 = (-x*cosp1 - y*sinp1)/cosi1;        // Y in plane of galaxy.
    r = sqrt(x1*x1 + y1*y1);                // Distance from centre.
  
    cost1 = x1 / r;                         // Azimutal angle (theta).
    sint1 = y1 / r;                         

    float vrad = vs + (vc*cost1 + vr*sint1) * sini1;
        
    return vrad;    
}


void derv (float *c, float *p, float *d, int npar) {
    
  /// The function derv calculates the partial derivatives 
  /// with respect to the parameters.
  ///
  /// \param  c     Grid position in plane of galaxy. Dimension = 2.
  /// \param  p     List of parameters of ring.
  /// \param  d     Partial derivatives, returned to user.
  /// \param  npar  Number of parameters.

    float   vc, vr;                         // Parameters of velocity field.
    float   x, y;                           // Sky coordinates.
    float   cost1, cost2, sint1, sint2;
    float   x1, y1, r;
    
    static float phi= 0.0, inc = 0.0;       // Saved parameters (static type).
    static float cosp1 = 1.0;
    static float sinp1 = 0.0;
    static float cosi1 = 1.0, cosi2 = 1.0;
    static float sini1 = 0.0, sini2 = 0.0;
    
    const double F = M_PI/180.;

    vc = p[VROT];                           // Circular velocity.
    vr = p[VEXP];                           // Expansion velocity.
    
    if (p[PA] != phi) {                     // Position angle.
        phi = p[PA];                            
        cosp1 = cos( F * phi );                 
        sinp1 = sin ( F * phi );                
    }
    
    if (p[INC] != inc) {                    // Inclination.
        inc = p[INC];               
        cosi1 = cos( F * inc );         
        cosi2 = cosi1 * cosi1;          
        sini1 = sin ( F * inc );            
        sini2 = sini1 * sini1;          
    }
    x = c[0] - p[X0];                       // Calculate X.
    y = c[1] - p[Y0];                       // Calculate Y
    x1 = (-x*sinp1 + y*cosp1);              // X in plane of galaxy.
    y1 = (-x*cosp1 - y*sinp1)/cosi1;        // Y in plane of galaxy.
    r = sqrt(x1*x1 + y1*y1);                // Distance from centre.
    
    cost1 = x1 / r;                         // Azimutal angle (theta).
    sint1 = y1 / r;             
    cost2 = cost1 * cost1;          
    sint2 = sint1 * sint1;          

    d[VSYS] = 1.0;                          // Now calculate derivatives.
    d[VROT] = sini1 * cost1;            
    d[VEXP] = sini1 * sint1;            
                    
    d[PA]   = F * vc * (1.0 - sini2 * sint2) * sint1 * sini1 / cosi1 -
              F * vr * (1.0 - sini2 * cost2) * cost1 * sini1 / cosi1;

    d[INC]  = F * vc * (cosi2 - sini2 * sint2) * cost1 / cosi1 +
              F * vr * (cosi2 + sini2 * cost2) * sint1 / cosi1;

    d[X0]   = vc * (sint1 * sinp1 - cost1 * cosp1 / cosi1) * sint1 * sini1 / r -
              vr * (sint1 * sinp1 - cost1 * cosp1 / cosi1) * cost1 * sini1 / r;
        
    d[Y0]   = - vc * (sint1 * cosp1 + cost1 * sinp1 / cosi1) * sint1 * sini1 / r +
                vr * (sint1 * cosp1 + cost1 * sinp1 / cosi1) * cost1 * sini1 / r;
}

