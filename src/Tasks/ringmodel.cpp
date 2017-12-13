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
#include <iomanip>
#include <Tasks/ringmodel.hh>
#include <Arrays/cube.hh>
#include <Arrays/param.hh>
#include <Tasks/moment.hh>
#include <Utilities/lsqfit.hh>
#include <Utilities/utils.hh>
#include <Utilities/progressbar.hh>



#define MAXPIX  8192
#define MAXPAR  7

#define max(x,y)    (x>y ? x:y)
#define min(x,y)    (x<y ? x:y)
#define nint(x)     (x>0.0?(int)(x+0.5):(int)(x-0.5))


void Ringmodel::defaults() {
    
    tol = 0.001;
    cor[0] = -1;
    cor[1] = -1;    
    
    mask = new bool[7];
    for (int i=0; i<7; i++) mask[i] = true;
    
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
    int nrings = p.NRADII-1;
    float widths = p.RADSEP/(arcsconv(c->Head().Cunit(0))*c->Head().PixScale());
    float vsys = atof(p.VSYS.c_str());
    float vrot = atof(p.VROT.c_str());
    float vexp = p.VRAD=="-1" ? 0 : atof(p.VRAD.c_str());
    float posang = atof(p.PHI.c_str());
    float incl = atof(p.INC.c_str());
    float xcenter = atof(p.XPOS.c_str());
    float ycenter = atof(p.YPOS.c_str());
    float *radii = new float[nrings];

    for (int i=0; i<nrings; i++) radii[i]=(i+1)*widths-widths/2.;

    set(nrings, radii, widths, vsys, vrot, vexp, posang, incl, xcenter, ycenter);


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

    MomentMap<float> map;
    map.input(c);
    map.FirstMoment(true);
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
    
    for (int i=0; i<7; i++) mask[i] = maskpar[i];
    
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
   
    int iring;
    
    
    if (fieldAllocated && allAllocated) {
        
        ProgressBar bar;
        
        bar.init(nrad);
        
        for (nfit = 0, iring = 0; iring < nrad; iring++) {
            
            bar.update(iring+1);
            
            int n;
            float   e[MAXPAR];
            float   p[MAXPAR];
            float   ri, ro;
            float   q = 0.0;

            p[0] = vsysi[iring];
            p[1] = vroti[iring];
            p[2] = vexpi[iring];
            p[3] = posai[iring];
            p[4] = incli[iring];
            p[5] = xposi;
            p[6] = yposi;
            ri = rads[iring] - 0.5 * wids[iring];
            ro = rads[iring] + 0.5 * wids[iring];
            if ( ri < 0.0 ) ri = 0.0;
        
            if (rotfit(ri, ro, p, e, n, q)>0) {
                if (nfit < iring) {
                    rads[nfit] = rads[iring];
                    wids[nfit] = wids[iring];
                }
            
                vsysf[nfit] = p[0];
                if (e[0] < 999.99) vsyse[nfit] = e[0];
                else vsyse[nfit] = 999.99;
            
                vrotf[nfit] = p[1];
                if (e[1] < 999.99) vrote[nfit] = e[1];
                else vrote[nfit] = 999.99;
            
                vexpf[nfit] = p[2];
                if (e[2] < 999.99) vexpe[nfit] = e[2];
                else vexpe[nfit] = 999.99;
            
                posaf[nfit] = p[3];
                if (e[3] < 999.99) posae[nfit] = e[3];
                else posae[nfit] = 999.99;
            
                inclf[nfit] = p[4];
                if (e[4] < 999.99) incle[nfit] = e[4];
                else incle[nfit] = 999.99;
            
                xposf[nfit] = p[5];
                if (e[5] < 999.99) xpose[nfit] = e[5];
                else xpose[nfit] = 999.99;
            
                yposf[nfit] = p[6];
                if (e[6] < 999.99) ypose[nfit] = e[6];
                else ypose[nfit] = 999.99;
            
                elp[nfit][0] = elp4[0];
                elp[nfit][1] = elp4[1];
                elp[nfit][2] = elp4[2];
                elp[nfit][3] = elp4[3];
                npts[nfit] = n;
                chis[nfit] = q;
                nfit += 1;
            }   
        }
        
        bar.fillSpace(" Done.");
    
    }
    else {
        std::cout << "Tilted Ring Error: Arrays are not allocated!\n";
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
    
    float *x = new float [2*MAXPIX];    // (x,y) position.
    float *y = new float [MAXPIX];      // f(x,y).
    float *w = new float [MAXPIX];      // w(x,y).
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
        
        Lsqfit<float> lsqfit(x, xdim, y, w, n, pf, e, mask, npar, &func, &derv, tol, t, lab);
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
    
    delete [] x;
    delete [] y;
    delete [] w;
    delete [] pf;
    
    return ier;
}


int Ringmodel::getdat (float *x, float *y, float *w, float *p, 
                              float ri, float ro, float &q, int nfr) {
    
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
  
    int     n=0;                        // Return value.
    int     llo, lup, mlo, mup;         // Positions of corners.
    int     nlt;                    // Box sizes. 
   
    float   phi, inc, x0, y0;           // Elliptical parameters.
    float   free;                       // Free angle.
    float   cosp, cosi, sinp, sini;
    float   a, b;
    float   wi;

    const double F = M_PI/180.;
    
    q    = 0.0;                         // Definition of parameters.
    phi  = p[PA];       
    inc  = p[INC];              
    x0   = p[X0];           
    y0   = p[Y0];           
    free = fabs(sin(F*thetaf));
    sinp = sin(F*phi);      
    cosp = cos(F*phi);  
    sini = sin(F*inc);      
    cosi = cos(F*inc );     
    a    = sqrt(1.0-cosp*cosp*sini*sini);
    b    = sqrt(1.0-sinp*sinp*sini*sini); 
    llo  = max(blo[0], nint(x0-a*ro));
    lup  = min(bup[0], nint(x0+a*ro));
    mlo  = max(blo[1], nint(y0-b*ro));
    mup  = min(bup[1], nint(y0+b*ro));

    
    if ((llo > lup) || (mlo > mup)) {
        std::cout << "Ring is outside the map!"<<std::endl;
        q = FLT_MAX;
        return 0;
    }
   
    nlt = bup[0] - blo[0] + 1;                  // Number of pixels in X.
   
    for (int m=mlo; m<mup; m++) {
        float ry = m;                           // Y position in plane of galaxy.
        for (int l = llo; l < lup; l++) {
            float rx = l;                       // X position in plane of galaxy.
            float v;
            int   ip = (m-blo[1])*nlt+ (l-blo[0]);
            v = vfield[ip];                     // Radial velocity at this position.
            
            if (!(v!=v)) {
                float costh, r, theta, xr, yr;
                xr = (-(rx-x0)*sinp + (ry-y0)*cosp);
                yr = (-(rx-x0)*cosp - (ry-y0)*sinp)/cosi;
                r = sqrt(xr*xr+yr*yr);
                if (r<0.1)  theta = 0.0;            
                else theta = atan2(yr, xr)/F;   
                costh = fabs(cos(F*theta));
                    
                if (r>ri && r<ro && costh>free ) {      // If we are inside the ring.
                    float xx[2];
                    int use = 0;
                    wi = std::pow(double(costh), (double) wpow);        // Calculate weight of this point.
                    xx[0] = rx;             
                    xx[1] = ry;             
                    
                    switch (side) {                     // Which side of galaxy.
                        
                        case 1:                         //< Receding half.                              
                        use = (fabs(theta)<=90.0);      
                        break;
                        
                        case 2:                         //< Approaching half. 
                        use = (fabs(theta)>=90.0);
                        break;
                 
                        case 3:                         //< Both halves.
                        use = 1;
                        break;
                 
                        default: 
                        break;  
                    }
               
                    if (use) {                          // Load data point ? 
                        n += 1;
                        if (n<MAXPIX) {             
                            float s, vz;
                            int maxpar = MAXPAR;
                            vz = func (xx, p, maxpar);
                            s = v - vz;                 // Corrected difference.
                            x[2*n-2] = rx;              // Load X-coordinate.
                            x[2*n-1] = ry;              // Load Y-coordinate.
                            y[n-1] = v;                 // Load radial velocity.
                            w[n-1] = wi;                // Load weight.
                            q += s*s*wi;                // Calculate chi-squared.
                        }
                    }
                }
            }
        }   
    }

    
    if (n>MAXPIX) {
        std::cout << "Too many points in ring "<< n << ". Maximum is "<< MAXPIX; 
        n = MAXPIX;
    }
    
    if (n > nfr)                                        // Enough data points ? 
        q = sqrt(q/(float)(n-nfr));     
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
        
            for (int j=0; j<7; j++) 
                if (!mask[j]) {
                    
                    if (j==0) Stream << "  Vsys = " << vsysi[i];
                    if (j==1) Stream << "  Vrot = " << vroti[i];
                    if (j==2) Stream << "  Vexp = " << vexpi[i];
                    if (j==3) Stream << "  P.A. = " << posai[i];
                    if (j==4) Stream << "  Inclination = " << incli[i];
                    if (j==5) Stream << "  Xcenter = " << xposi;
                    if (j==6) Stream << "  Ycenter = " << yposi;
            }
            
            Stream <<"\n  Reduced chi-squared: " << chis[i];
            
            Stream  << endl << endl << endl << endl << setfill('-');
        }
    }
}

void Ringmodel::printfinal (std::ostream& Stream) {

    int m=10;
    Stream  << std::endl;
    Stream  << setw(m) << right << "#Radius"
            << setw(m) << right << "Radius"
            << setw(m) << right << "Vsys "
            << setw(m) << right << "Vrot "
            << setw(m) << right << "Vexp "
            << setw(m) << right << "P.A."
            << setw(m) << right << "Incl."
            << setw(m) << right << "Xcenter"
            << setw(m) << right << "Ycenter" << endl;

    Stream  << setw(m) << right  << "#[pix] "
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

