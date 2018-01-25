//----------------------------------------------------------
// ringmodel.hh: Definition of Ringmodel class.
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

#ifndef RINGMODEL_HH_
#define RINGMODEL_HH_

#include <vector>
#include <fstream>
#include <Arrays/cube.hh>

enum ALLPARS {VSYS, VROT, VEXP, PA, INC, X0, Y0};

class Ringmodel                         /// A class to make a least-square fitting
{                                       /// of velocity field with a tilted rings
public:                                 /// model.  
    
    Ringmodel();                        /// Default constructor.    
    Ringmodel(int nrings);              /// Alternative constructor.
    
    Ringmodel (Cube<float> *c);

    Ringmodel(int nrings, float *radii, float *widths, float *vsys, float *vrot,
              float *vexp, float *posang, float *incl, float xcenter, float ycenter);
              
    Ringmodel(int nrings, float *radii, float widths, float vsys, float vrot,
              float vexp, float posang, float incl, float xcenter, float ycenter);
    
    ~Ringmodel() {};                    /// Default destructor.
    
    void defaults();                    /// Default values for some variable.
    
    void set(int nrings, float *radii, float *widths, float *vsys, float *vrot,
             float *vexp, float *posang, float *incl, float xcenter, float ycenter);
              
    void set (int nrings, float *radii, float widths, float vsys, float vrot, 
              float vexp, float posang, float incl, float xcenter, float ycenter);
    
    /// Inline functions to access private members:
    float&   getRadius(int i) {return rads [i];};
    float&   getWidth (int i) {return wids [i];};
    float&   getVsysf (int i) {return vsysf[i];};
    float&   getVrotf (int i) {return vrotf[i];};
    float&   getVexpf (int i) {return vexpf[i];};
    float&   getPosaf (int i) {return posaf[i];};
    float&   getInclf (int i) {return inclf[i];};
    float&   getXposf (int i) {return xposf[i];};
    float&   getYposf (int i) {return yposf[i];};
    
    float&   getVsyse (int i) {return vsyse[i];};
    float&   getVrote (int i) {return vrote[i];};
    float&   getVexpe (int i) {return vexpe[i];};
    float&   getPosae (int i) {return posae[i];};
    float&   getIncle (int i) {return incle[i];};
    float&   getXpose (int i) {return xpose[i];};
    float&   getYpose (int i) {return ypose[i];};
    
    
    float&   getChisq (int i) {return chis [i];};
    int&     getNpts  (int i) {return npts [i];};
    int     getNradii () {return nrad;};
    
    float   getMatrix (int nring, int a) {return elp[nring][a];};
    float   **getMatrix () {return elp;};
        
    void    hold(const int i) {mask[i]=false;}            /// Hold a parameter fixed
    void    free(const int i) {mask[i]=true;}                              /// Release a fixed parameter
    
    /// Fitting specific functions:
    
    void    setoption (bool *maskpar, int hside, int wfunc, float freeangle);   
    void    setfield (float *Array, int xsize, int ysize, int *boxup, int *boxlow);
    
    void    ringfit();
    int     rotfit (float ri, float ro, float *p, float *e, int &n, float &q);
    int     getdat (std::vector<float> &x, std::vector<float> &y, std::vector<float> &w, float *p, float ri, float ro, float &q, int nfr);
    
    void    print (std::ostream& Stream);
    void    printfinal (std::ostream& Stream);
    void    writeModel (std::string fname);

    friend std::ostream& operator<< (std::ostream& Stream, Ringmodel& r);
    
protected:

    long    nrad;           ///< Number of rings.
    
    float   *rads;          ///< Radii of rings.
    float   *wids;          ///< Widths of rings.
    
    float   *vsysi;         ///< Initial systemic velocity.  
    float   *vsysf;         ///< Fitted systemic velocity.  
    float   *vsyse;         ///< Error in systemic velocity.  
    
    float   *vroti;         ///< Initial rotation velocity.
    float   *vrotf;         ///< Fitted rotation velocity.
    float   *vrote;         ///< Error rotation velocity.  
    
    float   *vexpi;         ///< Initial expansion velocity.    
    float   *vexpf;         ///< Fitted expansion velocity.
    float   *vexpe;         ///< Error in expansion velocity.
    
    float   *posai;         ///< Initial position angle (anticlockwise from north).
    float   *posaf;         ///< Fitted position angle.
    float   *posae;         ///< Error in position angle. 

    float   *incli;         ///< Initial inclination angle.
    float   *inclf;         ///< Fitted inclination angle.  
    float   *incle;         ///< Error in inclination angle.  
    
    float   xposi;          ///< Initial center position in X. 
    float   *xposf;         ///< Fitted center position in X.
    float   *xpose;         ///< Error in center position in X.
    
    float   yposi;          ///< Initial center position in Y. 
    float   *yposf;         ///< Fitted center position in Y. 
    float   *ypose;         ///< Error in center position in Y.  
    
    bool    *mask;          ///< Parameter fit mask
    int     *npts;          ///< Number of points in each ring.
    float   *chis;          ///< Chi-squared for each ring.
    
    float   *vfield;        ///< Velocity field array.
    float   **elp;          ///< Matrices of coefficients.
    
    bool    allAllocated;   ///< Have all array been allocated?
    bool    fieldAllocated; ///< Has the fitting field been constructed?

private:

    Cube<float> *in;
    int     blo[2];         ///< Lower edge of box.
    int     bup[2];         ///< Upper edge of box.
    int     nfit;           ///< Number of fits.
    int     side;           ///< Which half of velocity field:
    int     wpow;           ///< Weighting power (uniform, cos, cos^2);
    float   thetaf;         ///< Free angle.
    float   tol;            ///< Tolerance of fit.
    float   elp4[4];        ///< Matrix of coefficients.
    int     cor[2];         ///< Correlation ellipses.

};
  
  
  ///
  /// Definition of fitting function:
  ///
  /// 
  /// V(x,y) = VSYS + VROT*cos(theta)*sin(INCL) + VEXP*sin(theta)*sin(INCL)
  ///        
  ///
  ///                         - (x-XPOS) * sin(PA) + (y-YPOS) * cos(PA)
  /// with:     cos(theta) = ---------------------------------------------
  ///                                             r
  ///
  ///                          - (x-XPOS) * cos(PA) - (y-YPOS) * sin(PA)
  /// and:      sin(theta) = ---------------------------------------------
  ///                                            r * cos(INCL)
  ///
    
    float func (float *c, float *p, int npar);
    void  derv (float *c, float *p, float *d, int npar);


#endif
