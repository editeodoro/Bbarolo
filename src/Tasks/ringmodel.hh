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
#include <Arrays/rings.hh>

enum ALLPARS {VSYS, VROT, VEXP, PA, INC, X0, Y0, MAXPAR};

template <class T>
class Ringmodel
///    A class to make a least-square fitting
///    of velocity field with a tilted-ring model
{                                       
public:
    Ringmodel();                        /// Default constructor.    
    Ringmodel(int nrings);              /// Alternative constructors.
    
    Ringmodel (Cube<T> *c);

    Ringmodel(int nrings, T *radii, T *widths, T *vsys, T *vrot,
              T *vexp, T *posang, T *incl, T xcenter, T ycenter);
              
    Ringmodel(int nrings, T *radii, T widths, T vsys, T vrot,
              T vexp, T posang, T incl, T xcenter, T ycenter);
    
    ~Ringmodel();                       /// Default destructor.
    
    void defaults();                    /// Default values for some variable.
    
    void set(int nrings, T *radii, T *widths, T *vsys, T *vrot,
             T *vexp, T *posang, T *incl, T xcenter, T ycenter);
              
    void set (int nrings, T *radii, T widths, T vsys, T vrot, 
              T vexp, T posang, T incl, T xcenter, T ycenter);
    
    void setfromCube (Cube<T> *c, Rings<T> *r);
    
    /// Inline functions to access private members:
    T&   getRadius(int i) {return rads [i];}
    T&   getWidth (int i) {return wids [i];}
    T&   getVsysf (int i) {return vsysf[i];}
    T&   getVrotf (int i) {return vrotf[i];}
    T&   getVexpf (int i) {return vexpf[i];}
    T&   getPosaf (int i) {return posaf[i];}
    T&   getInclf (int i) {return inclf[i];}
    T&   getXposf (int i) {return xposf[i];}
    T&   getYposf (int i) {return yposf[i];}
    
    T&   getVsyse (int i) {return vsyse[i];}
    T&   getVrote (int i) {return vrote[i];}
    T&   getVexpe (int i) {return vexpe[i];}
    T&   getPosae (int i) {return posae[i];}
    T&   getIncle (int i) {return incle[i];}
    T&   getXpose (int i) {return xpose[i];}
    T&   getYpose (int i) {return ypose[i];}
    
    T&   getChisq (int i) {return chis [i];}
    int& getNpts  (int i) {return npts [i];}
    int  getNradii () {return nrad;}
    
    float getMatrix (int nring, int a) {return elp[nring][a];}
    float **getMatrix () {return elp;}
        
    void hold(const int i) {mask[i]=false;}  /// Hold a parameter fixed
    void free(const int i) {mask[i]=true;}   /// Release a fixed parameter
    
    /// Fitting specific functions:
    
    void setoption (bool *maskpar, int hside, int wfunc, float freeangle);   
    void setfield (T *Array, int xsize, int ysize, int *boxup, int *boxlow);
    
    void ringfit();
    int  rotfit (T ri, T ro, T *p, T *e, int &n, T &q);
    int  getdat (std::vector<T> &x, std::vector<T> &y, std::vector<T> &w, T *p, T ri, T ro, T &q, int nfr);
    
    void print (std::ostream& Stream);
    void printfinal (std::ostream& Stream);
    void writeModel (std::string fname);

    template <class M>
    friend std::ostream& operator<< (std::ostream& Stream, Ringmodel<M> &r);
    
protected:

    long nrad;           ///< Number of rings.
    
    T    *rads;          ///< Radii of rings.
    T    *wids;          ///< Widths of rings.
    
    T    *vsysi;         ///< Initial systemic velocity.  
    T    *vsysf;         ///< Fitted systemic velocity.  
    T    *vsyse;         ///< Error in systemic velocity.  
    
    T    *vroti;         ///< Initial rotation velocity.
    T    *vrotf;         ///< Fitted rotation velocity.
    T    *vrote;         ///< Error rotation velocity.  
    
    T    *vexpi;         ///< Initial expansion velocity.    
    T    *vexpf;         ///< Fitted expansion velocity.
    T    *vexpe;         ///< Error in expansion velocity.
    
    T    *posai;         ///< Initial position angle (anticlockwise from north).
    T    *posaf;         ///< Fitted position angle.
    T    *posae;         ///< Error in position angle. 

    T    *incli;         ///< Initial inclination angle.
    T    *inclf;         ///< Fitted inclination angle.  
    T    *incle;         ///< Error in inclination angle.  
    
    T    xposi;          ///< Initial center position in X. 
    T    *xposf;         ///< Fitted center position in X.
    T    *xpose;         ///< Error in center position in X.
    
    T    yposi;          ///< Initial center position in Y. 
    T    *yposf;         ///< Fitted center position in Y. 
    T    *ypose;         ///< Error in center position in Y.  
    
    bool *mask;          ///< Parameter fit mask
    int  *npts;          ///< Number of points in each ring.
    T    *chis;          ///< Chi-squared for each ring.
    
    T     *vfield;        ///< Velocity field array.
    float **elp;          ///< Matrices of coefficients.
    
    bool  allAllocated;   ///< Have all array been allocated?
    bool  fieldAllocated; ///< Has the fitting field been constructed?

private:
    Cube<T> *in;
    int     blo[2];         ///< Lower edge of box.
    int     bup[2];         ///< Upper edge of box.
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
    
template <class T> T func (T *c, T *p, int npar);
template <class T> void  derv (T *c, T *p, T *d, int npar);


#endif
