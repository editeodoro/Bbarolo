//----------------------------------------------------------------
// lsqfit.hh: Definition of Lsqfit class.
//----------------------------------------------------------------

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

#ifndef LSQFIT_HH_
#define LSQFIT_HH_

template <class T>
class Lsqfit
{
public:

    Lsqfit (T *x, int xd, T *y, T *w, int n, T *par, T *errpar,
            bool *maskpar, int numpar, T (*funk)(T*,T*, int), 
            void (*deriv)(T*, T*, T*, int), double tol=1.E-03, 
            int numiter=1000, double lab=1.E-03);
    
    ~Lsqfit ();
    
    double getChi(){return chi;}
    void hold(const int i, const T val) {mpar[i]=false; fpar[i]=val;}      /// Hold a parameter fixed
    void free(const int i) {mpar[i]=true;}                                 /// Release a fixed parameter

    int fit ();                                 /// Make the fit.
    int getvec ();                              /// Get correction vector.
    void getmat ();                             /// Calculates matrix of coefficients.
    int invmat ();                              /// Invert the matrix of coefficients.  


private:
    T       *xdat;                  ///< X values.
    int      xdim;                  ///< How many dimensions for X?
    T       *ydat;                  ///< Y values.
    T       *wdat;                  ///< Weights of data.
    int      ndat;                  ///< Number of data.
    T       *fpar;                  ///< Fitting parameters.
    T       *epar;                  ///< Parameter errors.
    bool    *mpar;                  ///< Mask for fixed parameters.
    int      npar;                  ///< Number of parameters.
    int      its;                   ///< Number of iterations.
    double   oldchi;                ///< Old reduced chi-squared.
    double   chi;                   ///< New reduced chi-squared.
    double   labda;                 ///< Mixing parameter.
    double   tolerance;             ///< Accuracy.
    int      nfree;                 ///< Number of free parameters. 
    int     *parptr;                ///< Parameter pointer.
    double  *vector;                ///< Correction vector.
    double  **matrix1;              ///< Matrix of coefficients.
    double  **matrix2;              ///< Inverse of matrix1.
    bool     allAllocated;          ///< Has memory been allocated?
    
    T       (*func)(T *c, T *p, int numpar);                ///< Function to fit.
    void    (*derv)(T *c, T *p, T *d, int numpar);      ///< Parameters derivatives.
    
};


// Some fitting functions
template <class T>
T coreExp (T *c, T *p, int npar) {
    return p[0]*(p[1]+1)/(p[1]+exp(c[0]/p[2]));
}


template <class T>
void coreExpd (T *c, T *p, T *d, int npar) {
    T expn = exp(c[0]/p[2]);
    T denm = p[1]+expn;
    d[0] = (p[1]+1)/denm;
    d[1] = p[0]*(expn-1)/(denm*denm);
    d[2] = p[0]*c[0]*(p[1]+1)*expn/(p[2]*p[2]*denm*denm);
}


template <class T>
T polyn (T *c, T *p, int npar) {
    T value=0;
    for (int i=0; i<npar; i++) value += p[i]*std::pow(double(c[0]),double(i));
    return value;
}

template <class T>
T dpolyn_dx (T *c, T *p, int npar) {
    T value=0;
    for (int i=1; i<npar; i++) value += p[i]*i*std::pow(double(c[0]),double(i-1));
    return value;
}

template <class T>
void polynd (T *c, T *p, T *d, int npar) {
    for (int i=0; i<npar; i++) d[i]=std::pow(double(c[0]),double(i));
}


template <class T>
void fpolyn (T x, T *p, int npar, T &y, T *dydp) {

    T value=0;
    for (int i=0; i<npar; i++) {
        value += p[i]*std::pow(double(x),double(i));
        dydp[i]=std::pow(double(x),double(i));
    }
    y = value;
}

#endif
