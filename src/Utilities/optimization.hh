//------------------------------------------------------------------------
// optimization.hpp: Classes for optimization purposes.
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
    Internet email: enrico.diteodoro@gmail.com
-----------------------------------------------------------------------*/

#ifndef OPTIMIZATION_HH_
#define OPTIMIZATION_HH_

#include <iostream>
#include <vector>
#include <cmath>
#include <Utilities/allocator.hpp>

using namespace std;

class NelderMead
/// A class to perform multidimensional minimization with the
/// Downhill-Simplex method of Nelder and Mead.
/// Implementation is derived from Numerical Recipes in C++ 2007
{
public:
    NelderMead() {}
    NelderMead(const double ftoll) : ftol(ftoll) {}
    ~NelderMead();

    template <class T>
    double* minimize(vector<double> &point, double del, T &func);

    template <class T>
    double* minimize(vector<double> &point, vector<double> &dels, T &func);

    template <class T>
    double* minimize(double **pp, int Ndim, T &func);

private:
     const double ftol = 1E-05; // Fractional convergence toleranc
     int     nfunc = 0;         // The number of function evaluations.
     int     ndim;              // Number of parameters to optimize.
     double  fmin;              // Function value at the minimum.
     double  *y = nullptr;      // Function values at the vertices of the simplex.
     double  **p = nullptr;     // Current simplex.
     bool    pAllocated = false;// Has been p allocated.

     void get_psum(vector<double> &psum);

     template <class T>
     double functry(vector<double> &psum, int ihi, double fac, T &func);
};

NelderMead::~NelderMead(){
    if (y!=nullptr) delete [] y;
    if (pAllocated) deallocate_2D<double>(p,ndim+1);
}

template <class T>
double* NelderMead::minimize(vector<double> &point, double del, T &func) {
    // Multidimensional minimization of the function func(x), where x[0..ndim-1]
    // is a vector in ndim dimensions. The initial simplex is specified by a
    // point[0..ndim-1] and a constant displacement del along each coordinate
    // direction. Returned is the location of the minimum.

    vector<double> dels(point.size(),del);
    return minimize(point,dels,func);
}


template <class T>
double* NelderMead::minimize(vector<double> &point, vector<double> &dels, T &func) {

    // Alternative interface that takes different displacements dels[0..ndim-1]
    // in different directions for the initial simplex.

    if (pAllocated) deallocate_2D<double>(p,ndim+1);
    ndim = point.size();
    p = allocate_2D<double>(ndim+1,ndim);
    pAllocated = true;

    for (int i=0;i<ndim+1; i++) {
        for (int j=0; j<ndim; j++) p[i][j] = point[j];
        if (i!=0) p[i][i-1] += dels[i-1];
    }

    return minimize(p,ndim,func);
}


template <class T>
double* NelderMead::minimize(double **pp, int Ndim, T &func) {

    // Most general interface: initial simplex specified by the matrix
    // pp[0..ndim][0..ndim-1]. Its ndim+1 rows are ndim-dimensional vectors
    // that are the vertices of the starting simplex. {

    const int NMAX=5000;        // Maximum allowed number of evaluations.
    const double TINY=1.0e-10;
    int ihi,ilo,inhi;
    ndim = Ndim;
    p = pp;

    vector<double> psum(ndim), x(ndim);
    double *pmin = new double[ndim];

    if (y!=nullptr) delete [] y;
    y =  new double[ndim+1];

    for (int i=0;i<ndim+1; i++) {
        for (int j=0; j<ndim; j++) x[j]=p[i][j];
        y[i]=func(x);
    }

    nfunc=0;
    get_psum(psum);

    for (;;) {
        ilo=0;
        //First we must determine which point is the highest (worst),
        // next-highest, and lowest (best), by looping over the points in the simplex.
        ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
        for (int i=0; i<ndim+1; i++) {
            if (y[i]<=y[ilo]) ilo=i;
            if (y[i]>y[ihi]) {
                inhi=ihi;
                ihi=i;
            }
            else if (y[i]>y[inhi] && i!=ihi) inhi=i;
        }

        //Compute the fractional range from highest to lowest and return if satisfactory.
        double rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);

        if (rtol<ftol) {
            //If returning, put best point and value in slot 0.
            std::swap(y[0],y[ilo]);
            for (int i=0; i<ndim; i++) {
                std::swap(p[0][i],p[ilo][i]);
                pmin[i]=p[0][i];
            }
            fmin = y[0];
            return pmin;
        }

        if (nfunc>=NMAX) throw("NELDER-MEAD Optimizer: Maximum number of iterations exceeded");

        nfunc += 2;

        // Begin a new iteration. First extrapolate by a factor 􏰘1 through the face of
        // the simplex across from the high point, i.e., reflect the simplex from the high point.
        double ytry = functry(psum,ihi,-1.0,func);
        if (ytry<=y[ilo]) {
            // If it gives a value better than the best point, try an
            // additional extrapolation by a factor 2.
            ytry = functry(psum,ihi,2.0,func);
        }
        else if (ytry>=y[inhi]) {
            // Otherwise, if the reflected point is worse than the second
            // highest, look for an intermediate lower point, i.e. do a
            // on dimensional contraction.
            double ysave = y[ihi];
            ytry = functry(psum,ihi,0.5,func);
            if (ytry>=ysave) {
                // Can't seem to get rid of that high point. Better contract
                // around the lowest (best) point.
                for (int i=0; i<ndim+1; i++) {
                    if (i!=ilo) {
                        for (int j=0; j<ndim; j++)
                            p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                        y[i]=func(psum);
                    }
                }
                nfunc += ndim;
                get_psum(psum);
            }
        }
        else --nfunc;
    }
}


inline void NelderMead::get_psum(vector<double> &psum) {
    for (int j=0; j<ndim; j++) {
        double sum=0.0;
        for (int i=0;i<ndim+1;i++) sum += p[i][j];
        psum[j]=sum;
    }
}


template <class T>
double NelderMead::functry(vector<double> &psum, int ihi, double fac, T &func) {
    // Extrapolates by a factor fac through the face of the simplex across
    // from the high point, tries it, and replaces the high point if the new
    // point is better.

    vector<double> ptry(ndim);
    double fac1=(1.0-fac)/ndim;
    double fac2=fac1-fac;

    for (int j=0; j<ndim; j++)
        ptry[j] = psum[j]*fac1-p[ihi][j]*fac2;

    double ytry=func(ptry);
    if (ytry<y[ihi]) {
        // Evaluate the function at the trial point.
        // If it’s better than the highest, then replace the highest.
        y[ihi]=ytry;
        for (int j=0; j<ndim; j++) {
               psum[j] += ptry[j]-p[ihi][j];
               p[ihi][j]=ptry[j];
        }
    }
    return ytry;
}


#endif
