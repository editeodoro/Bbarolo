//----------------------------------------------------------
// lsqfit.cpp: Member function for Lsqfit class.
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
#include <cmath>
#include <cfloat>
#include <Utilities/lsqfit.hh>
#include <Utilities/utils.hh>

#define LABMIN      1.0e-10         // Minimum value for labda.
#define LABFAC      10.0            // Labda step factor.
#define LABMAX      1.0e+10         // Maximum value for labda.
#define MAXPARAM    32              // Maximum number of parameters.


template <class T>
Lsqfit<T>::Lsqfit (T *x, int xd, T *y, T *w, int n, T *par, T *errpar,
              bool *maskpar, int numpar,  T (*funk)(T *, T *, int), 
              void (*deriv)(T *, T *, T *, int), double tol, int numiter, double lab) :
                
              xdat(x), xdim(xd), ydat(y), wdat(w), ndat(n), fpar(par), epar(errpar),
              mpar(maskpar), npar(numpar), its(numiter), func(funk), derv(deriv) 
{
                    
                    if (tol<(FLT_EPSILON*10.0)) tolerance = FLT_EPSILON * 10.0;
                    else tolerance = tol;           
                    
                    labda = fabs(lab)*LABFAC;
                    
                    matrix1 = allocate_2D<double>(npar, npar);
                    matrix2 = allocate_2D<double>(npar, npar);              
                    vector = new double [npar];     
                    parptr = new int [npar];    
                    allAllocated = true;                
}


template <class T>
Lsqfit<T>::~Lsqfit () { 
    
    if (allAllocated) { 
        deallocate_2D<double> (matrix1, npar);
        deallocate_2D<double> (matrix2, npar);
        delete [] vector;
        delete [] parptr;
        
    }
}


template <class T>
int Lsqfit<T>::fit () {
   
  /// This function make the fit.
  ///
  /// \return           The number of interation (>0) or an error code:
  ///                   -1 = Too many free parameters.
  ///                   -2 = No free parameters.
  ///                   -3 = Not enough degrees of freedom.
  ///                   -4 = Reached max number of iterations.
  ///                   -5 = Diagonal of matrix contains zeroes.
  ///                   -6 = Determinant of the coefficient matrix is zero.
  ///                   -7 = Square root of negative number.
  
    int r;
    static int itc, found, nuse;
    
    itc = 0;                
    found = 0;              
    nfree = 0;              
    nuse = 0;   
    
    for (int i=0; i<npar; i++) {
        if (mpar[i]) {
            if (nfree > MAXPARAM) return -1 ;       
            parptr[nfree++] = i;        
        }
    }
   
    if (nfree == 0) return -2;                  
   
    for (int n=0; n<ndat; n++) {
        if (wdat[n] > 0.0) nuse++;              // Legal weight.
    }
    
    if (nfree>=nuse) return -3;
    
    if (labda==0.0) {                           // Linear fit.
    
        for (int i=0; i<nfree; fpar[parptr[i++]]=0.0) ;
        getmat();
        r = getvec();
        
        if (r) return r;
      
        for (int i = 0; i<npar; i++) {
            fpar[i] = epar[i];                  // Save new parameters
            epar[i] = 0.0;                      // and set errors to zero.
        }
        
        oldchi = sqrt(oldchi/(double)(nuse-nfree));
        
        for (int i=0; i<nfree; i++) {
            if ((matrix1[i][i] <= 0.0) || (matrix2[i][i] <= 0.0)) return -7;
            epar[parptr[i]] = oldchi * sqrt( matrix2[i][i] ) / sqrt( matrix1[i][i]);
        }
    } 
    else {                                      // Non-linear fit.
        while (!found) {            
            if (itc++ == its) return -4 ;   
            getmat(); 
            if (labda>LABMIN) labda /= LABFAC;  
            r = getvec();
            if (r) return r;    
         
            while (oldchi >= chi) {      
                if (labda>LABMAX) break;
                labda *= LABFAC;    
                r = getvec();
                if (r) return r;
            }
            if (labda <= LABMAX) {      
                for (int i=0; i<npar; i++) fpar[i] = epar[i];
            }
            if (fabs(chi-oldchi)<=(tolerance*oldchi) || (labda>LABMAX)) {
                labda = 0.0;        
                getmat();
                r = getvec();
                if (r) return r;        
                
                for (int i=0; i<npar; i++) epar[i] = 0.0;
            
                chi = sqrt(chi/(double)(nuse - nfree));
                
                for (int i=0; i<nfree; i++) {
                    if ((matrix1[i][i]<=0.0) || (matrix2[i][i]<=0.0)) return -7;
                    epar[parptr[i]] = chi*sqrt(matrix2[i][i])/sqrt(matrix1[i][i]);
                }
                found = 1;
            }
        }
    }
    
    return itc; 
}


template <class T>
void Lsqfit<T>::getmat () {

  /// The function builds the matrix of coefficients.

    double wd, wn, yd;

    for (int j=0; j<nfree; j++) {
        vector[j] = 0.0;                
        for (int i=0; i<=j; i++) {      
            matrix1[j][i] = 0.0;            
        }
    }
    chi = 0.0;                  
   
    for (int n=0; n<ndat; n++) {        
        wn = wdat[n];
        if (wn > 0.0) {             
            yd = ydat[n] - func(&xdat[(xdim)*n], fpar, npar); 
            derv (&xdat[(xdim)*n], fpar, epar, npar);   
            chi += yd * yd * wn;            
            for (int j=0; j<nfree; j++) {
                wd = epar[parptr[j]] * wn;      
                vector[j] += yd * wd;       
                for (int i=0; i<=j; i++) {      
                    matrix1[j][i] += epar[parptr[i]] * wd;
                }
            }
        }
    }
}


template <class T>
int Lsqfit<T>::getvec () { 
    
  /// The function calculates the correction vector. The matrix has been built by
  /// getmat, we only have to rescale it for the current value for labda.
  /// The matrix is rescaled so that the diagonal gets the value 1 + labda.
  /// Next we calculate the inverse of the matrix and then the correction
  /// vector.
  
    double dj, dy, mii, mji, mjj, wn;
    int    r;

    for (int j=0; j<nfree; j++) {       
        mjj = matrix1[j][j];            
        if (mjj <= 0.0) return -5;       
        mjj = sqrt(mjj);
        for (int i=0; i<j; i++) {   
            mji = matrix1[j][i] / mjj / sqrt( matrix1[i][i] );
            matrix2[i][j] = matrix2[j][i] = mji;
        }
        matrix2[j][j] = 1.0 + labda;
    }
    
    r=invmat();
    if (r) return r;
    for (int i=0; i<npar; i++) epar[i] = fpar[i];
    
    for (int j=0; j<nfree; j++) {   
        dj = 0.0;                   
        mjj = matrix1[j][j];
        if (mjj <= 0.0) return -7;
        mjj = sqrt( mjj );
        for (int i = 0; i < nfree; i++) {
            mii = matrix1[i][i];
            if (mii <= 0.0) return -7;
            mii = sqrt(mii);
            dj += vector[i]*matrix2[j][i] / mjj / mii;
        }
        epar[parptr[j]] += dj;  
    }
    
    oldchi = 0.0;               
    for (int n=0; n<ndat; n++) {        
        wn = wdat[n];               
        if (wn > 0.0) {             
            dy = ydat[n] - func(&xdat[xdim*n], epar, npar);
            oldchi += wdat[n] * dy * dy;
        }
    }
    
    return 0;
}


template <class T>
int Lsqfit<T>::invmat () {

  /// The function calculates the inverse of matrix2. 
  /// The algorithm used is the Gauss-Jordan algorithm.
 
   double even, mjk, rowmax;
   double hv[npar];
   int   evin;
   int   per[npar];
   int   row;

   for (int i=0; i<nfree; i++) per[i] = i;
   for (int j=0; j<nfree; j++) {        
        rowmax = fabs( matrix2[j][j] );     
        row = j;                    
        for (int i=j+1; i<nfree; i++) {
            if (fabs(matrix2[i][j]) > rowmax) {
                rowmax = fabs( matrix2[i][j] );
                row = i;
            }
        }
        if (matrix2[row][j] == 0.0) return -6;
        if (row > j) {          
            for (int k=0; k<nfree; k++) {       
                even = matrix2[j][k];       
                matrix2[j][k] = matrix2[row][k];
                matrix2[row][k] = even;
            }
            evin = per[j];              
            per[j] = per[row];
            per[row] = evin;
        }
        even = 1.0 / matrix2[j][j]; 
        for (int i=0; i<nfree; i++) matrix2[i][j] *= even;
        matrix2[j][j] = even;
        for (int k=0; k<j; k++) {
            mjk = matrix2[j][k];
            for (int i=0; i<j; i++) matrix2[i][k] -= matrix2[i][j] * mjk;
            for (int i=j+1; i<nfree; i++) matrix2[i][k] -= matrix2[i][j] * mjk;
            matrix2[j][k] = -even * mjk;
        }
        for (int k=j+1; k<nfree; k++) {
            mjk = matrix2[j][k];
            for (int i=0; i<j; i++) matrix2[i][k] -= matrix2[i][j] * mjk;
            for (int i=j+1; i<nfree; i++) matrix2[i][k] -= matrix2[i][j] * mjk;
            matrix2[j][k] = -even * mjk;
        }
    }
    for (int i=0; i<nfree; i++) {       
        for (int k=0; k<nfree; k++) hv[per[k]] = matrix2[i][k]; 
        for (int k=0; k<nfree; k++) matrix2[i][k] = hv[k];
    }
   
    return 0;
}
