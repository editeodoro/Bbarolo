//-------------------------------------------------------------------------
// interpolation.cpp: Some functions for interpolation and fitting
//-------------------------------------------------------------------------

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

#include <vector>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <Utilities/utils.hh>


template <class Type>
double *spline(Type *x, Type *y, int n, double yp1, double ypn) {

  /// Given arrays x and y containing a function, and given optionally values  
  /// yp1 and ypn for the first derivative of the interpolating function
  /// at points 1 and n, this routine returns an array that contains the
  /// second derivatives of the interpolating function at the tabulated points.
  /// If yp1 and/or ypn are larger than 1e90, the routine sets the corresponding 
  /// boundary condition for a natural spline, with zero second derivative on 
  /// that boundary.

    double p, qn, sig, un;
    
    double *y2 = new double [n];
    double *u  = new double[n-1];
    
    if (yp1 > 1.e90)
        y2[0]=u[0]=0.0;
    else {
        y2[0] = -0.5;
        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    
    for (int i=1; i<n-1; i++) {
        sig   = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p     = sig*y2[i-1]+2.0;
        y2[i] = (sig-1.0)/p;
        u[i]  = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]  = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    
    if (ypn > 1.e90)
        qn = un = 0.0;
    else {
        qn = 0.5;
        un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    
    y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    
    for (int k=n-2; k>=0; k--)
        y2[k] = y2[k]*y2[k+1]+u[k];
    
    delete [] u;
    
    return y2;
}
template double* spline (int*, int*, int, double, double);
template double* spline (long*, long*, int, double, double);
template double* spline (float*, float*, int, double, double);
template double* spline (double*, double*, int, double, double);


template <class Type>
double splint(Type *xa, Type *ya, Type *y2a, int n, Type x) {
    
  /// Given the arrays xa and ya, which tabulate a function,
  /// and given the array y2a, containing the second derivatives, 
  /// and given a value of x, this routine returns a cubic-spline 
  /// interpolated value y.
    
    int klo,khi,k;
    double h,b,a;
    double y;

    klo=0;
    khi=n-1;
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xa[k] > x) khi=k;
        else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) throw("Bad xa input to routine splint");
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;

    return y;
}
template double splint (int*, int*, int*, int, int);
template double splint (long*, long*, long*, int, long);
template double splint (float*, float*, float*, int, float);
template double splint (double*, double*, double*, int, double);


template <class T> 
int linear_reg(int num, T *x, T *y, float &m, float &merr, float &b, float &berr, float &r, int ilow, int ihigh) {
  
  /// Computes the linear best fit to data y= mx + b, where x and y
  /// are arrays of size num, m is the slope and b the y-intercept.
  /// The values used in the arrays are those from ilow to ihigh,
  /// ie. if the full arrays are being used, then ilow=0 and high=num-1.
  /// Returns the values of slope & intercept (with errors) as well as
  /// r, the regression coefficient.
  /// 
  /// \param num        Size of the x & y arrays.
  /// \param x          Array of abscissae.
  /// \param y          Array of ordinates.
  /// \param ilow       Minimum index of the arrays to be used 
  ///                   (ilow=0 means   start at the beginning).
  /// \param ihigh      Maximum index of the arrays to be used 
  ///                   (ihigh=num-1 means finish at the end).
  /// \param m          Returns value of the slope of the best fit line.
  /// \param merr       Returns value of the estimated error in the
  ///                   slope value.
  /// \param b          Returns value of the y-intercept of the best
  ///                   fit line.
  /// \param berr       Returns value of the estimated error in the
  ///                   value of the y-intercept.
  /// \param r          Returns the value of the regression coefficient.
  /// \return           If everything works, returns 0. If slope 
  ///                   is infinite (eg, all points have same x value), 
  ///                   returns 1.
 

    if (ilow>ihigh) {
        std::cout << "Error! linear_regression.cc :: ilow (" << ilow 
                  << ") > ihigh (" << ihigh << ")!!\n";
        return 1;
    }
    if (ihigh>num-1) {
        std::cout << "Error! linear_regression.cc :: ihigh (" <<ihigh
                  << ") out of bounds of array (>" << num-1 << ")!!\n";
        return 1;
    }
    if(ilow<0){
        std::cout << "Error! linear_regression.cc :: ilow (" << ilow
                  << ") < 0. !!\n";
        return 1;
    }

    double sumx,sumy,sumxx,sumxy,sumyy;
    sumx=0.;
    sumy=0.;
    sumxx=0.;
    sumxy=0.;
    sumyy=0.;
    int count=0;
    for (int i=ilow;i<=ihigh;i++){
        count++;
        sumx = sumx + x[i];
        sumy = sumy + y[i];
        sumxx = sumxx + x[i]*x[i];
        sumxy = sumxy + x[i]*y[i];
        sumyy = sumyy + y[i]*y[i];
    }

    const float SMALLTHING=1.e-6;
    if(fabs(count*sumxx-sumx*sumx)<SMALLTHING) return 1;
    else{

        m = (count*sumxy - sumx*sumy)/(count*sumxx - sumx*sumx);
        merr = count / (count*sumxx - sumx*sumx);

        b = (sumy*sumxx - sumxy*sumx)/(count*sumxx - sumx*sumx);
        berr = sumxx / (count*sumxx - sumx*sumx);
    
        r = (count*sumxy - sumx*sumy) /
        (sqrt(count*sumxx-sumx*sumx) * sqrt(count*sumyy-sumy*sumy) );

        return 0;

    }
}
template int linear_reg (int, int*, int*, float&, float&, float&, float&, float&, int, int);
template int linear_reg (int, long*, long*, float&, float&, float&, float&, float&, int, int);
template int linear_reg (int, float*, float*, float&, float&, float&, float&, float&, int, int);
template int linear_reg (int, double*, double*, float&, float&, float&, float&, float&, int, int);


template <class T>
int linear_reg(int num, std::vector<T> &x, std::vector<T> &y, float *p, float *err, float &r, int ilow, int ihigh) {
  
    if (ilow>ihigh) {
        std::cout << "Error! linear_regression.cc :: ilow (" << ilow 
                  << ") > ihigh (" << ihigh << ")!!\n";
        return 1;
    }
    if (ihigh>num-1) {
        std::cout << "Error! linear_regression.cc :: ihigh (" <<ihigh
                  << ") out of bounds of array (>" << num-1 << ")!!\n";
        return 1;
    }
    if(ilow<0){
        std::cout << "Error! linear_regression.cc :: ilow (" << ilow
                  << ") < 0. !!\n";
        return 1;
    }

    double sumx,sumy,sumxx,sumxy,sumyy;
    sumx=0.;
    sumy=0.;
    sumxx=0.;
    sumxy=0.;
    sumyy=0.;
    int count=0;
    for (int i=ilow;i<=ihigh;i++){
        count++;
        sumx = sumx + x[i];
        sumy = sumy + y[i];
        sumxx = sumxx + x[i]*x[i];
        sumxy = sumxy + x[i]*y[i];
        sumyy = sumyy + y[i]*y[i];
    }

    const float SMALLTHING=1.e-6;
    if(fabs(count*sumxx-sumx*sumx)<SMALLTHING) return 1;
    else{

        p[0] = (count*sumxy - sumx*sumy)/(count*sumxx - sumx*sumx);
        err[0] = count / (count*sumxx - sumx*sumx);

        p[1] = (sumy*sumxx - sumxy*sumx)/(count*sumxx - sumx*sumx);
        err[1] = sumxx / (count*sumxx - sumx*sumx);
    
        r = (count*sumxy - sumx*sumy) /
        (sqrt(count*sumxx-sumx*sumx) * sqrt(count*sumyy-sumy*sumy) );
        
        return 0;

    }
}
template int linear_reg (int, std::vector<int>&, std::vector<int>&, float*, float*, float&, int, int);
template int linear_reg (int, std::vector<long>&, std::vector<long>&, float*, float*, float&, int, int);
template int linear_reg (int, std::vector<float>&, std::vector<float>&, float*, float*, float&, int, int);
template int linear_reg (int, std::vector<double>&, std::vector<double>&, float*, float*, float&, int, int);


template <class Type> 
Type func_gauss (Type *x, Type *p, int numpar) {
    
    Type ex, arg;
    
    Type y = 0.0;
    
    for (int i=0; i < numpar; i+=3) {
        arg = (*x-p[i+1])/p[i+2];
        ex = exp (-arg*arg/2.);
        y += p[i]*ex;
    }   
    
    return y;
}
template int func_gauss (int*, int*, int);
template long func_gauss (long*, long*, int);
template float func_gauss (float*, float*, int);
template double func_gauss (double*, double*, int);


template <class Type> 
void derv_gauss (Type *x, Type *p, Type *dyda, int numpar) {
    
    Type fac, ex, arg;
    
    for (int i=0; i < numpar; i+=3) {
        arg = (*x-p[i+1])/p[i+2];
        ex = exp (-arg*arg);
        fac = p[i]*ex*2.*arg;
        dyda[i] = ex;
        dyda[i+1] = fac/p[i+2];
        dyda[i+2] = fac*arg/p[i+2];     
    }   
}
template void derv_gauss (int*, int*, int*, int);
template void derv_gauss (long*, long*, long*, int);
template void derv_gauss (float*, float*, float*, int);
template void derv_gauss (double*, double*, double*, int);


template <class T>
double** RotMatrices (T alpha, T beta, T gamma) {

    /// Return three rotation matrices along the axis.
    /// input are the three angles along the three directions.

    double **matr = allocate_2D<double>(3,9);

    double angle[3] = {alpha*M_PI/180.,
                        beta*M_PI/180.,
                        gamma*M_PI/180.};

    int e[3] = {0,0,0};

    for (int i=3; i--;) {
        if (i==0) {e[0]=1;e[1]=0;e[2]=0;}
        else if (i==1) {e[0]=0; e[1]=1;e[2]=0;}
        else if (i==2) {e[0]=0;e[1]=0;e[2]=1;}
        
        double a=angle[i];
        double cosa = cos(a);
        double sina = sin(a);
        
        matr[i][0]=e[0]*e[0]+(1-e[0]*e[0])*cosa;
        matr[i][1]=(1-cosa)*e[0]*e[1]-sina*e[2];
        matr[i][2]=(1-cosa)*e[0]*e[2]+sina*e[1];
        matr[i][3]=(1-cosa)*e[1]*e[0]+sina*e[2];
        matr[i][4]=e[1]*e[1]+(1-e[1]*e[1])*cosa;
        matr[i][5]=(1-cosa)*e[1]*e[2]-sina*e[0];
        matr[i][6]=(1-cosa)*e[2]*e[0]-sina*e[1];
        matr[i][7]=(1-cosa)*e[2]*e[1]+sina*e[0];
        matr[i][8]=e[2]*e[2]+(1-e[2]*e[2])*cosa;
    }

    return matr;

}
template double** RotMatrices (float,float,float);
template double** RotMatrices (double,double,double);


template <class T> 
double* MatrixProduct (T *M1, int *size1, T *M2, int *size2) {
    
    assert (size1[1]==size2[0]);
    
    int size[2] = {size1[0],size2[1]};
    double *matrprod = new double[size[0]*size[1]];
    int n=size1[1]; 
    
    for (int y=0; y<size[1]; y++) {
        for (int x=0; x<size[0];x++) {
            matrprod[x+y*size[0]]=0.;
            for (int i=0; i<n; i++) {
                matrprod[x+y*size[1]]+=M1[x+i*size1[0]]*M2[i+y*size2[0]];
            }
        }
    }
    
    return matrprod;

}
template double* MatrixProduct (float*,int*, float*, int*);
template double* MatrixProduct (double*,int*, double*, int*);


template <class T>
void bezier_interp(std::vector<T> x_in,  std::vector<T> y_in,
                   std::vector<T> &x_out, std::vector<T> &y_out,
                   int fp, int np, int ns) {

    /// This function interpolates with a bezier function of degree n,
    /// where n is the number of datapoints.
    ///
    /// x_in,y_in           Input datapoints
    /// x_out,y_out         The output datapoints for Bezier function
    /// fp                  First datapoints (default=0)
    /// np                  Number of datapoints (default=x_in.size())
    /// ns                  Number of samples (default=np)

    if (x_in.size()!=y_in.size()) {
        std::cout << "BEZIER Error: x_in & y_in have different dimensions.\n";
        exit(EXIT_FAILURE);
    }

    if (np==-1) np=x_in.size();
    if (ns==-1) ns=np;

    if (np>x_in.size()-fp) {
        std::cout << "BEZIER Error: number of datapoints is greater than the size of input vectors.\n";
        exit(EXIT_FAILURE);
    }

    if (x_out.size()!=ns || y_out.size()!=ns) {
        std::cout << "BEZIER Error: The dimensions of output bezier vector must be = n_samples.\n";
        exit(EXIT_FAILURE);
    }

    double *bc = cp_binomial(np);

    for (int i=0; i<ns; i++) {
           double sr = (double)i/(double)(ns-1);
           unsigned int n = np-1;
           float px=0., py=0.;

           if (sr == 0.0) {
               px = x_in[fp];
               py = y_in[fp];
           }
           else if (sr == 1.0) {
               px = x_in[fp+n];
               py = y_in[fp+n];
           }
           else {
               double log_dsr_to_the_n = n*log(1-sr);
               double log_sr_over_dsr  = log(sr)-log(1-sr);

               for (unsigned j=0; j<=n; j++) {
                   double u = exp(bc[j]+log_dsr_to_the_n+j*log_sr_over_dsr);
                   px += x_in[j+fp]*u;
                   py += y_in[j+fp]*u;
               }
           }

           x_out[i]=px;
           y_out[i]=py;
       }

    delete [] bc;
}
template void bezier_interp(std::vector<float>,std::vector<float>,std::vector<float>&, std::vector<float>&,int,int,int);
template void bezier_interp(std::vector<double>,std::vector<double>,std::vector<double>&, std::vector<double>&,int,int,int);


double *cp_binomial(int points) {

    int e = points;
    int n = points-1;
    double *coeff = new double[e];
    e = n/2;
    coeff[0] = 0.0;
    for (int k = 0; k<e; k++)
        coeff[k+1] = coeff[k] + log(((double)(n-k))/((double)(k+1)));

    for (int k=n; k>=e; k--) coeff[k] = coeff[n-k];

    return (coeff);
}

