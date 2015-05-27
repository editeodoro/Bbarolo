//----------------------------------------------------------------
// lsqfit.hh: Definition of Lsqfit class.
//----------------------------------------------------------------

/*-----------------------------------------------------------------------
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 Bbarp;p is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with Bbarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning Bbarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
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
	
	void hold(const int i, const T val) {mpar[i]=false; fpar[i]=val;};		/// Hold a parameter fixed
	void free(const int i) {mpar[i]=true;};									/// Release a fixed parameter

	int fit ();									/// Make the fit.
	int getvec ();								/// Get correction vector.
	void getmat ();								/// Calculates matrix of coefficients.
	int invmat ();								/// Invert the matrix of coefficients.	


private:

	T	 	*xdat;					///< X values.
	int		 xdim;					///< How many dimensions for X?
	T		*ydat;					///< Y values.
	T		*wdat;					///< Weights of data.
	int		 ndat;					///< Number of data.
	T		*fpar;					///< Fitting parameters.
	T		*epar;					///< Parameter errors.
	bool	*mpar;					///< Mask for fixed parameters.
	int		 npar;					///< Number of parameters.
	int		 its;					///< Number of iterations.
	double	 oldchi;				///< Old reduced chi-squared.
	double	 chi;					///< New reduced chi-squared.
	double	 labda;					///< Mixing parameter.
	double	 tolerance;				///< Accuracy.
	int		 nfree;					///< Number of free parameters.	
	int		*parptr;				///< Parameter pointer.
	double	*vector;				///< Correction vector.
	double	**matrix1;				///< Matrix of coefficients.
	double	**matrix2;				///< Inverse of matrix1.
	bool  	 allAllocated;			///< Has memory been allocated?
	
	T 		(*func)(T *c, T *p, int numpar);				///< Function to fit.
	void  	(*derv)(T *c, T *p, T *d, int numpar);		///< Parameters derivatives.
	
};

#include "lsqfit.cpp"

#endif
