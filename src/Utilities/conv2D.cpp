//--------------------------------------------------------------------
// conv2D.hh: Structure for 2d convolution using FFTW.
//--------------------------------------------------------------------

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

#include "conv2D.hh"
#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif
 
int FACTORS[7] = {13,11,7,5,3,2,0}; // end with zero to detect the end of the array

void init_Conv2D (Conv2D &ws, CONV_MODE mode, int w_src, int h_src, int w_kernel, int h_kernel) {

	ws.h_src = h_src;
	ws.w_src = w_src;
	ws.h_ker = h_kernel;
	ws.w_ker = w_kernel;
	ws.mode = mode;
 
	switch(mode) {
		case LINEAR_FULL:
		// Full Linear convolution
			ws.h_fftw = find_closest_factor(h_src+h_kernel-1, FACTORS);
	 		ws.w_fftw = find_closest_factor(w_src+w_kernel-1, FACTORS);
			ws.h_dst = h_src+h_kernel-1;
			ws.w_dst = w_src+w_kernel-1;
			break;
			
		case LINEAR_SAME_UNPADDED:
		// Same Linear convolution
			ws.h_fftw = h_src+int(h_kernel/2.0);
			ws.w_fftw = w_src+int(w_kernel/2.0);
			ws.h_dst = h_src;
			ws.w_dst = w_src;
			break;
		
		case LINEAR_SAME:
		// Same Linear convolution
			ws.h_fftw = find_closest_factor(h_src+int(h_kernel/2.0),FACTORS);
			ws.w_fftw = find_closest_factor(w_src+int(w_kernel/2.0),FACTORS);
			ws.h_dst = h_src;
			ws.w_dst = w_src;
			break;

		case LINEAR_VALID:
		// Valid Linear convolution
			if(ws.h_ker>ws.h_src || ws.w_ker>ws.w_src) {
				std::cout << "CONV2D warning: The 'valid' convolution results is an empty matrix\n";
				ws.h_fftw = 0;
				ws.w_fftw = 0;
				ws.h_dst = 0;
				ws.w_dst = 0;
			}
			else {
				ws.h_fftw = find_closest_factor(h_src, FACTORS);
				ws.w_fftw = find_closest_factor(w_src, FACTORS);
				ws.h_dst = h_src-h_kernel+1;
				ws.w_dst = w_src-w_kernel+1;
			}
			break;
			 
		case CIRCULAR_SAME:
		// Circular convolution, modulo N
		// We don't padd with zeros because if we do, we need to padd with at least h_kernel/2; w_kernel/2 elements
		// plus the wrapp around which in facts leads to too much computations compared to the gain obtained with the optimal size
			ws.h_fftw = h_src;
			ws.w_fftw = w_src;
			ws.h_dst = h_src;
			ws.w_dst = w_src;
			break;

		case CIRCULAR_SAME_SHIFTED:
		// Circular convolution, modulo N, shifted by M/2
		// We don't padd with zeros because if we do, we need to padd with at least h_kernel/2; w_kernel/2 elements
		// plus the wrapp around which in facts leads to too much computations compared to the gain obtained with the optimal size
			ws.h_fftw = h_src;
			ws.w_fftw = w_src;
			ws.h_dst = h_src;
			ws.w_dst = w_src;
			break;
			
		case CIRCULAR_CUSTOM:
		// We here want to compute a circular convolution modulo h_dst, w_dst
		// These two variables must have been set before calling init_workscape !!
			ws.h_fftw = ws.h_dst;
			ws.w_fftw = ws.w_dst;
			break;
 
  		default:
			std::cout << "CONV2D error: Unrecognized convolution mode, possible modes are:\n"
					  << "   - LINEAR_FULL \n"
					  << "   - LINEAR_SAME \n"
					  << "   - LINEAR_SAME_UNPADDED\n"
					  << "   - LINEAR_VALID \n"
					  << "   - CIRCULAR_SAME \n"
					  << "   - CIRCULAR_SHIFTED\n"
					  << "   - CIRCULAR_CUSTOM\n";
            std::terminate();
	}

	ws.in_src  = new double[ws.h_fftw*ws.w_fftw];
	ws.out_src = (double*) fftw_malloc(sizeof(fftw_complex)*ws.h_fftw*(ws.w_fftw/2+1));
	ws.in_ker  = new double[ws.h_fftw*ws.w_fftw];
	ws.out_ker = (double*) fftw_malloc(sizeof(fftw_complex)*ws.h_fftw*(ws.w_fftw/2+1));
	ws.dst_fft = new double[ws.h_fftw*ws.w_fftw];
	ws.dst     = new double[ws.h_dst*ws.w_dst];

//#ifdef _OPENMP	
//	fftw_init_threads();
//	fftw_plan_with_nthreads(omp_get_max_threads());
//#endif

	// Initialization of the plans
	ws.p_forw_src = fftw_plan_dft_r2c_2d(ws.h_fftw, ws.w_fftw, ws.in_src, (fftw_complex*)ws.out_src, FFTW_ESTIMATE);
	ws.p_forw_ker = fftw_plan_dft_r2c_2d(ws.h_fftw, ws.w_fftw, ws.in_ker, (fftw_complex*)ws.out_ker, FFTW_ESTIMATE);
 
	// The backward FFT takes ws.out_kernel as input !!
	ws.p_back = fftw_plan_dft_c2r_2d(ws.h_fftw, ws.w_fftw, (fftw_complex*)ws.out_ker, ws.dst_fft, FFTW_ESTIMATE);
}


void clear_Conv2D (Conv2D &ws) {
     
	delete[] ws.in_src;
	fftw_free((fftw_complex*)ws.out_src);    
	delete[] ws.in_ker;
	fftw_free((fftw_complex*)ws.out_ker);

	delete[] ws.dst_fft;
	delete[] ws.dst;
 
	// Destroy the plans
	fftw_destroy_plan(ws.p_forw_src);
	fftw_destroy_plan(ws.p_forw_ker);
	fftw_destroy_plan(ws.p_back);
	//fftw_cleanup_threads();
}
 
// Compute the circular convolution of src and kernel modulo ws.h_fftw, ws.w_fftw
// using the Fast Fourier Transform. The result is in ws.dst
void circular_convolution(Conv2D &ws, double *src, double *kernel) {

	double *ptr, *ptr_end, *ptr2;
 
	// Reset the content of ws.in_src
	for(ptr=ws.in_src, ptr_end=ws.in_src+ws.h_fftw*ws.w_fftw; ptr!=ptr_end; ++ptr) *ptr = 0.0;
	for(ptr=ws.in_ker, ptr_end=ws.in_ker+ws.h_fftw*ws.w_fftw; ptr!=ptr_end; ++ptr) *ptr = 0.0;

	// Then we build our periodic signals
	//ptr = src;
	for(int i=0; i<ws.h_src ; ++i)
		for(int j=0; j<ws.w_src; ++j, ++ptr)
			ws.in_src[(i%ws.h_fftw)*ws.w_fftw+(j%ws.w_fftw)] += src[i*ws.w_src+j];
			//ptr = kernel;
	for(int i=0; i < ws.h_ker; ++i)
		for(int j=0 ; j<ws.w_ker; ++j, ++ptr)
			ws.in_ker[(i%ws.h_fftw)*ws.w_fftw+(j%ws.w_fftw)] += kernel[i*ws.w_ker+j];
	
	// And we compute their packed FFT
	fftw_execute(ws.p_forw_src);
	fftw_execute(ws.p_forw_ker);
	
	// Compute the element-wise product on the packed terms
	// Let's put the element wise products in ws.in_kernel
	double re_s, im_s, re_k, im_k;
	for(ptr=ws.out_src, ptr2=ws.out_ker, ptr_end=ws.out_src+2*ws.h_fftw*(ws.w_fftw/2+1); ptr!=ptr_end; ++ptr, ++ptr2) {
		re_s = *ptr;
		im_s = *(++ptr);
		re_k = *ptr2;
		im_k = *(++ptr2);
		*(ptr2-1) = re_s*re_k-im_s*im_k;
		*ptr2 = re_s*im_k+im_s*re_k;
	}

	// Compute the backward FFT
	// Carefull, The backward FFT does not preserve the output
	fftw_execute(ws.p_back);
	
	// Scale the transform
	for(ptr=ws.dst_fft, ptr_end=ws.dst_fft+ws.w_fftw*ws.h_fftw; ptr!=ptr_end; ++ptr) 
		*ptr /= double(ws.h_fftw*ws.w_fftw);

}


void convolve(Conv2D &ws, double *src, double *kernel) {
	
	if(ws.h_fftw<=0 || ws.w_fftw<=0) return;
 
	// Compute the circular convolution
	circular_convolution(ws, src, kernel);
 
	// Depending on the type of convolution one is looking for, we extract the appropriate part of the result from out_src
	int h_offset, w_offset;

	switch(ws.mode) {
		case LINEAR_FULL:
		// Full Linear convolution
		// Here we just keep the first [0:h_dst-1 ; 0:w_dst-1] real part elements of out_src
			for(int i=0; i<ws.h_dst; ++i)
				memcpy(&ws.dst[i*ws.w_dst], &ws.dst_fft[i*ws.w_fftw], ws.w_dst*sizeof(double));
				break;
   
	  case LINEAR_SAME_UNPADDED:
		// Same linear convolution
		// Here we just keep the first [h_filt/2:h_filt/2+h_dst-1 ; w_filt/2:w_filt/2+w_dst-1] real part elements of out_src
			h_offset = int(ws.h_ker/2.0);
			w_offset = int(ws.w_ker/2.0);
			for(int i = 0 ; i < ws.h_dst ; ++i)
				memcpy(&ws.dst[i*ws.w_dst], &ws.dst_fft[(i+h_offset)*ws.w_fftw+w_offset], ws.w_dst*sizeof(double));
		   //for(int j = 0 ; j < ws.w_dst ; ++j)
		   //   ws.dst[i*ws.w_dst + j] = ws.out_src[2*((i+h_offset)*ws.w_fftw+j+w_offset)+0]/double(ws.w_fftw * ws.h_fftw);
			break;

		case LINEAR_SAME:
			// Same linear convolution
			// Here we just keep the first [h_filt/2:h_filt/2+h_dst-1 ; w_filt/2:w_filt/2+w_dst-1] real part elements of out_src
			h_offset = int(ws.h_ker/2.0);
			w_offset = int(ws.w_ker/2.0);
			for(int i = 0 ; i < ws.h_dst ; ++i)
				memcpy(&ws.dst[i*ws.w_dst], &ws.dst_fft[(i+h_offset)*ws.w_fftw+w_offset], ws.w_dst*sizeof(double));
			break;
			
		case LINEAR_VALID:
		// Valid linear convolution
		// Here we just take [h_dst x w_dst] elements starting at [h_kernel-1;w_kernel-1]
			h_offset = ws.h_ker-1;
			w_offset = ws.w_ker-1;
			for(int i=0; i<ws.h_dst; ++i)
			 	memcpy(&ws.dst[i*ws.w_dst], &ws.dst_fft[(i+h_offset)*ws.w_fftw+w_offset], ws.w_dst*sizeof(double));
			break;
			
		case CIRCULAR_SAME:
		// Circular convolution
		// We copy the first [0:h_dst-1 ; 0:w_dst-1] real part elements of out_src
			for(int i=0; i<ws.h_dst; ++i)
				memcpy(&ws.dst[i*ws.w_dst], &ws.dst_fft[i*ws.w_fftw], ws.w_dst*sizeof(double));
			break;

		case CIRCULAR_SAME_SHIFTED:
		// Shifted Circular convolution
		// We copy the [h_ker/2:h_ker/2+h_dst-1 ; w_ker/2:w_ker/2+w_dst-1] real part elements of out_src
			for(int i=0; i<ws.h_dst; ++i)
				for(int j=0; j<ws.w_dst; ++j)
					ws.dst[i*ws.w_dst+j] = ws.dst_fft[((i+int(ws.h_ker/2.0))%ws.h_fftw)*ws.w_fftw+(j+int(ws.w_ker/2.0))%ws.w_fftw];
			break;
			
		case CIRCULAR_CUSTOM:
			for(int i=0; i<ws.h_dst; ++i)
				memcpy(&ws.dst[i*ws.w_dst], &ws.dst_fft[i*ws.w_fftw], ws.w_dst*sizeof(double) );	
			break;
 
    	default:
			std::cout << "CONV2D error: Unrecognized convolution mode, possible modes are:\n"
					  << "   - LINEAR_FULL \n"
					  << "   - LINEAR_SAME \n"
					  << "   - LINEAR_SAME_UNPADDED\n"
					  << "   - LINEAR_VALID \n"
					  << "   - CIRCULAR_SAME \n"
					  << "   - CIRCULAR_SHIFTED\n"
					  << "   - CIRCULAR_CUSTOM\n";
            std::terminate();
	}
}
 
 

// ******************** Begin of factorization code ***********************//
// A piece of code to determine if a number "n" can be written as products of 
// only the integers given in implem_fact

void factorize (const int n, int *n_fact, int fact[],int *implem_fact) {
	
	int nf = 0;
	int ntest = n;
	int factor;
	int i = 0;
 
	if (n == 0) {
		printf("Length n must be positive integer\n");
		return;
	}
 
	if (n == 1) {
		fact[0] = 1;
		*n_fact = 1;
		return ;	
	}
	
	/* deal with the implemented factors */
	
	while (implem_fact[i] && ntest!=1) {
		factor = implem_fact[i];
		while ((ntest % factor) == 0) {
			ntest = ntest/factor;
			fact[nf] = factor;
			nf++;
		}
		i++;
	}
 
	// Ok that's it
	if(ntest!=1) {
		fact[nf] = ntest;
		nf++;
	}
 
	/* check that the factorization is correct */
	int product = 1;

	for (i=0; i<nf; i++) product *= fact[i];
  		if (product != n) std::cout << "Factorization failed";
	*n_fact = nf;
}


bool is_optimal(int n, int *implem_fact) {
	
	int nf=0;
	int fact[64];
	int i = 0;
	factorize(n, &nf, fact, implem_fact);
 
	// We just have to check if the last factor belongs to GSL_FACTORS
	while(implem_fact[i]) {
		if(fact[nf-1] == implem_fact[i]) return true;
		++i;
	}
	return false;
}
 
 
int find_closest_factor(int n, int *implem_fact) {

	int j;
	if(is_optimal(n,implem_fact)) return n;
	else {
		j = n+1;
		while(!is_optimal(j,implem_fact)) ++j;
		return j;
	}
}

