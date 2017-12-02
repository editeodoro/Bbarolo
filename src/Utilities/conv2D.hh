//--------------------------------------------------------------------
// conv2D.hh: Structure for 2d convolution using FFTW.
//--------------------------------------------------------------------

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

#ifndef CONV2D_HH
#define CONV2D_HH
 

#include <fftw3.h>
 
void factorize (const int n, int *n_fact, int fact[],int *implem_fact);
bool is_optimal(int n, int *implem_fact);
int find_closest_factor(int n, int *implem_fact);

typedef enum
{
    LINEAR_FULL,
    LINEAR_SAME_UNPADDED,
    LINEAR_SAME,
    LINEAR_VALID,
    CIRCULAR_SAME,
    CIRCULAR_SAME_SHIFTED,
    CIRCULAR_CUSTOM
} CONV_MODE;


struct Conv2D
{
    double  *in_src;            
    double  *out_src;
    double  *in_ker;
    double  *out_ker;
    int     h_src, w_src;
    int     h_ker, w_ker;
    double  *dst_fft;
    double  *dst;                       // The array containing the result.
    int     h_dst, w_dst;               // Height and width of output array. 
    int     w_fftw, h_fftw;
    CONV_MODE mode;
    fftw_plan p_forw_src;
    fftw_plan p_forw_ker;
    fftw_plan p_back;
};

void init_Conv2D (Conv2D &ws, CONV_MODE mode, int w_src, int h_src, int w_kernel, int h_kernel);
void clear_Conv2D (Conv2D &ws);
void circular_convolution(Conv2D &ws, double *src, double *kernel);
void convolve(Conv2D &ws, double *src, double *kernel);

#endif 


