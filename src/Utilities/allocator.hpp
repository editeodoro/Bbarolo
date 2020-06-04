// -----------------------------------------------------------------------
// allocator.hpp: Routines for allocate and deallocate memory in arrays
// -----------------------------------------------------------------------

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

#ifndef ALLOCATOR_HPP_
#define ALLOCATOR_HPP_

#include <iostream>


template <class Type>
inline Type **allocate_2D (int xdim, int ydim) {
    
    Type **array2d = new Type*[xdim];
    for (int i=0; i<xdim; i++)
        array2d[i] = new Type[ydim];
    return array2d;
}


template <class Type>
inline void deallocate_2D (Type **array2d, int xdim) {
    
    for (int i=0; i<xdim; i++) delete [] array2d[i];
    delete [] array2d;
}


template <class Type>
inline Type ***allocate_3D (int xdim, int ydim, int zdim) {

    Type ***array3d = new Type**[xdim];
    
    for (int i=0; i<xdim; i++) {
        array3d[i] = new Type*[ydim];
        for (int j=0; j<ydim; j++) {
            array3d[i][j] = new Type[zdim];
            for (int k=0; k<zdim; k++) {
                array3d[i][j][k] = 0;
            }
        }
    }
    
    return array3d;
}


template <class Type>
inline void deallocate_3D (Type ***array3d, int xdim, int ydim) {

    for (int i=0; i<xdim; i++) {
        for (int j=0; j<ydim; j++) {
            delete [] array3d[i][j];
        }
        delete [] array3d[i];
    }
    
    delete [] array3d;
}

#endif
