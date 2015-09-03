// -----------------------------------------------------------------------
// voxel.hh: Definition of the Voxel class, storing a single 3D voxel
//           plus an associated flux.
// -----------------------------------------------------------------------

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
 along with BBarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning BBarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#ifndef VOXEL_H
#define VOXEL_H

#include <iostream>

namespace PixelInfo
{
  
  template <class T>
  class Voxel							/// A class to describe 3-dimensional pixel, 
  {										/// with x,y,z position + flux.
  public:
     
    Voxel(){};								 /// Default constructor.
    Voxel(long x, long y, long z, T f); 	 /// Default constructor.										 
    Voxel(long x, long y, long z);		 /// Specific constructor with flux=0. 
   
    Voxel(const Voxel& v);					 /// Copy constructor. 
    Voxel& operator= (const Voxel& v);		 /// Assignment operator. 
    virtual ~Voxel(){};						 /// Default destructor.

	
	/// Overloaded operators: print out << and equality ==.
	template <class Type>	
    friend std::ostream& operator<<( std::ostream& theStream, Voxel<T>& vox);
    template <class Type>
    friend bool operator== (Voxel<T> lhs, Voxel<T> rhs);	 
    
    
    /// Accessor functions to private class members.
    
    void   setX(long x){itsX = x;};
    void   setY(long y){itsY = y;};
    void   setZ(long z){itsZ = z;};
    void   setF(T f){itsF = f;}; 
    void   setXY(long x, long y){itsX = x; itsY = y;}; 
    void   setXYZ(long x, long y, long z){itsX = x; itsY = y; itsZ = z;}; 
    void   setXYF(long x, long y, T f){itsX = x; itsY = y; itsF = f;};
    void   setXYZF(long x, long y, long z, T f){itsX = x; itsY = y; itsZ = z; itsF = f;};
    long   getX(){return itsX;};
    long   getY(){return itsY;};
    long   getZ(){return itsZ;};
    T	   getF(){return itsF;};
    									
    long   arrayIndex(long *dim);					/// Return an index value for an array 	
	bool   match(Voxel other);						/// Function to test for equality of positions only.
	  

  protected:
    long  itsX;         ///< x-position of pixel
    long  itsY;         ///< y-position of pixel
    long  itsZ;         ///< z-position of pixel
    T	  itsF;         ///< flux of pixel
  };

  //==========================================================================

  template <class T>
  class Pixel : public Voxel<T>		/// A derived class of Voxel class:
  {										///  a 2-dimensional pixel, with just x & y position + flux
  public:
    Pixel(){this->itsZ=0;};						/// Default constructor.
    Pixel(long x, long y, T f);			/// Specific constructor.
    Pixel(const Pixel& p);					/// Copy constructor
    Pixel& operator= (const Pixel& p);		/// Assignement operator =.
    virtual ~Pixel(){};						/// Default destructor.
    
    /// Accessor functions
    
    void  setXY(long x, long y){this->itsX = x; this->itsY = y;};
    void  setXYF(long x, long y, T f){this->itsX = x; this->itsY = y; this->itsF = f;};

  };

}

#include "voxel.cpp"

#endif 
