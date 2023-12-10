// -----------------------------------------------------------------------
// voxel.cpp: Member functions for the Voxel class and the Pixel class
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

#include <iostream>
#include <iomanip>
#include <Map/voxel.hh>

namespace PixelInfo
{
  
  template <class T>
  Voxel<T>::Voxel(long x, long y, long z, T f) :
    itsX(x), itsY(y), itsZ(z), itsF(f) {}

  
  template <class T>
  Voxel<T>::Voxel(long x, long y, long z) :
    itsX(x), itsY(y), itsZ(z) { 
    
    this->itsF=0.;
  }


  template <class T>
  Voxel<T>::Voxel(const Voxel& v) {
    
    operator=(v);
  }
  
  
  template <class T>
  Voxel<T>& Voxel<T>::operator= (const Voxel<T>& v) {
      
    if(this == &v) return *this;
    this->itsX=v.itsX; 
    this->itsY=v.itsY; 
    this->itsZ=v.itsZ; 
    this->itsF=v.itsF;
    return *this;
  }
  
 
  template <class T>
  std::ostream& operator<< (std::ostream& theStream, Voxel<T>& vox) {
      
  /// A convenient way of printing the coordinate and flux values of
  /// a voxel. They are all printed to a single line (with no
  /// carriage-return), with the flux to precision of 4.

    theStream << std::setw(4) << vox.itsX ;
    theStream << " " << std::setw(4) << vox.itsY;
    theStream << " " << std::setw(4) << vox.itsZ;
    theStream << std::setprecision(4);
    theStream << "  " << vox.itsF;
    return theStream;

  }
 

  template <class T>
  bool operator== (Voxel<T> lhs, Voxel<T> rhs) {
      
  /// For two voxels to be equal, all four parameters must be equal.

    return (lhs.itsX == rhs.itsX) &&
      (lhs.itsY == rhs.itsY) &&
      (lhs.itsZ == rhs.itsZ) &&
      (lhs.itsF == rhs.itsF);
  }
  
  
  template <class T>
  bool Voxel<T>::match(Voxel<T> other) {
      
  /// This function just tests for equality of position. The flux is ignored.

    return (this->itsX == other.itsX) &&
           (this->itsY == other.itsY) &&
           (this->itsZ == other.itsZ);
  }

  
  template <class T>
  long Voxel<T>::arrayIndex(long *dim) {
      
    ///  Return the index value corresponding to the Voxel 
    ///  for an array with dimensions given by dim.

    long ind = itsX + dim[0]*itsY + dim[0]*dim[1]*itsZ;
    return ind;

  }



//--------------------------------------------------------------------
//  Member functions for the Pixel class.
//--------------------------------------------------------------------

  template <class T>
  Pixel<T>::Pixel(long x, long y, T f)  { 
    
    this->itsX=x; 
    this->itsY=y; 
    this->itsF=f;
  }
  
  
  template <class T>
  Pixel<T>::Pixel(const Pixel<T>& p) : Voxel<T>(p) {
   
    operator=(p);
  }
  
  
  template <class T>
  Pixel<T>& Pixel<T>::operator= (const Pixel<T>& p) {
      
    if(this == &p) return *this;
    this->itsX=p.itsX; 
    this->itsY=p.itsY; 
    this->itsF=p.itsF;
    return *this;
  }


// Explicit instantiation of the class
  template class Pixel<short>;
  template class Pixel<int>;
  template class Pixel<long>;
  template class Pixel<float>;
  template class Pixel<double>;
  
  template class Voxel<short>;
  template class Voxel<int>;
  template class Voxel<long>;
  template class Voxel<float>;
  template class Voxel<double>;

}
