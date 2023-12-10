// -----------------------------------------------------------------------
// Scan.hh: Definition of the Scan class, used to store row
//          information as part of a 2D object.

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

#ifndef SCAN_H
#define SCAN_H

#include <iostream>
#include <vector>

/// This namespace will store all the classes and functions necessary
/// to encode shapes and objects in 1-, 2- and 3-dimensions.
namespace PixelInfo
{

  /// A class to store the basic unit of pixel information, a scan
  /// encoded by an (x,y) starting point, and a length (in the
  /// x-direction).
  /// 
  /// This class is used by other classes to store objects in 2- and
  /// 3-dimensions.
  template <class T>
  class Scan
  {
  public:
    Scan();
    Scan(long y, long x, long xl);
    Scan(const Scan& s);
    Scan& operator= (const Scan& s);
    virtual ~Scan(){}
    
    /// Overloaded operators:
    template <class Type>
    friend std::ostream& operator<< ( std::ostream& theStream, Scan<Type>& scan);
    
    template <class Type>
    friend bool operator< (Scan<Type> lhs, Scan<Type> rhs); 
    
    template <class Type>
    friend bool operator== (Scan<Type> lhs, Scan<Type> rhs);

    /// Inline functions to access private members:
    
    long getY() {return itsY;}
    void setY(long l) {itsY=l;}
    long getX() {return itsX;}
    void setX(long l) {itsX=l;}
    long getXlen() {return itsXLen;}
    void setXlen(long l) {itsXLen=l;}
    void define(long y, long x, long xl){itsY=y; itsX=x; itsXLen=xl;}
    void clear(){itsY=-1;itsX=-1;itsXLen=0;}
    long getXmax(){return itsX+itsXLen-1;}
    void setXmax(long l){itsXLen = l-itsX+1;}
    void growLeft(){itsX--;itsXLen++;}
    void growRight(){itsXLen++;}
    void addOffsets(long xoff, long yoff){itsY+=yoff; itsX+=xoff;}
    
    /// Other utility functions:
    
    bool isNull();                      /// Whether a scan is null.
    bool addScan(const Scan &other);    /// Combine scans. 
    bool isInScan(long x, long y);      /// Tests whether a given (x,y) point is in the scan.
    bool touches(const Scan &other);
    bool overlaps(const Scan &other);
    bool isAdjacentTo(const Scan &other);
   

    template <class Type> friend class Object2D; 
    template <class Type> friend class Object3D;
  
  protected:
    long itsY;                          ///< The y-value of each point in the scan.
    long itsX;                          ///< The x-value of the start (left-hand end) of the scan.
    long itsXLen;                       ///< The length of the scan (number of pixels in the scan).

  };

  
  /// Other functions that use the class above:
  
  template <class T> Scan<T> unite(Scan<T> &s1, Scan<T> &s2);               /// Combine two scans into one. 
  template <class T> Scan<T> intersect(Scan<T> &s1, Scan<T> &s2);           /// Keep only the pixels in both the two scans.  
  template <class T> bool touching(Scan<T> &s1, Scan<T> &s2);           /// Test whether two scans either overlap or are adjacent.  
  template <class T> bool overlap(Scan<T> &s1, Scan<T> &s2);                /// Test whether two scans have pixels in common.   
  template <class T> bool adjacent(Scan<T> &scan1, Scan<T> &scan2);     /// Test whether two scans lie adjacent to each other. 
  template <class T> Scan<T> nullScan();                                /// Return the null scan, y=-1, x=-1, xlen=0.
  template <class T> void mergeList(std::vector<Scan<T> > scanlist);    /// Merge Scans that are touching. 
  template <class T> float minSep(Scan<T> &s1, Scan<T> &s2);                /// Get the minimum separation, in pixels, between two scans. 

}

#endif 
