// -----------------------------------------------------------------------
// scan.cpp: Member functions for the Scan class & other functions that	
//			 use the Scan class.
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
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#include <iostream>
#include <math.h>
#include "scan.hh"

namespace PixelInfo
{
  template <class T>
  Scan<T>::Scan() {
    
    itsY=-1;
    itsX=-1;
    itsXLen=0;
  }

  template <class T>
  Scan<T>::Scan(long y, long x, long xl):
    itsY(y), itsX(x),itsXLen(xl) {}

  template <class T>
  Scan<T>::Scan(const Scan<T>& s) {
    
    operator=(s);
  }
  

  template <class T>  
  Scan<T>& Scan<T>::operator= (const Scan<T>& s) {  
	    
    if(this == &s) return *this;
    this->itsY=s.itsY; 
    this->itsX=s.itsX; 
    this->itsXLen=s.itsXLen;
    return *this;
  }
 

  template <class T>
  bool Scan<T>::addScan(const Scan<T> &other) {
	  
    bool altered=touches(other);
    if(altered){
		long x = std::min(itsX,other.itsX);
		long xmax = std::max(itsX+itsXLen-1,other.itsX+other.itsXLen-1);
		itsX=x;
		itsXLen=xmax-x+1;
    }
    return altered;
  }


  template <class T> 
  Scan<T> nullScan() {
	  
    /// A simple way of returning a scan with zero length.
    Scan<T> null(-1,-1,0); 
    return null;
  }


  template <class T>
  bool Scan<T>::isNull() {
    
    return (itsY==-1 && itsX==-1 && itsXLen==0);
  }


  template <class T>
  Scan<T> unite(Scan<T> &scan1, Scan<T> &scan2) {
	  
  /// Return a scan that includes all pixels from both scans, but
  /// only if they overlap. If they do not, return the null scan.

    Scan<T> joined;
    if(!touching(scan1,scan2)){
		joined = nullScan<T>();
    }
    else {
		long y = scan1.getY();
		long x = std::min(scan1.getX(),scan2.getX());
		long xmax = std::max(scan1.getXmax(),scan2.getXmax());
		joined.define(y,x,xmax-x+1);
    }
    return joined;
  }


  template <class T>
  Scan<T> intersect(Scan<T> &scan1, Scan<T> &scan2) {
	  
  /// Return a scan that includes all pixels that lie in both scans. 
  /// If they do not overlap, return the null scan.

    Scan<T> intersection;
    if(!scan1.overlaps(scan2)){
		intersection = nullScan<T>();
    }
    else {
		long y = scan1.getY();
		long x = std::max(scan1.getX(),scan2.getX());
		long xmax = std::min(scan1.getXmax(),scan2.getXmax());
		intersection.define(y,x,xmax-x+1);
    }
    return intersection;
  }
  

  template <class T>
  bool Scan<T>::touches(const Scan<T> &other) {
	  
    return overlaps(other) || isAdjacentTo(other);
  }

 
  template <class T>
  bool Scan<T>::overlaps(const Scan<T> &other) {
	  
    if(itsY != other.itsY) return false;
    else if(itsX <= other.itsX){
		return (other.itsX < (itsX+itsXLen));
    }
    else {
		return (itsX < (other.itsX+other.itsXLen));
    }
  }


  template <class T>  
  bool Scan<T>::isAdjacentTo(const Scan<T> &other) {
	  
    if(itsY != other.itsY) return false;
    else if(itsX <= other.itsX){
		return (other.itsX == (itsX+itsXLen));
    }
    else {
		return (itsX == (other.itsX+other.itsXLen));
    }
  }


  template <class T>  
  bool Scan<T>::isInScan(long x, long y) {
	  
    return (y == itsY) && ((x>= itsX) && (x < (itsXLen+itsX)));
  }
  

  template <class T>  
  bool operator< (Scan<T> lhs, Scan<T> rhs) {
	  
  /// Test for less-than first on the y-values, and if they are
  /// equal, test on the starting x-value, and then finally on the
  /// length.

    if(lhs.itsY != rhs.itsY)      return (lhs.itsY    < rhs.itsY);
    else if(lhs.itsX != rhs.itsX) return (lhs.itsX    < rhs.itsX);
    else                          return (lhs.itsXLen < rhs.itsXLen);
  }
 
 
  template <class T>
  bool operator== (Scan<T> lhs, Scan<T> rhs) {
    
  /// For two scans to be equal, all three parameters must be equal.

    return (lhs.itsY == rhs.itsY) &&
      (lhs.itsX == rhs.itsX) &&
      (lhs.itsXLen == rhs.itsXLen);
  }
  

//=============================================================== 

  template <class T>
  bool touching(Scan<T> &scan1, Scan<T> &scan2) {
	  
  ///  Test whether two scans either overlap, or lie adjacent
  ///  (ie. there are no pixels lying between the two scans).

    return overlap(scan1,scan2) || adjacent(scan1,scan2);
  }


  template <class T>
  bool overlap(Scan<T> &scan1, Scan<T> &scan2) {
	  
  ///  Test whether two scans overlap, ie. they have pixels in
  ///  common.

    if(scan1.getY()!=scan2.getY()) return false;
    else if(scan1.getX() <= scan2.getX())
		return (scan2.getX() <= scan1.getXmax());
    else
		return (scan1.getX() <= scan2.getXmax());
  
  }
 
 
  template <class T>
  bool adjacent(Scan<T> &scan1, Scan<T> &scan2) {
	  
  /// Test whether two scans lie adjacent (ie. there are no pixels
  /// lying between the two scans).  If they overlap, return false.
  
    if(scan1.getY()!=scan2.getY()) return false;
    else if(scan1.getX() <= scan2.getX())
		return (scan2.getX() == scan1.getXmax()+1);
    else
		return (scan1.getX() == scan2.getXmax()+1);
  }


  template <class T>
  std::ostream& operator<< ( std::ostream& theStream, Scan<T>& scan) {
    
  /// Output the three key parameters of the scan.

    if(scan.isNull()) theStream << "NULL";
    else{
		theStream << scan.itsX;
		theStream << "-" << scan.getXmax();
		theStream << ", " << scan.itsY;
    }
    return theStream;
  }
 

  template <class T>
  float minSep(Scan<T> &s1, Scan<T> &s2) {
 
    if(s1.getX() > s2.getXmax()) return hypot(s1.getX()-s2.getXmax(),s1.getY()-s2.getY());
    else if(s2.getX() > s1.getXmax()) return hypot(s2.getX()-s1.getXmax(),s1.getY()-s2.getY());
    else return float(labs(s1.getY()-s2.getY()));
   
  }
  

  template <class T>  
  void mergeList(std::vector<Scan<T> > scanlist) {
	  
    typename std::vector<Scan<T> >::iterator iter;
    unsigned int counter=0,compCounter;
    while(counter<(scanlist.size()-1)) { 
		compCounter = counter+1;
		do {	
			if(touching(scanlist[counter],scanlist[compCounter])){
				Scan<T> temp = unite(scanlist[counter],scanlist[compCounter]);
				iter = scanlist.begin()+compCounter;
				scanlist.erase(iter);
				iter = scanlist.begin()+counter;
				scanlist.erase(iter);
				scanlist.push_back(temp);
			}
			else compCounter ++;
		} while(compCounter < scanlist.size());
		counter++;
    }

  }


}
