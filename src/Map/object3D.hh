// -----------------------------------------------------------------------
// Object3D.hh: Definition of the Object3D class that holds pixel
//              information for a three-dimensional object.
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

#ifndef OBJECT3D_H
#define OBJECT3D_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include "voxel.hh"
#include "scan.hh"
#include "object2D.hh"


namespace PixelInfo
{

  /// A set of pixels in 3D.  
  /// 
  /// This stores the pixels in a STL map, connecting a
  /// channel number (z-value) to a unique Object2D. Also recorded
  /// are the average x-, y- and z-values (via their sums), as well as
  /// their extrema.
  
  template <class T>
  class Object3D
  {
  public:
    Object3D();									/// Default constructor.
    Object3D(const Object3D& o);				/// Copy constructor.
    Object3D& operator= (const Object3D& o);  	/// Assignement operator.
    virtual ~Object3D(){};						/// Default destructor.
	
	template <class Type>
	friend std::ostream& operator<< (std::ostream& theStream, Object3D<Type>& obj);	/// Overloaded <<.
    template <class Type>
    friend Object3D<Type> operator+ (Object3D<Type> lhs, Object3D<Type> rhs);						/// Overloaded +.
    
	/// Inline functions to access private member:

	long getXmin(){return xmin;}; 
    long getYmin(){return ymin;};
    long getZmin(){return zmin;};
    long getXmax(){return xmax;};
    long getYmax(){return ymax;};
    long getZmax(){return zmax;};
    long getNumChanMap(){return chanlist.size();};
    unsigned int getSize(){return numVox;};

    /// Utility functions:
     
    bool isInObject(long x, long y, long z);			/// Is a 3-D voxel in the Object? 
    bool isInObject(Voxel<T> v){return isInObject(v.getX(),v.getY(),v.getZ());};
   
    virtual void addPixel(long x, long y, long z);		/// Add a single 3-D voxel to the Object. 
    virtual void addPixel(Voxel<T> v){addPixel(v.getX(),v.getY(),v.getZ());};
   
    void addScan(Scan<T> s, long z);						/// Add a scan to the object.  
    void addChannel(const long &z, Object2D<T> &obj);		/// Add a full channel map to the Object. 

    void calcParams();						/// Calculate the averages and extrema of the three coordinates.     
    float getXaverage();					/// Return the average x-value.
    float getYaverage();					/// Return the average y-value.
    float getZaverage();					/// Return the average z-value.
    unsigned long getSpatialSize();		 	/// Return the number of spatial pixels.   
    std::vector<long> getChannelList();		/// Return a vector list of the channel numbers.
    int getMaxAdjacentChannels();	 		/// Return the largest number of adjacent channels.
    Object2D<T> getChanMap(long z);			/// Get the channel map for channel z.     
    Object2D<T> getSpatialMap();				/// Return an Object2D with spatial (x,y) distribution of voxels.  
    void print(std::ostream& theStream);	/// Class function to print to a stream.
    
    std::vector<Voxel<T> > getPixelSet();		/// Return a vector set of all voxels in the Object. 
    std::vector<Voxel<T> > getPixelSet(T *array, int *dim);  /// Return set of voxel + flux of array. 
    
    virtual void addOffsets(long xoff, long yoff, long zoff);

  protected:
    std::map<long,Object2D<T> > chanlist;  		///< The list of 2D channel maps
    unsigned long           numVox;    		///< How many voxels in the Object?
    float                   xSum;      		///< Sum of the x-values
    float                   ySum;      		///< Sum of the y-values
    float                   zSum;      		///< Sum of the z-values
    long                    xmin,xmax; 		///< min and max x-values of object
    long                    ymin,ymax; 		///< min and max y-values of object
    long                    zmin,zmax; 		///< min and max z-values of object
  }; 

}

//#include "object3D.cpp"

#endif
