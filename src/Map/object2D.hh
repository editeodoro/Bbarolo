// -----------------------------------------------------------------------
// Object2D.hh: Definition of Object2D, a class to hold pixel
//              information for a 2-dimensional object.
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

#ifndef OBJECT2D_H
#define OBJECT2D_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
#include <Map/scan.hh>
#include <Map/voxel.hh>

namespace PixelInfo
{

  /// @brief A class to store the pixel information for an object
  /// defined in a 2-dimensional plane.
  /// 
  /// @details The object consists of a set of Scans (stored as a
  /// std::vector), together with basic information on the maximum,
  /// minimum, and average values (actually stored as xSum, ySum) for
  /// the x and y pixel values.
    
  class Object2D
  {
  public:
    Object2D() {};
    Object2D(const Object2D& o);
    Object2D& operator= (const Object2D& o);  
    virtual ~Object2D(){}
  
    /// @brief Clear the Object and set the number of pixels to zero.
    void  clear(){scanlist.clear(); numPix=0;}

    /// @brief Make sure there are no overlapping or adjacent scans in the Scan list. 
    void  cleanup();

    /// @brief Add a pixel to the Object, making sure no scans overlap afterwards. 
    void  addPixel(long &x, long &y);

    /// @brief Add a full Scan to the Object, making sure there are no overlapping scans afterwards. 
    void  addScan(Scan &scan);

    /// @brief Test whether a pixel (x,y) is in the Object. 
    bool  isInObject(long x, long y);
    /// @brief Test whether the (x,y) part of a Voxel is in the Object.
    template <class T>
    bool  isInObject(Voxel<T> v){return this->isInObject(v.getX(),v.getY());}
    /// @brief Test whether a Pixel is in the Object.
    template <class T>
    bool  isInObject(Pixel<T> p){return this->isInObject(p.getX(),p.getY());}

    /// @brief Test whether a given Scan overlaps with any pixels in the Object. 
    bool  scanOverlaps(Scan &scan);

    bool canMerge(Object2D &other, float threshS, bool flagAdj);
    bool isNear(Object2D &other, long gap);
    bool isClose(Object2D &other, float threshS, bool flagAdj);   

    /// @brief Return the number of pixels in the Object. 
    unsigned long  getSize(){return numPix;}

    /// @brief Return the number of Scans in the Object. 
    long  getNumScan(){return scanlist.size();}

    /// @brief Order the Scans in the list, using the < operator for Scans. 
    void  order(){std::stable_sort(scanlist.begin(),scanlist.end());}

    /// @brief Add values to the x- and y-axes, making sure to add the offsets to the sums and min/max values. 
    void  addOffsets(long xoff, long yoff);

    /// @brief Calculate the sums, mins, maxs for x&y -- should not be necessary as these are done when pixels & Scans are added. 
    void  calcParams();

    /// @brief Return the average x-value. 
    float getXaverage(){return xSum/float(numPix);}

    /// @brief Return the average y-value. 
    float getYaverage(){return ySum/float(numPix);}

    long  getXmin(){return xmin;}
    long  getYmin(){return ymin;}
    long  getXmax(){return xmax;}
    long  getYmax(){return ymax;}

    /// @brief Return the number of distinct y-values.
    long  getNumDistinctY();
    /// @brief Return the number of distinct x-values.
    long  getNumDistinctX();
  
    /// @brief Return a measurement of the primary position angle
    double getPositionAngle();
    /// @brief Return the lengths of the principal axes.
    double *getPrincipalAxes();

    /// @brief A stream output operator. 
    friend std::ostream& operator<< (std::ostream& theStream, Object2D& obj);

    /// @brief Adding two Objects together.
    friend Object2D operator+ (Object2D lhs, Object2D rhs);
    
    friend class Object3D; 
    template <class T> friend class Detection;

  protected:
    std::vector<Scan> scanlist;       ///< The list of Scans
    unsigned long     numPix = 0 ;    ///< Number of pixels in the Object
    float             xSum;           ///< Sum of x values
    float             ySum;           ///< Sum of y values
    long              xmin,xmax;      ///< min and max x-values of object
    long              ymin,ymax;      ///< min and max y-values of object
  };

}

#endif //OBJECT2D_H
