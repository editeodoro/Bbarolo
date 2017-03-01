// -----------------------------------------------------------------------
// objectgrower.hh: Implementation of the object growing functions
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

#ifndef OBJECT_GROWER_H
#define OBJECT_GROWER_H

#include <iostream>
#include <Map/detection.hh>
#include <Map/voxel.hh>
#include <Arrays/cube.hh>
#include <Arrays/stats.hh>

/// @brief A class to manage the growing of objects to a secondary
/// threshold
/// @details This class provides a mechanism for handling the
/// growing of objects. By keeping track of the state of each pixel,
/// through an array of flags indicating whether a pixel is
/// available or in an object, it is able to efficiently grow the
/// objects from pixels on the edge, rather than spending time
/// examining pixels that are completely surrounded by other object
/// pixels.

template <class T>
class ObjectGrower
{
public:
	ObjectGrower();
    virtual ~ObjectGrower(){};
    ObjectGrower(ObjectGrower &o);
    ObjectGrower& operator=(const ObjectGrower &o);

    /// Set up the class with parameters & pointers from the cube
    void define(Cube<T> *theCube);
    /// Update a Cube's detectMap based on the flag array
    void updateDetectMap(short *map);
    /// Grow an object
    virtual void grow(Detection<T> *theObject);
    /// Grow out from a single voxel, returning the list of new voxels.
    std::vector<Voxel<T> > growFromPixel(Voxel<T> &vox);

protected:
    std::vector<STATE> itsFlagArray;                   	///< The array of pixel flags
    std::vector<size_t> itsArrayDim;                    ///< The dimensions of the array
    Statistics::Stats<T> itsGrowthStats;  				///< The statistics used to determine membership of an object
    int itsSpatialThresh;                              	///< The spatial threshold for merging
    int itsVelocityThresh;                             	///< The spectral threshold for merging
    T* itsFluxArray;                               		///< The location of the pixel values
};



//#include "objectgrower.cpp"
#endif
