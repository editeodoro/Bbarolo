//--------------------------------------------------------------------
// cube.hh: Definition of Search class
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

#ifndef SEARCH_HH_
#define SEARCH_HH_

#include <iostream>
#include <string>
#include <vector>
#include <Arrays/stats.hh>
#include <Arrays/param.hh>
#include <Map/detection.hh>

using namespace Statistics;
using namespace std;

// Some convenient aliases for detection classes
template <typename T> using DetVec   = vector<Detection<T> >;
template <typename T> using ScanVec  = vector<PixelInfo::Scan<T> >;
template <typename T> using Obj2DVec = vector<Object2D<T> >;


/////////////////////////////////////////////////////////////////////////////////////
/// A class to search for sources in a datacube or in an image
/////////////////////////////////////////////////////////////////////////////////////
template <class T>
class Search
{
/// Search is a class to find sources a 3D or a 2D array. A detection is found
/// based on a threshold and sources are reconstructed based on proximity criteria.
/// A 3D array can be searched along 2D maps (Lutz+80) or along 1D spectra.
/// The algorithms contained in this class are based on Duchamp source-finder
/// (see Whiting+12)
///
/// The correct way to use the Search class is the following:
///
/// 1) call a constructor to set SEARCH paramters:
///
/// 2) call search(...) front-end functions, providing the array
///    to be searched
///
/// A list of detections is written in objecList and a 2D map of detections
/// in detectMap.
///
public:
    // Default constructor.
    Search(string searchtype, float snrCut, float threshold, bool adjacent, int threshSpatial,
           int threshVelocity, int minPixels, int minChannels, int minVoxels, int maxChannels,
           float maxAngSize, bool flagGrowth, float growthCut, float growthThreshold, bool RejectBefore, bool TwoStage);
    // Simplified constructor with SEARCH_PAR structure
    Search(SEARCH_PAR &p) : par(p) {}

    virtual ~Search();                                      // Destructor.
    Search(const Search &c);                                // Copy constructor.
    Search& operator=(const Search &c);                     // Copy operator.

    // Inline functions to access private members
    long          getNumObj(){return objectList->size();}
    short*        DetectMap () {return detectMap;}
    short&        DetectMap (size_t npix) {return detectMap[npix];}
    short&        DetectMap (size_t x, size_t y) {return detectMap[x+y*xSize];}

    Detection<T>  getObject(size_t i){return objectList->at(i);}
    Detection<T>* pObject(size_t i){return &(objectList->at(i));}
    DetVec<T>     getObjectList(){return *objectList;}
    DetVec<T>     *pObjectList(){return objectList;}
    DetVec<T>     &ObjectList(){DetVec<T> &rlist=*objectList; return rlist;}

    // Front-end functions to search in an array.
    void search(T *Array, Stats<T> &stat, size_t xsize, size_t ysize=1, size_t zsize=1,
                bool useRobust=true, int nthreads=1, bool Verbose=true, bool Showbar=true);
    void search(T *Array, size_t xsize, size_t ysize=1, size_t zsize=1,
                bool useRobust=true, int nthreads=1, bool Verbose=true, bool Showbar=true);
    // A function to select the largest detection.
    Detection<T>* LargestDetection ();

private:
    SEARCH_PAR par;                           //< Paramters for SEARCH task
    DetVec<T> *objectList = new DetVec<T>;    //< The list of detected objects.
    short     *detectMap;                     //< X,Y locations of detected pixels.
    bool      mapAllocated = false;           //< Has the detection map been allocated?
    T*        array;                          //< A pointer to the array to be searched.
    Stats<T>  stats;                          //< Statistics for array.
    int       xSize = 1;                      //< X-size of array.
    int       ySize = 1;                      //< Y-size of array.
    int       zSize = 1;                      //< Z-size of array.
    bool      verbose = true;                 //< Whether to print output messages.
    bool      showbar = true;                 //< Whether to use a progressbar.
    int       nthreads = 1;                   //< Number of threads for searching.

    // Functions to perform the search
    DetVec<T>   search3DArray();                                     // Switch functions for spectral or spatial.
    DetVec<T>   search3DArraySpectral();                             // Research objects along 1-D spectra.
    DetVec<T>   search3DArraySpatial();                              // Research objects along 2-D channels maps.
    ScanVec<T>  findSources1D(T *spectrum, int minSize);             // Front-end function to spectrumDetect.
    ScanVec<T>  spectrumDetect(vector<bool> &arraybool, int minSize);// Find sources in a 1-D spectrum.
    Obj2DVec<T> findSources2D(T *image, int minSize);                // Front-end function to imageDetect.
    Obj2DVec<T> imageDetect(vector<bool> &arraybool, int minSize);   // Find sources in a 2-D map.
    void    ObjectMerger();                                          // Front-end function to mergeList & finaliseList.
    void    ObjectMergerSimple();                                    // Front-end function to mergeList & finaliseList.
    void    mergeList(DetVec<T> &objList);                           // Merge a list of pixel in a single object.
    void    finaliseList(DetVec<T> &objList);                        // Verify if a detection can be considered an object.
    void    rejectObjects(DetVec<T> &objList);                       // Verify if a detection can be considered an object.
    void    mergeIntoList(Detection<T> &obj, DetVec<T> &objList);    // Add an object in a detection list.
    void    updateDetectMap();                                       // Update the map of detected pixels.
    void    updateDetectMap(Detection<T> obj);                       // Update the map of detected pixels for a Detection.

};


/// Some enumeration types and the FoundObject class

enum STATUS { NONOBJECT,    //< Pixel not above the threshold.
              OBJECT,       //< Pixel above the threshold.
              COMPLETE,     //< Object is complete
              INCOMPLETE    //< Object not yet complete
};  
    
enum ROW {PRIOR=0, CURRENT};
    
enum NULLS { NULLSTART=-1,  //< Default start/end value, obviously outside valid range.
             NULLMARKER=45  //< ASCII 45 = '-', which eases printing for debugging purposes
}; 


enum STATE {AVAILABLE, DETECTED, BLANK, MW};

template <class T>
class FoundObject           // Keeps a track of a detection, as well as the start and finish
{                           // locations of the detection on the current row.
public:
    FoundObject(){start=NULLSTART; end=NULLSTART;}
    int start;              //< Pixel on the current row where the detection starts.
    int end;                //< Pixel on the current row where the detection finishes.
    Object2D<T> info;       //< Collection of detected pixels.
};
//*/
#endif

