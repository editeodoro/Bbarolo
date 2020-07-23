//------------------------------------------------------------------------
//  search.cpp Search-related member functions for the Cube class.
//------------------------------------------------------------------------

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
#include <Tasks/search.hh>
#include <Arrays/stats.hh>
#include <Map/objectgrower.hh>
#include <Map/detection.hh>
#include <Utilities/progressbar.hh>

#ifdef _OPENMP
#include <omp.h>
#endif

template <class T>
Search<T>::Search(std::string searchtype, float snrCut, float threshold, bool adjacent, int threshSpatial,
                  int threshVelocity, int minPixels, int minChannels, int minVoxels, int maxChannels, float maxAngSize,
                  bool flagGrowth, float growthCut, float growthThreshold, bool RejectBefore, bool TwoStage) {

    par.searchType        = searchtype;
    par.snrCut            = snrCut;
    par.threshold         = threshold;
    par.flagAdjacent      = adjacent;
    par.threshSpatial     = threshSpatial;
    par.threshVelocity    = threshVelocity;
    par.minPix            = minPixels;
    par.minChannels       = minChannels;
    par.minVoxels         = minVoxels;
    par.maxChannels       = maxChannels;
    par.maxAngSize        = maxAngSize;
    par.flagGrowth        = flagGrowth;
    par.growthCut         = growthCut;
    par.growthThreshold   = growthThreshold;
    par.RejectBeforeMerge = RejectBefore;
    par.TwoStageMerging   = TwoStage;
    if (par.threshold!=0) par.UserThreshold = true;
    if (par.growthThreshold!=0) par.flagUserGrowthT = true;
}


template <class T>
Search<T>::Search(const Search<T> &s) {

    this->operator=(s);

}


template <class T>
Search<T>& Search<T>::operator=(const Search<T> &s) {

    if(this==&s) return *this;

    this->par       = s.par;
    this->array     = s.array;
    this->xSize     = s.xSize;
    this->ySize     = s.ySize;
    this->zSize     = s.zSize;
    this->stats     = s.stats;
    this->verbose   = s.verbose;
    this->showbar   = s.showbar;
    this->nthreads  = s.nthreads;

    std::cout << "QUI !" << std::endl;
    this->mapAllocated = s.mapAllocated;
    if(this->mapAllocated) {
        this->detectMap = new short[this->xSize*this->ySize];
        for(int i=0;i<(this->xSize*this->ySize);i++) this->detectMap[i] = s.detectMap[i];
    }

    this->objectList = new DetVec<T>(s.objectList->size());
    for (int i=0; i<s.objectList->size(); i++)
        this->objectList->at(i) = s.objectList->at(i);

    return *this;
}


template <class T>
Search<T>::~Search(){
    if (mapAllocated) delete [] detectMap;
    delete objectList;
}



template <class T>
void Search<T>::search(T *Array, Stats<T> &stat, size_t xsize, size_t ysize, size_t zsize, bool useRobust, int Nthreads, bool Verbose, bool Showbar) {

    xSize = xsize;
    ySize = ysize;
    zSize = zsize;
    stats = stat;
    array = Array;
    verbose = Verbose;
    showbar = Showbar;
    nthreads = Nthreads;

    if (mapAllocated) delete [] detectMap;
    detectMap = new short[xsize*ysize];
    mapAllocated = true;

    if (par.UserThreshold) stats.setThreshold(par.threshold);

    if(verbose)
        std::cout << "\n\nStarting research for possible sources in the cube... " << std::endl;

    *objectList = search3DArray();

    if(verbose) std::cout << "  Updating detection map... " << std::flush;
    updateDetectMap();
    if(verbose) std::cout << "Done.\n";

    if (verbose) {
        std::cout << "  Intermediate list has " << getNumObj();
        if(getNumObj()==1) std::cout << " object.\n";
        else std::cout << " objects.\n";
    }

    if(getNumObj()>0){
        if (verbose) std::cout << "  Merging and Rejecting...  "<< std::flush;
        ObjectMerger();
        if (verbose) std::cout << "Done.                      " << std::endl;
     }
     if (verbose) std::cout << "  ... All done.\n\nFinal object count = " << getNumObj() << std::endl << std::endl;

}


template <class T>
void Search<T>::search(T *Array, size_t xsize, size_t ysize, size_t zsize, bool useRobust, int Nthreads, bool Verbose, bool Showbar) {

    // Computing statistics for the input array
    Stats<T> stat;
    stat.setRobust(useRobust);
    size_t numPix = xsize*ysize*zsize;
    bool *blanks = new bool[numPix];
    for (size_t i=0; i<numPix; i++) blanks[i] = isBlank(Array[i]) ? false : true;
    stat.calculate(array,numPix,blanks);
    stat.setThresholdSNR(par.snrCut);
    delete [] blanks;

    search(Array,stat,xsize,ysize,zsize,useRobust,Nthreads,Verbose,Showbar);
}


template <class T>
DetVec<T> Search<T>::search3DArray() {

    /// A simple front end function to choose spatial
    /// or spectral research. If [searchType] parameter
    /// is not "spatial" or "spectral", return an empty
    /// list of object detected.

    std::string stype = par.searchType;
    if (stype=="spatial" || zSize==1)
        return search3DArraySpatial();
    else if(stype=="spectral" )
        return search3DArraySpectral();
    else {
        std::cout << "Unknown search type : " << stype << std::endl;
        return DetVec<T>(0);
    }
}


template <class T>
DetVec<T> Search<T>::search3DArraySpectral() {

  /// The function searches for detections in just
  /// the 1D spectra.
  /// Returns a vector list of Detections.
  ///
  /// \return       Vector of detected objects.

    DetVec<T> outputList;
    int num = 0;

    if(zSize>1){

        ProgressBar bar(false,verbose,showbar);
        bar.init("Searching in progress... ",ySize);

        T *spectrum = new T[zSize];

        for(int y=0; y<ySize; y++){
            bar.update(y+1);
            for(int x=0; x<xSize; x++){
                int npix = y*xSize + x;
                for(int z=0; z<zSize; z++) spectrum[z] = array[z*xSize*ySize+npix];

                ScanVec<T> objlist = findSources1D(spectrum,1);

                typename ScanVec<T>::iterator obj;
                num += objlist.size();
                for(obj=objlist.begin();obj<objlist.end();obj++){
                    Detection<T> newObject;
                    for(int z=obj->getX();z<=obj->getXmax();z++){
                        newObject.addPixel(x,y,z);
                    }
                    newObject.setOffsets();
                    if(par.TwoStageMerging) mergeIntoList(newObject,outputList);
                    else outputList.push_back(newObject);
                }
            }
        }

        delete [] spectrum;

        bar.fillSpace("Found "+to_string(num)+" items.\n");
    }

    return outputList;
}


template <class T>
DetVec<T> Search<T>::search3DArraySpatial() {

  ///  The function searches for detections just in the channel maps.
  ///  Returns a vector list of Detections.
  ///
  ///  \return A std::vector of detected objects.

    DetVec<T> outputList;
    int num = 0;

    bool useBar = (zSize>1);

    ProgressBar bar(false,verbose,useBar&&showbar);

#pragma omp parallel num_threads(nthreads)
{
    bar.init("Searching in progress... ",zSize);
#pragma omp for schedule(dynamic) reduction (+:num)
    for(int z=0; z<zSize; z++) {
        bar.update(z+1);

        Obj2DVec<T> objlist = findSources2D(&array[z*xSize*ySize],1);
        typename Obj2DVec<T>::iterator obj;
        num += objlist.size();

        for(obj=objlist.begin();obj!=objlist.end();obj++){
            Detection<T> newObject;
            newObject.addChannel(z,*obj);
            newObject.setOffsets();
#pragma omp critical
{
            if(par.TwoStageMerging) mergeIntoList(newObject,outputList);
            else outputList.push_back(newObject);
}
        }
    }
}

    bar.fillSpace("Found "+to_string(num)+" items.\n");

    return outputList;

}


template <class T>
ScanVec<T> Search<T>::findSources1D(T *spectrum, int minSize) {

    std::vector<bool> thresholdedArray(zSize);
    for(int z=0; z<zSize; z++)
        thresholdedArray[z] = stats.isDetection(spectrum[z]);
    return spectrumDetect(thresholdedArray,minSize);
}


template <class T>
Obj2DVec<T> Search<T>::findSources2D(T *image, int minSize) {

    std::vector<bool> thresholdedArray(xSize*ySize);
    for (int i=0; i<xSize*ySize; i++)
        thresholdedArray[i] = stats.isDetection(image[i]);
    return imageDetect(thresholdedArray,minSize);
}


template <class T>
ScanVec<T> Search<T>::spectrumDetect(std::vector<bool> &arraybool, int minSize) {

  ///  A detection algorithm that searches in a single 1-D spectrum. It
  ///  simply scans along the spectrum, storing connected sets of
  ///  detected pixels as Scans, where "detected" means according to the
  ///  Image::isDetection(long,long) function.
  ///
  ///  When finished a vector of the detected scans is returned.

    enum STATUS {NONOBJECT, OBJECT};
    STATUS status;
    PixelInfo::Scan<T> obj;
    ScanVec<T> outputlist;
    bool isObject;
    long dim = zSize;
    status = NONOBJECT;
    for(int pos=0;pos<(dim+1);pos++){
        if(pos<dim){
            isObject = arraybool[pos];
        }
        else isObject=false;

        if(isObject){
            if(status != OBJECT){
                status = OBJECT;
                obj.define(0, pos, 1);
            }
            else obj.growRight();
        }
        else {
            if(status == OBJECT){
                if(obj.getXlen() >= int(minSize)){
                    outputlist.push_back(obj);
                }
            obj.clear();
            }
            status = NONOBJECT;
        }

    }
    return outputlist;
}


template <class T>
Obj2DVec<T> Search<T>::imageDetect(std::vector<bool> &arraybool, int minSize) {

    ///  A detection algorithm for 2-dimensional images based on that of
    ///  Lutz (1980).
    ///
    ///  The image is raster-scanned, and searched row-by-row. Objects
    ///  detected in each row are compared to objects in subsequent rows,
    ///  and combined if they are connected (in an 8-fold sense).
    ///

    int xdim = xSize;
    int ydim = ySize;
    Obj2DVec<T> outputlist;
    STATUS *status  = new STATUS[2];
    Object2D<T> *store = new Object2D<T>[xdim+1];
    char *marker    = new char[xdim+1];
    for(int i=0; i<(xdim+1); i++) marker[i] = NULLMARKER;
    std::vector<FoundObject<T> > oS;
    std::vector<STATUS>      psS;

    Pixel<T> pix;
    size_t loc=0;

    for (long posY=0; posY<(ydim+1); posY++) {          // Loop over each row.
        // Consider rows one at a time.
        status[PRIOR] = COMPLETE;
        status[CURRENT] = NONOBJECT;

        for(long posX=0; posX<(xdim+1); posX++) {       // Now the loop for a given row,
            // looking at each column individually.
            char currentMarker = marker[posX];
            marker[posX] = NULLMARKER;
            bool isObject;
            if ((posX<xdim) && (posY<ydim))             // If we are in the original image.
                isObject = arraybool[loc++];            // Else we're in the padding row/col and
            else isObject = false;                      // isObject=FALSE;

            // -------------------------START SEGMENT ------------------------------
            // If the current pixel is object and the previous pixel is not, then
            // start a new segment.
            // If the pixel touches an object on the prior row, the marker is either
            // an S or an s, depending on whether the object has started yet.
            // If it doesn't touch a prior object, this is the start of a completly
            // new object on this row.

            if ((isObject) && (status[CURRENT] != OBJECT)) {

                status[CURRENT]=OBJECT;
                if(status[PRIOR]==OBJECT) {
                    if(oS.back().start==NULLSTART){
                        marker[posX] = 'S';
                        oS.back().start = posX;
                    }
                    else  marker[posX] = 's';
                }
                else {
                    psS.push_back(status[PRIOR]);
                    marker[posX] = 'S';
                    oS.resize(oS.size()+1);
                    oS.back().start = posX;
                    status[PRIOR] = COMPLETE;
                }
            }

            // ------------------------ PROCESS MARKER -----------------------------
            // If the current marker is not blank, then we need to deal with it.
            // Four cases:
            //   S --> start of object on prior row. Push priorStatus onto PSSTACK
            //         and set priorStatus to OBJECT
            //   s --> start of a secondary segment of object on prior row.
            //         If current object joined, pop PSSTACK and join the objects.
            //         Set priorStatus to OBJECT.
            //   f --> end of a secondary segment of object on prior row.
            //         Set priorStatus to INCOMPLETE.
            //   F --> end of object on prior row. If no more of the object is to
            //         come (priorStatus=COMPLETE), then finish it and output data.
            //         Add to list, but only if it has more than the minimum number
            //         of pixels.

            if(currentMarker != NULLMARKER) {
                if(currentMarker == 'S') {
                    psS.push_back(status[PRIOR]);
                    if(status[CURRENT] == NONOBJECT) {
                        psS.push_back(COMPLETE);
                        oS.resize(oS.size()+1);
                        oS.back().info = store[posX];
                    }
                    else oS.back().info = oS.back().info + store[posX];
                    status[PRIOR] = OBJECT;
                }

                if(currentMarker == 's'){
                    if((status[CURRENT]==OBJECT) && (status[PRIOR]==COMPLETE)){
                        status[PRIOR] = psS.back();
                        psS.pop_back();
                        oS[oS.size()-2].info = oS[oS.size()-2].info + oS.back().info;
                        if(oS[oS.size()-2].start == NULLSTART)
                            oS[oS.size()-2].start = oS.back().start;
                        else marker[oS.back().start] = 's';
                        oS.pop_back();
                    }
                    status[PRIOR] = OBJECT;
                }

                if(currentMarker == 'f') status[PRIOR] = INCOMPLETE;

                if(currentMarker == 'F') {
                    status[PRIOR] = psS.back();
                    psS.pop_back();
                    if((status[CURRENT]==NONOBJECT) && (status[PRIOR]==COMPLETE)){
                        if(oS.back().start == NULLSTART){               // The object is completed.
                            if(oS.back().info.getSize() >= minSize)     // If it is big enough,
                                outputlist.push_back(oS.back().info);   // add to the end of the
                        }                                               // output list.
                        else {
                            marker[oS.back().end] = 'F';
                            store[oS.back().start] = oS.back().info;
                        }

                        oS.pop_back();
                        status[PRIOR] = psS.back();
                        psS.pop_back();
                    }
                }

            }

            if (isObject) {
                oS.back().info.addPixel(posX,posY);
            }
            else {

                // ----------------------------- END SEGMENT -------------------------
                // If the current pixel is background and the previous pixel was an
                // object, then finish the segment.
                // If the prior status is COMPLETE, it's the end of the final segment
                // of the object section.
                // If not, it's end of the segment, but not necessarily the section.

                if (status[CURRENT]==OBJECT) {
                    status[CURRENT] = NONOBJECT;
                    if(status[PRIOR] != COMPLETE){
                        marker[posX] = 'f';
                        oS.back().end = posX;
                    }
                    else {
                        status[PRIOR] = psS.back();
                        psS.pop_back();
                        marker[posX] = 'F';
                        store[oS.back().start] = oS.back().info;
                        oS.pop_back();
                    }
                }
            }
        }
    }


    delete [] marker;
    delete [] store;
    delete [] status;

    return outputlist;

}


template <class T>
void Search<T>::ObjectMerger() {

    /// A Function that takes a Cube's list of Detections and
    /// combines those that are close
    /// It also excludes those that do not make the minimum
    /// number of channels requirement.
    /// A front end to simpler functions mergeList and finaliseList,
    /// with code to cover the option of growing objects.

    int startSize = objectList->size();

    if(startSize>0){

        DetVec<T> currentList(startSize);
        for(int i=0;i<startSize;i++) currentList[i] = objectList->at(i);
        objectList->clear();

        if(par.RejectBeforeMerge) finaliseList(currentList);

        mergeList(currentList);

        // Do growth stuff
        if(par.flagGrowth) {
            ObjectGrower<T> grower;
            grower.define(stats,array,xSize,ySize,zSize,objectList,par);
#pragma omp parallel for num_threads(nthreads)
            for(size_t i=0;i<currentList.size();i++){
                if(verbose && showbar){
#ifdef _OPENMP
                    int tid = omp_get_thread_num();
                    if(tid==0) {
#endif
                        std::cout.setf(std::ios::right);
                        std::cout << "Growing: " << std::setw(6) << i+1 << "/";
                        std::cout.unsetf(std::ios::right);
                        std::cout.setf(std::ios::left);
                        std::cout << std::setw(6) << currentList.size() << std::flush;
                        printBackSpace(22);
                        std::cout << std::flush;
#ifdef _OPENMP
                    }
#endif
                }
                grower.grow(&currentList[i]);
            }
            grower.updateDetectMap(detectMap);
            std::cout.unsetf(std::ios::left);
            mergeList(currentList);
        }

        if(!par.RejectBeforeMerge) finaliseList(currentList);
        if(par.maxChannels!=-1 || par.maxAngSize!=-1) rejectObjects(currentList);

        objectList->resize(currentList.size());
        for(size_t i=0;i<currentList.size();i++)
            objectList->at(i) = currentList[i];

      currentList.clear();

    }
}


template <class T>
void Search<T>::ObjectMergerSimple() {

    ///   A simple front-end to the mergeList() and finaliseList() functions,
    ///    so that if you want to merge a single list, it will
    ///    do both the merging and the cleaning up afterwards.

    mergeList(*objectList);
    finaliseList(*objectList);
}


template <class T>
void Search<T>::mergeList(DetVec<T> &objList) {

    /// A function that merges any objects in the list of
    /// Detections that are within stated threshold distances.
    /// Determination of whether objects are close is done by
    /// the function areClose.

    if(objList.size() > 0){

        typename DetVec<T>::iterator iter;
        std::vector<bool> isValid(objList.size(),true);
        int numRemoved=0;
        size_t counter=0, compCounter,goodCounter=0;

        while(counter < (objList.size()-1)){
            if(verbose && showbar){
                std::cout.setf(std::ios::right);
                std::cout << "Merging objects: ";
                std::cout << std::setw(6) << goodCounter+1 << "/" ;
                std::cout.unsetf(std::ios::right);
                std::cout.setf(std::ios::left);
                std::cout << std::setw(6) << objList.size()-numRemoved;
                printBackSpace(30);
                std::cout << std::flush;
                std::cout.unsetf(std::ios::left);
            }
            if(isValid[counter]){
                compCounter = counter + 1;
                do {
                    if(isValid[compCounter]){
                        bool close = objList[counter].canMerge(objList[compCounter], par);
                        if(close){
                            objList[counter].addDetection(objList[compCounter]);
                            isValid[compCounter]=false;
                            numRemoved++;
                            if(verbose && showbar){
                                std::cout.setf(std::ios::right);
                                std::cout << "Merging objects: ";
                                std::cout << std::setw(6) << goodCounter+1 << "/";
                                std::cout.unsetf(std::ios::right);
                                std::cout.setf(std::ios::left);
                                std::cout << std::setw(6) << objList.size()-numRemoved;
                                printBackSpace(30);
                                std::cout << std::flush;
                                std::cout.unsetf(std::ios::left);

                            }
                            compCounter = counter + 1;
                        }
                        else compCounter++;
                    }
                    else compCounter++;

                } while( (compCounter<objList.size()) );
            }
            counter++;
            if(isValid[counter]) goodCounter++;
        }

        DetVec<T> newlist(objList.size()-numRemoved);
        size_t ct=0;
        for(size_t i=0;i<objList.size();i++){
            if(isValid[i]) newlist[ct++]=objList[i];
        }
        objList.clear();
        objList=newlist;

    }
}


template <class T>
void Search<T>::mergeIntoList(Detection<T> &object, DetVec<T> &objList) {

    /// A function to add a detection to a list of detections, checking
    /// first to see if it can be combined with existing members of the
    /// list.
    ///
    /// The areClose testing and combining is done with the
    /// parameters as given by the Param set.
    ///
    /// \param object The Detection to be merged into the list.
    /// \param objList The vector list of Detections.
    /// \param par The Param set, used for testing if merging needs to be done.

    bool haveMerged = false;
    typename DetVec<T>::iterator iter;

    for(iter=objList.begin(); (!haveMerged && iter<objList.end()); iter++) {
        if(iter->canMerge(object, par)){
            iter->addDetection(object);
            haveMerged = true;
        }
    }

    if(!haveMerged) objList.push_back(object);

}


template <class T>
void Search<T>::rejectObjects(DetVec<T> &objList) {

    if(verbose && showbar){
        std::cout << "Rejecting:" << std::setw(6) << objList.size();
        printSpace(6);
        printBackSpace(22);
        std::cout << std::flush;
    }

    int maxchan = par.maxChannels;
    float maxsize = par.maxAngSize;

    DetVec<T> newlist;

    typename DetVec<T>::iterator obj = objList.begin();
    int numRej=0;
    for(;obj<objList.end();obj++){
        obj->setOffsets();

        int nchan = obj->getNumChannels();
        float mSize = std::max(fabs(obj->getXmax()-obj->getXmin()),fabs(obj->getYmax()-obj->getYmin()));

        if((maxchan>0 && maxchan<nchan) || (maxsize>0 && maxsize<mSize)) {
            numRej++;
            if(verbose && showbar){
                std::cout << "Rejecting:" << std::setw(6) << objList.size()-numRej;
                printSpace(6);
                printBackSpace(22);
                std::cout << std::flush;
            }
        }
        else newlist.push_back(*obj);
    }

    objList.clear();
    objList = newlist;

}



template <class T>
void Search<T>::finaliseList(DetVec<T> &objList) {

    ///  A function that looks at each object in the Detection vector
    ///  and determines whether is passes the requirements for the
    ///  minimum number of channels and spatial pixels, as provided by
    ///  the Param set par.
    ///  If it does not pass, it is removed from the list.
    ///  In the process, the offsets are set.

    if(verbose && showbar){
        std::cout << "Rejecting:" << std::setw(6) << objList.size();
        printSpace(6);
        printBackSpace(22);
        std::cout << std::flush;
    }

    DetVec<T> newlist;

    typename DetVec<T>::iterator obj = objList.begin();
    int numRej=0;
    for(;obj<objList.end();obj++){
        obj->setOffsets();

        if((obj->hasEnoughChannels(par.minChannels)
            && (int(obj->getSpatialSize()) >= par.minPix)
            && (int(obj->getSize()) >= par.minVoxels))){

                newlist.push_back(*obj);

        }
        else{
            numRej++;
            if(verbose && showbar){
                std::cout << "Rejecting:" << std::setw(6) << objList.size()-numRej;
                printSpace(6);
                printBackSpace(22);
                std::cout << std::flush;
            }

        }
    }

    objList.clear();
    objList = newlist;

}


template <class T>
void Search<T>::updateDetectMap() {

    ///  A function that, for each detected object in the
    ///  cube's list, increments the cube's detection map by the
    ///  required amount at each pixel.

    typename DetVec<T>::iterator obj;
    for(obj=objectList->begin();obj<objectList->end();obj++)
        updateDetectMap(*obj);

}


template <class T>
void Search<T>::updateDetectMap(Detection<T> obj) {

    ///  A function that, for the given object, increments the cube's
    ///  detection map by the required amount at each pixel.
    ///
    ///  \param obj     A Detection object that is being
    ///                 incorporated into the map.

    std::vector<Voxel<T> > vlist = obj.getPixelSet();
    typename std::vector<Voxel<T> >::iterator vox;
    for(vox=vlist.begin();vox<vlist.end();vox++)
        detectMap[vox->getX()+vox->getY()*xSize]++;

}


template <class T>
Detection<T>* Search<T>::LargestDetection () {

    int numObj = objectList->size();
    if (numObj==0) return NULL;
    uint n=0, size=0;
    for (int i=0; i<numObj; i++) {
        Detection<T> *obj = pObject(i);
        if (obj->getSize()>size) {n=i;size=obj->getSize();}
    }
    return pObject(n);
}


// Explicit instantiation of the class
template class Search<short>;
template class Search<int>;
template class Search<long>;
template class Search<float>;
template class Search<double>;


