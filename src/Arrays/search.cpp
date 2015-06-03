//------------------------------------------------------------------------
//	search.cpp Search-related member functions for the Cube class.
//------------------------------------------------------------------------

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
 along with Bbarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning Bbarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <climits>
#include <cfloat>
#include <iomanip>
#include <algorithm>
#include "cube.hh"
#include "image.hh"
#include "../Map/objectgrower.hh"
#include "../Map/detection.hh"
#include "../Utilities/gnuplot.hh"
#include "../Utilities/progressbar.hh"
#include "../Utilities/utils.hh"



template <class T>
void Cube<T>::Search() {
	
	Header &h = Head();	
	Param  &p = pars();	 
	if (!statsDefined) setCubeStats();
	if (h.BeamArea()==0) {
		cout << "\n Beam information is not available in the header: assuming a "
			 << p.getBeamFWHM()*3600 << " arcsec beam.\n";
		h.setBmaj(p.getBeamFWHM());
		h.setBmin(p.getBeamFWHM());
		h.calcArea();
	}
	float PixScale = (fabs(h.Cdelt(0))+fabs(h.Cdelt(1)))/2.;
	int thresS  = p.getThreshS()!=-1 	 ? p.getThreshS() 	  : ceil(h.Bmaj()/PixScale);
	int thresV  = p.getThreshV()!=-1 	 ? p.getThreshV()     : 3;
	int minchan = p.getMinChannels()!=-1 ? p.getMinChannels() : 2;
	int minpix  = p.getMinPix()!=-1 	 ? p.getMinPix() 	  : ceil(h.BeamArea());
	int minvox  = p.getMinVoxels()!=-1   ? p.getMinVoxels()   : minchan*minpix;
	p.setThreshS(thresS);
	p.setThreshV(thresV);
	p.setMinChannels(minchan);
	p.setMinPix(minpix);
	p.setMinVoxels(minvox);			
	
	CubicSearch();
	
    std::cout << "  Intermediate list has " << getNumObj();
    if(getNumObj()==1) std::cout << " object.\n";
    else std::cout << " objects.\n";

    if(getNumObj() > 0){
		std::cout << "  Merging and Rejecting...  "<< std::flush;
		ObjectMerger();	
		std::cout << "Done.                      " << std::endl;
    }
    std::cout << "  ... All done.\n\nFinal object count = " << getNumObj() << std::endl;
    std::cout << std::endl; 

    isSearched = true;
}
template void Cube<short>::Search();
template void Cube<int>::Search();
template void Cube<long>::Search();
template void Cube<float>::Search();
template void Cube<double>::Search();


template <class T>
void Cube<T>::CubicSearch() {

  ///  A front end to the cubic searching routine that does not
  ///  involve any wavelet reconstruction. 
  ///  The statistics of the cube are calculated first of all.
  ///  If baseline-removal is required that is done prior to searching.
  ///  Once searching is complete, the detection map is updated and
  ///  the intermediate detections are logged in the log file.

	
	std::cout<<"\n\nStarting research for possible sources in the cube..."<<std::endl;
	if(par.isVerbose()) std::cout << "  ";
	
	if (par.getFlagUserThreshold()) stats.setThreshold(par.getThreshold());
	else { 
		if (!statsDefined) setCubeStats(); 
	}
	
	*this->objectList = search3DArray();

	if(par.isVerbose()) std::cout << "  Updating detection map... " << std::flush;
 
	updateDetectMap();
	if(par.isVerbose()) std::cout << "Done.\n";
	
}
template void Cube<short>::CubicSearch();
template void Cube<int>::CubicSearch();
template void Cube<long>::CubicSearch();
template void Cube<float>::CubicSearch();
template void Cube<double>::CubicSearch();


template <class T> 
std::vector <Detection<T> > Cube<T>::search3DArray() {
	
  /// A simple front end functions to choose spatial
  /// or spectral research. If [searchType] parameter
  /// is not "spatial" or "spectral", return an empty
  /// list of object detected.

	if(par.getSearchType()=="spectral")
		return search3DArraySpectral();
	else if(par.getSearchType()=="spatial")
		return search3DArraySpatial();
	else {
		std::cout << "Unknown search type : " << par.getSearchType()<< std::endl;
    return std::vector<Detection<T> >(0);
	}
}
template std::vector <Detection<short> > Cube<short>::search3DArray();
template std::vector <Detection<int> > Cube<int>::search3DArray();
template std::vector <Detection<long> > Cube<long>::search3DArray();
template std::vector <Detection<float> > Cube<float>::search3DArray();
template std::vector <Detection<double> > Cube<double>::search3DArray();


template <class T>
std::vector <Detection<T> > Cube<T>::search3DArraySpectral() {
	
  /// The function searches for detections in just 
  /// the 1D spectra.
  /// Returns a vector list of Detections.
  ///
  /// \return 		Vector of detected objects.

	std::vector <Detection<T> > outputList;
	long zdim = axisDim[2];
	long xySize = axisDim[0] * axisDim[1];
	int num = 0;

	if(zdim>1){
		
        ProgressBar bar("Searching in progress... ");
        bar.setShowbar(par.getShowbar());
        if(par.isVerbose()) bar.init(xySize);
		
		int *specdim = new int[2];
		specdim[0] = zdim; specdim[1]=1;
		Image2D<T> *spectrum = new Image2D<T> (specdim);
		delete [] specdim;
		spectrum->saveStats(stats);
		spectrum->saveParam(par);
		spectrum->setMinSize(1);

		for(int y=0; y<axisDim[1]; y++){
			for(int x=0; x<axisDim[0]; x++){

				int npix = y*axisDim[0] + x;
                if(par.isVerbose()) bar.update(npix+1);
				spectrum->extractSpectrum(array,axisDim,npix);
				std::vector<Scan<T> > objlist = spectrum->findSources1D();
				typename std::vector<Scan<T> >::iterator obj;
				num += objlist.size();
				
				for(obj=objlist.begin();obj<objlist.end();obj++){
					Detection<T> newObject;
					for(int z=obj->getX();z<=obj->getXmax();z++) {
						newObject.addPixel(x,y,z);
					}
					newObject.setOffsets();
					if(par.getTwoStageMerging()) mergeIntoList(newObject,outputList);
					else outputList.push_back(newObject);
				}
			}
		}

		delete spectrum;
  
        if(par.isVerbose()){
			bar.remove();
			std::cout << "Found " << num << " items.\n";
		}
		
	}
	
	return outputList;
}
template std::vector <Detection<short> > Cube<short>::search3DArraySpectral();
template std::vector <Detection<int> > Cube<int>::search3DArraySpectral();
template std::vector <Detection<long> > Cube<long>::search3DArraySpectral();
template std::vector <Detection<float> > Cube<float>::search3DArraySpectral();
template std::vector <Detection<double> > Cube<double>::search3DArraySpectral();


template <class T>
std::vector <Detection<T> > Cube<T>::search3DArraySpatial() {
	
  ///  The function searches for detections just in the channel maps.
  ///  Returns a vector list of Detections.
  ///
  ///  \return A std::vector of detected objects.

	std::vector <Detection<T> > outputList;
	long zdim = axisDim[2];
	int num = 0;

    ProgressBar bar("Searching in progress... ");
    bool useBar = (zdim>1);
    bar.setShowbar(par.getShowbar());
    if(useBar && par.isVerbose()) bar.init(zdim);

	int *imdim = new int[2];
	imdim[0] = axisDim[0]; imdim[1] = axisDim[1];
	Image2D<T> *channelImage = new Image2D<T>(imdim);
	delete [] imdim;
	channelImage->saveParam(par);
	channelImage->saveStats(stats);
	channelImage->setMinSize(1);

	for(int z=0; z<zdim; z++){
		
        if(par.isVerbose() && useBar) bar.update(z+1);
		channelImage->extractImage(array,axisDim,z);
		std::vector<Object2D<T> > objlist = channelImage->findSources2D();
		typename std::vector<Object2D<T> >::iterator obj;
		num += objlist.size();
		for(obj=objlist.begin();obj!=objlist.end();obj++){
			Detection<T> newObject;
			newObject.addChannel(z,*obj);
			newObject.setOffsets();
			if(par.getTwoStageMerging()) mergeIntoList(newObject,outputList);
			else outputList.push_back(newObject);
		}
	}

	delete channelImage;

	if(par.isVerbose()){
        if(useBar) bar.remove();
        std::cout << "Found " << num << " items.\n";
	}

	return outputList;

}
template std::vector <Detection<short> > Cube<short>::search3DArraySpatial();
template std::vector <Detection<int> > Cube<int>::search3DArraySpatial();
template std::vector <Detection<long> > Cube<long>::search3DArraySpatial();
template std::vector <Detection<float> > Cube<float>::search3DArraySpatial();
template std::vector <Detection<double> > Cube<double>::search3DArraySpatial();


template <class T>
void Cube<T>::updateDetectMap() {
	
  ///  A function that, for each detected object in the
  ///  cube's list, increments the cube's detection map by the
  ///  required amount at each pixel.

    typename std::vector<Detection<T> >::iterator obj;
    for(obj=objectList->begin();obj<objectList->end();obj++)
		updateDetectMap(*obj);
    

}
template void Cube<short>::updateDetectMap();
template void Cube<int>::updateDetectMap();
template void Cube<long>::updateDetectMap();
template void Cube<float>::updateDetectMap();
template void Cube<double>::updateDetectMap();


template <class T>
void Cube<T>::updateDetectMap(Detection<T> obj) {
    
    ///  A function that, for the given object, increments the cube's
    ///  detection map by the required amount at each pixel.
    /// 
    ///  \param obj 	A Detection object that is being 
    ///					incorporated into the map.

    std::vector<Voxel<T> > vlist = obj.getPixelSet();
    typename std::vector<Voxel<T> >::iterator vox;
    for(vox=vlist.begin();vox<vlist.end();vox++) 
		detectMap[vox->getX()+vox->getY()*axisDim[0]]++;

}
template void Cube<short>::updateDetectMap(Detection<short>);
template void Cube<int>::updateDetectMap(Detection<int>);
template void Cube<long>::updateDetectMap(Detection<long>);
template void Cube<float>::updateDetectMap(Detection<float>);
template void Cube<double>::updateDetectMap(Detection<double>);


template <class T>
void Cube<T>::ObjectMerger() {

    /// A Function that takes a Cube's list of Detections and
    /// combines those that are close 
    /// It also excludes those that do not make the minimum
    /// number of channels requirement.
    /// A front end to simpler functions mergeList and finaliseList,
    /// with code to cover the option of growing objects.

    int startSize = objectList->size();

    if(startSize > 0){
		
		std::vector <Detection<T> > currentList(startSize);
		for(int i=0;i<startSize;i++) currentList[i] = objectList->at(i);
		objectList->clear();

        if(par.getRejectBeforeMerge()) finaliseList(currentList);

		mergeList(currentList);

		// Do growth stuff
		if(par.getFlagGrowth()) {
			ObjectGrower<T> grower;
			grower.define(this);
			for(size_t i=0;i<currentList.size();i++){
                if(par.isVerbose() && par.getShowbar()){
					std::cout.setf(std::ios::right);
                    std::cout << "Growing: " << std::setw(6) << i+1 << "/";
                    std::cout.unsetf(std::ios::right);
                    std::cout.setf(std::ios::left);
                    std::cout << std::setw(6) << currentList.size() << std::flush;
                    printBackSpace(22);
                    std::cout << std::flush;

				}
				grower.grow(&currentList[i]);
			}
			grower.updateDetectMap(detectMap);
			std::cout.unsetf(std::ios::left);
			mergeList(currentList);
		}	

		if(!par.getRejectBeforeMerge()) finaliseList(currentList);
        if(par.getMaxChannels()!=-1 || par.getMaxAngSize()!=-1) rejectObjects(currentList);

		objectList->resize(currentList.size());
		for(size_t i=0;i<currentList.size();i++)
			objectList->at(i) = currentList[i];
    
      currentList.clear();

    }
}
template void Cube<short>::ObjectMerger();
template void Cube<int>::ObjectMerger();
template void Cube<long>::ObjectMerger();
template void Cube<float>::ObjectMerger();
template void Cube<double>::ObjectMerger();


template <class T>
void Cube<T>::ObjectMergerSimple() {
    
    ///   A simple front-end to the mergeList() and finaliseList() functions,
    ///    so that if you want to merge a single list, it will
    ///    do both the merging and the cleaning up afterwards.

    mergeList(*objectList);
    finaliseList(*objectList);
}
template void Cube<short>::ObjectMergerSimple();
template void Cube<int>::ObjectMergerSimple();
template void Cube<long>::ObjectMergerSimple();
template void Cube<float>::ObjectMergerSimple();
template void Cube<double>::ObjectMergerSimple();


template <class T>
void Cube<T>::mergeList(std::vector<Detection<T> > &objList) {

    /// A function that merges any objects in the list of 
    /// Detections that are within stated threshold distances.
    /// Determination of whether objects are close is done by
    /// the function areClose. 

    if(objList.size() > 0){      
      
      bool isVerb = par.isVerbose();
      typename std::vector <Detection<T> >::iterator iter;
      std::vector<bool> isValid(objList.size(),true);
      int numRemoved=0;
      size_t counter=0, compCounter,goodCounter=0;

      while(counter < (objList.size()-1)){
          if(isVerb && par.getShowbar()){
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
                        if(isVerb && par.getShowbar()){
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

    std::vector<Detection<T> > newlist(objList.size()-numRemoved);
    size_t ct=0;
    for(size_t i=0;i<objList.size();i++){
		if(isValid[i]) newlist[ct++]=objList[i];
    }
    objList.clear();
    objList=newlist;

    }
}
template void Cube<short>::mergeList(std::vector<Detection<short> >&);
template void Cube<int>::mergeList(std::vector<Detection<int> >&);
template void Cube<long>::mergeList(std::vector<Detection<long> >&);
template void Cube<float>::mergeList(std::vector<Detection<float> >&);
template void Cube<double>::mergeList(std::vector<Detection<double> >&);


template <class T>
void Cube<T>::finaliseList(std::vector<Detection<T> > &objList) {

    ///  A function that looks at each object in the Detection vector
    ///  and determines whether is passes the requirements for the
    ///  minimum number of channels and spatial pixels, as provided by
    ///  the Param set par.
    ///  If it does not pass, it is removed from the list.
    ///  In the process, the offsets are set.

    if(par.isVerbose() && par.getShowbar()){
        std::cout << "Rejecting:" << std::setw(6) << objList.size();
        printSpace(6);
        printBackSpace(22);
        std::cout << std::flush;
    }

    std::vector<Detection<T> > newlist;
    
    typename std::vector<Detection<T> >::iterator obj = objList.begin();
    int numRej=0;
    for(;obj<objList.end();obj++){
		obj->setOffsets();
      
		if((obj->hasEnoughChannels(par.getMinChannels())
			&& (int(obj->getSpatialSize()) >= par.getMinPix())
			&& (int(obj->getSize()) >= par.getMinVoxels()))){

				newlist.push_back(*obj);

		}      
		else{
			numRej++;
            if(par.isVerbose() && par.getShowbar()){
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
template void Cube<short>::finaliseList(std::vector<Detection<short> >&);
template void Cube<int>::finaliseList(std::vector<Detection<int> >&);
template void Cube<long>::finaliseList(std::vector<Detection<long> >&);
template void Cube<float>::finaliseList(std::vector<Detection<float> >&);
template void Cube<double>::finaliseList(std::vector<Detection<double> >&);



template <class T>
void Cube<T>::rejectObjects(std::vector<Detection<T> > &objList) {

    if(par.isVerbose() && par.getShowbar()){
        std::cout << "Rejecting:" << std::setw(6) << objList.size();
        printSpace(6);
        printBackSpace(22);
        std::cout << std::flush;
    }

    int maxchan = par.getMaxChannels();
    float maxsize = par.getMaxAngSize()/(head.PixScale()*arcsconv(head.Cunit(0))/60.);

    std::vector<Detection<T> > newlist;

    typename std::vector<Detection<T> >::iterator obj = objList.begin();
    int numRej=0;
    for(;obj<objList.end();obj++){
        obj->setOffsets();

        int nchan = obj->getNumChannels();
        float mSize = std::max(fabs(obj->getXmax()-obj->getXmin()),fabs(obj->getYmax()-obj->getYmin()));

        if((maxchan>0 && maxchan<nchan) || (maxsize>0 && maxsize<mSize)) {
            numRej++;
            if(par.isVerbose() && par.getShowbar()){
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
template void Cube<short>::rejectObjects(std::vector<Detection<short> >&);
template void Cube<int>::rejectObjects(std::vector<Detection<int> >&);
template void Cube<long>::rejectObjects(std::vector<Detection<long> >&);
template void Cube<float>::rejectObjects(std::vector<Detection<float> >&);
template void Cube<double>::rejectObjects(std::vector<Detection<double> >&);



template <class T>
void Cube<T>::mergeIntoList(Detection<T> &object, std::vector <Detection<T> > &objList) {

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
    typename std::vector<Detection<T> >::iterator iter;
	
    for(iter=objList.begin(); (!haveMerged && iter<objList.end()); iter++) {
		if(iter->canMerge(object, par)){
			iter->addDetection(object);
			haveMerged = true;
		}  
    }
  
    if(!haveMerged) objList.push_back(object);

}
template void Cube<short>::mergeIntoList(Detection<short>&,std::vector<Detection<short> >&);
template void Cube<int>::mergeIntoList(Detection<int>&,std::vector<Detection<int> >&);
template void Cube<long>::mergeIntoList(Detection<long>&,std::vector<Detection<long> >&);
template void Cube<float>::mergeIntoList(Detection<float>&,std::vector<Detection<float> >&);
template void Cube<double>::mergeIntoList(Detection<double>&,std::vector<Detection<double> >&);


template <class T>
void Cube<T>::printDetections (std::ostream& Stream) {
	
	using namespace std;
	
	float RA=-1, DEC=-1, VEL=-1;
	int numObj = getNumObj();
	int m = 10;
    int k=29;
    string str;

	Stream 	<< showpoint << fixed;
    Stream 	<< endl << endl;
	
	if (headDefined) {		
        Stream 	<< "  Detections for " << head.Name() << " " << endl;
        Stream 	<< setw(126) << setfill('_') << " " << endl << endl;
	}
	else { 
        Stream 	<< "  Detections for " << par.getImageFile() << " " << endl;
				
        Stream 	<< setw(94) << setfill('_') << " " << endl << endl;
	}
	
	Stream 	<< setfill(' ');
	
	Stream 	<< setw(m-2) << left << "  Source" 
            << setw(m+6)  << right << "Center    ";
	
    if (headDefined) Stream << setw(k) << right << "Center (WCS)       ";
	
    Stream 	<< setw(m-3) << right << "Xwidth"
            << setw(m-3) << right << "Ywidth"
            << setw(m-3) << right << "Zwidth"
            << setw(m) << right << "NumPix"
            << setw(m) << right << "Vsys "
            << setw(m) << right << "W20  "
            << setw(m) << right << "Flux "
            << setw(m) << right << "PeakSNR";

	Stream 	<< endl;
	
	Stream 	<< setw(m-2)<< left  << "    [#]" 
            << setw(m+6)  << right << "[pix,pix,chan]";
	
    if (headDefined) {
        str = "["+head.Cunit(0)+","+head.Cunit(1)+",KM/S]      ";
        Stream 	<< setw(k) << right << str;
    }
	
    Stream 	<< setw(m-3) << right << "[pix] "
            << setw(m-3) << right << "[pix] "
            << setw(m-3) << right << "[chan]"
            << setw(m) << right << "[pix] "
            << setw(m) << right << "[km/s]"
            << setw(m) << right << "[km/s]"
            << setw(m) << right << "[JY]  ";
	
	Stream 	<< endl;
	
    if (headDefined) Stream << setw(126) << setfill('_') << " " << endl << endl;
    else Stream << setw(94) << setfill('_') << " " << endl << endl;
	
	Stream 	<< setfill(' ');
		
	for (int i=0; i<numObj; i++){
		Detection<T> *obj = new Detection<T>;
        *obj = objectList->at(i);

        obj->calcFluxes(obj->getPixelSet(array, axisDim));
        obj->calcWCSparams(head);
        obj->calcIntegFlux(DimZ(), obj->getPixelSet(array, axisDim), head);

		float Xcenter = obj->getXcentre();
		float Ycenter = obj->getYcentre();
		float Zcenter = obj->getZcentre();
		int Xmin = obj->getXmin();
		int Xmax = obj->getXmax();
		int Ymin = obj->getYmin();
		int Ymax = obj->getYmax();
		int Zmin = obj->getZmin();
		int Zmax = obj->getZmax();
				
		if (headDefined) {
            double pix[3] = {Xcenter,Ycenter,Zcenter};
            double world[3];
            pixToWCSSingle(head.WCS(), pix, world);
            RA  = world[0];
            DEC = world[1];
            VEL = AlltoVel (world[2],head);
            //float zval = ((Zcenter+1-head.Crpix(2))*head.Cdelt(2)+head.Crval(2));
            //VEL = AlltoVel(zval, head);
            if (RA<0) RA += 360;
            else if (RA>360) RA -= 360;
		}
		
        //str = to_string(Xcenter,0)+"  "+to_string(Ycenter,0)+"  "+to_string(Zcenter,0);
        Stream 	<< "     " << setw(m-7) << left << i+1
                << setw(6) << right << to_string(Xcenter,0)
                << setw(5) << right << to_string(Ycenter,0)
                << setw(5) << right << to_string(Zcenter,0);

				  
		if (headDefined) {
           // str = to_string(RA,3)+"  "+to_string(DEC,3)+"  "+to_string(VEL,1);
            //Stream 	<< setw(k) << right << str;
            Stream << setw(11) << right << to_string(RA,3)
                   << setw(9) << right << to_string(DEC,3)
                   << setw(8) << right << to_string(VEL,1);
		}
        //std::string Xint = to_string<int>(Xmin)+"-"+to_string<int>(Xmax);
        //std::string Yint = to_string<int>(Ymin)+"-"+to_string<int>(Ymax);
        //std::string Zint = to_string<int>(Zmin)+"-"+to_string<int>(Zmax);

        double pix[3] = {Xcenter,Ycenter,Zcenter};
        double world[3];
        pixToWCSSingle(head.WCS(), pix, world);

        std::string Xint = to_string<int>(fabs(Xmax-Xmin)+1);
        std::string Yint = to_string<int>(fabs(Ymax-Ymin)+1);
        std::string Zint = to_string<int>(fabs(Zmax-Zmin)+1);

        Stream 	<< setw(m-3) << right << Xint
                << setw(m-3) << right << Yint
                << setw(m-3) << right << Zint
				<< setw(m) << right << obj->getSize();

        Stream << right << setw(m)  << setprecision(1) << obj->getVsys()
               << right << setw(m)  << setprecision(1) << obj->getW20()
               << right << setw(m) << setprecision(3)  << obj->getIntegFlux()
               << right << setw(m) << setprecision(1)  << obj->getPeakFlux()/stats.getSpread();


		Stream 	<< endl;
		
		delete obj;
	}
	
	Stream 	<< endl;
	
}
template void Cube<short>::printDetections(std::ostream&);
template void Cube<int>::printDetections(std::ostream&);
template void Cube<long>::printDetections(std::ostream&);
template void Cube<float>::printDetections(std::ostream&);
template void Cube<double>::printDetections(std::ostream&);


template <class T>
void Cube<T>::plotDetections() {
	
	int numObj = objectList->size();
    /*
    static int a=0;
	static std::ofstream fileout;
	if (a==0) fileout.open("satellites.txt");
	a++;
	if (numObj>1) {
		float highflux=0;
		int numO = 0;
		for (int i=0; i<numObj; i++) {
			Detection<T> *obj = pObject(i);
			obj->calcFluxes(obj->getPixelSet(array, axisDim));
			obj->calcWCSparams(head);
			obj->calcIntegFlux(DimZ(), obj->getPixelSet(array, axisDim), head);
			if (obj->getIntegFlux()>highflux) {highflux=obj->getIntegFlux(); numO=i;}
		}
		//std::cout << pObject(numO)->getIntegFlux() << std::endl; 
		fileout << head.Name() << std::endl << fixed;
		
		int num=1;
		for (int i=0; i<numObj; i++) {
			if (i!=numO) {
				Detection<T> *obj = pObject(i); 
				float ramain = (pObject(numO)->getXcentre()+1-head.Crpix(0))*head.Cdelt(0)+head.Crval(0);
				float ra  = (pObject(i)->getXcentre()+1-head.Crpix(0))*head.Cdelt(0)+head.Crval(0);
				float demain = (pObject(numO)->getYcentre()+1-head.Crpix(1))*head.Cdelt(1)+head.Crval(1);
				float de = (pObject(i)->getYcentre()+1-head.Crpix(1))*head.Cdelt(1)+head.Crval(1);
				std::string center = to_string(obj->getXcentre(),0)+","+to_string(obj->getYcentre(),0)+","+to_string(obj->getZcentre(),0);
				fileout << "#" << setw(4) << num << "    " << setprecision(0)
						<< left << setw(14) << center 
						<< right << setw(16) << decToDMS(getXphys(obj->getXcentre()), head.Ctype(0),0) 
						<< right << setw(16) << decToDMS(getYphys(obj->getYcentre()), head.Ctype(1),1) 
						<< right << setw(8)  << setprecision(0) << obj->getVsys() 
						<< right << setw(8)  << setprecision(1) << obj->getW20()
						<< right << setw(14) << setprecision(5) << obj->getIntegFlux()  
						<< right << setw(10) << angularSeparation(ra,de,ramain,demain) << std::endl;
				num++;
				//fileout << pObject(numO)->getVsys()-obj->getVsys() << endl;
			}
		}
		fileout << std::endl << std::endl;
	}
    */
	
	
	Image2D<T> *HImap = new Image2D<T>(axisDim);
	Image2D<T> *Vemap = new Image2D<T>(axisDim);
	Image2D<T> *Dimap = new Image2D<T>(axisDim);

	ProgressBar bar (" Extracting maps... ", true);
    bar.setShowbar(par.getShowbar());
	if (par.isVerbose()) bar.init(axisDim[1]);
	
	std::vector<bool> isObj(numPix,false);
	
	for (int i=0; i<numObj; i++) {
		typename std::vector<Voxel<T> > voxelList = objectList->at(i).getPixelSet(array, axisDim);
		typename std::vector<Voxel<T> >::iterator vox;
		for(vox=voxelList.begin();vox<voxelList.end();vox++){
			//if(objectList->at(i).isInObject(*vox)){
				long pos = vox->getX()+vox->getY()*axisDim[0]+vox->getZ()*axisDim[0]*axisDim[1];
				isObj[pos] = true;
			//}
		}
	}

	/// Write a datacube with just the detected objects
    Cube<T> *det = new Cube<T>(axisDim);
    for (int i=0; i<det->NumPix(); i++) det->Array()[i] = array[i]*isObj[i];
    det->saveHead(head);
    det->saveParam(par);
    det->Head().setMinMax(0,0);
    det->fitswrite_3d((par.getOutfolder()+"detections.fits").c_str());
    delete det;


	for (int y=0; y<axisDim[1]; y++) {
		if (par.isVerbose()) bar.update(y+1);
		for (int x=0; x<axisDim[0]; x++) {
			float fluxint = 0;
			float fluxsum = 0;
			long mappix = x+y*axisDim[0];
			HImap->Array()[mappix] = 0;
			Vemap->Array()[mappix] = 0;
			Dimap->Array()[mappix] = 0;
			for (int z=0; z<axisDim[2]; z++) {
				long pix = mappix+z*axisDim[1]*axisDim[0];
				if (isObj[pix]) {
					float flux = FluxtoJy(array[pix],head);
					//Pbcor<float>(x,y,z,flux,head);
					fluxsum += flux;
					float zval = ((z+1-head.Crpix(2))*head.Cdelt(2)+head.Crval(2));
					fluxint += flux*zval;
				}
			}
			HImap->Array()[mappix] = fluxsum*fabs(DeltaVel<T>(head));
            Vemap->Array()[mappix] = AlltoVel(fluxint/fluxsum, head);
			
			T num=0;
			for (int z=0; z<axisDim[2]; z++) {
				long pix = mappix+z*axisDim[1]*axisDim[0];
				if (isObj[pix]) {
					float flux = FluxtoJy(array[pix],head);
					//Pbcor<float>(x,y,z,flux,head);
                    float vel = AlltoVel(getZphys(z), head);
					T firstmoment = Vemap->Array(mappix);
					num += flux*(vel-firstmoment)*(vel-firstmoment);		
				}				
			}
			num *= fabs(DeltaVel<T>(head));
			Dimap->Array()[mappix]=sqrt(num/HImap->Array(mappix));	
		}
	}


	std::string name = par.getOutfolder()+head.Obname()+"_mom0th.fits";
	HImap->copyHeader(head);
	HImap->Head().setMinMax(0,0);
	HImap->Head().setBtype("intensity");
	HImap->Head().setBunit("JY * KM/S");	
	HImap->fitswrite_2d(name.c_str());

    Vemap->copyHeader(head);
    Vemap->Head().setMinMax(0,0);
    Vemap->Head().setBtype("velocity");
    Vemap->Head().setBunit("KM/S");
    name = par.getOutfolder()+head.Obname()+"_mom1th.fits";
    Vemap->fitswrite_2d(name.c_str());
	
    Dimap->copyHeader(head);
    Dimap->Head().setMinMax(0,0);
    Dimap->Head().setBtype("vdisp");
    Dimap->Head().setBunit("KM/S");
    name = par.getOutfolder()+head.Obname()+"_mom2nd.fits";
    Dimap->fitswrite_2d(name.c_str());
	
	Image2D<int> *DetMap = new Image2D<int>(axisDim);
	for (int i=0; i<axisDim[0]*axisDim[1];i++) DetMap->Array()[i] = detectMap[i];
	DetMap->copyHeader(head);
	DetMap->Head().setMinMax(0,0);
	DetMap->Head().setBtype("detected_chan");
	name = par.getOutfolder()+"DetectMap.fits";
	DetMap->fitswrite_2d(name.c_str());
	
	delete DetMap;
	delete HImap;
	delete Vemap;
	delete Dimap;
	
	if (par.isVerbose()) bar.fillSpace(" OK.\n");
	
	std::ofstream fileo;

	for (int i=0; i<numObj; i++) {
		Detection<T> *obj = pObject(i);
		float *intSpec = new float[axisDim[2]];
		obj->calcWCSparams(head);
		for(int z=0; z<axisDim[2]; z++) intSpec[z]=0;      
		
		for (int z=obj->getZmin(); z<=obj->getZmax(); z++) {
			for (int x=obj->getXmin(); x<=obj->getXmax(); x++) {
				for (int y=obj->getYmin(); y<=obj->getYmax(); y++) {
					if (obj->isInObject(x,y,z)) {
						long pix = x+y*axisDim[0]+z*axisDim[1]*axisDim[0];
						double flux = array[pix];
						Pbcor<double>(x,y,z,flux,head);
						intSpec[z] += flux;
					}
				}
			}
		
		}
		
		std::string name = par.getOutfolder()+"spectrum"+to_string(i)+".dat";
		fileo.open(name.c_str());
		for (int z=0; z<axisDim[2]; z++) {
			double vel = ((z+1-head.Crpix(2))*head.Cdelt(2)+head.Crval(2));
            vel = AlltoVel(vel, head);
			intSpec[z] = FluxtoJy(intSpec[z],head);
			fileo << vel << "  " << intSpec[z] << endl;
		}
		
		delete [] intSpec;
		fileo.close();
		
	}

#ifdef HAVE_GNUPLOT		
	Gnuplot gp;	
	
	gp.begin();	
	gp.commandln("set terminal postscript eps color");
	gp.commandln("unset key");
	gp.commandln("set xlabel 'Velocity [km/s]'");
	gp.commandln("set ylabel 'Flux density [JY]");
	gp.commandln("set size square");
	
	std::string globalcmd;
	for (int i=0; i<numObj; i++) {
		float vmin = objectList->at(i).getVelMin();
		float vmax = objectList->at(i).getVelMax();
		float delta = 2*fabs(DeltaVel<float>(head));
		if (vmin>vmax) std::swap(vmin, vmax);
		vmin -= delta;
		vmax += delta;		
		std::string minn = to_string<int>(vmin);
		std::string maxx = to_string<int>(vmax);
		gp.commandln(("set xrange ["+minn+":"+maxx+"]").c_str());
		std::string titlecmd = "set title 'Individual spectrum for " +to_string(i)+"'";
		std::string outcmd = "set output '"+par.getOutfolder()+"spectrum #"+to_string(i)+".eps'";
		std::string plotcmd = "plot '"+par.getOutfolder()+"spectrum"+to_string(i)+".dat"+"' with lp ls 3 lt 7";
		gp.commandln(titlecmd.c_str());
		gp.commandln(outcmd.c_str());
		gp.commandln((plotcmd).c_str());
	}
	
	gp.end();
#endif
	
	for (int i=0; i<numObj; i++)
		remove((par.getOutfolder()+"spectrum"+to_string(i)+".dat").c_str());
	
	
	
}
template void Cube<short>::plotDetections();
template void Cube<int>::plotDetections();
template void Cube<long>::plotDetections();
template void Cube<float>::plotDetections();
template void Cube<double>::plotDetections();

template <class Type>
bool Cube<Type>::isDetection(long x, long y, long z) {

  /// Is a given voxel at position (x,y,z) a detection, based on the statistics
  /// in the Cube's Stats? 
  /// If the pixel lies outside the valid range for the data array, return false.
  ///
  /// \param x 	X-value of the Cube's voxel to be tested.
  /// \param y 	Y-value of the Cube's voxel to be tested.
  /// \param z 	Z-value of the Cube's voxel to be tested.

    long voxel = z*axisDim[0]*axisDim[1] + y*axisDim[0] + x;
    return stats.isDetection(array[voxel]);
}
template bool Cube<short>::isDetection(long,long,long);
template bool Cube<int>::isDetection(long,long,long);
template bool Cube<long>::isDetection(long,long,long);
template bool Cube<float>::isDetection(long,long,long);
template bool Cube<double>::isDetection(long,long,long);


template <class Type>
Detection<Type>* Cube<Type>::LargestDetection () {
    int numObj = objectList->size();
    if (numObj==0) return NULL;
    uint n=0, size=0;
    for (int i=0; i<numObj; i++) {
        Detection<Type> *obj = pObject(i);
        if (obj->getSize()>size) {n=i;size=obj->getSize();}
    }
    return pObject(n);
}
template Detection<short>* Cube<short>::LargestDetection ();
template Detection<int>* Cube<int>::LargestDetection ();
template Detection<long>* Cube<long>::LargestDetection ();
template Detection<float>* Cube<float>::LargestDetection ();
template Detection<double>* Cube<double>::LargestDetection ();
