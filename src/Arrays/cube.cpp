//--------------------------------------------------------------------
// cube.cpp: Member functions for the Cube class
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
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include "cube.hh"
#include "stats.hh"
#include "../Map/detection.hh"
#include "../Utilities/utils.hh"
#include "../Utilities/progressbar.hh"
#include "../Utilities/smooth3D.hh"


template <class T>
void Cube<T>::defaults() {
	
	numAxes = 3;
	arrayAllocated = false;
	headDefined = false; 
	axisDimAllocated = false; 
	statsDefined = false;
	maskAllocated = false;
    mapAllocated  = false;
    isSearched = false;
	objectList = new std::vector<Detection<T> >; 
	
}
template void Cube<short>::defaults(); 
template void Cube<int>::defaults(); 
template void Cube<long>::defaults(); 
template void Cube<float>::defaults(); 
template void Cube<double>::defaults(); 



template <class T>
Cube<T>::Cube() {	
	
	defaults();
	
}
template Cube<short>::Cube();
template Cube<int>::Cube();
template Cube<long>::Cube();
template Cube<float>::Cube();
template Cube<double>::Cube();



template <class T>
Cube<T>::~Cube () {
	
	if (arrayAllocated) delete [] array;
	arrayAllocated=false;
	if (maskAllocated)  delete [] mask;
	maskAllocated=false;
	if (axisDimAllocated) delete [] axisDim;
	axisDimAllocated=false;
	if (mapAllocated) delete [] detectMap;
	mapAllocated=false;
	delete objectList;
} 
template Cube<short>::~Cube();
template Cube<int>::~Cube();
template Cube<long>::~Cube();
template Cube<float>::~Cube();
template Cube<double>::~Cube();
 


template <class T>
Cube<T>::Cube(std::string fname) {
	
	defaults();
	
	par.setImageFile(fname);
	numAxes = 3;
	
	head.header_read(par.getImageFile());
	headDefined = true;
	axisDim = new int [numAxes];
	axisDimAllocated = true; 
	for (int i=0; i<numAxes; i++) axisDim[i] = head.DimAx(i);
	numPix = axisDim[0]*axisDim[1]*axisDim[2];
	fitsread_3d();	
	detectMap = new short [axisDim[0]*axisDim[1]];
	mapAllocated = true;

}
template Cube<short>::Cube(std::string);
template Cube<int>::Cube(std::string);
template Cube<long>::Cube(std::string);
template Cube<float>::Cube(std::string);
template Cube<double>::Cube(std::string);



template <class T>
Cube<T>::Cube(int *dimensions) {
    
    defaults();
    
    int size   = dimensions[0]*dimensions[1]*dimensions[2];
    int imsize = dimensions[0]*dimensions[1];
    axisDimAllocated = false;
    arrayAllocated = false;
    numPix = numAxes = 0;
    if((size<0) || (imsize<0) ) {
        std::cout << "Error [Cube(dimArray)]: Negative size -- could not define Cube "<< std::endl;
        std::terminate();
    }
    else {
		numPix = size;
		if(size>0){
			array = new T[size];
			arrayAllocated = true;
			detectMap = new short[imsize];
			mapAllocated = true;
		}
		numAxes  = 3;
		axisDim = new int[numAxes];
		axisDimAllocated = true;
		for(int i=0; i<numAxes; i++) axisDim[i] = dimensions[i];   
		for(int i=0; i<imsize; i++) detectMap[i] = 0; 
    }
    
}
template Cube<short>::Cube(int *);
template Cube<int>::Cube(int*);
template Cube<long>::Cube(int*);
template Cube<float>::Cube(int*);
template Cube<double>::Cube(int*);



template <class T>
Cube<T>::Cube(const Cube<T> &c) {
  
	this->operator=(c);

}
template Cube<short>::Cube(const Cube<short> &);
template Cube<int>::Cube(const Cube<int> &);
template Cube<long>::Cube(const Cube<long> &);
template Cube<float>::Cube(const Cube<float> &);
template Cube<double>::Cube(const Cube<double> &);


template <class T>
Cube<T>& Cube<T>::operator=(const Cube<T> &c) {
    
    if(this==&c) return *this;
    
    if(this->arrayAllocated) delete [] array;
	if(this->axisDimAllocated) delete [] axisDim;
	if(this->maskAllocated) delete [] mask;
	
    this->numPix	= c.numPix;
    this->numAxes	= c.numAxes;
    this->par		= c.par;
    
    this->axisDimAllocated = c.axisDimAllocated;
    if (axisDimAllocated) {
		this->axisDim = new int [numAxes];
		for (int i=0; i<numAxes;  i++) this->axisDim[i] = c.axisDim[i];
	}
	
    this->arrayAllocated = c.arrayAllocated;
    if(this->arrayAllocated) {
		this->array = new T[this->numPix];
		for(int i=0; i<this->numPix; i++) this->array[i] = c.array[i];
    }
    
    this->maskAllocated = c.maskAllocated; 
    if(this->maskAllocated) {
		this->mask = new bool[numPix];
		for(int i=0;i<this->numPix;i++) this->mask[i] = c.mask[i];
    }
    
    this->mapAllocated = c.mapAllocated; 
    if(this->mapAllocated) {
		this->detectMap = new short[this->axisDim[0]*this->axisDim[1]];
		for(int i=0;i<(this->axisDim[0]*this->axisDim[1]);i++) this->detectMap[i] = c.detectMap[i];
    }
    
    this->headDefined = c.headDefined;
	if (this->headDefined) this->head = c.head;
	this->statsDefined = c.statsDefined;
	if (this->statsDefined) this->stats = c.stats;
    this->isSearched = c.isSearched;
	
    return *this;
}
template Cube<short>& Cube<short>::operator=(const Cube<short>&);
template Cube<int>& Cube<int>::operator=(const Cube<int>&);
template Cube<long>& Cube<long>::operator=(const Cube<long>&);
template Cube<float>& Cube<float>::operator=(const Cube<float>&);
template Cube<double>& Cube<double>::operator=(const Cube<double>&);



/**=====================================================================================*/ 
/** FUNCTIONS FOR INPUT AND OUTPUT */

template <class T>
void Cube<T>::setCube (T *input, int *dim) {

	if (arrayAllocated) delete [] array;
	if (axisDimAllocated) delete [] axisDim;
	numAxes = 3;
	axisDim = new int [numAxes];
	for (int i=0; i<numAxes; i++) axisDim[i] = dim[i];
	axisDimAllocated = true;
	numPix = axisDim[0]*axisDim[1]*axisDim[2];
	array = new T [numPix];
	arrayAllocated=true;
	for (int i=0; i<axisDim[2]; i++)
		for (int j=0; j<axisDim[1]; j++) 
			for (int k=0; k<axisDim[0]; k++) {
				long nPix = k+j*axisDim[0]+i*axisDim[0]*axisDim[1];
				array[nPix]=input[nPix]; 
			}
	if (mapAllocated) delete [] detectMap;
	detectMap = new short [axisDim[0]*axisDim[1]];
	mapAllocated=true;
	
}
template void Cube<short>::setCube (short*, int*);
template void Cube<int>::setCube (int*, int*);
template void Cube<long>::setCube (long*, int*);
template void Cube<float>::setCube (float*, int*);
template void Cube<double>::setCube (double*, int*);



template <class T>
bool Cube<T>::readCube (std::string fname) {
	
	par.setImageFile(fname);
	numAxes = 3;

    if(!head.header_read(par.getImageFile())) return false;

	headDefined = true;
        // I do not like the two folliwing 2 lines, I should think something better
        if (par.getRedshift()!=-1) head.setRedshift(par.getRedshift());
        if (par.getRestwave()!=-1) head.setWave0(par.getRestwave());
	axisDim = new int [numAxes];
	axisDimAllocated = true; 
	for (int i=0; i<numAxes; i++) axisDim[i] = head.DimAx(i);
    if (head.NumAx()<3) axisDim[2] = 1;
    if (head.NumAx()<2) axisDim[1] = 1;
	numPix = axisDim[0]*axisDim[1]*axisDim[2];
	if (!fitsread_3d()) return false;	
	objectList = new std::vector<Detection<T> >; 
	detectMap = new short [axisDim[0]*axisDim[1]];
	mapAllocated = true;
	for (int i=0; i<axisDim[0]*axisDim[1]; i++) detectMap[i] = 0;
	return true;
}
template bool Cube<short>::readCube (std::string);
template bool Cube<int>::readCube (std::string);
template bool Cube<long>::readCube (std::string);
template bool Cube<float>::readCube (std::string);
template bool Cube<double>::readCube (std::string);



template <class T>
bool Cube<T>::fitsread_3d() {

	fitsfile *fptr3;
	int status, anynul, fpixel;
	float nulval;

	if (par.isVerbose()) { 
        std::cout << "\nOpening file "<< par.getImageFile() << std::endl;
		std::cout << "Reading "<<axisDim[0]<<" x "<<axisDim[1]<<" x "<<axisDim[2]
				  << " pixels FITS file... ";
    }

	// Open the FITS file
	status = 0;
    if(fits_open_file(&fptr3, par.getImageFile().c_str(), READONLY, &status) ){
      fits_report_error(stderr, status);
      return false;
    }

	// Read elements from the FITS data array    
    if (!arrayAllocated) array = new T[numPix];						
    arrayAllocated = true;
    fpixel=1;

    status=0;
    if (fits_read_img(fptr3, selectDatatype<T>(), fpixel, numPix, &nulval, array, &anynul, &status)){
		fits_report_error(stderr, status);
		return false;
    }  

    // Close the FITS File
    if (fits_close_file(fptr3, &status)){
		fits_report_error(stderr, status);
	}

	if (par.isVerbose()) std::cout << "Done.\n" << std::endl;

	return true;
}
template bool Cube<short>::fitsread_3d();
template bool Cube<int>::fitsread_3d();
template bool Cube<long>::fitsread_3d();
template bool Cube<float>::fitsread_3d();
template bool Cube<double>::fitsread_3d();



template <class T>
bool Cube<T>::fitswrite_3d(const char *outfile, bool fullHead) {
	
	fitsfile *fptr;      
	long  fpixel = 1;
	long dnaxes[3] = {axisDim[0], axisDim[1], axisDim[2]};
	int status=0;
  
	remove(outfile);             
	 
	if (fits_create_file(&fptr, outfile, &status)) {
		fits_report_error(stderr, status); 
		return false; 
	}		
	
	status=0; 
	if (fits_create_img(fptr, selectBitpix<T>(), 3, dnaxes, &status)) {
		fits_report_error(stderr, status); 
		return false; 
	}
    
    if (headDefined) {
        if (head.NumAx()==2) head.headwrite_2d(fptr,fullHead);
        else head.headwrite_3d (fptr, fullHead);
    }
	
	status=0;
	if (fits_write_img(fptr, selectDatatype<T>(), fpixel, numPix, array, &status)) {
		fits_report_error(stderr, status);
		return false;
	}

	if (fits_close_file(fptr, &status)) {
		fits_report_error(stderr, status);
	}      
	
	return true;
}
template bool Cube<short>::fitswrite_3d(const char*,bool);
template bool Cube<int>::fitswrite_3d(const char*,bool);
template bool Cube<long>::fitswrite_3d(const char*,bool);
template bool Cube<float>::fitswrite_3d(const char*,bool);
template bool Cube<double>::fitswrite_3d(const char*,bool);


/**=====================================================================================*/ 
/** STATISTICAL FUNCTIONS */

template <class T>
void Cube<T>::setCubeStats() {

  /// Calculates the full statistics for the cube: mean, rms, median, madfm.
  /// Also work out the threshold and store it in the stats set.

    if(par.isVerbose()) std::cout << " Calculating statistics for the cube... "<<std::flush;
    
    stats.setRobust(par.getFlagRobustStats());
    bool *blanks = new bool[numPix];
    for (int i=0; i<numPix; i++) blanks[i] = isBlank(array[i]) ? false : true;

    stats.calculate(array,numPix,blanks);
	stats.setThresholdSNR(par.getCut());

    if(par.isVerbose()) {
		std::cout << "Using flux threshold of: ";
		T thresh;
		if (par.getFlagUserThreshold()) thresh = par.getThreshold();
        else thresh = stats.getThreshold();
        if (thresh<1E-04) std::cout << std::scientific;
        else std::cout << std::fixed;
        std::cout << std::setprecision(5) << thresh << " " << head.Bunit() << std::endl;
		std::cout << std::setw(52) << std::right << "(middle = " 
                  << stats.getMiddle() << ", spread = " << stats.getSpread() << ")\n" << std::fixed;
	}

    delete [] blanks;
	
	statsDefined = true;
    
}
template void Cube<short>::setCubeStats();
template void Cube<int>::setCubeStats();
template void Cube<long>::setCubeStats();
template void Cube<float>::setCubeStats();
template void Cube<double>::setCubeStats();



template <class T>
void Cube<T>::BlankCube (T *Array, long size) {
	
  /// A function for blanking an array with Cube mask data.
  ///
  /// \param Array		The array to blank.
  /// \param size		The size of array. It must be equal
  ///					to the size of Cube-object array.

	if (size!=numPix) 
        std::cout << "Error blanking cube: array size is different from cube size" << std::endl;
	else {
		if (!maskAllocated) BlankMask();
		for (int i=0; i<size; i++) Array[i] *= mask[i];
	}
}
template void Cube<short>::BlankCube (short*, long);
template void Cube<int>::BlankCube (int*, long);
template void Cube<long>::BlankCube (long*, long);
template void Cube<float>::BlankCube (float*, long);
template void Cube<double>::BlankCube (double*, long);



template <class T>
void Cube<T>::BlankMask (float *channel_noise){
	
     /*/////////////////////////////////////////////////////////////////////////////////
     * This function builds a mask for the cube. The type of mask depends on the
     * parameter MASK:
     *
     *  - SEARCH:      Uses the source finding algorith and masks the largest object.
     *                 All search-related parameters can be used.
     *  - THRESHOLD:   Applies a simple threshold cut given by THRESHOLD parameter.
     *  - SMOOTH:      Smooths the cube by a factor FACTOR and applies a S/N threshold
     *                 on the smoothed cube given by BLANKCUT. Default are FACTOR=2
     *                 and BLANKCUT=3. If BMAJ and BMIN parameters are, it smooths to
     *                 these values.
     *  - NEGATIVE:    Calculates the noise statitistics just on the negative pixels
     *                 and builds the mask based on the S/N threshold BLANKCUT.
     *  - FILE(Name):  User-provided mask. 'Name' is a fitsfile with same size of the
     *                 cube and filled with 0(false) or 1(true).
     *  - NONE:        No mask, just pixels > 0.
     *
     /////////////////////////////////////////////////////////////////////////////////*/

    if (maskAllocated) delete [] mask;
    mask = new bool[numPix];

    for (int i=0; i<numPix; i++) mask[i]=0;

    bool verb = par.isVerbose();
    if (verb) {
        std::cout << " Creating mask (" << par.getMASK() << ") ..." << std::flush;
        par.setVerbosity(false);
    }

    Statistics::Stats<T> *st = new Statistics::Stats<T>;
    st->setRobust(par.getFlagRobustStats());

    if (par.getMASK()=="SEARCH") {
        // Masking using the search algorithm and mask the largest of object.
        if (!isSearched) Search();
        Detection<T> *larg = LargestDetection();
        if (larg==NULL) {
            std::cout << "3DFIT error: No sources detected in the datacube. Cannot build mask!!! \n";
            std::terminate();
        }
        std::vector<Voxel<T> > voxlist = larg->getPixelSet();
        typename std::vector<Voxel<T> >::iterator vox;
        for(vox=voxlist.begin();vox<voxlist.end();vox++) {
            mask[nPix(vox->getX(),vox->getY(),vox->getZ())]=1;
        }
    }
    else if (par.getMASK()=="THRESHOLD") {
        float thresh = par.getThreshold();
        for (uint i=numPix; i--;) {
            if (array[i]>thresh) mask[i] = 1;
        }
    }
    else if (par.getMASK()=="SMOOTH") {
        double bmaj  = head.Bmaj()*arcsconv(head.Cunit(0));
        double bmin  = head.Bmin()*arcsconv(head.Cunit(0));
        double bpa   = head.Bpa();
        float factor = par.getFactor()==-1 ? 2 : par.getFactor();
        double nbmaj = par.getBmaj()==-1 ? factor*bmaj : par.getBmaj();
        double nbmin = par.getBmin()==-1 ? factor*bmin : par.getBmin();
        double nbpa  = par.getBpa()==-1  ? bpa    : par.getBpa();
        Beam oldbeam = {bmaj,bmin,bpa};
        Beam newbeam = {nbmaj,nbmin,nbpa};

        Smooth3D<T> *sm = new Smooth3D<T>;
        sm->smooth(this, oldbeam, newbeam);
        bool *blanks = new bool[numPix];
        for (int i=0; i<numPix; i++) blanks[i] = isBlank(sm->Array(i)) ? false : true;
        st->calculate(sm->Array(),numPix,blanks);
        st->setThresholdSNR(par.getBlankCut());

        ///* Without three consecutive channels requirement
        for (size_t i=0; i<numPix; i++) {
            if (sm->Array(i)>st->getThreshold()) mask[i] = 1;
        }
        //*/

        /* With three consecutive channels requirement
        T thr = st->getThreshold();
        T *Array = sm->Array();
        for (int z=1; z<axisDim[2]-1; z++) {
            for (int y=0; y<axisDim[1]; y++) {
                for (int x=0; x<axisDim[0]; x++) {
                    long npix = nPix(x,y,z);
                    long nchan = nPix(x,y,z+1);
                    long pchan = nPix(x,y,z-1);
                    mask[npix] = Array[npix]>thr && Array[pchan]>thr && Array[nchan]>thr;
                }
            }
        }
        for (int y=0; y<axisDim[1]; y++) {
            for (int x=0; x<axisDim[0]; x++) {
                mask[nPix(x,y,0)]=Array[nPix(x,y,0)]>thr && Array[nPix(x,y,1)]>thr && Array[nPix(x,y,2)]>thr;
                int l = axisDim[2]-1;
                mask[nPix(x,y,l)]=Array[nPix(x,y,l)]>thr && Array[nPix(x,y,l-1)]>thr && Array[nPix(x,y,l-2)]>thr;
            }
        }
        */

        delete sm;
        delete [] blanks;
    }
    else if (par.getMASK()=="NEGATIVE") {
        for (int z=0; z<DimZ(); z++) {
            std::vector<T> onlyneg;
            for (int i=0; i<DimX()*DimY(); i++)  {
                if (array[i+z*DimX()*DimY()]<0) {
                    onlyneg.push_back(array[i+z*DimX()*DimY()]);
                    onlyneg.push_back(-array[i+z*DimX()*DimY()]);
                }
            }
            st->calculate(&onlyneg[0], onlyneg.size());
            st->setThresholdSNR(par.getBlankCut());

            for (int i=0; i<DimX()*DimY(); i++)  {
                if (array[i+z*DimX()*DimY()]>st->getThreshold()) mask[i+z*DimX()*DimY()]=1;
            }
            if (channel_noise!=NULL) channel_noise[z]=st->getSpread();
        }
    }
    else if (par.getMASK().find("FILE(")!=-1) {
        std::string str = par.getMASK();
        size_t first = str.find_first_of("(");
        size_t last = str.find_last_of(")");
        if (first==-1 || last==-1) {
            std::cerr << "\n  ERROR: MASK parameter is not correct. Provide file(Maskfitsfile)\n";
            std::terminate();
        }
        std::string filename = str.substr (first+1,last-first-1);
        Cube<short> *ma = new Cube<short>;

        if (!fexists(filename) || !ma->readCube(filename)) {
            std::cerr << "\n ERROR: Mask " << filename
                      << " is not a readable FITS image! Exiting ...\n";
            std::terminate();
        }

        if (ma->NumPix()!=numPix || ma->DimX()!=axisDim[0] ||
            ma->DimY()!=axisDim[1] || ma->DimZ()!=axisDim[2]) {
            std::cerr << "\n ERROR: Mask file and data file have different dimensions." << filename
                      << " Exiting ...\n";
            std::terminate();
        }

        for (size_t i=0; i<numPix; i++) mask[i] = ma->Array(i);

        delete ma;

    }
    else if (par.getMASK()=="NONE") {
        for (uint i=numPix; i--;) {
            if (array[i]>0) mask[i] = 1;
        }
    }

    delete st;


    if (verb) {
        std::cout << " Done." << std::endl;
        par.setVerbosity(true);
    }

    maskAllocated = true;
	
	/*
	//	OLD VERSION
	 
	/// This function is used for blanking a cube.
	/// It creates a boolean mask: the mask takes into account only
	/// those regions which show emission in three consecutive
	/// channels above a set level (3σ) and set bool value to 1.
	/// Else set the mask value to 0.
	
	
	if (!statsDefined) {
		if (!par.isVerbose()) setCubeStats();
		else {
			par.setVerbosity(false);
			setCubeStats();
			par.setVerbosity(true);
		}
	}
	
	
	T thresh, middle, spread;
	T snrCut = par.getBlankCut(); 
	if (par.getFlagRobustStats()) {
		middle = stats.getMedian();
		spread = stats.getMadfm()/0.6744888 ;
	}
	else {
		middle = stats.getMean();
		spread = stats.getStddev();
	}
	
	thresh = middle + snrCut*spread;
	
	mask = new bool [numPix];
	maskAllocated = true;
	
	if (par.isVerbose()) std::cout << " Blanking the data cube..." << std::flush;
	ProgressBar bar;
	if (par.isVerbose()) bar.init(axisDim[0]);
	
	for (int x=0; x<axisDim[0]; x++) {
		if (par.isVerbose()) bar.update(x+1);
		for (int y=0; y<axisDim[1]; y++) {
			for (int z=0; z<axisDim[2]-2; z++) {
				long npix = x+y*axisDim[0]+z*axisDim[0]*axisDim[1];
				long nextchan = x+y*axisDim[0]+(z+1)*axisDim[0]*axisDim[1];
				long nnchan = x+y*axisDim[0]+(z+2)*axisDim[0]*axisDim[1];
				bool isGood = array[npix]>thresh && array[nextchan]>thresh && array[nnchan]>thresh;
				if (isGood) {
					mask[npix] = 1;
					continue;
				}
				if (z>0) {
					long prechan = x+y*axisDim[0]+(z-1)*axisDim[0]*axisDim[1];
					isGood = array[npix]>thresh && array[nextchan]>thresh && array[prechan]>thresh;
					if (isGood) {
						mask[npix] = 1;
						continue;
					}
				}
				if (z>1) {
					long prechan = x+y*axisDim[0]+(z-1)*axisDim[0]*axisDim[1];
					long ppchan = x+y*axisDim[0]+(z-2)*axisDim[0]*axisDim[1];
					isGood = array[npix]>thresh && array[prechan]>thresh && array[ppchan]>thresh;
					if (isGood) {
						mask[npix] = 1;
						continue;
					}
				}
				
				mask[npix] = 0;			
			}
			for (int z=axisDim[2]-1; z>axisDim[2]-3; z--) {
				long npix = x+y*axisDim[0]+z*axisDim[0]*axisDim[1];
				long prechan = x+y*axisDim[0]+(z-1)*axisDim[0]*axisDim[1];
				long ppchan = x+y*axisDim[0]+(z-2)*axisDim[0]*axisDim[1];
				bool isGood = array[npix]>thresh && array[prechan]>thresh && array[ppchan]>thresh;
				if (isGood) mask[npix] = 1;
				else mask[npix] = 0;
			}
		}		
	}
	if (par.isVerbose()) bar.fillSpace(" Done.\n");
	*/
}
template void Cube<short>::BlankMask (float*);
template void Cube<int>::BlankMask (float*);
template void Cube<long>::BlankMask (float*);
template void Cube<float>::BlankMask (float*);
template void Cube<double>::BlankMask (float*);


/**=====================================================================================*/ 



template <class T>
Cube<T>* Cube<T>::Reduce (int fac) {
	
	int dim[3] = {axisDim[0]/fac,axisDim[1]/fac,axisDim[2]};
	
	
	if (par.isVerbose()) std::cout << " Reducing..." << std::flush;

	Cube<T> *reduced = new Cube<T>(dim);
	reduced->saveParam(par);
	reduced->saveHead(head);
    reduced->Head().setCdelt(0, fac*head.Cdelt(0));
	reduced->Head().setCdelt(1, fac*head.Cdelt(1));
	reduced->Head().setCrpix(0, lround(head.Crpix(0)/double(fac)));
	reduced->Head().setCrpix(1, lround(head.Crpix(1)/double(fac)));
	reduced->Head().calcArea();
	
	T *Array = reduced->Array();
	
	std::string obeamsize = "  Old beam size: "+to_string(Head().BeamArea())+" pixels";
	std::string nbeamsize = "  New beam size: "+to_string(reduced->Head().BeamArea())+" pixels";
	
	int xx, yy;
	for (int z=0; z<dim[2]; z++) {
		xx = yy = 0;
		for (int y=0; y<dim[1]; y++) {
			for (int x=0; x<dim[0]; x++) {
				long Arraypix = x+y*dim[0]+z*dim[0]*dim[1];
				T sum = 0;
				for (yy=fac*y; yy<fac*y+fac; yy++) {
					for (xx=fac*x; xx<fac*x+fac; xx++) {
						long arraypix = xx+yy*axisDim[0]+z*axisDim[0]*axisDim[1];
						sum += array[arraypix];
					}
				}
				Array[Arraypix] = sum/double(fac*fac);
			}
		}
	}
		
	if (par.isVerbose()) std::cout << "OK\n" << obeamsize << std::endl << nbeamsize << std::endl;
	
	T minn,maxx;
	findMinMax<T>(reduced->Array(), reduced->NumPix(), minn, maxx);
	reduced->Head().setDataMax(double(maxx));
	reduced->Head().setDataMin(double(minn));
	reduced->Head().Keys().push_back("HISTORY BBAROLO RESAMPLING: "+nbeamsize);
	reduced->Head().Keys().push_back("HISTORY BBAROLO RESAMPLING: "+obeamsize);
	
	return reduced;
	
}
template Cube<short>* Cube<short>::Reduce (int);
template Cube<int>* Cube<int>::Reduce (int);
template Cube<long>* Cube<long>::Reduce (int);
template Cube<float>* Cube<float>::Reduce (int);
template Cube<double>* Cube<double>::Reduce (int);	
	
	
	
template <class T>
void Cube<T>::CheckChannels () {
	
	int xySize = axisDim[0]*axisDim[1];
	int zdim = axisDim[2];
	T *mapchan = new T[xySize];
	
	if (!statsDefined) setCubeStats();
	
	std::cout << "\n\n";
	
	ProgressBar bar(" Checking for bad channels... ", true);
    bar.setShowbar(par.getShowbar());
	if (par.isVerbose()) bar.init(zdim);
	
	std::vector<int> bad;
	
	for (int z=0; z<axisDim[2]; z++) {
		if(par.isVerbose()) bar.update(z+1);		
		for (int pix=0; pix<xySize; pix++) mapchan[pix] = array[z*xySize+pix];
		T median = findMedian<T>(mapchan, xySize, false);
		T MADFM = findMADFM<T>(mapchan, xySize, median, false);
		if (MADFM > 1.5*stats.getMadfm()) {
			bad.push_back(z);
			for (int pix=0; pix<xySize; pix++) array[z*xySize+pix] = 0.;
		} 
	}
	
	
	if (bad.size()==0) {
		if(par.isVerbose()){
			bar.fillSpace(" All channels are good.");
			std::cout << "\n\n";
			delete [] mapchan;
			return;
		}
	}
	
	int count=0;
	if (par.isVerbose()) {
		bar.fillSpace(" Following channels have been erased: ");
		for (unsigned int i=0; i<bad.size(); i++) {
			if (bad[i]==int(i)) count++;
			std::cout << bad[i]+1 << " ";
		}
		std::cout << "\n\n";
	}
	
	
	float med = stats.getMedian();
	bool oldv = par.isVerbose();
	if (par.isVerbose()) std::cout << " Re-calculating statistics for new cube ..." << std::flush;
	par.setVerbosity(false);
	setCubeStats();
	stats.setMedian(med);
	par.setVerbosity(oldv);
	if (par.isVerbose()) {
		std::cout << " Done. New threshold is " << stats.getThreshold() << " " << head.Bunit() 
				  << std::endl << std::endl;
	}
	
	string name = par.getImageFile();	
	name = par.getImageFile();
	int found = name.find(".fits");
	name.insert(found, "ean");
	
	int NdatZ = axisDim[2]-bad.size(); 
	int ax[3] = {axisDim[0], axisDim[1], NdatZ};
	Cube<T> *out = new Cube<T>(ax);
	for (int i=0; i<out->NumPix(); i++) 
		out->Array()[i] = array[i+count*xySize];
 	out->saveHead(head);
 	out->saveParam(par);
 	out->Head().setMinMax(0,0);
 	out->Head().setDimAx(0, axisDim[0]);
 	out->Head().setDimAx(1, axisDim[1]);
 	out->Head().setDimAx(2, NdatZ);
 	out->Head().setCrpix(0, head.Crpix(0));
 	out->Head().setCrpix(1, head.Crpix(1));
 	out->Head().setCrpix(2, head.Crpix(2)-count);
	//out->fitswrite_3d(name.c_str(),true);
	
	/// MODO ALTERNATIVO - COPIA DELL'HEADER ORIGINALE
	fitsfile *infptr, *outfptr;
	int status=0;
	char com[100];
	remove (name.c_str());
	fits_open_file(&infptr, par.getImageFile().c_str(), READONLY, &status);
	fits_create_file(&outfptr, name.c_str(), &status); 	
	fits_copy_header(infptr, outfptr, &status);
	fits_update_key_lng(outfptr, "BITPIX", Head().Bitpix(), com, &status);
	fits_update_key_lng(outfptr, "NAXIS3", NdatZ, com, &status);
	fits_update_key_dbl(outfptr, "CRPIX3", head.Crpix(2)-count, 12, com, &status);
	if (head.Bmaj()!=0) fits_update_key_dbl(outfptr, "BMMAJ", head.Bmaj(), 12, com, &status);	
	if (head.Bmin()!=0) fits_update_key_dbl(outfptr, "BMMIN", head.Bmin(), 12, com, &status);	
	fits_write_img(outfptr, TFLOAT, 1, xySize*NdatZ, out->Array(), &status);
	fits_close_file(infptr, &status);
	fits_close_file(outfptr,&status);    
	fits_report_error(stderr, status);
	///  -------------------------------------------------------
	
	
	delete out;
	delete [] mapchan;
	
}
template void Cube<short>::CheckChannels ();
template void Cube<int>::CheckChannels ();
template void Cube<long>::CheckChannels ();
template void Cube<float>::CheckChannels ();
template void Cube<double>::CheckChannels ();
