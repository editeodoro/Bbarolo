//--------------------------------------------------------------------
// image.cpp: Members functions of the Image class.
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

#include <iostream>
#include <fstream>
#include <cstring>
#include <fitsio.h>
#include <Arrays/stats.hh>
#include <Arrays/header.hh>
#include <Arrays/image.hh>
#include <Map/object2D.hh>
#include <Map/voxel.hh>
#include <Map/scan.hh>
#include <Utilities/utils.hh>

template <class Type>
void Image2D<Type>::defaults() {
    
    numPix  = 0;
    numAxes = 2;
    arrayAllocated = false;
    headDefined = false;
    statsDefined= false;
}
template void Image2D<short>::defaults();
template void Image2D<int>::defaults();
template void Image2D<long>::defaults();
template void Image2D<float>::defaults();
template void Image2D<double>::defaults();


template <class Type> 
Image2D<Type>::Image2D () {
    
    defaults();
    
}
template Image2D<short>::Image2D();
template Image2D<int>::Image2D();
template Image2D<long>::Image2D();
template Image2D<float>::Image2D();
template Image2D<double>::Image2D();


template <class Type> 
Image2D<Type>::~Image2D () {
    
    if (arrayAllocated) delete [] array;
}
template Image2D<short>::~Image2D();
template Image2D<int>::~Image2D();
template Image2D<long>::~Image2D();
template Image2D<float>::~Image2D();
template Image2D<double>::~Image2D();


template <class Type>
Image2D<Type>::Image2D(int *dim)  {
    
    defaults();
    int size = dim[0]*dim[1];
    arrayAllocated = false;
    numPix = numAxes = 0;
    if(size<0)
        std::cout << "Error [Image2D(dimArray)]: Negative size -- could not define Image"<<std::endl;
    else {
        numPix = size;
        if(size>0){
            array = new Type[size];
            arrayAllocated = true;
        }
        numAxes  = 2;
        for(int i=0; i<numAxes; i++) axisDim[i] = dim[i];

    }
}
template Image2D<short>::Image2D(int*);
template Image2D<int>::Image2D(int*);
template Image2D<long>::Image2D(int*);
template Image2D<float>::Image2D(int*);
template Image2D<double>::Image2D(int*);


template <class Type>
Image2D<Type>::Image2D(const Image2D<Type> &i) {
    
    this->operator=(i);
}
template Image2D<short>::Image2D(const Image2D<short>&);
template Image2D<int>::Image2D(const Image2D<int>&);
template Image2D<long>::Image2D(const Image2D<long>&);
template Image2D<float>::Image2D(const Image2D<float>&);
template Image2D<double>::Image2D(const Image2D<double>&);


template <class Type>
Image2D<Type>& Image2D<Type>::operator=(const Image2D<Type> &i) {

    if(this==&i) return *this;
   
    this->numPix    = i.numPix;
    this->numAxes   = i.numAxes;
    this->par       = i.par;
   
    for (int j=0; j<numAxes; j++) this->axisDim[j] = i.axisDim[j];
    
    if(this->arrayAllocated) delete [] array; 
    this->arrayAllocated = i.arrayAllocated;
    if(this->arrayAllocated) {
        this->array = new Type[this->numPix];
        for(int j=0; j<numPix; j++) this->array[j] = i.array[j];
    }
    
    this->headDefined = i.headDefined;
    if (this->headDefined) this->head = i.head;
    this->statsDefined = i.statsDefined;
    if (this->statsDefined) this->stats = i.stats;
   
    return *this;
}
template Image2D<short>& Image2D<short>::operator=(const Image2D<short>&);
template Image2D<int>& Image2D<int>::operator=(const Image2D<int>&);
template Image2D<long>& Image2D<long>::operator=(const Image2D<long>&);
template Image2D<float>& Image2D<float>::operator=(const Image2D<float>&);
template Image2D<double>& Image2D<double>::operator=(const Image2D<double>&);

/**==============================================================================*/

template <class Type>
void Image2D<Type>::setImage(Type *input, int *dim) {

    if (arrayAllocated) delete [] array;
    numAxes = 2;
    for (int i=0; i<numAxes; i++) axisDim[i] = dim[i];
    numPix = axisDim[0]*axisDim[1];
    array = new Type [numPix];
    arrayAllocated=true;
    for (int y=0; y<axisDim[1]; y++) 
        for (int x=0; x<axisDim[0]; x++) {
                long nPix = x+y*axisDim[0];
                array[nPix]=input[nPix]; 
            }
    
}
template void Image2D<short>::setImage(short*,int*);
template void Image2D<int>::setImage(int*,int*);
template void Image2D<long>::setImage(long*,int*);
template void Image2D<float>::setImage(float*,int*);
template void Image2D<double>::setImage(double*,int*);


template <class Type>
bool Image2D<Type>::readImage(std::string fname) {
    
    par.setImageFile(fname);
    numAxes = 2;
    
    if(!head.header_read(par.getImageFile())) return false;
    headDefined = true;
    for (int i=0; i<numAxes; i++) axisDim[i] = head.DimAx(i);
    numPix = axisDim[0]*axisDim[1];
    if (!fitsread_2d()) return false;
    return true;
}
template bool Image2D<short>::readImage(std::string);
template bool Image2D<int>::readImage(std::string);
template bool Image2D<long>::readImage(std::string);
template bool Image2D<float>::readImage(std::string);
template bool Image2D<double>::readImage(std::string);



template <class Type>
bool Image2D<Type>::fitsread_2d () {

    fitsfile *fptr;
    int status, anynul, fpixel;
    float nulval;
    
    std::cout << "\nOpening file "<< par.getImageFile() << std::endl; 
    if (par.isVerbose()) { 
        std::cout << "Reading "<<axisDim[0]<<" x "<<axisDim[1]  << " pixels FITS file... ";
    }
    
    // Open the FITS file
    char *filename = new char [100];
    strcpy(filename, (par.getImageFile()).c_str());
    status = 0;
    if(fits_open_file(&fptr, filename, READONLY, &status) ){
      fits_report_error(stderr, status);
      return false;
    }
    
    // Read elements from the FITS data array    
    if (!arrayAllocated) array = new Type[numPix];                      
    arrayAllocated = true;
    
    fpixel=1;    
    status=0;
    if (fits_read_img(fptr, selectDatatype<Type>(), fpixel, numPix, &nulval, array, &anynul, &status)){
        fits_report_error(stderr, status);
        return false;
    }  
     
    // Close the FITS File
    if (fits_close_file(fptr, &status)){
        fits_report_error(stderr, status);
    }
    
    delete [] filename;
    
    if (par.isVerbose()) std::cout << "Done.\n" << std::endl;
    
    return true;
    
}
template bool Image2D<short>::fitsread_2d ();
template bool Image2D<int>::fitsread_2d ();
template bool Image2D<long>::fitsread_2d ();
template bool Image2D<float>::fitsread_2d ();
template bool Image2D<double>::fitsread_2d ();


template <class Type>
bool Image2D<Type>::fitswrite_2d (const char *outname) {
    
    fitsfile *fptr;
    int status;
    long firstPix = 1;
    long dnaxes[2] = {axisDim[0], axisDim[1]};
    
    remove (outname);
    
    //Create the new file
    status = 0;
    if (fits_create_file (&fptr, outname, &status)) {
        fits_report_error (stderr, status);
        return false;
        }
    
    //Create the primary array Image2D  
    if (fits_create_img (fptr, selectBitpix<Type>(), numAxes, dnaxes, &status)){
        fits_report_error (stderr, status);
        return false;
    }
    
    if (headDefined) head.headwrite_2d (fptr,false);
    
    if (fits_write_img (fptr, selectDatatype<Type>(), firstPix, numPix, array, &status)){
        fits_report_error (stderr, status);
        return false;
    }
    
    // Close the FITS File
    if (fits_close_file(fptr, &status)){
        fits_report_error(stderr, status);
    }
    
    return true;
    
}
template bool Image2D<short>::fitswrite_2d (const char*);
template bool Image2D<int>::fitswrite_2d (const char*);
template bool Image2D<long>::fitswrite_2d (const char*);
template bool Image2D<float>::fitswrite_2d (const char*);
template bool Image2D<double>::fitswrite_2d (const char*);



template <class Type>
void Image2D<Type>::copyHeader (Header &c) {

    headDefined = true;
    if (headDefined) {
        head.setNumAx(2);
        for (int i=0; i<numAxes; i++) {
            head.setDimAx(i, c.DimAx(i));
            head.setCdelt(i, c.Cdelt(i));
            head.setCrpix(i, c.Crpix(i));
            head.setCrval(i, c.Crval(i));
            head.setCunit(i, c.Cunit(i));
            head.setCtype(i, c.Ctype(i));
        }
        head.setBunit(c.Bunit());
        head.setBmaj(c.Bmaj());
        head.setBmin(c.Bmin());
        head.setBpa(c.Bpa());
        head.setEpoch(c.Epoch());
        head.setName(c.Name());
    }
}
template void Image2D<short>::copyHeader(Header&);
template void Image2D<int>::copyHeader(Header&);
template void Image2D<long>::copyHeader(Header&);
template void Image2D<float>::copyHeader(Header&);
template void Image2D<double>::copyHeader(Header&);


template <class Type>
void Image2D<Type>::setImageStats() {

  /// Calculates the full statistics for the Image2D: mean, rms, median, madfm.
  /// Also work out the threshold and store it in the stats set.

    if(par.isVerbose()) 
        std::cout << " Calculating statistics for the image... "<<std::flush;
    
    stats.setRobust(par.getFlagRobustStats()); 

    stats.calculate(array,numPix);
    
    stats.setUseFDR(false);

    stats.setThresholdSNR(par.getCut());
      
    if(par.isVerbose()) {
        std::cout << "Using ";
        if(stats.getUseFDR()) std::cout << "effective ";
        std::cout << "flux threshold of: ";
        Type thresh = stats.getThreshold();
        std::cout << thresh << " " << head.Bunit() << std::endl;
    }
    
    statsDefined = true;
    
}
template void Image2D<short>::setImageStats();
template void Image2D<int>::setImageStats();
template void Image2D<long>::setImageStats();
template void Image2D<float>::setImageStats();
template void Image2D<double>::setImageStats();


template <class Type>
void Image2D<Type>::extractSpectrum(Type *Array, int *dim, long pixel) {
    
  ///  A function to extract a 1-D spectrum from a 3-D array.
  ///  The array is assumed to be 3-D with the third dimension the spectral one.
  ///  The spectrum extracted is the one lying in the spatial pixel referenced
  ///  by the third argument.
  ///  The extracted spectrum is stored in the pixel array Image::array.
  ///
  ///  \param Array         The array containing the pixel values, from which
  ///                       the spectrum is extracted.
  ///  \param dim           The array of dimension values.
  ///  \param pixel         The spatial pixel that contains the desired spectrum.

    if((pixel<0)||(pixel>=dim[0]*dim[1]))
        std::cout << "Image::extractSpectrum: Requested spatial pixel outside allowed range. Cannot save."<<std::endl;
    else if(dim[2] != numPix)
        std::cout << "Image::extractSpectrum: Input array different size to existing array. Cannot save."<<std::endl;
    else {
        if(numPix>0 && arrayAllocated) delete [] array;
        numPix = dim[2];
        if(numPix>0){
            array = new Type[dim[2]];
            arrayAllocated = true;
            for(int z=0;z<dim[2];z++) array[z] = Array[z*dim[0]*dim[1] + pixel];
        }
    }
}
template void Image2D<short>::extractSpectrum(short*,int*,long);
template void Image2D<int>::extractSpectrum(int*,int*,long);
template void Image2D<long>::extractSpectrum(long*,int*,long);
template void Image2D<float>::extractSpectrum(float*,int*,long);
template void Image2D<double>::extractSpectrum(double*,int*,long);


template <class Type>
void Image2D<Type>::extractSpectrum(Cube<Type> &cube, long pixel)  {
    
  ///  A function to extract a 1-D spectrum from a Cube class
  ///  The spectrum extracted is the one lying in the spatial pixel referenced
  ///  by the second argument.
  ///  The extracted spectrum is stored in the pixel array Image::array.
  ///
  ///  \param cube      The Cube containing the pixel values, 
  ///                   from which the spectrum is extracted.
  ///  \param pixel     The spatial pixel that contains the desired spectrum.

    long zdim = cube.DimZ();
    long spatSize = cube.DimX()*cube.DimY();
    if((pixel<0)||(pixel>=spatSize))
        std::cout << "Image::extractSpectrum: Requested spatial pixel outside allowed range. Cannot save."<< std::endl;
    else if(zdim != numPix)
        std::cout << "Image::extractSpectrum: Input array different size to existing array. Cannot save."<< std::endl;
    else {
        if(numPix>0 && arrayAllocated) delete [] array;
        numPix = zdim;
        if(numPix>0){
            array = new Type[zdim];
            arrayAllocated = true;
        for(int z=0;z<zdim;z++) 
            array[z] = cube.Array(z*spatSize + pixel);
        }
    }
}
template void Image2D<short>::extractSpectrum(Cube<short>&,long);
template void Image2D<int>::extractSpectrum(Cube<int>&,long);
template void Image2D<long>::extractSpectrum(Cube<long>&,long);
template void Image2D<float>::extractSpectrum(Cube<float>&,long);
template void Image2D<double>::extractSpectrum(Cube<double>&,long);

template <class Type>
void Image2D<Type>::extractGlobalSpectrum(Cube<Type> *cube)  {
    
  ///  A function to extract a 1-D spectrum global spectrum from a Cube class
  ///
  ///  \param cube      The Cube containing the pixel values, 
  ///                   from which the spectrum is extracted.
  ///  \param pixel     The spatial pixel that contains the desired spectrum.

    long zdim = cube->DimZ();
    if(numPix>0 && arrayAllocated) delete [] array;
    numPix = zdim;
    if(numPix>0){
        array = new Type[zdim];
        arrayAllocated = true;
        for(int z=0;z<zdim;z++) { 
            array[z] = 0;
            for(int y=0;y<cube->DimY();y++) 
                for(int x=0;x<cube->DimX();x++) array[z] += cube->Array(cube->nPix(x,y,z));
        }
    }
}
template void Image2D<short>::extractGlobalSpectrum(Cube<short>*);
template void Image2D<int>::extractGlobalSpectrum(Cube<int>*);
template void Image2D<long>::extractGlobalSpectrum(Cube<long>*);
template void Image2D<float>::extractGlobalSpectrum(Cube<float>*);
template void Image2D<double>::extractGlobalSpectrum(Cube<double>*);

template <class Type>
void Image2D<Type>::extractImage(Type *Array, int *dim, long channel) {
    
  ///  A function to extract a 2-D image from a 3-D array.
  ///  The array is assumed to be 3-D with the third dimension the spectral one.
  ///  The dimensions of the array are in the dim[] array.
  ///  The image extracted is the one lying in the channel referenced
  ///  by the third argument.
  ///  The extracted image is stored in the pixel array Image::array.
  ///
  ///  \param           Array The array containing the pixel values, 
  ///                   from which the image is extracted.
  ///  \param dim       The array of dimension values.
  ///  \param channel   The spectral channel that contains the desired image.

    long spatSize = dim[0]*dim[1];
    if((channel<0)||(channel>=dim[2]))
        std::cout<<"Image::extractImage: Requested channel outside allowed range. Cannot save."<<std::endl;
    else if(spatSize != numPix)
        std::cout<<"Image::extractImage: Input array different size to existing array. Cannot save."<<std::endl;
    else {
        if(numPix>0 && arrayAllocated) delete [] array;
        numPix = spatSize;
        if(numPix>0){
            array = new Type[spatSize];
            arrayAllocated = true;
            for(int npix=0; npix<spatSize; npix++)
                array[npix] = Array[channel*spatSize + npix];
        }
    }
}
template void Image2D<short>::extractImage(short*,int*,long);
template void Image2D<int>::extractImage(int*,int*,long);
template void Image2D<long>::extractImage(long*,int*,long);
template void Image2D<float>::extractImage(float*,int*,long);
template void Image2D<double>::extractImage(double*,int*,long);


template <class Type>
void Image2D<Type>::extractImage(Cube<Type> &cube, long channel) {
    
  ///  A function to extract a 2-D image from Cube class.
  ///  The image extracted is the one lying in the channel referenced
  ///  by the second argument.
  ///  The extracted image is stored in the pixel array Image::array.
  ///   
  ///  \param cube      The Cube containing the pixel values, 
  ///                   from which the image is extracted.
  ///  \param channel   The spectral channel that contains the desired image.
    
    long spatSize = cube.DimX()*cube.DimY();
    if((channel<0)||(channel>=cube.DimZ()))
        std::cout<<"Image::extractImage: Requested channel outside allowed range. Cannot save."<<std::endl;
    else if(spatSize != numPix)
        std::cout<<"Image::extractImage: Input array different size to existing array. Cannot save."<<std::endl;
    else {
        if(numPix>0 && arrayAllocated) delete [] array;
        numPix = spatSize;
        if(numPix>0){
            array = new Type[spatSize];
            arrayAllocated = true;
        for(int npix=0; npix<spatSize; npix++) 
            array[npix] = cube.Array(channel*spatSize + npix);
        }
    }
}
template void Image2D<short>::extractImage(Cube<short>&,long);
template void Image2D<int>::extractImage(Cube<int>&,long);
template void Image2D<long>::extractImage(Cube<long>&,long);
template void Image2D<float>::extractImage(Cube<float>&,long);
template void Image2D<double>::extractImage(Cube<double>&,long);


template <class Type>
std::vector<PixelInfo::Scan<Type> > Image2D<Type>::findSources1D() {
    
    std::vector<bool> thresholdedArray(axisDim[0]);
        for(int posX=0;posX<axisDim[0];posX++){
            thresholdedArray[posX] = isDetection(posX,0);
        }
    return spectrumDetect(thresholdedArray);
}
template std::vector<PixelInfo::Scan<short> > Image2D<short>::findSources1D();
template std::vector<PixelInfo::Scan<int> > Image2D<int>::findSources1D();
template std::vector<PixelInfo::Scan<long> > Image2D<long>::findSources1D();
template std::vector<PixelInfo::Scan<float> > Image2D<float>::findSources1D();
template std::vector<PixelInfo::Scan<double> > Image2D<double>::findSources1D();


template <class Type>
std::vector<Object2D<Type> > Image2D<Type>::findSources2D() {
    
    std::vector<bool> thresholdedArray(axisDim[0]*axisDim[1]);
    for(int posY=0;posY<axisDim[1];posY++) {
        for(int posX=0;posX<axisDim[0];posX++) {
            int loc = posX + axisDim[0]*posY;
                thresholdedArray[loc] = isDetection(posX,posY);
        }
    }
    return imageDetect(thresholdedArray);
}
template std::vector<Object2D<short> > Image2D<short>::findSources2D();
template std::vector<Object2D<int> > Image2D<int>::findSources2D();
template std::vector<Object2D<long> > Image2D<long>::findSources2D();
template std::vector<Object2D<float> > Image2D<float>::findSources2D();
template std::vector<Object2D<double> > Image2D<double>::findSources2D(); 


template <class Type>  
bool Image2D<Type>::isDetection(long x, long y) {
      
      /// Test whether a pixel (x,y) is a statistically
      /// significant detection, according to the set of statistics in
      /// the local Stats object.

      long voxel = y*axisDim[0] + x;
      return stats.isDetection(array[voxel]);
}
template bool Image2D<short>::isDetection(long,long);
template bool Image2D<int>::isDetection(long,long);
template bool Image2D<long>::isDetection(long,long);
template bool Image2D<float>::isDetection(long,long);
template bool Image2D<double>::isDetection(long,long); 


template <class Type>
std::vector<PixelInfo::Scan<Type> > Image2D<Type>::spectrumDetect(std::vector<bool> &arraybool) {
    
  ///  A detection algorithm that searches in a single 1-D spectrum.  It
  ///  simply scans along the spectrum, storing connected sets of
  ///  detected pixels as Scans, where "detected" means according to the
  ///  Image::isDetection(long,long) function.
  /// 
  ///  When finished a vector of the detected scans is returned.
    
    enum STATUS { NONOBJECT, OBJECT };
    STATUS status;
    PixelInfo::Scan<Type> obj;
    std::vector<PixelInfo::Scan<Type> > outputlist;
    bool isObject;
    long dim = axisDim[0];
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
template std::vector<PixelInfo::Scan<short> > Image2D<short>::spectrumDetect(std::vector<bool>&);
template std::vector<PixelInfo::Scan<int> > Image2D<int>::spectrumDetect(std::vector<bool>&);
template std::vector<PixelInfo::Scan<long> > Image2D<long>::spectrumDetect(std::vector<bool>&);
template std::vector<PixelInfo::Scan<float> > Image2D<float>::spectrumDetect(std::vector<bool>&);
template std::vector<PixelInfo::Scan<double> > Image2D<double>::spectrumDetect(std::vector<bool>&);


template <class Type>
std::vector<Object2D<Type> > Image2D<Type>::imageDetect(std::vector<bool> &arraybool) {

    ///  A detection algorithm for 2-dimensional images based on that of
    ///  Lutz (1980).
    ///  
    ///  The image is raster-scanned, and searched row-by-row. Objects
    ///  detected in each row are compared to objects in subsequent rows,
    ///  and combined if they are connected (in an 8-fold sense).
    /// 
    
    
    
    long xdim = axisDim[0];
    long ydim = axisDim[1];
    std::vector<Object2D<Type> > outputlist;
    STATUS *status  = new STATUS[2];
    Object2D<Type> *store = new Object2D<Type>[xdim+1];
    char *marker    = new char[xdim+1];
    for(int i=0; i<(xdim+1); i++) marker[i] = NULLMARKER;
    std::vector<FoundObject<Type> > oS;
    std::vector<STATUS>      psS;
    

    Pixel<Type> pix;
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
template std::vector<Object2D<short> > Image2D<short>::imageDetect(std::vector<bool>&);
template std::vector<Object2D<int> > Image2D<int>::imageDetect(std::vector<bool>&);
template std::vector<Object2D<long> > Image2D<long>::imageDetect(std::vector<bool>&);
template std::vector<Object2D<float> > Image2D<float>::imageDetect(std::vector<bool>&);
template std::vector<Object2D<double> > Image2D<double>::imageDetect(std::vector<bool>&); 


