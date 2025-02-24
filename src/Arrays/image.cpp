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


template <class Type> 
Image2D<Type>::~Image2D () {
    
    if (arrayAllocated) delete [] array;
}


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


template <class Type>
Image2D<Type>::Image2D(const Image2D<Type> &i) {
    
    this->operator=(i);
}


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

/**==============================================================================*/

template <class Type>
void Image2D<Type>::setImage(int *dim) {

    if (arrayAllocated) delete [] array;
    numAxes = 2;
    for (int i=0; i<numAxes; i++) axisDim[i] = dim[i];
    numPix = axisDim[0]*axisDim[1];
    array = new Type [numPix];
    arrayAllocated=true;
    head.setNumAx(numAxes);
}


template <class Type>
void Image2D<Type>::setImage(Type *input, int *dim) {

    setImage(dim);
    for (int i=axisDim[0]*axisDim[1]; i--;) array[i]=input[i];

}


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


template <class Type>
bool Image2D<Type>::fitswrite_2d (const char *outname, bool fullHead) {
    
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
    
    if (headDefined) head.headwrite(fptr,2,fullHead);
    
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


template <class Type>
void Image2D<Type>::copyHeader (Header &h) {

    headDefined = true;
    if (headDefined) {
        head.setNumAx(2);
        for (int i=0; i<numAxes; i++) {
            head.setDimAx(i, h.DimAx(i));
            head.setCdelt(i, h.Cdelt(i));
            head.setCrpix(i, h.Crpix(i));
            head.setCrval(i, h.Crval(i));
            head.setCunit(i, h.Cunit(i));
            head.setCtype(i, h.Ctype(i));
        }
        head.setBunit(h.Bunit());
        head.setBmaj(h.Bmaj());
        head.setBmin(h.Bmin());
        head.setBpa(h.Bpa());
        head.setEpoch(h.Epoch());
        head.setName(h.Name());
        head.setTelesc(h.Telesc());
        head.setFreq0(h.Freq0());
        head.setRaDeSys(h.RaDeSys());
        head.setSpecSys(h.SpecSys());
    }
}


template <class Type>
void Image2D<Type>::setImageStats() {

  /// Calculates the full statistics for the Image2D: mean, rms, median, madfm.
  /// Also work out the threshold and store it in the stats set.

    if(par.isVerbose()) 
        std::cout << " Calculating statistics for the image... "<<std::flush;
    
    stats.setRobust(par.getFlagRobustStats()); 
    stats.calculate(array,numPix);
    stats.setUseFDR(false);
    stats.setThresholdSNR(par.getParSE().snrCut);
      
    if(par.isVerbose()) {
        std::cout << "Using ";
        if(stats.getUseFDR()) std::cout << "effective ";
        std::cout << "flux threshold of: ";
        Type thresh = stats.getThreshold();
        std::cout << thresh << " " << head.Bunit() << std::endl;
    }
    
    statsDefined = true;
    
}


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


// Explicit instantiation of the class
template class Image2D<short>;
template class Image2D<int>;
template class Image2D<long>;
template class Image2D<float>;
template class Image2D<double>;

