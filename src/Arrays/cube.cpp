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
    Internet email: enrico.diteodoro@gmail.com
-----------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <Arrays/cube.hh>
#include <Arrays/image.hh>
#include <Arrays/stats.hh>
#include <Map/detection.hh>
#include <Utilities/utils.hh>
#include <Utilities/progressbar.hh>
#include <Utilities/gnuplot.hh>
#include <Tasks/smooth3D.hh>
#include <Tasks/moment.hh>


template <class T>
void Cube<T>::defaults() {
    
    numAxes = 3;
    arrayAllocated = false;
    headDefined = false; 
    axisDimAllocated = false; 
    statsDefined = false;
    maskAllocated = false;
    isSearched = false;
    
}


template <class T>
Cube<T>::~Cube () {
    
    if (arrayAllocated) delete [] array;
    arrayAllocated=false;
    if (maskAllocated)  delete [] mask;
    maskAllocated=false;
    if (axisDimAllocated) delete [] axisDim;
    axisDimAllocated=false;
    if (isSearched) delete sources;
} 


template <class T>
Cube<T>::Cube(std::string fname) {
    
    defaults();
    par.setImageFile(fname);
    this->readCube(fname);
}


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
        }
        numAxes  = 3;
        axisDim = new int[numAxes];
        axisDimAllocated = true;
        for(int i=0; i<numAxes; i++) axisDim[i] = dimensions[i];   
    }
    
}


template <class T>
Cube<T>::Cube(const Cube<T> &c) {
  
    defaults();
    this->operator=(c);

}


template <class T>
Cube<T>& Cube<T>::operator=(const Cube<T> &c) {
    
    if(this==&c) return *this;
    
    if(this->arrayAllocated) delete [] array;
    if(this->axisDimAllocated) delete [] axisDim;
    if(this->maskAllocated) delete [] mask;
    
    this->numPix    = c.numPix;
    this->numAxes   = c.numAxes;
    this->par       = c.par;
    
    this->axisDimAllocated = c.axisDimAllocated;
    if (axisDimAllocated) {
        this->axisDim = new int [numAxes];
        for (short i=0; i<numAxes;  i++) this->axisDim[i] = c.axisDim[i];
    }
    
    this->arrayAllocated = c.arrayAllocated;
    if(this->arrayAllocated) {
        this->array = new T[this->numPix];
        for(size_t i=0; i<this->numPix; i++) this->array[i] = c.array[i];
    }
    
    this->maskAllocated = c.maskAllocated; 
    if(this->maskAllocated) {
        this->mask = new bool[numPix];
        for(size_t i=0;i<this->numPix;i++) this->mask[i] = c.mask[i];
    }
    
    this->headDefined = c.headDefined;
    if (this->headDefined) this->head = c.head;
    this->statsDefined = c.statsDefined;
    if (this->statsDefined) this->stats = c.stats;
    this->isSearched = c.isSearched;
    if (this->isSearched)
        *this->sources = *c.sources;

    return *this;
}


template <class T>
void Cube<T>::checkBeam() {
    
    // If beam size is not defined in the header
    if (head.BeamArea()==0) {
        // Try if BMAJ, BMIN, BPA parameters are defined (in arcs)
        float bmaj = par.getBmaj()/3600.;
        float bmin = par.getBmin()/3600.;
        float bpa  = par.getBpa();
        
        if (bmaj>0 && bmin>0);
        else if (bmaj>0 && bmin<0) bmin = bmaj;
        else bmaj = bmin = par.getBeamFWHM();   // Try BEAMFWHM (is already in deg)
        
        std::cout << "\n WARNING: Beam not available in the header: using a "
                  << bmaj*3600. << "x" << bmin*3600. << " arcsec (BPA=" << bpa 
                  << "). You can set the beam with BMAJ/BMIN/BPA or BeamFWHM params (in arcsec).\n\n";
        head.setBmaj(bmaj);
        head.setBmin(bmin);
        head.setBpa(bpa);
        head.calcArea();
    }
}


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
    for (size_t i=0; i<numPix; i++) array[i]=input[i]; 

}


template <class T>
bool Cube<T>::readCube (std::string fname) {
    
    par.setImageFile(fname);
    numAxes = 3;

    if(!head.header_read(par.getImageFile())) return false;

    headDefined = true;
    // I do not like the following 3 lines, I should think something better
    head.setRedshift(par.getRedshift());
    head.setWave0(par.getRestwave());
    // The FREQ0 can be in the header. Override if given by the user.
    if (par.getRestfreq()!=-1) head.setFreq0(par.getRestfreq());
    axisDim = new int [numAxes];
    axisDimAllocated = true; 
    for (short i=0; i<numAxes; i++) axisDim[i] = head.DimAx(i);
    if (head.NumAx()<3) axisDim[2] = 1;
    if (head.NumAx()<2) axisDim[1] = 1;
    numPix = axisDim[0]*axisDim[1]*axisDim[2];
    if (!fitsread_3d()) return false;
    return true;
}


template <class T>
bool Cube<T>::fitsread_3d() {

    fitsfile *fptr3;
    int status, anynul, fpixel;

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
    if (fits_read_img(fptr3, selectDatatype<T>(), fpixel, numPix, NULL, array, &anynul, &status)){
        fits_report_error(stderr, status);
        return false;
    } 
    
    // Close the FITS File
    if (fits_close_file(fptr3, &status)){
        fits_report_error(stderr, status);
    }

    ///////// I need to solve the problem with NANs in a better way ///////////////
    for (size_t i=numPix; i--;)
        if (isNaN(array[i])) array[i] = 0;
    ///////////////////////////////////////////////////////////////////////////////

    if (par.isVerbose()) std::cout << "Done.\n" << std::endl;

    return true;
}


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


/**=====================================================================================*/ 
/** STATISTICAL FUNCTIONS */

template <class T>
void Cube<T>::setCubeStats() {

  /// Calculates the full statistics for the cube: mean, rms, median, madfm.
  /// Also work out the threshold and store it in the stats set.

    if(par.isVerbose()) std::cout << " Calculating statistics for the cube... "<<std::flush;
    
    stats.setRobust(par.getFlagRobustStats());
    bool *blanks = new bool[numPix];
    for (size_t i=0; i<numPix; i++) blanks[i] = isBlank(array[i]) ? false : true;

    stats.calculate(array,numPix,blanks);
    stats.setThresholdSNR(par.getParSE().snrCut);

    if(par.isVerbose()) {
        std::cout << "Using flux threshold of: ";
        T thresh;
        if (par.getParSE().UserThreshold) thresh = par.getParSE().threshold;
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


template <class T>
void Cube<T>::BlankCube (T *Array, size_t size) {
    
  /// A function for blanking an array with Cube mask data.
  ///
  /// \param Array      The array to blank.
  /// \param size       The size of array. It must be equal
  ///                   to the size of Cube-object array.

    if (size!=numPix) 
        std::cout << "Error blanking cube: array size is different from cube size" << std::endl;
    else {
        if (!maskAllocated) BlankMask();
        for (size_t i=0; i<size; i++) Array[i] *= mask[i];
    }
}


template <class T>
void Cube<T>::BlankMask (float *channel_noise, bool onlyLargest){
    
     ///////////////////////////////////////////////////////////////////////////////////
     /// This function builds a mask for the cube. The type of mask depends on the
     /// parameter MASK:
     ///
     /// - SEARCH:      Uses the source finding algorith and masks the largest object.
     ///                All search-related parameters can be used.
     /// - THRESHOLD:   Applies a simple threshold cut given by THRESHOLD parameter.
     /// - SMOOTH:      Smooths the cube by a factor FACTOR and applies a S/N threshold
     ///                on the smoothed cube given by BLANKCUT. Default are FACTOR=2
     ///                and BLANKCUT=3. If BMAJ and BMIN parameters are, it smooths to
     ///                these values.
     /// - NEGATIVE:    Calculates the noise statitistics just on the negative pixels
     ///                and builds the mask based on the S/N threshold BLANKCUT.
     /// - FILE(Name):  User-provided mask. 'Name' is a fitsfile with same size of the
     ///                cube and filled with 0(false) or 1(true).
     /// - NONE:        No mask, just pixels > 0.
     ///
     ///////////////////////////////////////////////////////////////////////////////////

    if (maskAllocated) delete [] mask;
    mask = new bool[numPix];

    for (size_t i=0; i<numPix; i++) mask[i]=0;

    bool verb = par.isVerbose();
    if (verb) {
        std::cout << " Creating mask (" << par.getMASK() << ") ..." << std::flush;
        par.setVerbosity(false);
    }

    Statistics::Stats<T> *st = new Statistics::Stats<T>;
    st->setRobust(par.getFlagRobustStats());
    
    if (par.getMASK().find("SMOOTH&SEARCH")!=std::string::npos) {
        // Smoothing first and searching for the largest object
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
        
        Cube<T> *smoothed = new Cube<T>();
        smoothed->setCube(sm->Array(),axisDim);
        smoothed->saveHead(head);
        smoothed->saveParam(par);
        smoothed->Head().setBmaj(nbmaj/3600.);
        smoothed->Head().setBmin(nbmin/3600.);
        smoothed->Head().setBpa(nbpa);
        smoothed->Head().calcArea();
        smoothed->setCubeStats();
        smoothed->search();
        size_t numObj = smoothed->getNumObj();
        if (numObj==0) {
            std::cout << "MASKING error: No sources detected in the datacube. Cannot build mask!!! \n";
            std::terminate();
        }
        
        if (onlyLargest || numObj==1) {
            Detection<T> *larg = smoothed->getSources()->LargestDetection();
            std::vector<Voxel<T> > voxlist = larg->getPixelSet();
            typename std::vector<Voxel<T> >::iterator vox;
            for(vox=voxlist.begin();vox<voxlist.end();vox++) {
                mask[nPix(vox->getX(),vox->getY(),vox->getZ())]=1;
            }
        }
        else {
            for (size_t i=0; i<numObj; i++) {
                Detection<T> *obj = smoothed->pObject(i);
                std::vector<Voxel<T> > voxlist = obj->getPixelSet();
                typename std::vector<Voxel<T> >::iterator vox;
                for(vox=voxlist.begin();vox<voxlist.end();vox++) {
                    mask[nPix(vox->getX(),vox->getY(),vox->getZ())]=1;
                }
            }
            
        }
        delete smoothed;
        delete sm;
    }
    else if (par.getMASK().find("SEARCH")!=std::string::npos) {
        // Masking using the search algorithm.
        if (!isSearched) search();
        
        uint numObj = getNumObj();
        if (numObj==0) {
            std::cout << "MASKING error: No sources detected in the datacube. Cannot build mask!!! \n";
            std::terminate();
        }
        
        if (onlyLargest || numObj==1) {
            Detection<T> *larg = sources->LargestDetection();
            std::vector<Voxel<T> > voxlist = larg->getPixelSet();
            typename std::vector<Voxel<T> >::iterator vox;
            for(vox=voxlist.begin();vox<voxlist.end();vox++) {
                mask[nPix(vox->getX(),vox->getY(),vox->getZ())]=1;
            }
        }
        else {
            for (size_t i=0; i<numObj; i++) {
                Detection<T> *obj = pObject(i);
                std::vector<Voxel<T> > voxlist = obj->getPixelSet();
                typename std::vector<Voxel<T> >::iterator vox;
                for(vox=voxlist.begin();vox<voxlist.end();vox++) {
                    mask[nPix(vox->getX(),vox->getY(),vox->getZ())]=1;
                }
            }
            
        }
    }
    else if (par.getMASK()=="SMOOTH") {
        // Smooth and cut
        double bmaj  = head.Bmaj()*arcsconv(head.Cunit(0));
        double bmin  = head.Bmin()*arcsconv(head.Cunit(0));
        double bpa   = head.Bpa();
        float factor = par.getFactor()==-1 ? 2 : par.getFactor();
        double nbmaj = par.getBmaj()==-1 ? factor*bmaj : par.getBmaj();
        double nbmin = par.getBmin()==-1 ? factor*bmin : par.getBmin();
        double nbpa  = par.getBpa()==-1  ? bpa    : par.getBpa();   
        //if (nbmaj/bmaj<1.1) nbmaj = factor*bmaj;
        //if (nbmin/bmin<1.1) nbmin = factor*bmin;
        Beam oldbeam = {bmaj,bmin,bpa};
        Beam newbeam = {nbmaj,nbmin,nbpa};
                
        Smooth3D<T> *sm = new Smooth3D<T>;
        sm->smooth(this, oldbeam, newbeam);
        bool *blanks = new bool[numPix];
        for (size_t i=0; i<numPix; i++) blanks[i] = isBlank(sm->Array(i)) ? false : true;
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
    else if (par.getMASK()=="THRESHOLD") {
        // Simple cut
        float thresh = par.getParSE().threshold;
        for (uint i=numPix; i--;) {
            if (array[i]>thresh) mask[i] = 1;
        }
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
    else if (par.getMASK().find("FILE(")!=std::string::npos) {
        std::string str = par.getMASK();
        size_t first = str.find_first_of("(");
        size_t last = str.find_last_of(")");
        if (first==std::string::npos || last==std::string::npos) {
            std::cerr << "\n  ERROR: MASK parameter is not correct. Provide file(Maskfitsfile)\n";
            std::terminate();
        }
        std::string filename = str.substr (first+1,last-first-1);
        Cube<float> *ma = new Cube<float>;

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
            //if (array[i]>0) 
                mask[i] = 1;
        }
    }

    delete st;
    
    // Writing mask to FITS file
    Cube<short> *m = new Cube<short>(axisDim);
    m->saveHead(head);
    m->saveParam(par);
    m->Head().setMinMax(0.,0);
    for (size_t i=numPix; i--;) m->Array(i) = short(mask[i]);
    m->fitswrite_3d((par.getOutfolder()+"mask.fits").c_str());
    delete m;

    if (verb) {
        std::cout << " Done." << std::endl;
        par.setVerbosity(true);
    }

    maskAllocated = true;
    
}


/**=====================================================================================*/ 


template <class T>
Cube<T>* Cube<T>::Reduce (int fac, std::string rtype) {
    
    /// This function reduces the size of a cube by averaging pixels/channels
    ///
    /// \param fac      Reduction factor. 
    /// \param rtype    "spectral" or "spatial" averaging
    
    if (par.isVerbose()) std::cout << " Reducing..." << std::flush;
    
    // Defining dimensions of the output cube
    int dim[3];
    if (rtype=="spectral") {
        dim[0] = axisDim[0];
        dim[1] = axisDim[1];
        dim[2] = axisDim[2]/fac;
    }
    else {
        dim[0] = axisDim[0]/fac;
        dim[1] = axisDim[1]/fac;
        dim[2] = axisDim[2];
    }

    // Initializing the new datacube 
    Cube<T> *reduced = new Cube<T>(dim);
    reduced->saveParam(par);
    reduced->saveHead(head);
    
    if (rtype=="spectral") {                // Spectral reduction
        for (int i=0; i<dim[0]*dim[1]; i++) {
            for (int z=0; z<dim[2]; z++) {
                T sum = 0;
                for (int zz=fac*z; zz<fac*z+fac; zz++)
                    sum += array[i+zz*axisDim[0]*axisDim[1]];
                reduced->Array(i+z*dim[0]*dim[1]) = sum/double(fac);
            }
        }
        
        reduced->Head().setCdelt(2, fac*head.Cdelt(2));
        reduced->Head().setCrpix(2, lround(head.Crpix(2)/double(fac)));
        
        std::string ochsize = "  Old channel width: "+to_string(Head().Cdelt(2))+" "+Head().Cunit(2);
        std::string nchsize = "  New channel width: "+to_string(reduced->Head().Cdelt(2))+" "+Head().Cunit(2);
        reduced->Head().addKey("HISTORY BBAROLO SPECTRAL AVERAGING: "+nchsize);
        reduced->Head().addKey("HISTORY BBAROLO SPECTRAL AVERAGING: "+ochsize);
        if (par.isVerbose()) std::cout << "OK\n" << ochsize << std::endl << nchsize << std::endl;
        
    }
    else {                                  // Spatial reduction
        for (int z=0; z<dim[2]; z++) {
            for (int y=0; y<dim[1]; y++) {
                for (int x=0; x<dim[0]; x++) {
                    T sum = 0;
                    for (int yy=fac*y; yy<fac*y+fac; yy++)
                        for (int xx=fac*x; xx<fac*x+fac; xx++) 
                            sum += array[xx+yy*axisDim[0]+z*axisDim[0]*axisDim[1]];
                    reduced->Array(x,y,z) = sum/double(fac*fac);
                }
            }
        }
    
        reduced->Head().setCdelt(0, fac*head.Cdelt(0));
        reduced->Head().setCdelt(1, fac*head.Cdelt(1));
        reduced->Head().setCrpix(0, lround(head.Crpix(0)/double(fac)));
        reduced->Head().setCrpix(1, lround(head.Crpix(1)/double(fac)));
        reduced->Head().calcArea();
    
        std::string obeamsize = "  Old beam size: "+to_string(Head().BeamArea())+" pixels";
        std::string nbeamsize = "  New beam size: "+to_string(reduced->Head().BeamArea())+" pixels";
        reduced->Head().addKey("HISTORY BBAROLO SPATIAL AVERAGING: "+nbeamsize);
        reduced->Head().addKey("HISTORY BBAROLO SPATIAL AVERAGING: "+obeamsize);
        if (par.isVerbose()) std::cout << "OK\n" << obeamsize << std::endl << nbeamsize << std::endl;
    }
    
    T minn,maxx;
    findMinMax<T>(reduced->Array(), reduced->NumPix(), minn, maxx);
    reduced->Head().setDataMax(double(maxx));
    reduced->Head().setDataMin(double(minn));
    
    return reduced;
    
}


template <class T>
void Cube<T>::CheckChannels () {
    
    int xySize = axisDim[0]*axisDim[1];
    int zdim = axisDim[2];
    
    if (!statsDefined) setCubeStats();
    
    std::cout << "\n\n";
    
    ProgressBar bar(true,par.isVerbose(),par.getShowbar());
    bar.init(" Checking for bad channels... ",zdim);
    
    std::vector<int> bad;
    
    for (int z=0; z<axisDim[2]; z++) {
        bar.update(z+1);
        T median = findMedian<T>(&array[z*xySize], xySize, false);
        T MADFM = findMADFM<T>(&array[z*xySize], xySize, median, false);
        if (MADFM > 1.5*stats.getMadfm()) {
            bad.push_back(z);
            for (int pix=0; pix<xySize; pix++) array[z*xySize+pix] = 0.;
        } 
    }
    
    
    if (bad.size()==0) {
        bar.fillSpace(" All channels are good.\n\n");
        return;
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
        std::cout << " Done. New threshold is " << setprecision(5) << stats.getThreshold()
                  << " " << head.Bunit() << std::endl << std::endl;
    }
    
    string name = par.getImageFile();
    int found = name.find(".fits");
    name.insert(found, "ean");
    
    int NdatZ = axisDim[2]-bad.size(); 
    int ax[3] = {axisDim[0], axisDim[1], NdatZ};
    Cube<T> *out = new Cube<T>(ax);
    for (size_t i=0; i<out->NumPix(); i++) 
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
    
}


///=====================================================================================
/// SOURCE-FINDING RELATED FUNCTIONS
///=====================================================================================
template <class T>
void Cube<T>::search() {

    Header &h = Head();
    SEARCH_PAR p = par.getParSE();

    if (!statsDefined) setCubeStats();
    checkBeam();

    float PixScale = (fabs(h.Cdelt(0))+fabs(h.Cdelt(1)))/2.;
    int thresS  = p.threshSpatial!=-1  ? p.threshSpatial     : ceil(h.Bmaj()/PixScale);
    int thresV  = p.threshVelocity!=-1 ? p.threshVelocity   : 3;
    int minchan = p.minChannels!=-1 ? p.minChannels : 2;
    int minpix  = p.minPix!=-1      ? p.minPix      : ceil(h.BeamArea());
    int minvox  = p.minVoxels!=-1   ? p.minVoxels   : minchan*minpix;
    p.threshSpatial = thresS;
    p.threshVelocity = thresV;
    p.minChannels = minchan;
    p.minPix = minpix;
    p.minVoxels = minvox;
    p.maxAngSize = p.maxAngSize/(head.PixScale()*arcsconv(head.Cunit(1))/60.);

    sources = new Search<T>(p);
    sources->search(array,stats,axisDim[0],axisDim[1],axisDim[2],
                    par.getFlagRobustStats(),par.getThreads(),par.isVerbose(),par.getShowbar());

    isSearched = true;

}


template <class T>
void Cube<T>::search(std::string searchtype, float snrCut, float threshold, bool adjacent, int threshSpatial,
                     int threshVelocity, int minPixels, int minChannels, int minVoxels, int maxChannels,
                     float maxAngSize, bool flagGrowth, float growthCut, float growthThreshold, bool RejectBefore,
                     bool TwoStage, int NTHREADS) {

    SEARCH_PAR &p       = par.getParSE();
    p.searchType        = searchtype;
    p.snrCut            = snrCut;
    p.threshold         = threshold;
    p.flagAdjacent      = adjacent;
    p.threshSpatial     = threshSpatial;
    p.threshVelocity    = threshVelocity;
    p.minPix            = minPixels;
    p.minChannels       = minChannels;
    p.minVoxels         = minVoxels;
    p.maxChannels       = maxChannels;
    p.maxAngSize        = maxAngSize;
    p.flagGrowth        = flagGrowth;
    p.growthCut         = growthCut;
    p.growthThreshold   = growthThreshold;
    p.RejectBeforeMerge = RejectBefore;
    p.TwoStageMerging   = TwoStage;
    if (p.threshold!=0) p.UserThreshold = true;
    if (p.growthThreshold!=0) p.flagUserGrowthT = true;
    par.setThreads(NTHREADS);

    setCubeStats();
    search();

}


template <class T>
void Cube<T>::printDetections (std::ostream& Stream) {

    using namespace std;

    float RA=-1, DEC=-1, VEL=-1;
    int numObj = getNumObj();
    int m = 10;
    int k=29;
    string str;

    Stream  << showpoint << fixed;
    Stream  << "#" << endl;

    if (headDefined) {
        Stream  << "# Detections for " << head.Name() << " " << endl;
        Stream  << "#" << setw(150) << setfill('_') << " " << endl << "#" << endl;
    }
    else {
        Stream  << "# Detections for " << par.getImageFile() << " " << endl;

        Stream  << "#" << setw(104) << setfill('_') << " " << endl << endl;
    }

    Stream  << setfill(' ');

    Stream  << "#" << setw(m-2) << left << "  Source"
            << setw(m+6)  << right << "Center    ";

    if (headDefined) Stream << setw(k) << right << "Center (WCS)       ";

    Stream  << setw(m-3) << right << "Xwidth"
            << setw(m-3) << right << "Ywidth"
            << setw(m-3) << right << "Zwidth"
            << setw(m) << right << "NumPix"
            << setw(m) << right << "NumVox"
            << setw(m) << right << "Vsys "
            << setw(m) << right << "W20  "
            << setw(m+5) << right << "Flux "
            << setw(m-2) << right << "PeakSNR"
            << setw(m) << right << "PeakFlux";

    Stream  << endl;

    Stream  << "#" << setw(m-2)<< left  << "    [#]"
            << setw(m+6)  << right << "[pix,pix,chan]";

    if (headDefined) {
        str = "["+head.Cunit(0)+","+head.Cunit(1)+",KM/S]      ";
        Stream  << setw(k) << right << str;
    }

    Stream  << setw(m-3) << right << "[pix] "
            << setw(m-3) << right << "[pix] "
            << setw(m-3) << right << "[chan]"
            << setw(m) << right << "[pix] "
            << setw(m) << right << "[pix] "
            << setw(m) << right << "[km/s]"
            << setw(m) << right << "[km/s]"
            << setw(m+5) << right << "[JY*km/s]"
            << setw(m-2) << right << " "
            << setw(m) << right << head.Bunit();

    Stream  << endl;

    if (headDefined) Stream << "#" << setw(150) << setfill('_') << " " << endl <<  "#" << endl;
    else Stream << "#" << setw(104) << setfill('_') << " " << endl << "#" << endl;

    Stream  << setfill(' ');

    for (int i=0; i<numObj; i++){
        Detection<T> *obj = sources->pObject(i);

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
        Stream  << "      " << setw(m-7) << left << i+1
                << setw(6) << right << to_string(Xcenter,0)
                << setw(5) << right << to_string(Ycenter,0)
                << setw(5) << right << to_string(Zcenter,0);


        if (headDefined) {
           // str = to_string(RA,3)+"  "+to_string(DEC,3)+"  "+to_string(VEL,1);
            //Stream    << setw(k) << right << str;
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

        Stream  << setw(m-3) << right << Xint
                << setw(m-3) << right << Yint
                << setw(m-3) << right << Zint
                << setw(m) << right << obj->getSpatialSize()
                << setw(m) << right << obj->getSize();

        Stream << right << setw(m)  << setprecision(1) << obj->getVsys()
               << right << setw(m)  << setprecision(1) << obj->getW20()
               << right << setw(m+5) << setprecision(3)  << obj->getIntegFlux()
               << right << setw(m-2) << setprecision(1)  << obj->getPeakFlux()/stats.getSpread()
               << right << setw(m) << setprecision(3)  << obj->getPeakFlux();

        Stream  << endl;

    }

    Stream  << endl;

}


template <class T>
void Cube<T>::plotDetections() {

    int numObj = getNumObj();
    if (numObj==0) return;

    // Getting regions of detections
    bool *isObj = new bool[numPix];
    for (int j=0; j<numPix; j++) isObj[j]=false;
    for (int i=0; i<numObj; i++) {
        for(auto &vox : sources->pObject(i)->getPixelSet(array, axisDim)){
            long pos = vox.getX()+vox.getY()*axisDim[0]+vox.getZ()*axisDim[0]*axisDim[1];
            isObj[pos] = true;
        }
    }

    // Writing a datacube with just the detected objects
    Cube<T> *det = new Cube<T>(axisDim);
    for (size_t i=0; i<det->NumPix(); i++) det->Array()[i] = array[i]*isObj[i];
    det->saveHead(head);
    det->saveParam(par);
    det->Head().setMinMax(0,0);
    det->fitswrite_3d((par.getOutfolder()+"detections.fits.gz").c_str());
    delete det;

    // Writing detection map
    Image2D<int> *DetMap = new Image2D<int>(axisDim);
    for (int i=0; i<axisDim[0]*axisDim[1];i++) DetMap->Array()[i] = sources->DetectMap(i);
    DetMap->copyHeader(head);
    DetMap->Head().setMinMax(0,0);
    DetMap->Head().setBtype("detected_chan");
    DetMap->fitswrite_2d((par.getOutfolder()+"DetectMap.fits.gz").c_str());
    delete DetMap;

    // Writing kinematic maps
    std::vector< MomentMap<T> > allmaps = getAllMoments<T>(this,true,isObj,"MOMENT");
    allmaps[0].fitswrite_2d((par.getOutfolder()+head.Obname()+"_mom0th.fits").c_str());
    allmaps[1].fitswrite_2d((par.getOutfolder()+head.Obname()+"_mom1st.fits").c_str());
    allmaps[2].fitswrite_2d((par.getOutfolder()+head.Obname()+"_mom2nd.fits").c_str());

}



template <class T>
Cube<T>* Cube<T>::extractCubelet(Detection<T> *obj, int edges, int *starts) {

    /// Extract a sub-cube from a detection

    // Finding boundaries of the cubelet
    long objmin[3] =  {obj->getXmin(),obj->getYmin(),obj->getZmin()};
    long objmax[3] =  {obj->getXmax(),obj->getYmax(),obj->getZmax()};
    int stops[3], newdim[3];
    for (int i=0; i<3; i++) {
        starts[i] = objmin[i]-edges>=0 ? objmin[i]-edges : 0;
        stops[i]  = objmax[i]+edges<axisDim[i] ? objmax[i]+edges : axisDim[i]-1;
        newdim[i] = stops[i]-starts[i]+1;
    }

    // Creating the sub-cube and copying data
    Cube<T> *clet = new Cube<T>(newdim);
    clet->saveHead(head);
    clet->saveParam(par);

    for (int z=starts[2]; z<=stops[2]; z++)
        for (int y=starts[1]; y<=stops[1]; y++)
            for (int x=starts[0]; x<=stops[0]; x++)
                clet->Array(x-starts[0],y-starts[1],z-starts[2]) = Array(x,y,z);

    // Upgrading header. WCS is recentred
    double pix[3] = {round(obj->getXaverage()),round(obj->getYaverage()),0}, world[3];
    int status = head.pixToWCS(pix,world);
    for (int i=0; i<3; i++) {
        clet->Head().setDimAx(i,newdim[i]);
        if (status==0) {
            clet->Head().setCrpix(i,pix[i]-starts[i]+1);
            clet->Head().setCrval(i,world[i]);
        }
        else clet->Head().setCrpix(i,head.Crpix(i)-starts[i]);
        clet->Head().updateWCS();
    }

    clet->Head().setMinMax(0,0);
    clet->Head().Keys().clear();
    clet->Head().Keys().push_back("HISTORY BBAROLO SOURCE FINDER DATA PRODUCT");

    return clet;
}


template <class T>
void Cube<T>::writeCubelets() {

    /// Write FITS cubelets and separate kinematic maps of all detections

    if (getNumObj()==0) return;

    std::string outfold = par.getOutfolder()+"sources/";
    std::string object  = head.Name();
    mkdirp(outfold.c_str());


    ProgressBar bar(true,par.isVerbose(),par.getShowbar());
    bar.init(" Writing cubelets ...",getNumObj());

    for (int i=0; i<getNumObj(); i++) {
        bar.update(i+1);
        std::string srcname = object+"_"+to_string(i+1);
        mkdirp((outfold+srcname).c_str());

        Detection<T> *obj = sources->pObject(i);
        obj->calcWCSparams(head);

        // Writing sub-cube
        int starts[3];
        Cube<T> *c = extractCubelet(obj,par.getParSE().edges,starts);
        c->Head().setName(srcname);
        c->fitswrite_3d((outfold+srcname+"/"+srcname+".fits").c_str(),true);

        // Creating a submask for cubelet
        bool *isObj = new bool[c->NumPix()];
        for (int j=0; j<c->NumPix(); j++) isObj[j]=false;

        for(auto &v : obj->getPixelSet(array, axisDim)){
            long pos = c->nPix(v.getX()-starts[0],v.getY()-starts[1],v.getZ()-starts[2]);
            isObj[pos] = true;
        }

        // Writing sub-kinematic maps
        c->pars().setVerbosity(false);
        // Writing kinematic maps
        std::vector< MomentMap<T> > allmaps = getAllMoments<T>(c,true,isObj,"MOMENT");
        allmaps[0].fitswrite_2d((outfold+srcname+"/"+srcname+"_mom0.fits").c_str());
        allmaps[1].fitswrite_2d((outfold+srcname+"/"+srcname+"_mom1.fits").c_str());
        allmaps[2].fitswrite_2d((outfold+srcname+"/"+srcname+"_mom2.fits").c_str());

        delete [] isObj;

        // Writing total spectrum to file
        std::ofstream fileo((outfold+srcname+"/"+srcname+"_spectrum.dat").c_str());
        fileo << "#Velocity(km/s)     Flux" << std::endl;

        for (int z=0; z<c->DimZ(); z++) {
            float intSpec = 0;
            double vel = AlltoVel(c->getZphys(z),c->Head());
            for (int y=0; y<c->DimY(); y++) {
                for (int x=0; x<c->DimX(); x++) {
                    if (obj->isInObject(x+starts[0],y+starts[1],z+starts[2])) {
                        double flux = c->Array(x,y,z);
                        Pbcor<double>(x+starts[0],y+starts[1],z+starts[2],flux,head);
                        intSpec += flux;
                    }
                }
            }
            fileo << left << setw(18) << vel << "  " << FluxtoJy(intSpec,c->Head()) << endl;
        }

        fileo.close();
        delete c;

#ifdef HAVE_GNUPLOT
        Gnuplot gp;

        gp.begin();
        gp.commandln("set terminal postscript eps color");
        gp.commandln("unset key");
        gp.commandln("set xlabel 'Velocity [km/s]'");
        gp.commandln("set ylabel 'Flux density'");
        //gp.commandln("set size square");

        std::string titlecmd = "set title 'Individual spectrum for " +srcname+"'";
        std::string outcmd = "set output '"+outfold+srcname+"/"+srcname+"_spectrum.eps'";
        std::string plotcmd = "plot '"+outfold+srcname+"/"+srcname+"_spectrum.dat' with lp ls 3 lt 7";
        gp.commandln(titlecmd.c_str());
        gp.commandln(outcmd.c_str());
        gp.commandln((plotcmd).c_str());

    gp.end();
#endif

    }

    bar.fillSpace("Done.\n");

}



// Explicit instantiation of the class
template class Cube<short>;
template class Cube<int>;
template class Cube<long>;
template class Cube<float>;
template class Cube<double>;
