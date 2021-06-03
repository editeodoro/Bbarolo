//--------------------------------------------------------------------
// cube.hh: Definition of Cube class
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

#ifndef CUBE_HH_
#define CUBE_HH_

#include <string>
#include <Arrays/header.hh>
#include <Arrays/stats.hh>
#include <Arrays/param.hh>
#include <Tasks/search.hh>
#include <Map/detection.hh>

using namespace Statistics;

template <class T>
class Cube
{   
public:
    
    Cube() {defaults();}                                /// Default constructor.
    Cube(std::string fname);                            /// Basic constructor.
    Cube(int *dimensions);                              /// Alternative constructor.
    virtual ~Cube();                                    /// Destructor. 
    Cube(const Cube &c);                                /// Copy constructor.
    Cube& operator=(const Cube &c);                     /// Copy operator.
    void defaults();

    // Overloadad () operator for easy access the main array. No controls on the index.
    inline T& operator() (size_t x, size_t y, size_t z) {return *(array+nPix(x,y,z));}
    inline T& operator() (size_t i) {return *(array+i);}
    inline T& operator[] (size_t i) {return *(array+i);}
    
    /// Obvious inline functions to access a private member of class:
    int     NumAx () {return numAxes;}
    size_t  NumPix() {return numPix;}
    int*    AxisDim () {return axisDim;}
    int     AxesDim (int i) {return axisDim[i];}
    int     DimX(){return axisDim[0];}
    int     DimY(){return axisDim[1];}
    int     DimZ(){return axisDim[2];}
    size_t  nPix  (size_t x,size_t y,size_t z) {return x+y*axisDim[0]+z*axisDim[0]*axisDim[1];}
    T*      Array () {return array;}
    T&      Array (size_t npix) {return this->operator()(npix);}
    T&      Array (size_t x,size_t y,size_t z) {return this->operator()(x,y,z);}
    void    setArray (T *ar) {array = ar;}
    double  getZphys (double z) {return (z+1-head.Crpix(2))*head.Cdelt(2)+head.Crval(2);}
    double  getXphys (double x) {return (x+1-head.Crpix(0))*head.Cdelt(0)+head.Crval(0);}
    double  getYphys (double y) {return (y+1-head.Crpix(1))*head.Cdelt(1)+head.Crval(1);}
    double  getZgrid (double v) {return (v-head.Crval(2))/head.Cdelt(2)+head.Crpix(2)-1;}
    double* getXYZphys (double x, double y, double z) {double p[3] = {x,y,z}; double *w = new double[3]; head.pixToWCS(p,w); return w;}
    double* getXYphys (double x, double y) {return getXYZphys(x,y,0);}

    Header& Head     () {Header &h = head; return h;}
    void    setHeadDef (bool b) {headDefined = b;}
    bool    HeadDef (){return headDefined;}
    void    saveHead (Header &h) {head = h; headDefined=true;}
    void    setBeam  (float a, float b, float c) {head.setBeam(a,b,c);}
    float*  getBeam() {float *f=new float; f[0]=head.Bmaj(); f[1]=head.Bmin(),f[2]=head.Bpa(); return f;}
    
    bool*   Mask    () {return mask;}
    bool&   Mask    (size_t npix) {return mask[npix];}
    bool&   Mask    (size_t x,size_t y,size_t z) {return mask[x+y*axisDim[0]+z*axisDim[0]*axisDim[1]];}
    bool    MaskAll () {return maskAllocated;}
    
    void    printStats() {std::cout << stats << std::endl;}
    Stats<T>  getStats(){ return stats;}
    Stats<T>& stat(){Stats<T> &rstats = stats; return rstats;}
    void    saveStats(Stats<T> newStats){stats = newStats; statsDefined=true;}
    bool    StatsDef () {return statsDefined;}
    int     getopts(int argc, char **argv){return par.getopts(argc,argv);}
    Param   getParam(){return par;}
    Param&  pars(){ Param &rpar = par; return rpar;}
    void    showParam(std::ostream &stream){stream << par;}
    void    saveParam(Param &newpar){par = newpar;}
    long    getNumObj(){return sources->getNumObj();}
    Detection<T>* pObject(long n){return sources->pObject(n);}
    Search<T>* getSources () {return sources;}
    bool    getIsSearched () {return isSearched;}
    void    saveSources (Search<T> &s) {if (isSearched) delete sources; sources = new Search<T>(s); isSearched=true;}



    /// Functions for Fitsfile I/O:
    void    setCube  (T *input, int *dim);
    bool    readCube (std::string fname);                                   /// Front-end to read array from Fits.
    bool    fitsread_3d ();                                                 /// Read data array from Fits file.                                             
    bool    fitswrite_3d (const char *outfile, bool fullHead=false);        /// Write a Fits cube.                                      
    
    /// Statistics functions:
    void    setCubeStats();                              /// Calculate statistical parameters for cube. 
    bool    isDetection(long x, long y, long z) {return stats.isDetection(Array(x,y,z));}         /// Can be a voxel considered a detection ?

    /// Searching functions, defined in search.cpp.
    void    search();                                    /// Front-end function to search in a 3-D cube.
    void    search(std::string searchtype, float snrCut, float threshold, bool adjacent,
                   int threshSpatial, int threshVelocity, int minPixels, int minChannels,
                   int minVoxels, int maxChannels, float maxAngSize, bool flagGrowth,
                   float growthCut, float growthThreshold, bool RejectBefore, bool TwoStage,int NTHREADS);
    void    printDetections (std::ostream& Stream);      /// An easy way to print the detection list.
    Cube<T>* extractCubelet(Detection<T> *obj, int edges, int *starts);
    void    writeCubelets();
    void    writeDetections();
    int     plotDetections();




    /// Blanking and Maps functions.
    void    BlankCube (T *Array, size_t size);            /// Blank a input array using Cube::mask.
    void    BlankMask(float *channel_noise=NULL,bool onlyLargest=true);       /// Define Cube::mask;
    
    Cube<T>*    Reduce (int fac,std::string rtype="spatial");
    void    CheckChannels ();
    void    continuumSubtract();
    void    checkBeam();
    
protected:
    T           *array;                     ///< The cube data array.
    bool        arrayAllocated;             ///< Is array allocated?
    short       numAxes;                    ///< Number of axis.
    size_t      numPix;                     ///< Total number of pixel.
    int         *axisDim;                   ///< Array of axis dimensions of cube
    bool        axisDimAllocated;           ///< has axisDim been allocated?
    Header      head;                       ///< Fits header information.
    bool        headDefined;                ///< Has been an header defined?            
    Stats<T>    stats;                      ///< The statistics for the data array.
    bool        statsDefined;               ///< Have been statistics defined?
    Param       par;                        ///< A parameter list.
    

private:
    int         bitpix;                     ///< Data type code values for FITS images.
    int         datatype;                   ///< Data type when reading or writing data.
    bool        *mask;                      ///< A mask for blanked cube.
    bool        maskAllocated;              ///< Has mask been allocated?
    Search<T>   *sources;                   ///< A pointer to the source-finder.
    bool        isSearched;                 ///< Already searched?

};

#endif
