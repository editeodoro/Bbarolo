//--------------------------------------------------------------------
// image.hh: Definition of the Image class.
//--------------------------------------------------------------------

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
 along with BBarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning BBarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#ifndef IMAGE_HH_
#define IMAGE_HH_

#include <iostream>
#include <string>
#include "stats.hh"
#include "header.hh"
#include "cube.hh"

using namespace Statistics;

template <class Type>
class Image2D
{
public:
    Image2D();										/// Default constructor.
    Image2D(int *dimensions);						/// Alternative constructor.
    virtual ~Image2D();								/// Default destructor.
    Image2D(const Image2D &i);						/// Copy constructor.
    Image2D& operator=(const Image2D &i);			/// Assignement operator.
	
    // Overloadad () operator for easy access the main array. Any control on the index.
    inline Type& operator() (unsigned x, unsigned y) {return array[x+y*axisDim[0]];};
    inline Type  operator() (unsigned x, unsigned y) const {return array[x+y*axisDim[0]];};
    inline Type& operator() (unsigned i) {return array[i];};
    inline Type  operator() (unsigned i) const {return array[i];};

	void defaults();

    /// Inline functions to access the data:
	    
	long 	NumPix() {return numPix;};
	int     DimX() {return axisDim[0];}; 
    int     DimY() {return axisDim[1];};
	long 	nPix  (int x, int y) {return y*axisDim[0]+x;};
	Type*	Array () {return array;};
    Type 	Array (long npix) {return array[npix];};
    Type	Array (int x, int y) {return array[y*axisDim[0]+x];};
    Header&	Head	 () {Header &h = head; return h;};
    bool	HeadDef (){return headDefined;};
    void	setXsize (int i) {axisDim[0] = i;};
	void	setYsize (int i) {axisDim[1] = i;};
	void	setHeadDef (bool b) {headDefined = b;};
	void	saveHead (Header &h) {head = h; headDefined=true;};
	Type 	printStats() {std::cout << stats << std::endl;}; 
	int     getopts(int argc, char **argv){return par.getopts(argc,argv);};	
	void 	saveStats(Stats<Type> newStats){stats = newStats;};
	bool	StatsDef () {return statsDefined;};
	Param&  pars(){Param &rpar = par; return rpar;};
	void    showParam(std::ostream &stream){stream << par;};
	void    saveParam(Param &newpar){par = newpar;};
   	Stats<Type>  getStats() { return stats;};
	Stats<Type>& stat() {Stats<Type> &rstats = stats; return rstats;};
	void setMinSize(int i){minSize=i;};
    double	getYphys (double y) {return (y+1-head.Crpix(1))*head.Cdelt(1)+head.Crval(1);};
    double	getXphys (double x) {return (x+1-head.Crpix(0))*head.Cdelt(0)+head.Crval(0);};

    ///	I/O functions:

    void setImage(Type *input, int *dim);		
    bool readImage (std::string fname);				
	void copyHeader (Header &c);	
	void setImageStats();
	bool fitswrite_2d (const char *outname="./output/image.fits"); 

	void extractSpectrum(Type *Array, int *dim, long pixel);  /// Extract a spectrum from an array. 
    void extractImage(Type *Array, int *dim, long channel);	/// Extract an image from an array.																	
    void extractSpectrum(Cube<Type> &cube, long pixel);   			/// Extract a spectrum from a cube.
    void extractImage(Cube<Type> &cube, long channel);				/// Extract an image from a cube.
    void extractGlobalSpectrum(Cube<Type> *cube);
    
	/// Detection-related functions:
    
    bool isDetection(long x, long y);				 			/// A detection test for a pixel location. 
    std::vector<PixelInfo::Scan<Type> > findSources1D();	 			/// Front-end function for spectrumDetect.
    std::vector<PixelInfo::Object2D<Type> > findSources2D();			/// Front-end function for imageDetect.   
    
    std::vector<PixelInfo::Scan<Type> > spectrumDetect(std::vector<bool> &arraybool);   /// Detect objects in a 1-D spectrum. 
    std::vector<PixelInfo::Object2D<Type> > imageDetect(std::vector<bool> &arraybool);	 ///Detect objects in a 2-D image. 
	
protected:
	Type 		*array;						///< The cube data array.
	bool		arrayAllocated;				///< Is array allocated?
	int			numAxes;					///< Number of axis.
	long		numPix;						///< Total number of pixel.
	int         axisDim[2];   				///< Array of axis dimensions of cube
	Header		head;						///< Fits header information.
	bool		headDefined;				///< Has been an header defined?			
	Stats<Type> stats; 						///< The statistics for the data array.
	bool		statsDefined;				///< Have been statitistics defined?
	Param		par;				        ///< A parameter list.
    unsigned int minSize;    				///< The minimum number of pixels for a detection to be accepted.
    
	bool fitsread_2d ();

};

//#include "image.cpp"


#endif
