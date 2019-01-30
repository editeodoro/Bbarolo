//---------------------------------------------------------------
// smooth3d.hh: Definition of the Smooth3D class.
//---------------------------------------------------------------

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

// Smooth3D is a class for smoothing datacubes with a gaussian. 
// NB; All beams are in arcsec, position angles in degrees.
//
// The correct way to call the class is the following:
//
// 1) call smooth(...) functions: 
//      -Cube *c:       the Cube object to be smoothed.
//      -int *Bhi:      an array with upper coordinate box limits.
//      -int *Blo:      an array with lower coordinate box limits.
//      -Beam oldbeam:  a Beam object with {oldbmaj, oldbin, oldpa}. 
//      -Beam newbeam;  the desired beam after the smoothing process.
//
//  
// The smoothed array is written in the 'array' variable. 
//
//


#ifndef SMOOTH3D_HH_
#define SMOOTH3D_HH_

#include <iostream>
#include <Arrays/cube.hh>


struct Beam {
    double bmaj;                    //< Major axis of gaussian beam FWHM.
    double bmin;                    //< Minor axis of gaussian beam FWHM.
    double bpa;                     //< Position angle fo beam (N->E).
};



template <class T>
class Smooth3D 
{
public:
    Smooth3D();                                 //< Default constructor.
    virtual ~Smooth3D();                        //< Destructor.
    Smooth3D(const Smooth3D &g);                //< Copy constructor.
    Smooth3D& operator=(const Smooth3D &g);     //< Copy operator.
    void defaults();

    
    /// Obvious inline functions
    T       Array (int i) {return array[i];};
    T       *Array () {return array;};
    T       *Confie() {return confie;};                 
    Beam    Oldbeam() {return oldbeam;};                    
    Beam    Newbeam() {return newbeam;};                    
    Beam    Conbeam() {return conbeam;};                    
    double  Scalefac(){return scalefac;};
    void    setUseScalefac (bool ff) {usescalefac=ff;};
    void    setUseBlanks(bool b) {useBlanks=b;};
    
    void cubesmooth(Cube<T> *c);
    void smooth(Cube<T> *c, Beam Oldbeam, Beam Newbeam);
    void smooth(Cube<T> *c, int *Bhi, int *Blo, Beam Oldbeam, Beam Newbeam);
    bool smooth(Cube<T> *c, Beam Oldbeam, Beam NewBeam, T *OldArray, T *NewArray);
    void fitswrite();

private:
    
    Cube<T> *in;                        //< A pointer to the Cube to be smoothed.
    T       *array;                     //< The smoothed array.
    bool    arrayAllocated;             //< Have been array allocated?
    bool    *blanks;                    //< An array with blank pixels.
    bool    blanksAllocated;            //< Have been blanks allocated?
    bool    useBlanks;
    int     NdatX;                      //< X-box dimension.
    int     NdatY;                      //< Y-box dimension.
    int     NdatZ;                      //< Number of channels.
    double  *confie;                    //< The convolution field.
    bool    confieAllocated;            //< Have been confie allocated?
    int     NconX;                      //< X-dimension of confie.
    int     NconY;                      //< Y-dimension of confie.
    int     bhi[3], blo[3];             //< Box boundaries.
    int     fhi[3], flo[3];             //< Full cube boundaries.
    int     dimAxes[3];                 //< Full cube X-Y dimensions.
    double  gridspace[2];               //< Grid spacing of channel maps.
    Beam    oldbeam;                    //< Input beam.
    Beam    newbeam;                    //< Output beam.
    Beam    conbeam;                    //< Convolution beam.
    bool    beamDefined;                
    double  scalefac;                   //< Scale factor.
    bool    usescalefac;
    float   crota;                      //< Rotation angle of maps.
    double  cutoffratio;                //< Cutoff for gaussian kernel.
    size_t  maxconv;                    //< Max convolution size: 512*512.
    bool    fft;                        //< Convolution with FFT.   

    /// Pointer to the function to be convolved with (Gaussian or Moffat)
    typedef bool (Smooth3D<T>::*funcPtr) (Beam, Beam);
    funcPtr func_psf;
    
    bool defineBeam_Gaussian(Beam Oldbeam, Beam Newbeam);
    bool defineBeam_Moffat(Beam Oldbeam, Beam Newbeam);
    bool calculate(T *OldArray, T *NewArray);   
    bool calculatefft(T *OldArray, T *NewArray);
    bool Convpars();                
    bool Fillgauss2d(Beam varbeam, float ampl, bool norm, int &NconX, 
                     int &NconY, double *cfie);
    bool FillMoffat2d(Beam varbeam, float ampl, bool norm, int &NconX, 
                                      int &NconY, double *cfie);
    int Convolve(double *cfie, int ncx, int ncy, T *dat1, 
                 T *dat2, int ndx, int ndy);

};


/////////////////////////////////////////////////////////////////////////////////////
/// A class for Hanning smoothing a datacube
/////////////////////////////////////////////////////////////////////////////////////
template <class T>
class Hanning3D 
{
public:
    Hanning3D(size_t window_size) {window=window_size;}         //< Constructor.
    virtual ~Hanning3D() {if (arrayAllocated) delete [] array;} //< Destructor.

    /// Obvious inline functions
    T&   Array (int i) {return array[i];};
    T    *Array () {return array;};

    void compute(Cube<T> *c);
    void compute(T *inarray, size_t xsize, size_t ysize, size_t zsize);
    void fitswrite(Cube<T> *templ);
   
private: 
    T       *array;                     //< The smoothed array.
    bool    arrayAllocated;             //< Have been array allocated?
    size_t  window;                     //< Size of the Hanning window
};

#endif
