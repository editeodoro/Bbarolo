// -----------------------------------------------------------------------
// utils.hh: Prototypes for utility functions.
// -----------------------------------------------------------------------

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

#ifndef UTILS_HH_
#define UTILS_HH_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <time.h>
#include <fitsio.h>
#include <Arrays/param.hh>
#include <Arrays/header.hh>
#include <Arrays/rings.hh>
#include <Map/voxel.hh>
#include <sys/stat.h>


#define DEFAULT_MODE      S_IRWXU | S_IRGRP |  S_IXGRP | S_IROTH | S_IXOTH


/// Functions to allocate or deallocate multidimensional arrays.
/// Defined in "allocator.cpp".

template <class Type> Type **allocate_2D (int xdim, int ydim);
template <class Type> void deallocate_2D (Type **array2d, int xdim);
template <class Type> Type ***allocate_3D (int xdim, int ydim, int zdim);
template <class Type> void deallocate_3D (Type ***array3d, int xdim, int ydim);


/// Functions to calculate statistical parameters.
/// Defined in "statistics.cpp".

template <class T> T absval(T value);
template <class T> void findMinMax(const T *array, const size_t size, T &min, T &max);
template <class T> void findAllStats(T *array, size_t size, T &mean, T &stddev, T &median, T &madfm);
template <class T> void findAllStats(T *array, size_t size, T &mean, T &stddev, T &median, T &madfm, T &maxx, T &minn);
template <class T> void findAllStats(T *array, size_t size, bool *mask, T &mean, T &stddev, T &median, T &madfm, T &maxx, T &minn);
template <class T> T findMean(T *array, size_t size); 
template <class T> T findMean(T *array, bool *mask, size_t size);
template <class T> T findStddev(T *array, size_t size);
template <class T> T findStddev(T *array, bool *mask, size_t size);
template <class T> T findMedian(T *array, size_t size, bool changeArray=false);
template <class T> T findMedian(T *array, bool *mask, size_t size);
template <class T> T findMedian(std::vector<T> array, size_t size);
template <class T> T findMADFM(T *array, size_t size, T median, bool changeArray=false);
template <class T> T findMADFM(T *array, bool *mask, size_t size, T median);


/// Functions to manipulate strings.
/// Defined in string.cpp

std::string makelower( std::string s );
std::string makeupper( std::string s );
std::string stringize(bool b);
bool boolify(std::string s);
std::string readSval(std::stringstream& ss);
std::string readFilename(std::stringstream& ss);
bool readFlag(std::stringstream& ss);
float readFval(std::stringstream& ss);
double readDval(std::stringstream& ss);
int readIval(std::stringstream& ss);
template <class T> void readArray(std::stringstream& ss, T *val, int n);
template <class T> std::vector<T> readVec(std::stringstream& ss);
std::string removeLeadingBlanks(std::string s);
std::string deblank(std::string s);
std::string deblankAll (std::string s);
template <class T> std::string to_string (const T& t, int precision=-1);
void printBackSpace(std::ostream &stream, int num);
void printBackSpace(int num);
void printSpace(std::ostream &stream, int num);
void printSpace(int num);
void printHash(std::ostream &stream, int num);
void printHash(int num);
void checkHome(std::string &s); 
std::string randomAdjective (int type);
    
/// Fuctions to interpolation and fit.
/// Defined in interpolation.cpp

template <class T> double *spline(T *x, T *y, int n, double yp1=1.e91, double ypn=1.e91);
template <class T> double splint(T *xa, T *ya, T *y2a, int n, T x);
template <class T> int linear_reg(int num, T *x, T *y, float &m, float &merr, float &b, 
                                  float &berr, float &r, int ilow, int ihigh);
template <class T> int linear_reg(int num, std::vector<T> &x, std::vector<T> &y, float *p, 
                                  float *err, float &r, int ilow, int ihigh);

template <class Type> Type func_gauss (Type *x, Type *p, int numpar);
template <class Type> void derv_gauss (Type *x, Type *p, Type *dyda, int numpar);
template <class T> double** RotMatrices (T alpha, T beta, T gamma);
template <class T> double* MatrixProduct (T *M1, int *size1, T *M2, int *size2);
double *cp_binomial(int points);
template <class T> void bezier_interp(std::vector<T> x_in,  std::vector<T> y_in,
                                      std::vector<T> &x_out,std::vector<T> &y_out,
                                      int fp=0, int np=-1, int ns=-1);

/// Other functions:
/// Defined in utils.cpp

bool mkdirp(const char* path, mode_t mode = DEFAULT_MODE);
void FitsWrite_2D (const char *filename, float *image, long xsize, long ysize);
void FitsWrite_2D (const char *filename, double *image, long xsize, long ysize);
void FitsWrite_3D (const char *outfile, float *outcube, long *dimAxes);
void FitsWrite_3D (const char *outfile, short *outcube, long *dimAxes);

template <class T> T AlltoVel (T in, Header &h);
template <class T> T DeltaVel (Header &h);
template <class T> T FluxtoJy (T in, Header &h);
template <class T> T Pbcor (PixelInfo::Voxel<T> &v, Header &h);
template <class T> void Pbcor(long x, long y, long z, T &flux, Header &h);

// Function to obtain a linear scale from an angular scale. 
std::string get_currentpath();
// d is in Mpc and the returned value is in Kpc/arc. 
double KpcPerArc(double d);
// Function for the distance from systemic velocity by Hubble law.
// vsys is given in km/s, the returned value is in Mpc.
double VeltoDist(double vsys); 
double RedtoDist(double redshift);

template <class Type> bool isNaN (Type n);
template <class T> bool isBlank(T val) {return (val==0 || isNaN(val));}
bool fexists(std::string filename);

double arcsconv(std::string cunit);
double degconv(std::string cunit);
std::string decToDMS(const double dec, std::string type, int decPrecision=2);
double dmsToDec(std::string dms);
template <class T> T angularSeparation(T &ra1, T &dec1, T &ra2, T &dec2);


template <class T> T unifrand(T maxs, T mins) {
    //srand(time(NULL)); 
    return (rand()/(double)RAND_MAX)*(maxs-mins)+mins;
}

template <class T> bool getData (std::vector<std::vector<T> > &allData, std::string file, bool err_verbose=true);
template <class T> bool getDataColumn (std::vector<T> &data, std::string filestring);
template <class T> Rings<T> *readRings(GALFIT_PAR &par, Header &h);

double* getCenterCoordinates(std::string *pos, Header &h);

template <class T> T* RingRegion (Rings<T> *r, Header &h);
template <class T> T* SimulateNoise(double stddev, size_t size);
template <class T> T* HanningSmoothing(T *inarray, size_t npts, size_t hanningSize);

template <class T> int selectBitpix();
template <class T> int selectDatatype();



// WCS Utilities
int pixToWCSSingle(struct wcsprm *wcs, const double *pix, double *world);
int wcsToPixSingle(struct wcsprm *wcs, const double *world, double *pix);


#endif
