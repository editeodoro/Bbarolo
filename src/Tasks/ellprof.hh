////////////////////////////////////////////////////////////////////////
// ellprof.hh. A class for deriving density profiles from a density map
////////////////////////////////////////////////////////////////////////

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


// Ellprof derives a radial profile from a total intensity map.
// The class is derived from the ELLINT task from GIPSY software (see
// ELLINT documentation for further information about the parameters).
//
// Output parameters:
//      - Num:            Number of pixels per ring per segment
//      - Sum:            Sum of pixels per ring per segment
//      - Mean:.......... Mean of pixels per ring per segment
//      - Median:         Median of pixels per ring per segment
//      - Var:            Variance of pixels per ring per segment.
//                        Error RMS is sqrt(Var).
//      - Datamin:        Minimum value per ring per segment.
//      - Datamax:        Maximum value per ring per segment.
//      - Area:           Area covered per ring per segment.
//      - Surfdens:       Surface density per ring per segment.
//                        sd =  Sum / Area in map unit / arcsec2.
//      - Mass_Surfdens:  Mass surface density in Msun/pc2
//                        NEED TO SET Mass and Distance(see below)
//                        or it returns all NaN values!
//
//      - NumBlanks:       Number of blank pixels (NaN) per ring/seg.
//      - BlankArea:       Area covered by blank pixels (NaN).
//      - Surfdens_Bl:     as Surfdens but taking into account blanks.
//                         sd =  Sum / (Area+BlankArea).
//      -Mass_Surfdens_Bl: Mass surface density with blanks areas
//
//
//
// The correct way to call the class is the following:
//
// 1) call constructor function:
//
//      -MomentMap *image:  the image from which the profile is built.
//                          it should have NaN as blank pixels.
//      -Rings *rings:      a Rings object containing the input rings.
//                          In input, distances have to be in ARCSEC,
//                          angles in DEGREE. Required quantities are nr,
//                          radsep, radii, inc and phi.
//      -int nseg:          an array with upper coordinate limits.
//      -float *segments:   array of size 2*nseg with segment ranges.
//                          Profiles can be calculated in a part of a
//                          ring, a so called segment. Default segment
//                          is the entire ring (0,360). Required to de-
//                          fine a valid segment are two angles. The
//                          first pair of segment angles MUST be 0,360.
//                          Example: segments = [0,360,180,360]
//                                   --> full ring + half galaxy
//
// 2) Optional: call setOptions(...) to change default options:
//
//       - bool overlap:    Weight data in overlapping rings?
//                          If true, data are weighted in overlapping rings.
//                          If at certain position a pixel with image value
//                          X is encountered in N rings, its image value
//                          decreases to X/N. DEFAULT = TRUE
//       - float *range:    Range of intensities for a pixel to be used.
//                          It is an array of two values,
//                          such that range[0] < pixval < range[1].
//       - float *subp:     arryas of size two, giving number of 'sub'pixels
//                          in x and y per pixel. DEFAULT [2,2]
//                          Make a better approximation of the true area in a
//                          ring by dividing a pixel into subpixels.
//       - float *mass:     Mass of the gas enclosed in the last ring in Msun.
//                          NEEDED if you want the mass density profile!
//       - float *distance: Distance in Mpc.
//                          NEEDED if you want the mass density profile!
//
// 3) call RadialProfile()
//
// 4) Optional: call printProfile(std::ostream) to print output values on the
//              screen (std::cout) or in a file (std::ofstream).
//
//


#ifndef ELLPROF_HH_
#define ELLPROF_HH_

#include <iostream>
#include <vector>
#include <Arrays/cube.hh>
#include <Tasks/moment.hh>
#include <Tasks/galmod.hh>

namespace Tasks {
template <class T>
class Ellprof
{
public:
    Ellprof(Cube<T> *c);
    Ellprof(Cube<T> *c, Rings<T> *inR) {setFromCube(c,inR);}
    Ellprof(MomentMap<T> *image, size_t nrad, float width, float phi, float inc, float *pos, size_t nseg=1, float* segments=NULL);
    Ellprof(MomentMap<T> *image, Rings<T> *rings, size_t nseg=1, float* segments=NULL);
    Ellprof(const Ellprof& p);
    Ellprof& operator= (const Ellprof& p);
    virtual ~Ellprof() {deallocateArrays();}
    
    void   setFromCube(Cube<T> *c, Rings<T> *inR);
    void   setOptions (bool overlap, float *range, float *subp);
    void   setOptions (float mass, float distance);
    void   RadialProfile ();
    void   printProfile (ostream& theStream=std::cout, int seg=0);
    void   writeMap (std::string fname) {im->fitswrite_2d(fname.c_str());}


    // Inline functions to access class members.
    size_t getNrad () {return Nrad;}
    size_t getNseg () {return Nseg;}
    float  getRadius (int i) {return Radius[i];}
    float  getRadiusKpc (int i) {return Radius_kpc[i];}
    float  getWidth (int i) {return Width[i];}
    float  getPhi (int i) {return Phi[i];}
    float  getInc (int i) {return Inc[i];}
    float  getPosition (size_t i) {return Position[i];}

    double getSum(size_t i, size_t j=0) {return Sum[i][j];}
    long   getNum(size_t i, size_t j=0) {return Num[i][j];}
    long   getNumBlanks(size_t i, size_t j=0) {return Numblanks[i][j];}
    double getDatamin(size_t i, size_t j=0) {return Datamin[i][j];}
    double getDatamax(size_t i, size_t j=0) {return Datamax[i][j];}
    double getVariance(size_t i, size_t j=0) {return Var[i][j];}
    double getRMS(size_t i, size_t j=0) {return sqrt(fabs(Var[i][j]));}
    double getMean(size_t i, size_t j=0) {return Mean[i][j];}
    double getMedian(size_t i, size_t j=0) {return Median[i][j];}
    double getArea(size_t i, size_t j=0) {return Area[i][j];}
    double getBlankArea(size_t i, size_t j=0) {return Blankarea[i][j];}
    double getSurfDens(size_t i, size_t j=0) {return Surfdens[i][j];}
    double getSurfDensFaceOn(size_t i, size_t j=0) {return Cosinc[i]*Surfdens[i][j];}
    double getSurfDens_Bl(size_t i, size_t j=0) {return Surfdens_Bl[i][j];}
    double getSurfDensFaceOn_Bl(size_t i, size_t j=0) {return Cosinc[i]*Surfdens_Bl[i][j];}

private:

    // Input parameters
    MomentMap<T> *im;
    size_t  Nrad;
    size_t  Nseg;
    T       *Radius;        /* Radii */
    T       *Radius_kpc;    /* Radii in kpc (need distance)*/
    T       *Width;         /* Ring widths */
    T       *Phi;           /* Pos. angle of major axis */
    T       *Inc;           /* Inclination of object */
    T       **Annuli;       /* Inner and outer radius of ring */
    T       Position[2];    /* Central position of all ellipses */
    float   Mass;           /* Mass of galaxy in Msun */
    float   Distance;       /* Distance to galaxy in Mpc */
    float   Range[2];       /* Min/max value for a pixel to be accepted */
    float   *Segments;            /* An array of size 2*Nseg with angle intervals */

    // Output parameters
    double  **Sum;                /* Sums of pixels values */
    double  **Sumsqr;             /* Squared sums */
    long    **Num;                /* Number of pixels in each ring */
    long    **Numblanks;          /* Number of 'subpixel' BLANK hits */
    double  **Datamin;            /* Minimum value  per ring per segment */
    double  **Datamax;            /* Maximum value  per ring per segment */
    
    double  **Mean;               /* Means of pixel values */
    double  **Median;             /* Medians of pixel values */
    double  **Var;                /* Store the variances */
    double  **MAD;                /* Median absolute deviation from the median */
    double  **Area;               /* Areas covered by pixels */
    double  **Blankarea;          /* Areas covered by BLANK pixels */
    double  **Surfdens;           /* Surface density in BUNIT / arcsec2 */
    double  **Surfdens_Bl;        /* Same as above, but averaging on blank too */
    double  *Mass_Surfdens;       /* Mass surface density in Msun/pc2 knowing the total mass*/
    double  *Mass_Surfdens_Bl;    /* Same as above, but averaging on blank too */

    // Working parameters
    double *Cosphi;
    double *Sinphi;
    double *Cosinc;
    float   Dx;             /* Spacings in arcsec*/
    float   Dy;
    T       Rmax;
    bool    Overlap;        /* Weight data in overlapping regions? */
    int     subpix[2];
    float   stepxy[2];
    float   maprotation;
    long    **Contrib;      /* Number of different pixels in a ring/segment */
    std::vector<double> **medianArray;


    void   defaults();
    void   init(MomentMap<T> *image, Rings<T> *rings, size_t nseg, float* segments);
    void   allocateArrays (size_t nrad, size_t nseg);
    void   deallocateArrays ();
    void   processpixel(int x, int y, float imval);
    bool   IsInRange(float value, float *Range);
    bool   IsInRing( float  Xr, float  Yr, int radnr);
    float  gettheta(float X,float Y,float Phi,float Crota);
    bool   IsInSegment(float Angle, float Segm1, float Segm2);
    float  toangle(float Angle);
};
}

#endif
