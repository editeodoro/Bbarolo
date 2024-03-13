//-----------------------------------------------------------------------
// header.hh: Definition of Header class, a class to collect header info.
//-----------------------------------------------------------------------

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

#ifndef HEADER_HH_
#define HEADER_HH_

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fitsio.h> 
#include <wcslib/wcs.h>

//enum CUNITS {DEGREE,ASEC,AMIN,M_S,KM_S,HZ,KHZ,MHZ,GHZ,MUM};
//enum CTYPES {RA,DEC,LAT,LONG,VELO,FREQ,WAVE};

class Header
{
public: 
    Header();                                           /// Default constructor.
    Header(std::string f) : Header() {header_read(f);}  /// Alternative constructor.
    virtual ~Header();                                  /// Default destructor.
    Header(const Header& h);                            /// Copy constructor.
    Header& operator= (const Header& h);                /// Assignement operator.
     
    /// Obvious inline functions
    
    int    NumAx   () {return numAxes;}
    int    Bitpix  () {return bitpix;}
    long   *DimAx  () {return dimAxes;}
    long   DimAx (int i) {return dimAxes[i];}
    double *Crpix  () {return crpix;}
    double *Crval  () {return crval;}
    double *Cdelt  () {return cdelt;}
    double Bmin    () {return bmin;}
    double Bmaj    () {return bmaj;}
    double Bpa     () {return bpa;}
    float  BeamArea() {return beamArea;}
    float  Bzero   () {return bzero;}
    float  Bscale  () {return bscale;}
    float  Blank   () {return blank;}
    float  Epoch   () {return epoch;}
    double Freq0   () {return freq0;}
    double Crota   () {return crota;}
    double DataMax () {return datamax;}
    double DataMin () {return datamin;}
    double Cdelt   (int i) {return cdelt[i];}
    double Crval   (int i) {return crval[i];}
    double Crpix   (int i) {return crpix[i];}
    double Drval3  () {return drval3;}
    double PixScale () {return (fabs(cdelt[0])+fabs(cdelt[1]))/2.;}
    struct wcsprm *WCS () {return wcs;}
    double Wave0 () {return wave0;}
    double Redshift () {return redshift;}

    std::vector<std::string>& Keys () {std::vector<std::string> &k=keys; return k;}
    std::string Name () {return object;}
    std::string Bunit () {return bunit;}
    std::string Btype () {return btype;}
    std::string Ctype (int i) {return ctype[i];}
    std::string* Ctype () {return ctype;}
    std::string Cunit (int i) {return cunit[i];}
    std::string Dunit3 () {return dunit3;}
    std::string Obname () {return object;}
    std::string Telesc () {return telescope;}
    std::string VelDef () {return veldef;}
    std::string SpectralType() {return sptype;}
    
    void setBitpix (int i) {bitpix = i;}
    void setDimAx (int i, long val) {dimAxes[i] = val;}
    void setCrpix (int i, float val) {crpix[i]=val;}
    void setCrval (int i, float val) {crval[i]=val;}
    void setCdelt (int i, float val) {cdelt[i]=val;}
    void setDrval3 (double val) {drval3=val;}
    void setDunit3 (std::string s) {dunit3=s;}
    void setBmaj  (float val) {bmaj = val;}
    void setBmin  (float val) {bmin = val;}
    void setBpa   (float val) {bpa = val;}
    void setBeam  (float a, float b, float c) {bmaj=a; bmin=b; bpa=c; calcArea();}
    void setBzero (float val) {bzero = val;}
    void setBscale(float val) {bscale = val;}
    void setBlank (float val) {blank = val;}
    void setEpoch (float val) {epoch = val;}
    void setBunit (std::string ch) {bunit = ch;}
    void setDataMax (double val) {datamax=val;}
    void setDataMin (double val) {datamin=val;}
    void setMinMax (double minn, double maxx) {datamin=minn;datamax=maxx;}
    void setFreq0 (double val) {freq0=val;}
    void setCtype (int i, std::string s) {ctype[i] = s;}
    void setCunit (int i, std::string s) {cunit[i] = s;}
    void setBtype (std::string s) {btype = s;}
    void setName  (std::string s) {object = s;}
    void setTelesc(std::string s) {telescope = s;}
    void setFitsName(std::string s) {fitsname = s;}
    void setPointAllocated (bool b) {pointAllocated=b;}
    void setWarning (bool b) {warning=b;}
    void setWave0 (double w) {wave0=w;}
    void setRedshift (double r) {redshift=r;}
    void setVelDef(std::string s) {veldef=s;}
        
    void Warning(std::string s) {if (warning) std::cout << s << std::endl;}
    void addKey(std::string s) {keys.push_back(s);}


    /// Functions defined in header.cpp.
    
    void    setNumAx (int n);
    void    calcArea ();                                        /// Calculate beam area from bmaj & bmin.
    bool    header_read (std::string fname);                    /// Read from header of a Fits file.
    void    headwrite (fitsfile *fptr, short numDim, bool fullHead); /// Write header of a Fits cube.
    void    updateWCS();                                        /// Update WCS structure
    bool    checkHeader();                                      /// Check header is ok for BBarolo
    int     wcsToPix(const double *world, double *pix, size_t npts=1);
    int     pixToWCS(const double *pix, double *world, size_t npts=1);

    template <class T>                                          /// Read the request keyword and write on "key".
    bool read_keyword(std::string keyword, T &key, bool err=false); 
    
private:
    int     numAxes;                ///< Number of axes.
    int     bitpix;                 ///< Image type.
    long    *dimAxes;               ///< Dimensions of axes.
    double  *crpix;                 ///< Central pixels.
    double  *crval;                 ///< Values of central pixels.
    double  *cdelt;                 ///< Delta pixel.
    double  drval3;                 ///< Secondary reference value of third axis.
    bool    pointAllocated;         ///< Have been the pointers allocated?
    double  bmaj;                   ///< The major main-beam FWHM.
    double  bmin;                   ///< The minor main-beam FWHM.
    double  bpa;                    ///< The beam position angle.
    float   bzero;                  ///< Bias for real values.
    float   bscale;                 ///< Scale for physical values.
    float   blank;                  ///< Value for blank pixel.
    float   beamArea;               ///< The area of the beam.
    float   epoch;                  ///< Epoch for coordinates.
    double  freq0;                  ///< Frequency at rest.
    double  wave0;                  ///< Wavelength at rest
    double  crota;                  ///< Rotation angle.
    double  datamin;                ///< Minimum pixel value.
    double  datamax;                ///< Maximum data value.
    double  redshift;
    std::string fitsname;           ///< The name of the fitsfile.
    std::string btype;              ///< Beam type.
    std::string bunit;              ///< Units of pixel value.
    std::string object;             ///< The name of the object.
    std::string *ctype;             ///< Type of axis.
    std::string *cunit;             ///< Unity of axis.
    std::string dunit3;             ///< Secondary units of third axis.
    std::string telescope;          ///< Instrument.
    std::string veldef;             ///< Velocity definition.
    std::string sptype;             ///< Spectral type (frequency, wavelength, velocity radio or optical)
    std::string radesys;            ///< System in RA and DEC.
    std::vector<std::string> keys;  ///< Whole header as strings.

    struct wcsprm *wcs;             ///< The WCS parameters in a struct from the wcslib library.
    int    nwcs;                    ///< The number of WCS parameters
    bool   wcsIsGood;               ///< A flag indicating whether there is a valid WCS


    bool    warning;               ///< Write warning on std::cout.
};

#endif
