// -----------------------------------------------------------------------
// param.hh: Definition of the Param class
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

#ifndef PARAM_HH_
#define PARAM_HH_

#include <iostream>
#include <string>
#include <vector>

using namespace std;

// Container for input parameters for GALMOD (generates a 3D model galaxy)
struct GALMOD_PAR {
    bool   flagGALMOD = false;  ///< Enable task 
    int    NRADII     = -1;     ///< Number of rings.
    string RADII      = "-1";   ///< A list of radii.
    string XPOS       = "-1";   ///< X center of the galaxy (pixel).
    string YPOS       = "-1";   ///< Y center of the galaxy (pixel).
    double RADSEP     = -1;     ///< Separation between rings (arcs).
    string VSYS       = "-1";   ///< Systemic velocity (km/s).
    string VROT       = "-1";   ///< Circular velocity (km/s).
    string VRAD       = "-1";   ///< Radial velocity (km/s).
    string VVERT      = "-1";   ///< Vertical velocity (km/s)
    string VDISP      = "-1";   ///< Rotation velocity (km/s).
    string INC        = "-1";   ///< Inclination angle (degrees).
    string PHI        = "-1";   ///< Position angle from north anti-clockwise.
    string Z0         = "-1";   ///< Height scale of the disk (arcs).
    string DENS       = "-1";   ///< Column density of gas (atoms/cm2).
    string DVDZ       = "-1";   ///< Vertical rotational gradient (km/s/arcsec).
    string ZCYL       = "-1";   ///< Height where the rotational gradient starts.
    int    CDENS      = 10;     ///< Surface density of clouds in the plane of ring (1e20).
    int    CMODE      = 1;      ///< Mode for column density distribution.
    int    ISEED      = -1;     ///< Seed for random number generator
    int    LTYPE      = 1;      ///< Layer type along z.
    int    NV         = -1;     ///< Number of subclouds per profile.
    double REDSHIFT   = 0;      ///< Redshift of the galaxy.
    vector<double> RESTWAVE = {-1}; ///< Rest wavelengths.
    vector<double> RESTFREQ = {-1}; ///< Rest frequencies.
    vector<double> RELINT   = {1}; ///< Relative intensities of lines.
    bool   SM         = true;   /// If false, disable smoothing.
    double NOISERMS   = 0;      /// RMS noise to be added to the cube.
    string ringfile;            /// A BBarolo output ring file.
};

// Container for input parameters for GALFIT (fit a 3D galaxy model)
struct GALFIT_PAR : GALMOD_PAR {
    bool   flagGALFIT = false;    ///< Enable task 
    float  DELTAINC   = 5;        ///< Inclination angle variation (degrees).
    float  DELTAPHI   = 15;       ///< Position angle variation (degrees).
    float  DELTAVROT  = 1000.;    ///< Rotation velocity variation (km/s).
    float  MAXVDISP   = 1000.;    ///< Maximum velocity dispersion.
    float  MINVDISP   = 0.1;      ///< Minimum velocity dispersion.
    int    FTYPE      = 2;        ///< Type of function to be minimized;
    int    WFUNC      = 2;        ///< Weighting function.
    double TOL        = 1E-03;    ///< Tolerance for minimization.
    string NORM       = "AZIM";  ///< Normalization type: LOCAL, AZIM or NONE.
    string FREE       = "VROT VDISP INC PA"; ///< Free parameters.
    string SIDE       = "B";      ///< Approaching(A), Receding(R), Both(B), Single(S)
    bool   TWOSTAGE   = true;     ///< Whether fitting a second model after regularization.
    string REGTYPE    = "auto";   ///< Regularazation type for second stage
    float  BWEIGHT    = 1;        ///< Power of the weighting function for Blank pixels.
    int    STARTRAD   = 0;        ///< Starting radius
    bool   flagERRORS = false;    ///< Whether estimating errors.
    float  DISTANCE   = -1;       ///< Distance of the galaxy to convert arcs to kpc.
    bool   flagADRIFT = false;    ///< Whether correcting for asymmetric drift.
    int    ADRIFTPOL1 = 3;        ///< Degree of polynomial for fit of sigma in ADRIFT
    int    ADRIFTPOL2 = 3;        ///< Degree of polynomial for fit of log(density*disp2) in ADRIFT
    string REVERSE    = "false";  ///< Whether to use a reverse cumulative fitting.
    bool   NORMALCUBE = true;     ///< Whether to normalize the input flux values.
    bool   flagBADOUT = false;    ///< Whether to write bad rings (with flag) in output.
    bool   PLOTMASK   = false;    ///< Whether to show the mask in output plots
    float  PLOTMINCON = -1;       ///< Whether to show the mask in output plots
};

// Container for input parameters for GALWIND (generates a 3D biconical outflow)
struct GALWIND_PAR : GALMOD_PAR {
    bool   flagGALWIND =  false; ///< Enable task 
    string VWIND       = "-1";   ///< Outflow velocity (purely radial)
    string VROT        = "0";    ///< Wind can also rotate (km/s).
    string OPENANG     = "-1";   ///< Wind opening angle
    float  HTOT        = -1;     ///< Maximum height reached by the wind
    size_t NTOT        = 25;     ///< Number of layers/cylinders for each cone
    int    LTYPE       = 5;      ///< Layer type along z.    
    size_t DENSTYPE    = 1;      ///< Integer. How to distribute density in layers.
                                 ///   0=constant mass in layers     (Dens decreases)
                                 ///   1=constant density in layers  (Mass increases)
                                 ///   2=hollow cone with constant density
    short WTYPE        = 0;      ///< Type of wind: 0=cylindrical, 1=spherical
};


// Container for input parameters for the source finder 
struct SEARCH_PAR {
    bool   flagSearch        = false;     ///< Should search for sources in cube?
    string searchType        = "spatial"; ///< "Spectral" or "Spatial" search?
    float  snrCut            = 4.0;       ///< Signal to Noise for detection when sigma-clipping.
    bool   iternoise         = false;     ///< Whether to estimate noise iteratively.
    float  threshold         = 0;         ///< What the threshold is (when sigma-clipping).
    bool   UserThreshold     = false;     ///< Whether the user has defined a threshold of their own.
    bool   flagAdjacent      = true;      ///< Use the adjacent criterion for objects merger?
    int    threshSpatial     = -1;        ///< Maximum spatial separation between objects.
    int    threshVelocity    = 2;         ///< Maximum channels separation between objects.
    int    minChannels       = 2;         ///< Minimum channels to make an object.
    bool   RejectBeforeMerge = true;      ///< Whether to reject sources before merging.
    bool   TwoStageMerging   = true;      ///< Whether to do a partial merge during search.
    int    minVoxels         = -1;        ///< Minimum voxels required in an object.
    int    minPix            = -1;        ///< Minimum pixels required in an object.
    int    maxChannels       = -1;        ///< Maximum channels to accept an object.
    float  maxAngSize        = -1;        ///< Maximum angular size in the object in arcmin.
    bool   flagGrowth        = true;      ///< Are we growing objects once they are found?
    float  growthCut         = 2.5;       ///< The SNR that we are growing objects down to.
    bool   flagUserGrowthT   = false;     ///< Whether the user has manually defined a threshold
    float  growthThreshold   = 0;         ///< The threshold for growing objects down to
    bool   cubelets          = false;     ///< If true, produce a sub-cube of each detection.
    int    edges             = 20;        ///< Number of pixels at the edges of a cubelet.
    string sortsrcs          = "nvox";    ///< Type of sorting for detections.
    int    pbcorr            = 1;         ///< Whether to correct for primary beam attenuation.

};

// Container for input parameters for maps
struct MAPS_PAR {
    string maptype        = "MOMENT";       ///< How to extract kinematic map: GAUSSIAN OR MOMENT
    bool   globprof       = false;          ///< Whether to calculate the global profile.
    bool   massdensmap    = false;          ///< Whether to calculate the mass density HI map.
    bool   totalmap       = false;          ///< Whether to calculate the total map.
    bool   velocitymap    = false;          ///< Whether to calculate the velocity field.
    bool   dispersionmap  = false;          ///< Whether to calculate velocity dispersion field.
    bool   rmsmap         = false;          ///< Whether to calculate the RMS map.
    bool   SNmap          = false;          ///< Whether to calculate A S/N map of the moment 0 map. 
    vector<int> contChans = {1000000};      ///< Number of channels used for continuum subtraction.
    string taper          = "uniform";      ///< Whether input cube has been Hanning tapered.
    string veldef         = "auto";         ///< Velocity conversion when spectral axis not velocity.
};


class Param              /// Class to store general parameters for main program.
{
public:
     
    Param();                                /// Default constructor.
    Param(const Param& p);                  /// Copy constructor.
    Param& operator= (const Param& p);      /// Assignement operator =.
    virtual ~Param() {}                     /// Default destructor.
    void  defaultValues();                  /// Default values for constructor.
    
    /// Inline functions to access private member:
    
    string  getImageFile () {return imageFile;}    
    string  getImageList () {return imageList;}
    string  getImage (int i) {return images[i];}  
    int     getListSize () {return images.size();}
    void    setImageFile (std::string fname) {imageFile = fname;}     
    void    setImage (std::string fname) {images.push_back(fname);}
    string  getOutfolder () {return outFolder;}
    void    setOutfolder (std::string s) {if (s!="" && s[s.size()-1]!='/') s.append("/"); outFolder=s;}
    string  getOutPrefix () {return outPrefix;}
    void    setOutPrefix (std::string s) {outPrefix=s;}
    bool    getLogFile () {return logFile;}
    bool    isVerbose () {return verbose;}
    void    setVerbosity (bool f) {verbose=f;}
    bool    getShowbar () {return showbar;}
    void    setShowbar (bool s) {showbar = s;}
    int     getThreads () {return threads;}
    void    setThreads (int t) {threads=t;}
    bool    getFlagDebug() {return debug;}
    int     getFlagPlots() {return plots;}
    bool    getFlagAuto() {return AUTO;}
    bool    getFlagStats() {return flagStats;}
    bool    getFluxConvert() {return fluxConvert;}

    bool    getFlatContsub() {return contsub;}
    string  getExcludeWind() {return exclude_windows;}
    int     getContOrder() {return cont_order;}
    
    bool    getMakeMask() {return makeMask;}
    string  getMASK() {return MaskType;}
    void    setMASK(string s) {MaskType=s;}
    void    setBlankCut (float f) {blankCut=f;}
    float   getBlankCut() {return blankCut;}
    string  getVelDef() {return parMA.veldef;}

    float   getBeamFWHM() {return beamFWHM;}
    void    setBeamFWHM(float val) {beamFWHM=val;}
    int     getCheckCube () {return checkCube;}
    void    setCheckCube (bool a) {checkCube=a;}
    
    bool    getFlagRobustStats () {return flagRobustStats;}
    void    setFlagRobustStats (bool flag) {flagRobustStats=flag;}
    
    bool    getMaps() {return (parMA.globprof || parMA.totalmap || parMA.velocitymap || parMA.dispersionmap || parMA.massdensmap || parMA.rmsmap);}
    
    bool    getFlagRing () {return flagRing;}
    void    setFlagRing (bool b) {flagRing = b;}

    int     getBOX  (int i) {return BOX[i];}

    bool    getflagSearch () {return parSE.flagSearch;}
    bool    getflagGalFit () {return parGF.flagGALFIT;}
    bool    getflagGalMod () {return parGM.flagGALMOD;}
    
    // This should go on a different param struct
    double  getRedshift() {return parGF.REDSHIFT;}
    double  getRestwave() {return parGF.RESTWAVE[0];}
    double  getRestfreq() {return parGF.RESTFREQ[0];}
    double  getDistance () {return parGF.DISTANCE;}
    
    
    GALMOD_PAR&  getParGM() {return parGM;}
    GALFIT_PAR&  getParGF() {return parGF;}
    GALWIND_PAR& getParGW() {return parGW;}
    SEARCH_PAR&  getParSE() {return parSE;}
    MAPS_PAR&    getParMA() {return parMA;}
    
    bool    getflagSpace () {return flagSpace;}
    string  getP1 () {return P1;}
    string  getP2 () {return P2;}
    float   getP1p (int i) {return P1p[i];}
    float   getP2p (int i) {return P2p[i];}
    
    bool    getflagSmooth () {return flagSmooth;}
    double  getBmaj () {return bmaj;}
    double  getBmin () {return bmin;}
    double  getBpa  () {return bpa;}
    double  getOBmaj    () {return obmaj;}
    double  getOBmin    () {return obmin;}
    double  getOBpa () {return obpa;}
    bool    getflagFFT () {return flagFFT;}
    void    setflagFFT (bool f) {flagFFT=f;}
    float   getLinear () {return linear;}
    void    setLinear (float a) {linear=a;}
    float   getFactor () {return factor;}
    float   getScaleFactor () {return scalefactor;}
    void    setFactor (float f) {factor=f;}
    bool    getflagReduce() {return (reduce=="false" ? false : true);}
    void    setflagReduce(bool f) {f ? reduce="true" : reduce="false";}
    string  getReduce() {return reduce;}
    string  getSmoothOut () {return smo_out;}
    
    bool    getflagSmoothSpectral () {return flagSmoothSpectral;}
    string  getWindowType () {return window_type;}
    float   getWindowSize () {return window_size;}
    
    bool    getFlagSlitfit () {return flagSlitfit;}
    float   getSlitWidth () {return slitwidth;}

    bool    getFlagPV() {return flagPV;}
    string  getXPOS_PV() {return XPOS_PV;}
    string  getYPOS_PV() {return YPOS_PV;}
    float   getPA_PV() {return PA_PV;}
    float   getP1_PV (int i) {return P1_PV[i];}
    float   getP2_PV (int i) {return P2_PV[i];}
    float   getWIDTH_PV () {return WIDTH_PV;}
    float   getANTIALIAS () {return ANTIALIAS;}
    void    setANTIALIAS (float i) {ANTIALIAS=i;}
    
    bool    getFlagRend3D() {return flagRend3D;}
    float   getRendAngle() {return rendangle;}

    bool    getFlagEllProf() {return flagEllProf;}
    
    /// Utility functions:
    
    bool    getopts(int argc, char **argv);          /// Parse the command line parameters correctly. 
    bool    readParamFile(string paramfile);         /// Read in parameters from a disk file.
    bool    readParamString(string paramstring);     /// Read in parameters from a string.
    bool    readParamCL(std::string parstr);         /// Read a parameter from the command line.
    void    setParam(string &parstr);                /// Set a parameter value.
    bool    checkPars();                             /// Check the parameter list for inconsistencies. 
    void    printDefaults (ostream& theStream=cout, string wtask="ALL"); /// Print on screen the defaults values.
    void    createTemplate();                        /// Create a template file for 3DFIT.
    friend ostream& operator<< (ostream& theStream, Param& par);
    
    friend class Image;
  
private:
    string          imageFile;          ///< The image to be analysed.
    string          imageList;          ///< A file with list of images to be analized.
    vector<string>  images;             ///< A vector with single images in the list.
    string          outFolder;          ///< Folder where saving output files.
    string          outPrefix;          ///< Prefix for output files in the outFolder
    bool            logFile;            ///< A log file to redirect std::cout and std::cerr.
    bool            verbose;            ///< Is verbosity activated?
    bool            showbar;            ///< Show progress bar?
    int             checkCube;          ///< Checking for bad channels/rows/cols in the cube?
    float           beamFWHM;           ///< Beam to adopt if any information in header.
    bool            flagRobustStats;    ///< Whether to use robust statistics.
    int             plots;              ///< Whether producing output plots.
    bool            flagStats;          ///< Whether to calculate and return stats
    bool            fluxConvert;        ///< Whether to convert fluxes to Jy, when possible
    
    bool            contsub;            ///< Whether to subtract continuum from a cube. 
    string          exclude_windows;    ///< Exclude channels for continuum subtraction.
    int             cont_order;         ///< Order of polynomial fit for continuum.
    
    bool            makeMask;           ///< Whether to write a mask.
    string          MaskType;           ///< Type of mask: SEARCH,SMOOTH,THRESHOLD,NEGATIVE,SMOOTH&SEARCH or NONE.
    float           blankCut;           ///< SNR clipping cut for blanked area when mask==SMOOTH

    bool            flagRing;           ///< Do you want to fit a tilted ring model?

    int             BOX[6];             ///< A box in RA-DEC-VELO. Not used anymore!

    SEARCH_PAR      parSE;              ///< Input parameters for the SEARCH task
    GALMOD_PAR      parGM;              ///< Input parameters for the GALMOD task
    GALFIT_PAR      parGF;              ///< Input parameters for the GALFIT task
    GALWIND_PAR     parGW;              ///< Input parameters for the GALWIND task
    MAPS_PAR        parMA;              ///< Input parameter for MAPS
    
    bool            flagSpace;
    string          P1;
    string          P2;
    float           P1p[3];
    float           P2p[3];

    bool            flagSmooth;
    double          bmaj;               ///< Beam of the smoothed array (arcs).
    double          bmin;               ///< Beam of the smoothed array (arcs).
    double          bpa;                ///< Beam of the smoothed array (deg).
    double          obmaj;              ///< Beam of the original array.
    double          obmin;
    double          obpa;
    bool            flagFFT;            ///< Using FFT for convolution?
    float           linear;             ///< Linear resolution to be achieved.
    float           factor;             ///< The newbeam is a factor of the old.
    float           scalefactor;
    string          reduce;             ///< Whether to average pixel spatially/spectrally.
    string          smo_out;            ///< Output file.
    
    bool            flagSmoothSpectral; ///< Spectral smoothing the datacube?
    string          window_type;        ///< Type of smoothing window
    float           window_size;        ///< Type of smoothing window

    bool            flagSlitfit;        ///< Fitting a 3D model to a slit observation.
    float           slitwidth;          ///< Width of the slit in arcsec
    
    bool            flagPV;             ///< Extracting a position-velocity diagram.
    string          XPOS_PV;            ///< X-center of PV slice
    string          YPOS_PV;            ///< Y-center of PV slice
    float           PA_PV;              ///< Position angle of PV slice
    float           WIDTH_PV;           ///< Width of the slice 
    float           P1_PV[2];           ///< Point 1 defining the PV slice
    float           P2_PV[2];           ///< Point 2 defining the PV slice
    float           ANTIALIAS;          ///< Anti-aliasing to use
    
    bool            flagEllProf;
    int             threads;
    bool            debug;
    bool            AUTO;               ///< Whether it is an automated run.

    bool            flagRend3D;         ///< Whether to perform 3D rendering.
    float           rendangle;          ///< Azimuth angle for 3D rendering.

};

  /// Write out info on a parameter to e.g. the results file.
void recordParameters(ostream& theStream, string paramName, string paramDesc, string paramValue);

  /// A macro to handle streamed output to recordParameters.
#define recordParam(outstream,string1,string2,instream)          \
    do {                                                         \
        ostringstream oss;                                       \
        oss<<instream;                                           \
        recordParameters(outstream,string1,string2,oss.str());   \
    } while(0)

// Some utility functions
void listTasks (std::ostream& Str=cout);
void listFitsUtils(std::ostream& Str);
void helpscreen(ostream& Str=cout);
void versionInfo(ostream& Str, char ** argv);
void welcomeMessage(std::ostream& ostr=std::cout);
void printParams (ostream& Str, Param &p, bool defaults=false, string whichtask="ALL");


#endif
