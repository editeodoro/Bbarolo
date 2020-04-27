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
    int    LTYPE      = 1;      ///< Layer type along z.
    int    NV         = -1;     ///< Number of subclouds per profile.
    double REDSHIFT   = 0;      ///< Redshift of the galaxy.
    vector<double> RESTWAVE = {-1}; ///< Rest wavelengths.
    vector<double> RESTFREQ = {-1}; ///< Rest frequencies.
    vector<double> RELINT   = {1}; ///< Relative intensities of lines.
    bool   SM         = true;   /// If false, disable smoothing.
    double NOISERMS   = 0;      /// RMS noise to be added to the cube.
};

// Container for input parameters for GALFIT (fit a 3D galaxy model)
struct GALFIT_PAR : GALMOD_PAR {
    bool   flagGALFIT = false;    ///< Enable task 
    float  DELTAINC   = 5;        ///< Inclination angle variation (degrees).
    float  DELTAPHI   = 15;       ///< Position angle variation (degrees).
    float  DELTAVROT  = 1000;     ///< Rotation velocity variation (degrees).
    int    FTYPE      = 2;        ///< Type of function to be minimized;
    int    WFUNC      = 2;        ///< Weighting function.
    double TOL        = 1E-03;    ///< Tolerance for minimization.
    string NORM       = "LOCAL";  ///< Normalization type: LOCAL, AZIM or NONE.
    string FREE       = "VROT VDISP INC PA"; ///< Free parameters.
    string SIDE       = "B";      ///< Approaching(A), Receding(R), Both(B), Single(S)
    bool   TWOSTAGE   = true;     ///< Whether fitting a second model after regularization.
    string POLYN      = "bezier"; ///< Degree of polynomials fitting INC e PA.
    int    BWEIGHT    = 1;        ///< Power of the weighting function for Blank pixels.
    int    STARTRAD   = 0;        ///< Starting radius
    bool   flagERRORS = false;    ///< Whether estimating errors.
    float  DISTANCE   = -1;       ///< Distance of the galaxy to convert arcs to kpc.
    bool   flagADRIFT = false;    ///< Whether correcting for asymmetric drift.
    bool   PLOTMASK   = false;    ///< Whether to show the mask in output plots
    bool   CUMULATIVE = false;    ///< Whether to show the mask in output plots
    bool   NORMALCUBE = true;     ///< Whether to normalize the input flux values.
    
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
    float  snrCut            = 5;         ///< Signal to Noise for detection when sigma-clipping.
    float  threshold         = 0;         ///< What the threshold is (when sigma-clipping).
    bool   UserThreshold     = false;     ///< Whether the user has defined a threshold of their own.
    bool   flagAdjacent      = true;      ///< Use the adjacent criterion for objects merger?
    int    threshSpatial     = -1;        ///< Maximum spatial separation between objects.
    int    threshVelocity    = -1;        ///< Maximum channels separation between objects.
    int    minChannels       = -1;        ///< Minimum channels to make an object.
    bool   RejectBeforeMerge = true;      ///< Whether to reject sources before merging.
    bool   TwoStageMerging   = true;      ///< Whether to do a partial merge during search.
    int    minVoxels         = -1;        ///< Minimum voxels required in an object.
    int    minPix            = -1;        ///< Minimum pixels required in an object.
    int    maxChannels       = -1;        ///< Maximum channels to accept an object.
    float  maxAngSize        = -1;        ///< Maximum angular size in the object in arcmin.
    bool   flagGrowth        = true;      ///< Are we growing objects once they are found?
    float  growthCut         = 3;         ///< The SNR that we are growing objects down to.
    bool   flagUserGrowthT   = false;     ///< Whether the user has manually defined a threshold
    float  growthThreshold   = 0;         ///< The threshold for growing objects down to
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
    bool    isVerbose () {return verbose;}
    void    setVerbosity (bool f) {verbose=f;}
    bool    getShowbar () {return showbar;}
    int     getThreads () {return threads;}
    void    setThreads (int t) {threads=t;}
    bool    getFlagDebug() {return debug;}
    bool    getFlagPlots() {return plots;}
    
    bool    getMakeMask() {return makeMask;}
    string  getMASK() {return MaskType;}
    void    setMASK(string s) {MaskType=s;}
        
    float   getBeamFWHM() {return beamFWHM;}
    void    setBeamFWHM(float val) {beamFWHM=val;}
    bool    getCheckCh () {return checkChannels;}
    void    setCheckCh (bool flag) {checkChannels=flag;}
    
    bool    getFlagRobustStats () {return flagRobustStats;}
    void    setFlagRobustStats (bool flag) {flagRobustStats=flag;}
    
    void    setGlobProf (bool flag) {globprof = flag;}
    bool    getGlobProf () {return globprof;}
    void    setTotalMap (bool flag) {totalmap = flag;}
    bool    getTotalMap () {return totalmap;}
    bool    getMassDensMap () {return massdensmap;}
    void    setVelMap (bool flag) {velocitymap = flag;}
    bool    getVelMap () {return velocitymap;}
    void    setDispMap (bool flag) {dispersionmap = flag;}
    bool    getDispMap () {return dispersionmap;}
    bool    getRMSMap () {return rmsmap;}
    bool    getMaps() {return (globprof || totalmap || velocitymap || dispersionmap || massdensmap || rmsmap);}
    void    setBlankCut (float f) {blankCut=f;}
    float   getBlankCut() {return blankCut;}
    string  getMapType() {return maptype;}
    
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
    float   getFactor () {return factor;}
    float   getScaleFactor () {return scalefactor;}
    void    setFactor (float f) {factor=f;}
    bool    getflagReduce() {return flagReduce;}
    string  getSmoothOut () {return smo_out;}
    
    bool    getflagHanning () {return flagHanning;}
    size_t  getHanningWindow () {return hanning_window;}
    

    bool    getFlagSlitfit () {return flagSlitfit;}
    string  getWavefile () {return wavefile;}
    string  getIvarfile () {return ivarfile;}
    string  getLine() {return linetofit;}

    bool    getFlagPV() {return flagPV;}
    float   getXPOS_PV() {return XPOS_PV;}
    float   getYPOS_PV() {return YPOS_PV;}
    float   getPA_PV() {return PA_PV;}

    bool    getFlagEllProf() {return flagEllProf;}
    
    /// Utility functions:
    
    bool    getopts(int argc, char **argv);          /// Parse the command line parameters correctly. 
    int     readParams(string paramfile);            /// Read in parameters from a disk file.
    void    setParam(stringstream &ss);               /// Set a parameter value.
    bool    checkPars();                             /// Check the parameter list for inconsistencies. 
    void    printDefaults (ostream& theStream=cout, string wtask="ALL"); /// Print on screen the defaults values.
    void    createTemplate();                        /// Create a template file for 3DFIT.
    void    overrideParameter(std::string parstr);         /// Override a parameter in parameter file.
    friend ostream& operator<< (ostream& theStream, Param& par);
    
    friend class Image;
  
private:
    string          imageFile;          ///< The image to be analysed.
    string          imageList;          ///< A file with list of images to be analized.
    vector<string>  images;             ///< A vector with single images in the list.
    string          outFolder;          ///< Folder where saving output files.
    bool            verbose;            ///< Is verbosity activated?
    bool            showbar;            ///< Show progress bar?
    bool            checkChannels;      ///< Checking for bad channels in the cube?
    float           beamFWHM;           ///< Beam to adopt if any information in header.
    bool            flagRobustStats;    ///< Whether to use robust statistics.
    bool            plots;              ///< Whether producing output plots.
    
    bool            makeMask;           ///< Whether to write a mask.
    string          MaskType;           ///< Type of mask: SEARCH,SMOOTH,THRESHOLD,NEGATIVE, SMOOTH&SEARCH or NONE.
    
    bool            globprof;           ///< Whether the user wants the global profile.
    bool            massdensmap;        ///< Whether the user wants the mass density HI map.
    bool            totalmap;           ///< Whether the user wants the total HI map.
    bool            velocitymap;        ///< Whether the user wants the velocity field.
    bool            dispersionmap;      ///< Whether the user wants the velocity dispersion field.
    bool            rmsmap;             ///< Whether the user wants the RMS map.
    float           blankCut;           ///< SNR clipping cut for blanked area.
    string          maptype;            ///< How to extract kinematic map: GAUSSIAN OR MOMENT

    bool            flagRing;           ///< Do you want to fit a tilted ring model?

    int             BOX[6];             ///< A box in RA-DEC-VELO. Not used anymore!

    SEARCH_PAR      parSE;              ///< Input parameters for the SEARCH task
    GALMOD_PAR      parGM;              ///< Input parameters for the GALMOD task
    GALFIT_PAR      parGF;              ///< Input parameters for the GALFIT task
    GALWIND_PAR     parGW;              ///< Input parameters for the GALWIND task
    
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
    bool            flagReduce;
    string          smo_out;            ///< Output file.
    
    bool            flagHanning;        ///< Hanning smooth the datacube?
    size_t          hanning_window;     ///< Size of hanning window

    bool            flagSlitfit;        ///< Fitting a 3D model to a slit observation.
    string          wavefile;           ///< Fitsfile containing the wavelegths.
    string          ivarfile;           ///< Fitsfile containing the inverse-variance or a value.
    string          linetofit;          ///< Line to fit: Ha, Hb, OIII etc...
    
    bool            flagPV;                 ///< Extracting a position-velocity diagram.
    float           XPOS_PV;
    float           YPOS_PV;
    float           PA_PV;

    bool            flagEllProf;
    int             threads;
    bool            debug;
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
void helpscreen(ostream& Str=cout);
void versionInfo(ostream& Str, char ** argv);
void welcomeMessage(std::ostream& ostr=std::cout);
void printParams (ostream& Str, Param &p, bool defaults=false, string whichtask="ALL");


#endif
