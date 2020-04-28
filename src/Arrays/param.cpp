// -----------------------------------------------------------------------
// param.cpp: Member functions for the Param class.
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <thread>
#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include <Arrays/param.hh>
#include <Utilities/utils.hh>

#define BBVERSION "1.5.2"

// A list of all tasks available in BBarolo
vector<string> tasks = {
     "3DFIT","GALMOD","SEARCH","2DFIT","SMOOTH","HANNING",
     "ELLPROF","MAKEMASK","SPACEPAR","GALWIND","SLITFIT",
     "GLOBALPROFILE","TOTALMAP","VELOCITYMAP","DISPERSIONMAP","PV",
 };

// A brief description of tasks
vector<string> taskdescr = {
     "Fitting a 3D tilted-ring model to an emission-line datacube.",
     "Making a 3D model observation datacube of an emission line.",
     "Searching for sources in a 3D datacube or in a 2D image.",
     "Fitting a 2D tilted-ring model to a velocity field.",
     "Spatial smoothing of a datacube with a 2D Gaussian kernel.",
     "Spectral smoothing of a datacube with a Hanning window.",
     "Deriving radial profiles in elliptical rings from a map.",
     "Defining a mask for emission-line data.",
     "Computing the full 2D space for any pair of 3DFIT parameters.",
     "Making a 3D emission-line model of a bi-conical outflow.",
     "Fitting the kinematics of slit data (CURRENTLY BROKEN).",
     "Extracting a total spectral profile from a datacube.",
     "Computing a total intensity map from a datacube.",
     "Computing a velocity field from a datacube.",
     "Computing a velocity dispersion field from a datacube.",
     "Extracting a position-velocity slice from a datacube.",
};


Param::Param() {
    
    defaultValues();
}

void Param::defaultValues() {
    
    imageFile           = "";
    imageList           = "NONE";
    outFolder           = "";
    verbose             = true;
    showbar             = true;
    plots               = true;
    beamFWHM            = 30.;
    checkChannels       = false;
    flagRobustStats     = true;

    makeMask            = false;
    MaskType            = "SMOOTH";
    
    globprof            = false;
    massdensmap         = false;
    totalmap            = false;
    velocitymap         = false;
    dispersionmap       = false;
    rmsmap              = false;
    blankCut            = 3;
    maptype             = "MOMENT";
    
    flagRing            = false;

    flagSmooth          = false;
    flagFFT             = true;
    bmaj                = -1;
    bmin                = -1;
    bpa                 = 0;
    obmaj               = -1;
    obmin               = -1;
    obpa                = 0;
    linear              = -1;
    factor              = 2;
    scalefactor         = -1;
    flagReduce          = false;
    smo_out             = "NONE";
    for (int i=0; i<6; i++) BOX[i] = -1;

    flagHanning         = false;
    hanning_window      = 3;
    
    flagSlitfit         = false;
    wavefile            = "NONE";
    ivarfile            = "NONE";
    linetofit           = "Ha";

    flagSpace           = false;
    P1                  = "VROT";
    P2                  = "VDISP";
    P1p[0] = P2p[0]     = 0;
    P1p[1] = P2p[1]     = 100;
    P1p[2] = P2p[2]     = 2;
    
    flagPV              = false;
    XPOS_PV             = 0;
    YPOS_PV             = 0;
    PA_PV               = 0;
    
    flagEllProf         = false;

    debug               = false;    
    threads             = 1;
#ifdef _OPENMP
    threads = std::thread::hardware_concurrency();
    if (threads==0) threads = 1;
#endif
}

  
Param::Param (const Param& p) {
    
    operator=(p);
}

  
Param& Param::operator= (const Param& p) {
    
    if(this == &p) return *this;
    this->imageFile         = p.imageFile;
    this->imageList         = p.imageList;
    for (unsigned int i=0; i<images.size(); i++) 
        this->images[i]     = p.images[i];
    this->outFolder         = p.outFolder;
    this->beamFWHM          = p.beamFWHM;
    this->checkChannels     = p.checkChannels;
    this->verbose           = p.verbose; 
    this->showbar           = p.showbar;
    this->plots             = p.plots;
    this->flagRobustStats   = p.flagRobustStats;

    this->makeMask          = p.makeMask;
    this->MaskType          = p.MaskType;
    
    this->globprof          = p.globprof;
    this->velocitymap       = p.velocitymap;
    this->totalmap          = p.totalmap;
    this->massdensmap       = p.massdensmap;
    this->dispersionmap     = p.dispersionmap;
    this->rmsmap            = p.rmsmap;
    this->blankCut          = p.blankCut;
    this->maptype           = p.maptype;
         
    this->flagRing          = p.flagRing;
 
    for (int i=0; i<6; i++) this->BOX[i] = p.BOX[i];
    
    this->parSE             = p.parSE;
    this->parGM             = p.parGM;
    this->parGF             = p.parGF;
    this->parGW             = p.parGW;
    
    this->flagSpace     = p.flagSpace;
    this->P1                = p.P1;
    this->P2                = p.P2;
    for (int i=0; i<3; i++) {
        this->P1p[i]        = p.P1p[i];
        this->P2p[i]        = p.P2p[i];
    }   
        
    this->flagSmooth            = p.flagSmooth;
    this->flagFFT           = p.flagFFT;
    this->bmaj              = p.bmaj;                       
    this->bmin              = p.bmin;
    this->bpa               = p.bpa;
    this->obmaj             = p.obmaj;
    this->obmin             = p.obmin;
    this->obpa              = p.obpa;
    this->linear            = p.linear; 
    this->factor            = p.factor; 
    this->scalefactor       = p.scalefactor;
    this->flagReduce        = p.flagReduce;
    this->smo_out           = p.smo_out;

    this->flagHanning       = p.flagHanning;
    this->hanning_window    = p.hanning_window;
    
    this->flagSlitfit       = p.flagSlitfit;
    this->wavefile          = p.wavefile;
    this->ivarfile          = p.ivarfile;
    this->linetofit         = p.linetofit;

    this->flagPV            = p.flagPV;
    this->XPOS_PV           = p.XPOS_PV;
    this->YPOS_PV           = p.YPOS_PV;
    this->PA_PV             = p.PA_PV;
    
    this->flagEllProf       = p.flagEllProf;
    
    this->threads           = p.threads;
    this->debug             = p.debug;

    return *this;
}
 
 
bool Param::getopts(int argc, char ** argv) {
    
  /// A function that reads in the command-line options, in a manner 
  /// tailored for use with the main program.
  /// 
  /// \param argc The number of command line arguments.
  /// \param argv The array of command line arguments.

    bool returnValue = false;
    int status = 0;
    if(argc<2){
        helpscreen();
        returnValue = false;
    }
    else {
        
        defaultValues();
        
        // List of short and long options
        const char* const short_opts = "p:f:d:qvlth";
        const option long_opts[] = {
                {"paramfile", required_argument, nullptr, 'p'},
                {"fitsfile",  required_argument, nullptr, 'f'},
                {"defaults",  required_argument, nullptr, 'd'},
                {"list",      no_argument,       nullptr, 'l'},
                {"version",   no_argument,       nullptr, 'v'},
                {"template",  no_argument,       nullptr, 't'},
                {"help",      no_argument,       nullptr, 'h'},
                {"quote",     no_argument,       nullptr, 'q'},
                {0, 0, 0, 0}
        };
        
        string file;
        char c=0;
        // Turning off error messages from getopt_long
        opterr = 0;
        while((c = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1){
            
            switch(c) {
                
            case 'p':                    // Parameter file provided
                file = optarg;
                status = readParams(file);
                if(status==1)
                    cout << "Could not open parameter file " << file << ".\n";
                else if(status==2)
                    cout << "\nError opening list file: " << imageList << " doesn't exist.\n\n";
                else if(status==3) 
                    cout << "\nWrong input file. Exiting...\n\n";
                else returnValue = true;
                break;
                
            case 'f':                    // Fits file provided
                file = optarg;
                images.push_back(file);
                beamFWHM /= 3600.;
                parSE.flagSearch=true;
                parGF.flagGALFIT=true;
                returnValue = true;
                break;
                
            case 'd':                    // Print default paramaters
                file = makeupper(optarg);
                printDefaults(cout,file);
                break;
                
            case 't':                    // Write template paramater file
                createTemplate();
                break;
            
            case 'l':                    // List available tasks
                listTasks(cout);
                break;
                
            case 'v':                    // Print version info
                versionInfo(cout, argv);
                break;
            
            case 'q':                    // Random citation
                //cout << " \"BBarolo is such a " << randomAdjective(1)
                //     << " code!\" (cit. " << randomAdjective(3) << ")\n";
                cout << endl << " " << randomQuoting() << endl << endl;
                break;
                
            case '?':                     // Unrecognized option
                // If d is provided with no argument, print all defaults
                if (optopt == 'd') printDefaults(cout);
                break;
                
            case ':':                     // Missing option argument
            case 'h':                     // Print help screen
            default :
                helpscreen();
                break;
            }
        }
    
        // Check if any parameter has been changed from the command line
        if (optind<argc && returnValue) {
            while (optind<argc) overrideParameter(std::string(argv[optind++]));
            if (!checkPars()) returnValue = false;
        }
    
    }
    
    return returnValue;
}


int Param::readParams(std::string paramfile) {
    
  /// The parameters are read in from a file, on the assumption that each
  /// line of the file has the format "parameter value" (eg. alphafdr 0.1)
  /// 
  /// The case of the parameter name does not matter, nor does the
  /// formatting of the spaces (it can be any amount of whitespace or tabs). 
  /// 
  /// \param paramfile  A std::string containing the parameter filename.
  /// 
  /// \return           1 if the parameter file does not exist. SUCCESS if
  ///                   it is able to read it.


    std::ifstream fin(paramfile.c_str());
    if(!fin.is_open()) return 1;
    std::string line;
    while(!std::getline(fin,line,'\n').eof()){
        if(line[0]!='#' || line[0]!='/'){
            std::stringstream ss;
            ss.str(line);
            setParam(ss);
        }
    }

    if(imageList!="NONE") {
        std::ifstream file(imageList.c_str());
        if(!file) return 2;
        string s;
        while(file.good()) {
            getline(file, s);
            checkHome(s);
            if (s!="") images.push_back(s);
        }
        file.close();
    }
    else {
        checkHome(imageFile);
        images.push_back(imageFile);
    }

    beamFWHM /= 3600.;

    if (!checkPars()) return 3;

    return 0;

}


void Param::setParam(stringstream &ss) {
    
    string arg;
    ss >> arg;
    arg = makelower(arg);
    
    if(arg=="fitsfile")         imageFile = readFilename(ss);
    if(arg=="fitslist")         imageList = readFilename(ss);
    if(arg=="verbose")          verbose   = readFlag(ss);
    if(arg=="outfolder")        outFolder = readFilename(ss);
    if(arg=="threads")          threads   = readIval(ss);
    if(arg=="debug")            debug     = readFlag(ss);
    if(arg=="showbar")          showbar   = readFlag(ss);
    if(arg=="plots")            plots     = readFlag(ss);
    if(arg=="beamfwhm")         beamFWHM  = readFval(ss);
    if(arg=="checkchannels")    checkChannels = readFlag(ss);
    
    if(arg=="makemask")         makeMask  = readFlag(ss);
    if(arg=="mask")             MaskType  = readFilename(ss);
    
    if(arg=="flagrobuststats")  flagRobustStats = readFlag(ss); 
    
    if(arg=="snrcut")           parSE.snrCut = readFval(ss); 
    if(arg=="threshold"){
        parSE.threshold = readFval(ss);
        parSE.UserThreshold = true;
    }

    if(arg=="search")            parSE.flagSearch = readFlag(ss);
    if(arg=="searchtype")        parSE.searchType = readSval(ss);
    if(arg=="flagadjacent")      parSE.flagAdjacent = readFlag(ss);
    if(arg=="threshspatial")     parSE.threshSpatial = readIval(ss);
    if(arg=="threshvelocity")    parSE.threshVelocity = readIval(ss);
    if(arg=="minchannels")       parSE.minChannels = readIval(ss);
    if(arg=="minvoxels")         parSE.minVoxels = readIval(ss);
    if(arg=="minpix")            parSE.minPix = readIval(ss);
    if(arg=="maxchannels")       parSE.maxChannels = readIval(ss);
    if(arg=="maxangsize")        parSE.maxAngSize = readFval(ss);
    if(arg=="rejectbeforemerge") parSE.RejectBeforeMerge = readFlag(ss);
    if(arg=="twostagemerging")   parSE.TwoStageMerging = readFlag(ss);
    if(arg=="flaggrowth")        parSE.flagGrowth = readFlag(ss);
    if(arg=="growthcut")         parSE.growthCut = readFval(ss);
    if(arg=="growththreshold"){
        parSE.growthThreshold = readFval(ss);
        parSE.flagUserGrowthT  = true;
    }

    if(arg=="globalprofile")    globprof = readFlag(ss);
    if(arg=="totalmap")         totalmap = readFlag(ss);
    if(arg=="massdensmap")      massdensmap = readFlag(ss);
    if(arg=="velocitymap")      velocitymap = readFlag(ss);
    if(arg=="dispersionmap")    dispersionmap = readFlag(ss);
    if(arg=="rmsmap")           rmsmap = readFlag(ss);
    if(arg=="blankcut")         blankCut = readFval(ss);
    if(arg=="maptype")          maptype  = makeupper(readFilename(ss));
    

    if(arg=="2dfit")     flagRing = readFlag(ss);
    if (arg=="ellprof")  flagEllProf = readFlag(ss);

    // SHARED PARAMETERS BETWEEN GALMOD, GALFIT AND GALWIND
    if(arg=="galfit")    parGF.flagGALFIT  = readFlag(ss);
    if(arg=="3dfit")     parGF.flagGALFIT  = readFlag(ss);
    if(arg=="galmod")    parGM.flagGALMOD  = readFlag(ss);
    if(arg=="galwind")   parGW.flagGALWIND = readFlag(ss);
    if(arg=="spacepar")  flagSpace  = readFlag(ss);
    if(arg=="nradii")    parGM.NRADII = parGF.NRADII               = readIval(ss);
    if(arg=="radii")     parGM.RADII  = parGF.RADII = parGW.RADII  = readFilename(ss);
    if(arg=="radsep")    parGM.RADSEP = parGF.RADSEP               = readDval(ss);
    if(arg=="vrot")      parGM.VROT   = parGF.VROT = parGW.VROT    = readFilename(ss);
    if(arg=="vrad")      parGM.VRAD   = parGF.VRAD                 = readFilename(ss);
    if(arg=="z0")        parGM.Z0     = parGF.Z0                   = readFilename(ss);
    if(arg=="vvert")     parGM.VVERT  = parGF.VVERT                = readFilename(ss);
    if(arg=="dvdz")      parGM.DVDZ   = parGF.DVDZ                 = readFilename(ss);
    if(arg=="zcyl")      parGM.ZCYL   = parGF.ZCYL                 = readFilename(ss);
    if(arg=="vdisp")     parGM.VDISP  = parGF.VDISP  = parGW.VDISP = readFilename(ss);
    if(arg=="xpos")      parGM.XPOS   = parGF.XPOS   = parGW.XPOS  = makelower(readFilename(ss));
    if(arg=="ypos")      parGM.YPOS   = parGF.YPOS   = parGW.YPOS  = makelower(readFilename(ss));
    if(arg=="vsys")      parGM.VSYS   = parGF.VSYS   = parGW.VSYS  = readFilename(ss);
    if(arg=="inc")       parGM.INC    = parGF.INC    = parGW.INC   = readFilename(ss);
    if(arg=="pa")        parGM.PHI    = parGF.PHI    = parGW.PHI   = readFilename(ss);
    if(arg=="dens")      parGM.DENS   = parGF.DENS   = parGW.DENS  = readFilename(ss);
    if(arg=="cdens")     parGM.CDENS  = parGF.CDENS  = parGW.CDENS = readFval(ss);
    if(arg=="nv")        parGM.NV     = parGF.NV     = parGW.NV    = readIval(ss);
    if(arg=="sm")        parGM.SM     = parGF.SM     = parGW.SM    = readFlag(ss);
    if(arg=="ltype")     parGM.LTYPE  = parGF.LTYPE                = readIval(ss);
    if(arg=="redshift")  parGM.REDSHIFT = parGF.REDSHIFT           = readDval(ss);
    if(arg=="restwave")  parGM.RESTWAVE = parGF.RESTWAVE           = readVec<double>(ss);
    if(arg=="restfreq")  parGM.RESTFREQ = parGF.RESTFREQ           = readVec<double>(ss);
    if(arg=="relint")    parGM.RELINT = parGF.RELINT               = readVec<double>(ss);
    if(arg=="noiserms")  parGM.NOISERMS = parGF.NOISERMS           = readDval(ss);
    
    // GALFIT ONLY PARAMETERS
    if(arg=="deltainc")  parGF.DELTAINC   = readFval(ss);
    if(arg=="deltapa")   parGF.DELTAPHI   = readFval(ss);
    if(arg=="deltavrot") parGF.DELTAVROT  = readFval(ss);
    if(arg=="ftype")     parGF.FTYPE      = readIval(ss);
    if(arg=="wfunc")     parGF.WFUNC      = readIval(ss);
    if(arg=="tol")       parGF.TOL        = readDval(ss);
    if(arg=="free")      parGF.FREE       = readFilename(ss);
    if(arg=="side")      parGF.SIDE       = readFilename(ss);
    if(arg=="bweight")   parGF.BWEIGHT    = readIval(ss);
    if(arg=="twostage")  parGF.TWOSTAGE   = readFlag(ss);
    if(arg=="flagerrors")parGF.flagERRORS = readFlag(ss);
    if(arg=="norm")      parGF.NORM       = makeupper(readFilename(ss));
    if(arg=="polyn")     parGF.POLYN      = readFilename(ss);
    if(arg=="startrad")  parGF.STARTRAD   = readIval(ss);
    if(arg=="distance")  parGF.DISTANCE   = readFval(ss);
    if(arg=="adrift")    parGF.flagADRIFT = readFlag(ss);
    if(arg=="plotmask")  parGF.PLOTMASK   = readFlag(ss);
    if(arg=="cumulative")parGF.CUMULATIVE = readFlag(ss);
    if(arg=="normalcube")parGF.NORMALCUBE = readFlag(ss);

    // GALWIND ONLY PARAMETERS
    if(arg=="htot")      parGW.HTOT       = readFval(ss);
    if(arg=="openang")   parGW.OPENANG    = readFilename(ss);
    if(arg=="ntot")      parGW.NTOT       = readIval(ss);
    if(arg=="vwind")     parGW.VWIND      = readFilename(ss);
    if(arg=="denstype")  parGW.DENSTYPE   = readIval(ss);
    if(arg=="wtype")     parGW.WTYPE      = readIval(ss);

    if(arg=="p1")        P1 = readFilename(ss);
    if(arg=="p1par")     readArray<float>(ss,P1p,3);
    if(arg=="p2")        P2 = readFilename(ss);
    if(arg=="p2par")     readArray<float>(ss,P2p,3);
    
    if(arg=="smooth")    flagSmooth = readFlag(ss);
    if(arg=="box")       readArray<int>(ss,BOX,6);
    if(arg=="bmaj")      bmaj = readDval(ss); 
    if(arg=="bmin")      bmin = readDval(ss);
    if(arg=="bpa")       bpa = readDval(ss);
    if(arg=="obmaj")     obmaj = readDval(ss);
    if(arg=="obmin")     obmin = readDval(ss);
    if(arg=="obpa")      obpa = readDval(ss);
    if(arg=="linear")    linear = readDval(ss);
    if(arg=="factor")    factor = readDval(ss);
    if(arg=="scalefactor") scalefactor = readDval(ss);
    if(arg=="fft")       flagFFT = readFlag(ss);
    if(arg=="reduce")    flagReduce = readFlag(ss);
    if(arg=="smoothoutput") smo_out = readFilename(ss);

    if(arg=="hanning")    flagHanning = readFlag(ss);
    if(arg=="hanning_size") hanning_window = readIval(ss);

    if(arg=="slitfit")    flagSlitfit = readFlag(ss);
    if(arg=="wavefile")   wavefile = readFilename(ss);
    if(arg=="ivarfile")   ivarfile = readFilename(ss);
    if(arg=="linetofit")  linetofit = readFilename(ss);
    
    if (arg=="flagpv")    flagPV = readFlag(ss);
    if (arg=="pv")        flagPV = readFlag(ss);
    if (arg=="xpos_pv")   XPOS_PV = readFval(ss);
    if (arg=="ypos_pv")   YPOS_PV = readFval(ss);
    if (arg=="pa_pv")     PA_PV = readFval(ss);
    
}


bool Param::checkPars() {
    
    bool good = true;
    
    checkHome(imageFile);
    checkHome(outFolder);
    if (outFolder!="" && outFolder[outFolder.size()-1]!='/') outFolder.append("/");

    if (!verbose) showbar=false;

#ifndef _OPENMP
    if (threads!=1) {
        std::cerr << "\n WARNING: code has been compiled without OpenMP support. Using only 1 thread.";
        threads = 1;
    }
#endif

    // Checking MASK parameter
    std::string maskstr = makeupper(MaskType);
    if (maskstr=="NONE" || maskstr=="SMOOTH" || maskstr=="SEARCH" || maskstr=="SEARCHLARGEST" ||
        maskstr=="THRESHOLD" || maskstr=="NEGATIVE" ||  maskstr=="SMOOTH&SEARCH"  || 
        maskstr=="SMOOTH&SEARCHLARGEST") {
        MaskType = maskstr;
    }
    else if (maskstr.find("FILE(")!=std::string::npos) {
        size_t first = maskstr.find_first_of("(");
        std::string sub1 = maskstr.substr(0,first);
        std::string sub2 = MaskType.substr(first);
        MaskType = sub1+sub2;
    }
    else {
        cout << "\n ERROR: Unknown type of mask: " << MaskType << std::endl;
        cout << " Setting to SMOOTH" << std::endl;
        MaskType = "SMOOTH";
    }


    // Checking parameters for source finder
    if (parSE.flagSearch) {
        if(parSE.searchType != "spatial" && parSE.searchType != "spectral"){
            cout << "You have requested a search type of \""<<parSE.searchType<<"\".\n"
                      << "Only \"spectral\" and \"spatial\" are accepted. Setting to \"spectral\".\n";
            parSE.searchType = "spectral";
        }
    
        if(parSE.flagGrowth){
            if(parSE.UserThreshold && ((parSE.threshold<parSE.growthThreshold)||
              (parSE.snrCut<parSE.growthCut))) {
                  cout << "Your \"growthThreshold\" parameter" << parSE.growthThreshold
                            <<" is larger than your \"threshold\""  << parSE.threshold << std::endl;
                  good = false;
            }
            if(!parSE.UserThreshold && (parSE.snrCut<parSE.growthCut)) {
                cout << "Your \"growthCut\" parameter " << parSE.growthCut
                          << " is larger than your \"snrCut\"" << parSE.snrCut << std::endl;
                good = false;
            }
        
            if(!good) {
                cout << "The growth function is being turned off\n.";
                parSE.flagGrowth=false;
                good = true;
            }
        }
    }
    
    // Checking parameters for 3DFIT and GALMOD
    if (parGF.flagGALFIT || flagSpace || parGM.flagGALMOD || flagSlitfit) {
        std::string str =" is not an optional parameter. Please specify it in the input file or set flagSearch=true";
        
        if (parGF.flagGALFIT && parGM.flagGALMOD) {
            cout << "3DFIT warning: 3DFIT and GALMOD can not be run at the same time. Turning off GALMOD.";
            parGM.flagGALMOD = false;
        }
        
        if (parGF.NRADII==-1 && parGF.RADII=="-1")  {
            cout << "3DFIT error: NRADII or RADII" << str << std::endl;
            good = false;
        }

        if (parGM.flagGALMOD) {
            if (parGM.XPOS=="-1") {
                cout << "3DFIT error: XPOS" << str << std::endl;
                good = false;
            }
            if (parGM.YPOS=="-1") {
                cout << "3DFIT error: YPOS" << str << std::endl;
                good = false;
            }
            if (parGM.RADSEP==-1 && parGM.RADII=="-1")  {
                cout << "3DIT error: RADSEP" << str << std::endl;
                good = false;
            }
            if (parGM.VSYS=="-1") {
                cout << "3DFIT error: VSYS" << str << std::endl;
                good = false;
            }
            if (parGM.VROT=="-1") {
                cout << "3DFIT error: VROT" << str << std::endl;
                good = false;
            }
            if (parGM.VDISP=="-1")    {
                cout << "3DFIT error: VDISP" << str << std::endl;
                good = false;
            }
            if (parGM.INC=="-1")  {
                cout << "3DFIT error: INC" << str << std::endl;
                good = false;
            }
            if (parGM.PHI=="-1")  {
                cout << "3DFIT error: PHI" << str << std::endl;
                good = false;
            }
            if (parGM.Z0=="-1")   {
                cout << "3DFIT error: Z0" << str << std::endl;
                good = false;
            }
        }

        if (parGF.NORM!="NONE" && parGF.NORM!="AZIM" && parGF.NORM!="LOCAL") {
            if (!(parGM.flagGALMOD && parGF.NORM=="BOTH")) {
                cout << " ERROR: Unknown type of normalization: " << parGF.NORM << std::endl;
                cout << "Setting to LOCAL" << std::endl;
                parGF.NORM="LOCAL";
            }
        }

        if (parGF.flagGALFIT) {
            if (parGF.FREE=="") {
                cout << "3DFIT error: FREE" << str << std::endl;
                good = false;
            }

            if (parGF.FTYPE<1 || parGF.FTYPE>4) {
                cout << "3DFIT warning: ";
                cout << "Not valid argument for FTYPE parameter. ";
                cout << "Assuming 2 (|mod-obs|).\n";
                parGF.FTYPE = 2;
            }

            if (parGF.WFUNC<0 || parGF.WFUNC>2) {
                cout << "3DFIT warning: ";
                cout << "Not valid argument for WFUNC parameter. ";
                cout << "Assuming 1 (|cos(θ)| weighting function).\n";
                parGF.WFUNC = 1;
            }
            
            string SIDE = parGF.SIDE;
            if (SIDE=="A"||SIDE=="APP"||SIDE=="APPR"||SIDE=="APPROACHING") SIDE = "A";
            else if (SIDE=="R"||SIDE=="REC"||SIDE=="RECED"||SIDE=="RECEDING") SIDE = "R";
            else if (SIDE=="B"||SIDE=="BOTH") SIDE = "B";
            else if (SIDE=="BS"||SIDE=="S"||SIDE=="SINGLE"||SIDE=="BOTH SINGLE") SIDE = "S";
            else {
                cout << "3DFIT WARNING: ";
                cout << "Not valid argument for SIDE parameter. ";
                cout << "Assuming B (fitting both sides of the galaxy).\n";
                parGF.SIDE = "B";
            }

        }

        if (parGF.LTYPE<1 || parGF.LTYPE>5) {
            cout << "3DFIT warning: ";
            cout << "Not valid argument for LTYPE parameter. ";
            cout << "Assuming 1 (gaussian layer).\n";
            parGF.LTYPE = 1;
        }
        
        //if ((parGF.RESTWAVE!=-1 && parGF.REDSHIFT==-1) || (parGF.RESTWAVE==-1 && parGF.REDSHIFT!=-1)) {
         //   cout<< "3DFIT warning: Restwave and Redshift must be set both. Exiting...\n";
         //   std::terminate();
        //}
        
        if (flagSlitfit & parGF.flagGALFIT) {
            checkHome(wavefile);
            checkHome(ivarfile);
            cout<< "3DFIT WARNING: Galfit and Slitfit cannot be run at the same time. "
                     << "Switching off Slitfit. \n";
            flagSlitfit = false;
        }
        
        if (parGF.RESTWAVE.size()==0) parGF.RESTWAVE.push_back(-1);
        if (parGF.RESTFREQ.size()==0) parGF.RESTFREQ.push_back(-1);
        if (parGF.RELINT.size()==0) parGF.RELINT.push_back(-1);
        
        if (parGF.RESTWAVE[0]!=-1 && parGF.RELINT.size()!=parGF.RESTWAVE.size()) {
            std::cerr << "3DFIT ERROR: RESTWAVE and RELINT must have same size.";
            good = false;
        }
        if (parGF.RESTFREQ[0]!=-1 && parGF.RELINT.size()!=parGF.RESTFREQ.size()) {
            std::cerr << "3DFIT ERROR: RESTFREQ and RELINT must have same size.";
            good = false;
        }
    }
    
    // Check parameters for GALWIND task
    if (parGW.flagGALWIND) {
        std::string str =" is not an optional parameter. Please specify it in the input file.";
        if (parGW.XPOS=="-1") {
            cout << "GALWIND ERROR: XPOS" << str << std::endl;
            good = false;
        }
        if (parGW.YPOS=="-1") {
            cout << "GALWIND error: YPOS" << str << std::endl;
            good = false;
        }
        if (parGW.VSYS=="-1") {
            cout << "GALWIND error: VSYS" << str << std::endl;
            good = false;
        }
        if (parGW.VDISP=="-1")    {
            cout << "GALWIND error: VDISP" << str << std::endl;
            good = false;
        }
        if (parGW.INC=="-1")  {
            cout << "GALWIND error: INC" << str << std::endl;
            good = false;
        }
        if (parGW.PHI=="-1")  {
            cout << "GALWIND error: PHI" << str << std::endl;
            good = false;
        }
        if (parGW.DENS=="-1")   {
            cout << "GALWIND error: DENS" << str << std::endl;
            good = false;
        }
        if (parGW.OPENANG=="-1")   {
            cout << "GALWIND error: OPENANG" << str << std::endl;
            good = false;
        }
        if (parGW.HTOT==-1)   {
            cout << "GALWIND error: HTOT" << str << std::endl;
            good = false;
        }
        if (parGW.VWIND=="-1")   {
            cout << "GALWIND error: VWIND" << str << std::endl;
            good = false;
        }        
    }
    
    // Checking parameters for SMOOTH
    if (flagSmooth) {
        if (bmaj==-1 && bmin==-1 && linear==-1 && factor==-1) {
            cout << "SMOOTH error: "
                      << "you need to specify either the new beam (BMAJ, BMIN, BPA) "
                      << "or the parameters LINEAR and DISTANCE.\n";
            good = false;
        }
        if (linear!=-1 && parGF.DISTANCE==-1) {
            cout << "SMOOTH error: "
                      << "with LINEAR parameter you must specify also the DISTANCE. ";
            good = false;
        }
        if (bmaj!=-1 && bmin==-1) {
            bmin = bmaj;
        }
    }

    if (flagHanning) {
        if (hanning_window%2==0) {
            cout << "HANNING error: Hanning window must be an odd number\n";
            good = false;
        }
    }

    if (getMaps()) {
        if (maptype!="GAUSSIAN" && maptype!="MOMENT")
            cout << "MAP warning: MAPTYPE is either MOMENT or GAUSSIAN. Reverting to MOMENT.\n";
    }
    
    return good;
}


void Param::overrideParameter(std::string parstr){
    // Override a parameter in the parameter file from the commandline:
    // Example: BBarolo -p param.par INC=60
    int found = parstr.find("=");
    if (found!=std::string::npos) {
        std::stringstream ss(parstr.replace(found,1," "));
        setParam(ss);
    }
}


void Param::printDefaults (std::ostream& theStream, string wtask) {
    
    Param par;
    
    string whichtask = wtask;
    // Check that whichtask is ok
    bool isOK = false;
    for (auto i=0;i<tasks.size();i++) {
        isOK = whichtask==tasks[i];
        if (isOK) break;
    }
    if (!isOK) whichtask = "ALL";
    
    theStream.setf(std::ios::left);
    if (whichtask=="ALL") theStream  <<"\n---------------------- Default parameters ----------------------\n" << endl;
    else theStream  <<"\n----------------- Default parameters for " << whichtask << " -----------------\n" << endl;
    theStream  << std::setfill('.');
   
    printParams(theStream,par,true,whichtask);
    
}


void Param::createTemplate() {
    
    using namespace std;
    
    int m = 12;
    std::ofstream parf;
    parf.open("param.par");
    
    parf << "// This is a template input file for the 3DFIT task.\n";
    parf << "// Lines beginning with the double-slash or hash and blank \n"
         << "// lines are not read by the program.\n\n";
    
    parf << "// Name of the fitsfile to be modeled.\n";
    parf << setw(m) << left << "FITSFILE";
    parf << setw(m) << left << "/yourpath/yourfile.fits\n" << endl;
    
    parf << "// Using the 3DFIT utility? Must be true!!\n";
    parf << setw(m) << left << "3DFIT";
    parf << setw(m) << left << "true\n" << endl;
    
    parf << "// Number of radii to be modeled.\n";
    parf << setw(m) << left << "NRADII";
    parf << setw(m) << left << "10\n" << endl;
    
    parf << "// Separation between radii in arcsec.\n";
    parf << setw(m) << left << "RADSEP";
    parf << setw(m) << left << "30\n" << endl;
    
    parf << "// Systemic velocity of the galaxy (in km/s).\n";
    parf << setw(m) << left << "VSYS" << setw(m) << left << "132\n" << endl;
    
    parf << "// X-Y coordinates of the galaxy center (in pixel).\n";
    parf << setw(m) << left << "XPOS" << setw(m) << left << "256" << endl;
    parf << setw(m) << left << "YPOS" << setw(m) << left << "256\n" << endl;
    
    parf << "// Initial global values for parameters:\n";
    parf << "// Rotation and dispersion velocities (in km/s),\n";
    parf << "// inclination and position angles [measured\n"; 
    parf << "// anti-clockwise from north] (in degrees),\n";
    parf << "// height scale of the disk (in arcsec).\n";
    parf << setw(m) << left << "VROT"  << setw(m) << left << "100"<< endl;
    parf << setw(m) << left << "VRAD"  << setw(m) << left << "0"<< endl;
    parf << setw(m) << left << "VDISP" << setw(m) << left << "10" << endl;
    parf << setw(m) << left << "INC"   << setw(m) << left << "60" << endl;
    parf << setw(m) << left << "PA"    << setw(m) << left << "100"<< endl;
    parf << setw(m) << left << "Z0"    << setw(m) << left << "20\n" << endl;
    
    parf << "// Free parameters for the minimization.\n";
    parf << setw(m) << left << "FREE" << setw(m) << left << "VROT VDISP INC PA\n" << endl;
    
    parf << "// OPTIONAL: Function to be minimized (default is " << parGF.FTYPE << "):\n";
    parf << "// = 1: chi-squared.\n";
    parf << "// = 2: |mod-obs|.\n";
    parf << "// = 3: |mod-obs|/|mod+obs|.\n";
    parf << "// = 4: (mod-obs)^2.\n";
    parf << setw(m) << left << "FTYPE" << setw(m) << left << "1\n"<< endl;

    parf << "// OPTIONAL: Weighting function (default is " << parGF.WFUNC << "):\n";
    parf << "// = 0: uniform weight.\n";
    parf << "// = 1: |cos(θ)|.\n";
    parf << "// = 2: cos(θ)^2.\n";
    parf << "// θ is the azimuthal angle.\n";
    parf << setw(m) << left << "WFUNC" << setw(m) << left << "1\n"<< endl;

    parf << "// OPTIONAL: Layer type along z (default is " << parGF.LTYPE << "):\n";
    parf << "// = 1: gaussian layer.\n";
    parf << "// = 2: sech2 layer.\n";
    parf << "// = 3: exponential layer.\n";
    parf << "// = 4: Lorentzian layer.\n";
    parf << "// = 5: box layer.;\n";
    parf << setw(m) << left << "LTYPE" << setw(m) << left << "1\n"<< endl;

    parf << "// OPTIONAL: Number of subcloud in a velocity profile.\n";
    parf << "// (default is = total number of channels).\n";
    parf << setw(m) << left << "NV" << setw(m) << left << "60\n"<< endl;
    
    parf << "// OPTIONAL: Surface density of clouds in the plane of ring (1e20).\n";
    parf << "// (default is = " << parGF.CDENS << "):\n";
    parf << setw(m) << left << "CDENS" << setw(m) << left << "10\n"<< endl;

    parf << "// OPTIONAL: Tolerance for the minimization (default is " << parGF.TOL << "):\n";
    parf << setw(m) << left << "TOL" << setw(m) << left << "1E-03\n"<< endl;

    parf << "// OPTIONAL: Using a mask for the minimization (default is " << MaskType << "):\n";
    parf << setw(m) << left << "MASK" << setw(m) << left << "SMOOTH\n"<< endl;
    
    parf << "// OPTIONAL: Normalization type (default is " << parGF.NORM << "):\n";
    parf << setw(m) << left << "NORM" << setw(m) << left << "AZIM\n"<< endl;
    
    parf << "// OPTIONAL: Side of the galaxy to be fitted (default is " << parGF.SIDE << "):\n";
    parf << "// = A: Approaching.\n";
    parf << "// = R: Receding.\n";
    parf << "// = B: Both.\n";
//  parf << "// = S: Both but separated.\n";
    parf << setw(m) << left << "SIDE" << setw(m) << left << "B\n"<< endl;
    
    parf << "// OPTIONAL: Using a two stages minimization (default is " << stringize(parGF.TWOSTAGE) << "):\n";
    parf << setw(m) << left << "TWOSTAGE" << setw(m) << left << "false\n"<< endl;
    
    parf << "// OPTIONAL: Degree of polynomial fitting angles (default is " << parGF.POLYN << "):\n";
    parf << setw(m) << left << "POLYN" << setw(m) << left << "3\n"<< endl;
    
    parf << "// OPTIONAL: Enabling error estimation (default is " << stringize(parGF.flagERRORS) << "):\n";
    parf << setw(m) << left << "flagErrors" << setw(m) << left << "false\n"<< endl;
    
    parf.close();
    
    cout << "\n A template parameter input file (\"param.par\") for BBarolo " 
              << "has been generated.\n\n" ;
}


void recordParameters(std::ostream& theStream, std::string paramName, std::string paramDesc, std::string paramValue) {
    
    const int width = 60;
    int widthText = width - paramName.size();

    theStream << std::setw(widthText) << paramDesc
              << setiosflags(std::ios::right) << paramName
              << "  =  " << resetiosflags(std::ios::right) << paramValue 
              << std::endl;
}


std::string fileOption(bool flag, std::string file) {
    
    std::ostringstream ss;
    ss << stringize(flag);
    if(flag) ss << " --> " << file;
    return ss.str();
    
}

  
std::ostream& operator<< (std::ostream& theStream, Param& par) {
    
    /// Print out the parameter set in a formatted, easy to read style.
   
    theStream.setf(std::ios::left);
    theStream  <<"\n-------------------------- Parameters -------------------------\n"<<std::endl;
    theStream  << std::setfill('.');
    printParams(theStream,par,false);
    return theStream;    
}


void printParams (std::ostream& Str, Param &p, bool defaults, string whichtask) {
    
    bool isAll = whichtask=="ALL" && defaults;
    
    // BASIC PARAMETERS FOR ALL TASKS
    if (defaults) {
        recordParam(Str, "[FITSFILE]", "FITS file to be analysed", p.getImageFile());
        recordParam(Str, "[FITSLIST]", "List of FITS files to be analysed", p.getImageList());
    }
    else {
            if (p.getImageList()=="NONE") 
                recordParam(Str, "[FITSFILE]", "FITS file to be analysed", p.getImageFile());
            else recordParam(Str, "[FITSLIST]", "List of FITS files to be analysed", p.getImageList());
    }
    recordParam(Str, "[THREADS]", "Number of threads", p.getThreads());
    recordParam(Str, "[PLOTS]", "Producing output plots?", stringize(p.getFlagPlots()));
    recordParam(Str, "[SHOWBAR]", "Showing progress bars?", stringize(p.getShowbar()));
    recordParam(Str, "[VERBOSE]", "Printing output messages?", stringize(p.isVerbose()));
    if (p.getOutfolder()!="" || defaults)
        recordParam(Str, "[OUTFOLDER]", "Directory where outputs are written", p.getOutfolder());
    recordParam(Str, "[flagRobustStats]", "Using Robust statistics?", stringize(p.getFlagRobustStats()));
    
    Str  <<"-----------------"<<std::endl;
    
    // PARAMETERS FOR MAKEMASK TASK
    bool toPrint = isAll || p.getMakeMask() || (defaults && whichtask=="MAKEMASK");
    if (toPrint) {
        recordParam(Str, "[MAKEMASK]", "Writing the mask to a fitsfile", stringize(p.getMakeMask()));
        recordParam(Str, "[MASK]", "   Type of mask", p.getMASK());
    }
    
    // PARAMETERS FOR SEARCH TASK
    toPrint = isAll || p.getParSE().flagSearch || (defaults && whichtask=="SEARCH");
    if (toPrint) {
        recordParam(Str, "[SEARCH]", "Searching for sources in cube?", stringize(p.getParSE().flagSearch));
        recordParam(Str, "[searchType]", "   Type of searching performed", p.getParSE().searchType);
        recordParam(Str, "[minPix]", "   Minimum # Pixels in a detection", p.getParSE().minPix);
        recordParam(Str, "[minChannels]", "   Minimum # Channels in a detection", p.getParSE().minChannels);
        recordParam(Str, "[minVoxels]", "   Minimum # Voxels in a detection", p.getParSE().minVoxels);
        recordParam(Str, "[maxChannels]", "   Maximum # Channels in a detection", p.getParSE().maxChannels);
        recordParam(Str, "[maxAngsize]", "   Max angular size of a detection in arcmin", p.getParSE().maxAngSize);
        recordParam(Str, "[flagAdjacent]", "   Using Adjacent-pixel criterion?", stringize(p.getParSE().flagAdjacent));
        if(!p.getParSE().flagAdjacent || defaults)
            recordParam(Str, "[threshSpatial]", "   Max. spatial separation for merging", p.getParSE().threshSpatial);
        recordParam(Str, "[threshVelocity]", "   Max. velocity separation for merging", p.getParSE().threshVelocity);
        recordParam(Str, "[RejectBeforeMerge]", "   Reject objects before merging?", stringize(p.getParSE().RejectBeforeMerge));
        recordParam(Str, "[TwoStageMerging]", "   Merge objects in two stages?", stringize(p.getParSE().TwoStageMerging));
        if (defaults) {
            recordParam(Str, "[threshold]", "   Detection Threshold", p.getParSE().threshold);
            recordParam(Str, "[snrCut]", "   SNR Threshold (in sigma)", p.getParSE().snrCut);
        }
        else {
            if(p.getParSE().UserThreshold)
                recordParam(Str, "[threshold]", "   Detection Threshold", p.getParSE().threshold);
            else 
                recordParam(Str, "[snrCut]", "   SNR Threshold (in sigma)", p.getParSE().snrCut);
        }
        
        recordParam(Str, "[flagGrowth]", "   Growing objects after detection?", stringize(p.getParSE().flagGrowth));
        if(p.getParSE().flagGrowth || defaults) {
            if (defaults) {
                recordParam(Str, "[growthThreshold]", "   Threshold for growth", p.getParSE().growthThreshold);
                recordParam(Str, "[growthCut]", "   SNR Threshold for growth", p.getParSE().growthCut);
            }
            else {
                if(p.getParSE().flagUserGrowthT)
                    recordParam(Str, "[growthThreshold]", "   Threshold for growth", p.getParSE().growthThreshold);
                else
                    recordParam(Str, "[growthCut]", "   SNR Threshold for growth", p.getParSE().growthCut);
            }            
        }
    }
    
    // PARAMETERS FOR SMOOTH
    toPrint = isAll || p.getflagSmooth() || (defaults && whichtask=="SMOOTH");
    if (toPrint) {
        recordParam(Str, "[SMOOTH]", "Smoothing the datacube?", stringize(p.getflagSmooth()));
        std::string box;
        for (int i=0;i<6;i++) if (p.getBOX(i)!=-1) box += to_string<int>(p.getBOX(i))+" ";
        if (box=="") box = "NONE";
        recordParam(Str, "[BOX]", "   Sub-region to be used?", box);
        if ((p.getOBmaj()!=-1 && p.getOBmin()!=-1) || defaults) {
            recordParam(Str, "[OBMAJ]", "   Old major beam (arcsec)", p.getOBmaj());
            recordParam(Str, "[OBMIN]", "   Old minor beam (arcsec)", p.getOBmin());
            recordParam(Str, "[OBPA]", "   Old beam position angle (degree)", p.getOBpa()); 
        }
        if (defaults) {
            recordParam(Str, "[BMAJ]", "   New major beam (arcsec)", p.getBmaj());
            recordParam(Str, "[BMIN]", "   New minor beam (arcsec)", p.getBmin());
            recordParam(Str, "[BPA]", "   New beam position angle (degree)", p.getBpa()); 
            recordParam(Str, "[FACTOR]", "   New beam factor (times old beam)", p.getFactor());
            recordParam(Str, "[LINEAR]", "   New linear resolution (kpc)", p.getLinear());
        }
        else {
            if (p.getLinear()!=-1 || defaults) 
                recordParam(Str, "[LINEAR]", "   New linear resolution (kpc)", p.getLinear());
            else if (p.getBmaj()!=-1 && p.getBmin()!=-1){
                recordParam(Str, "[BMAJ]", "   New major beam (arcsec)", p.getBmaj());
                recordParam(Str, "[BMIN]", "   New minor beam (arcsec)", p.getBmin());
                recordParam(Str, "[BPA]", "   New beam position angle (degree)", p.getBpa());
            }
            else 
                recordParam(Str, "[FACTOR]", "   New beam factor (times old beam)", p.getFactor());
        }
        recordParam(Str, "[FFT]", "   Using FFT for convolution?", stringize(p.getflagFFT()));
        recordParam(Str, "[REDUCE]", "   Reducing datacube?", stringize(p.getflagReduce()));
    }
    
    // PARAMETERS FOR HANNING
    toPrint = isAll || p.getflagHanning() || (defaults && whichtask=="HANNING");
    if (toPrint) {
        recordParam(Str, "[HANNING]", "Hanning smoothing the datacube?", stringize(p.getflagHanning()));
        recordParam(Str, "[HANNING_SIZE]", "   Size of hanning window (channels)", p.getHanningWindow());
    }
    
    // GALMOD & 3DFIT parameters
    bool isGalfit   = isAll || p.getflagGalFit() || (defaults && whichtask=="3DFIT");
    bool isGalmod   = isAll || p.getflagGalMod() || (defaults && whichtask=="GALMOD");
    bool isSlitfit  = isAll || p.getFlagSlitfit() || (defaults && whichtask=="SLITFIT");
    bool isSpacepar = isAll || p.getflagSpace() || (defaults && whichtask=="SPACEPAR");
    toPrint = isGalfit || isGalmod || isSlitfit || isSpacepar;
    if (toPrint) {
        if (isGalmod)   recordParam(Str, "[GALMOD]", "Writing a 3D model?", stringize(p.getflagGalMod()));    
        if (isSlitfit)  recordParam(Str, "[SLITFIT]", "Fitting a model to slitdata?", stringize(p.getFlagSlitfit()));
        if (isGalfit)   recordParam(Str, "[3DFIT]", "Fitting a 3D model to the datacube?", stringize(p.getflagGalFit()));
        
        recordParam(Str, "[REDSHIFT]", "   Redshift of the galaxy?", p.getParGF().REDSHIFT);
        if (p.getParGF().RESTWAVE[0]!=-1 || defaults) {
            string rstr = "";
            for (unsigned i=0; i<p.getParGF().RESTWAVE.size(); i++) rstr = rstr + to_string<double>(p.getParGF().RESTWAVE[i],2) + " ";
            recordParam(Str, "[RESTWAVE]", "   Transition wavelength at rest?", rstr);
        }
        if (p.getParGF().RESTFREQ[0]!=-1 || defaults) {
            string rstr = "";
            for (unsigned i=0; i<p.getParGF().RESTFREQ.size(); i++) rstr = rstr + to_string<double>(p.getParGF().RESTFREQ[i],2) + " ";
            recordParam(Str, "[RESTFREQ]", "   Transition frequency at rest?", rstr);
        }
        if (p.getParGF().RELINT.size()>1 || defaults) {
            string rstr = "";
            for (unsigned i=0; i<p.getParGF().RELINT.size(); i++) rstr = rstr + to_string<double>(p.getParGF().RELINT[i],2) + " ";
            recordParam(Str, "[RELINT]", "   Relative intensities of lines?", rstr);
        }
        
        if (defaults) {
            recordParam(Str, "[RADII]", "   Radii for rings", p.getParGF().RADII);
            recordParam(Str, "[NRADII]", "   Number of radii", p.getParGF().NRADII);
            recordParam(Str, "[RADSEP]", "   Separation between radii (arcsec)", p.getParGF().RADSEP);
        }
        else {
            if (p.getParGF().RADII!="-1") 
                recordParam(Str, "[RADII]", "   Radii for rings", p.getParGF().RADII);
            else {
                recordParam(Str, "[NRADII]", "   Number of radii", p.getParGF().NRADII);
                recordParam(Str, "[RADSEP]", "   Separation between radii (arcsec)", p.getParGF().RADSEP);
            }
        }
        recordParam(Str, "[XPOS]", "   X center of the galaxy (pixel)", p.getParGF().XPOS);
        recordParam(Str, "[YPOS]", "   Y center of the galaxy (pixel)", p.getParGF().YPOS);
        recordParam(Str, "[VSYS]", "   Systemic velocity of the galaxy (km/s)", p.getParGF().VSYS);
        recordParam(Str, "[VROT]", "   Initial global rotation velocity (km/s)", p.getParGF().VROT);
        if (isGalfit) recordParam(Str, "[DELTAVROT]", "   Max rot. velocity variation from VROT (km/s)", p.getParGF().DELTAVROT);
        recordParam(Str, "[VRAD]", "   Initial global radial velocity (km/s)", p.getParGF().VRAD);
        recordParam(Str, "[VDISP]", "   Initial global velocity dispersion (km/s)", p.getParGF().VDISP);
        recordParam(Str, "[INC]", "   Initial global inclination (degrees)", p.getParGF().INC);
        if (isGalfit) recordParam(Str, "[DELTAINC]", "   Max inclination variation from INC (degrees)", p.getParGF().DELTAINC);
        recordParam(Str, "[PA]", "   Initial global position angle (degrees)", p.getParGF().PHI);
        if (isGalfit) recordParam(Str, "[DELTAPA]", "   Max position angle variation from PA (degrees)", p.getParGF().DELTAPHI);
        recordParam(Str, "[Z0]", "   Scale height of the disk (arcsec)", p.getParGF().Z0);
        recordParam(Str, "[DENS]", "   Global column density of gas (atoms/cm2)", p.getParGF().DENS);
        if (isGalfit) recordParam(Str, "[FREE]", "   Parameters to be minimized", p.getParGF().FREE);
        recordParam(Str, "[MASK]", "   Type of mask", p.getMASK());
        recordParam(Str, "[NORM]", "   Type of normalization", (p.getParGF().NORM));
        if (isGalfit) recordParam(Str, "[SIDE]", "   Side of the galaxy to be used", (p.getParGF().SIDE)); 
        if (p.getParGF().DISTANCE!=-1)
            recordParam(Str, "[DISTANCE]", "   Distance of the galaxy (Mpc)?", p.getParGF().DISTANCE);
        
        std::string typ = "";
        int t = p.getParGF().LTYPE;
        if (t==1)      typ = "gaussian";
        else if (t==2) typ = "sech2";
        else if (t==3) typ = "exponential";
        else if (t==4) typ = "Lorentzian";
        else if (t==5) typ = "box";
        recordParam(Str, "[LTYPE]", "   Layer type along z direction", typ);
        
        typ = "";
        t = p.getParGF().FTYPE;
        if (t==1)      typ = "chi-squared";
        else if (t==2) typ = "|m-o|";
        else if (t==3) typ = "|m-o|/|m+o|";
        else if (t==4) typ = "(m-o)^2";
        recordParam(Str, "[FTYPE]", "   Residuals to minimize", typ);
        
        typ = "";
        t = p.getParGF().WFUNC;
        if (t==0)      typ = "uniform";
        else if (t==1) typ = "|cos(θ)|";
        else if (t==2) typ = "cos(θ)^2";
        if (isGalfit) {
            recordParam(Str, "[WFUNC]", "   Weighting function", typ);
            recordParam(Str, "[BWEIGHT]", "   Weight for blank pixels", p.getParGF().BWEIGHT); 
            recordParam(Str, "[TOL]", "   Minimization tolerance", p.getParGF().TOL); 
            recordParam(Str, "[SIDE]", "   What side of the galaxy to be used", (p.getParGF().SIDE)); 
            recordParam(Str, "[TWOSTAGE]", "   Two stages minimization?", stringize(p.getParGF().TWOSTAGE));
            if (p.getParGF().TWOSTAGE || defaults)
                recordParam(Str, "[POLYN]", "     Degree of polynomial fitting angles?", p.getParGF().POLYN);
            recordParam(Str, "[FLAGERRORS]", "   Estimating errors?", stringize(p.getParGF().flagERRORS));
        
            recordParam(Str, "[ADRIFT]", "   Computing asymmetric drift correction?", stringize(p.getParGF().flagADRIFT));
            recordParam(Str, "[PLOTMASK]", "   Overlaying mask to output plots?", stringize(p.getParGF().PLOTMASK));
            recordParam(Str, "[CUMULATIVE]", "   Using cumulative rings during the fit?", stringize(p.getParGF().CUMULATIVE));
            recordParam(Str, "[NORMALCUBE]", "   Normalizing cube to help convergence?", stringize(p.getParGF().NORMALCUBE));
        }

        recordParam(Str, "[NOISERMS]", "   RMS noise to add to the model", stringize(p.getParGF().NOISERMS));
    
        // PARAMETERS FOR SPACEPAR

        if (isSpacepar) {
            recordParam(Str, "[SPACEPAR]", "Full parameter space for a pair of parameters", stringize(p.getflagSpace()));    
            recordParam(Str, "[P1]", "   First parameter to explore", p.getP1());    
            recordParam(Str, "[P2]", "   Second parameter to explore", p.getP2());   
            std::string pp = "";
            for (int i=0;i<3;i++) pp += to_string<float>(p.getP1p(i),2)+" ";
            recordParam(Str, "[P1PAR]", "   Range and step for P1", pp);
            pp = "";
            for (int i=0;i<3;i++) pp += to_string<float>(p.getP2p(i),2)+" ";
            recordParam(Str, "[P2PAR]", "   Range and step for P2", pp);
        }
    }
        
    // PARAMETERS FOR GALWIND
    toPrint = isAll || p.getParGW().flagGALWIND || (defaults && whichtask=="GALWIND");
    if (toPrint) {
        recordParam(Str, "[GALWIND]", "Generating a 3D datacube with a wind model?", stringize(p.getParGW().flagGALWIND));
        recordParam(Str, "[VWIND]",   "   Radial velocity of the wind (km/s)", p.getParGW().VWIND);
        recordParam(Str, "[VROT]",   "   Azimuthal velocity of the wind (km/s)", p.getParGW().VROT);
        recordParam(Str, "[OPENANG]", "   Wind opening angle (degrees)", p.getParGW().OPENANG);
        recordParam(Str, "[HTOT]",    "   Wind maximum distance (arcsec)", p.getParGW().HTOT);
        recordParam(Str, "[XPOS]",    "   X center of the galaxy (pixel)", p.getParGW().XPOS);
        recordParam(Str, "[YPOS]",    "   Y center of the galaxy (pixel)", p.getParGW().YPOS);
        recordParam(Str, "[VSYS]",    "   Systemic velocity of the galaxy (km/s)", p.getParGW().VSYS);
        recordParam(Str, "[VDISP]",   "   Global velocity dispersion (km/s)", p.getParGW().VDISP);
        recordParam(Str, "[INC]",     "   Global inclination (degrees)", p.getParGW().INC);
        recordParam(Str, "[PA]",      "   Global position angle (degrees)", p.getParGW().PHI);
        recordParam(Str, "[DENS]",    "   Global column density of gas (atoms/cm2)", p.getParGW().DENS);
        recordParam(Str, "[NTOT]",    "   Number of layers/cylinder for each cone", p.getParGW().NTOT);
        recordParam(Str, "[DENSTYPE]","   How to distribute density in layers", p.getParGW().DENSTYPE);
        recordParam(Str, "[WTYPE]","   Cylindrical or spherical model", p.getParGW().WTYPE);
        recordParam(Str, "[SM]","   Whether to smooth the model", stringize(p.getParGF().SM));
    }
    
    // PARAMETERS FOR 2DFIT
    toPrint = isAll || p.getFlagRing() || (defaults && whichtask=="2DFIT");
    if (toPrint) {
        recordParam(Str, "[2DFIT]", "Fitting velocity field with a ring model?", stringize(p.getFlagRing()));
        if (defaults) {
            recordParam(Str, "[RADII]", "   Radii for rings", p.getParGF().RADII);
            recordParam(Str, "[NRADII]", "   Number of radii", p.getParGF().NRADII);
            recordParam(Str, "[RADSEP]", "   Separation between radii (arcsec)", p.getParGF().RADSEP);
        }
        else {
            if (p.getParGF().RADII!="-1") 
                recordParam(Str, "[RADII]", "   Radii for rings", p.getParGF().RADII);

            else {
                recordParam(Str, "[NRADII]", "   Number of radii", p.getParGF().NRADII);
                recordParam(Str, "[RADSEP]", "   Separation between radii (arcsec)", p.getParGF().RADSEP);
            }
        }
        recordParam(Str, "[XPOS]", "   X center of the galaxy (pixel)", p.getParGF().XPOS);
        recordParam(Str, "[YPOS]", "   Y center of the galaxy (pixel)", p.getParGF().YPOS);
        recordParam(Str, "[VSYS]", "   Systemic velocity of the galaxy (km/s)", p.getParGF().VSYS);
        recordParam(Str, "[VROT]", "   Rotation velocity (km/s)", p.getParGF().VROT);
        recordParam(Str, "[VRAD]", "   Radial velocity (km/s)", p.getParGF().VRAD);
        recordParam(Str, "[INC]",  "   Inclination angle (degrees)", p.getParGF().INC);
        recordParam(Str, "[PA]",   "   Position angle (degrees)", p.getParGF().PHI);
        recordParam(Str, "[FREE]", "   Parameters to be fit", p.getParGF().FREE);
        recordParam(Str, "[MASK]", "   Type of mask for velocity map", p.getMASK());
        recordParam(Str, "[SIDE]", "   Side of the galaxy to be used", p.getParGF().SIDE); 
    }
    
    // PARAMETERS FOR ELLPROF
    toPrint = isAll || p.getFlagEllProf() || (defaults && whichtask=="ELLPROF");
    if (toPrint) {
        recordParam(Str, "[ELLPROF]", "Deriving radial intensity profile?", stringize(p.getFlagEllProf()));
        if (defaults) {
            recordParam(Str, "[RADII]", "   Radii for rings", p.getParGF().RADII);
            recordParam(Str, "[NRADII]", "   Number of radii", p.getParGF().NRADII);
            recordParam(Str, "[RADSEP]", "   Separation between radii (arcsec)", p.getParGF().RADSEP);
        }
        else {
            if (p.getParGF().RADII!="-1") 
                recordParam(Str, "[RADII]", "   Radii for rings", p.getParGF().RADII);

            else {
                recordParam(Str, "[NRADII]", "   Number of radii", p.getParGF().NRADII);
                recordParam(Str, "[RADSEP]", "   Separation between radii (arcsec)", p.getParGF().RADSEP);
            }
        }
        recordParam(Str, "[XPOS]", "   X center of the galaxy (pixel)", p.getParGF().XPOS);
        recordParam(Str, "[YPOS]", "   Y center of the galaxy (pixel)", p.getParGF().YPOS);
        recordParam(Str, "[INC]",  "   Inclination angle (degrees)", p.getParGF().INC);
        recordParam(Str, "[PA]",   "   Position angle (degrees)", p.getParGF().PHI);
        recordParam(Str, "[MASK]", "   Type of mask for intensity map", p.getMASK());
        recordParam(Str, "[SIDE]", "   Side of the galaxy to be used", p.getParGF().SIDE); 
    }
    
    // Print CHECKCHANNELS
    if (p.getCheckCh() || (defaults && isAll))
        recordParam(Str, "[checkChannels]", "Checking for bad channels in the cube", stringize(p.getCheckCh()));
    
    // PARAMETERS FOR MOMENT MAPS
    if (p.getGlobProf() || (defaults && isAll))
        recordParam(Str, "[globalProfile]", "Saving the global profile?", stringize(p.getGlobProf()));
    if (p.getTotalMap() || (defaults && isAll))
        recordParam(Str, "[totalMap]",      "Saving integrated intensity map to FITS file?", stringize(p.getTotalMap()));
    if (p.getMassDensMap() || (defaults && isAll))
        recordParam(Str, "[massdensMap]",   "Saving HI mass density map to FITS file?", stringize(p.getMassDensMap()));
    if (p.getVelMap() || (defaults && isAll))
        recordParam(Str, "[velocityMap]",   "Saving velocity map to FITS file?", stringize(p.getVelMap()));
    if (p.getDispMap() || (defaults && isAll))
        recordParam(Str, "[dispersionMap]", "Saving velocity dispersion map to FITS file?", stringize(p.getDispMap()));
    if (p.getRMSMap() || (defaults && isAll))
        recordParam(Str, "[rmsMap]",        "Saving RMS map to FITS file?", stringize(p.getRMSMap()));
    if (p.getMaps() || (defaults && isAll)) 
        recordParam(Str, "[MASK]",          "   Mask used for maps and profile?", p.getMASK());
    if (p.getTotalMap() || p.getVelMap() || p.getDispMap() || (defaults && isAll))
        recordParam(Str, "[MAPTYPE]",       "   How to extract the map (gaussian or moment)?", p.getMapType());
    
    // PARAMETERS FOR PV
    toPrint = isAll || p.getFlagPV() || (defaults && whichtask=="PV");
    if (toPrint) {
        recordParam(Str, "[PV]", "Extracting a position-velocity slice from a datacube?", stringize(p.getFlagPV()));
        recordParam(Str, "[XPOS_PV]", "   X center of the slice", p.getXPOS_PV());
        recordParam(Str, "[YPOS_PV]", "   Y center of the slice", p.getYPOS_PV());
        recordParam(Str, "[YPOS_PV]", "   Position angle of the slice", p.getPA_PV());
    }
    
    Str  << std::endl <<"-----------------------------";
    Str  << "------------------------------\n\n";
    Str  << std::setfill(' ');
    Str.unsetf(std::ios::left);
}


void listTasks (std::ostream& Str) {
    Str.setf(std::ios::left);
    
    Str << "     _ _                                                            \n"
        << ".-. | B |                                               .--.        \n"
        << "|A|_| A |                                         .-===-|==|---.    \n"
        << "|S| | R.|<\\            Available tasks            | C++ |  |___|-. \n"
        << "|T| | O | \\\\                                      |     |GR|---|=|\n"
        << "|R| | L |  \\\\                                     |     |  |___| |\n"
        << "|O| | O |   \\>                                    |-----|==|---|=| \n"
        << setfill('"') << setw(68) << right << "\n\n";

    for (auto i=0;i<tasks.size();i++) {
        Str << tasks[i] << ": " << endl;
        Str << "   "  << taskdescr[i] << endl << endl;
    }
    Str << setfill('"') << setw(68) << right << "\n\n";

}


void helpscreen(std::ostream& Str) {
    
    using std::endl;
    int m=20;

    Str << endl
        << "             ____                                _____        __           \n"
        << "   /\\      .'    /\\                             /.---.\\      (==)       \n"
        << "  |K -----;     |  |                           | ````` |     |~~|          \n"
        << "   \\/      '.____\\/                             \\     /      |  |       \n"
        << "                   °.                            `-.-'       |()|          \n" 
        << "        __          °.                             |        /`  `\\        \n" 
        << "   __  |_/         .°.  _________________        __|__     /      \\       \n"
        << "   \\_|\\\\ _          .°.                         `-----`   ;  ____  ;    \n" 
        << "      _\\(_)_           BBarolo quick guide                ||`    `||      \n"
        << "     (_)_)(_)_      _________________________     ___     ||BAROLO||       \n"
        << "    (_)(_)_)(_)                                  ||  \\\\   ||      ||     \n"
        << "     (_)(_))_)                                   ||__// _ || 2015 ||       \n"
        << "      (_(_(_)                                    ||  \\\\   | \\____/ |    \n"
        << "       (_)_)                                     ||__//   |        |       \n"
        << "        (_)                                               \\_.-\"\"-._/    \n "
        << "                                                          `\"\"\"\"\"\"`   \n\n\n "
        << "  Usage:           BBarolo <option> [arg] \n\n"
        << " Options: \n\n"
        << setw(m) << left << "   -f, --fitsfile" 
        << "BBarolo runs the 3DFIT task in automatic mode \n"
        << setw(m) << left << " " 
        << "using default parameters. [arg] is required   \n"
        << setw(m) << left << " " 
        << "and represents a FITS file to be analysed.    \n"
        << setw(m) << left << " " 
        << "BBarolo looks for sources, estimates initial  \n"
        << setw(m) << left << " " 
        << "parameters and fits a 3D tilted-ring model. \n\n" 
        << setw(m) << left << " "  
        << "Example:      BBarolo -f myfitsfile.fits      \n"
        << endl << endl
        << setw(m) << left << "   -p, --paramfile" 
        << "BBarolo runs with a parameter file. [arg] is \n"
        << setw(m) << left << " " 
        << "required and is a text file with a list of   \n"
        << setw(m) << left << " " 
        << "parameters. We recommend to run BBarolo this \n" 
        << setw(m) << left << " " 
        << "this way to achieve best results. \n\n"
        << setw(m) << left << " " 
        << "Example:      BBarolo -p param.par       \n"
        << endl << endl
        << setw(m) << left << "   -d, --defaults" 
        << "Print on the screen a list all available par-\n"
        << setw(m) << left << " " 
        << "ameters and their default values. A task name\n"
        << setw(m) << left << " " 
        << "[arg] can optionally be given to print only  \n"
        << setw(m) << left << " " 
        << "parameters relevant for that task. \n\n"
        << setw(m) << left << " " 
        << "Example:      BBarolo -d 3DFIT       \n"
        << endl << endl
        << setw(m) << left << "   -l, --list" 
        << "Print a list of available tasks.       \n"
        << endl << endl
        << setw(m) << left << "   -t, --template" 
        << "Create a template input parameter file named \n"
        << setw(m) << left << " " 
        << "param.par for the 3DFIT task.        \n"
        << endl << endl
        << setw(m) << left << "   -v, --version" 
        << "Version info    \n"
        << endl;
}
    

void versionInfo(std::ostream& ostr, char **argv) {
    // Print info on the code version and compiler
    
    std::string compiler = "\nBuilt with ";
#if defined(__clang__)
    /* Clang/LLVM. ---------------------------------------------- */
    compiler += "Clang/LLVM "+ std::string(__VERSION__);
#elif defined(__ICC) || defined(__INTEL_COMPILER)
    /* Intel ICC/ICPC. ------------------------------------------ */
    compiler += "Intel ICC/ICPC "
#elif defined(__GNUC__) || defined(__GNUG__) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    /* GNU GCC/G++. --------------------------------------------- */
    compiler += "GNU GCC ";
#else
    compiler += "Unknown compiler";
#endif
    
#if defined (__VERSION__)
    compiler += std::string(__VERSION__);
#endif
    
#if defined(__DATE__) && defined (__TIME__)
    compiler += " on "+std::string(__DATE__)+" "+std::string(__TIME__);
#endif

    std::string flags = "\nOptions: ";
#if defined(HAVE_GNUPLOT)
    flags += "GNUPLOT ";
#endif
#if defined(HAVE_PYTHON)
    flags += "PYTHON ";
#endif
#if defined(HAVE_FFTW3)
    flags += "FFTW3 ";
#endif
#if defined(MACOSX)
    flags += "MACOSX ";
#endif
#if defined(_OPENMP)
    flags += "OPENMP ";
#endif
#if defined(HAVE_MPI)
    flags += "MPI ";
#endif

    ostr << "\nThis is BBarolo version " << BBVERSION <<  std::endl;    
    ostr << "\nExecutable: " << argv[0];
        
    ostr << compiler << flags << std::endl << std::endl;  
}


void welcomeMessage(std::ostream& ostr) {

    std::string adj = randomAdjective(1);
    std::string welcomemsg = " Welcome to BBarolo v"+std::string(BBVERSION);
    
    int size = 70;
    int sizel = (size-10+welcomemsg.size())/2.;
    int sizer = size-sizel;
        
    ostr << "\n" << setfill('=') << setw(size) << right << "\n"
         << "====="  << setfill(' ') << setw(size-5) << "=====\n"
         << "=====" << setfill(' ') << setw(sizel) << right << welcomemsg << setw(sizer-5) << "=====\n"
         << "====="  << setfill(' ') << setw(size-5) << "=====\n"
         << setfill('=') << setw(size) << right << "\n"
         << setfill(' ') << endl << left;
    
}

