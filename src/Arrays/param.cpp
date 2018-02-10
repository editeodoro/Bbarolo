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
#include <unistd.h>
#include <Arrays/param.hh>
#include <Utilities/utils.hh>



Param::Param() {
    
    defaultValues();
}


void Param::defaultValues() {
    
    imageFile           = "";
    imageList           = "NONE";
    outFolder           = "";
    verbose             = true;
    showbar             = true;
    beamFWHM            = 30.;
    checkChannels       = false;
    flagRobustStats     = true;

    globprof            = false;
    massdensmap         = false;
    totalmap            = false;
    velocitymap         = false;
    dispersionmap       = false;
    rmsmap              = false;
    blankCut            = 3;

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

    flagSlitfit         = false;
    wavefile            = "NONE";
    ivarfile            = "NONE";
    linetofit           = "Ha";

    
    flagPV              = false;
    XPOS_PV             = 0;
    YPOS_PV             = 0;
    PA_PV               = 0;
    
    flagEllProf         = false;
    
    threads             = 1;
    debug               = false;
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
    this->flagRobustStats   = p.flagRobustStats;

    this->globprof          = p.globprof;
    this->velocitymap       = p.velocitymap;
    this->totalmap          = p.totalmap;
    this->massdensmap       = p.massdensmap;
    this->dispersionmap     = p.dispersionmap;
    this->rmsmap            = p.rmsmap;
    this->blankCut          = p.blankCut;
         
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
    if(argc<2 || argc>3){
        helpscreen();
        returnValue = false;
    }
    else {      
        std::string file;
        defaultValues();
        char c=0;
        while(( c = getopt(argc,argv,"p:f:dt")) != -1){
            switch(c) {
                case 'p':
                file = optarg;
                status = readParams(file);
                if(status==1){
                    std::cout << "Could not open parameter file " << file << ".\n";
                }
                else if(status==2){
                    std::cout << "\nError opening list file: " << imageList << " doesn't exist.\n\n";
                }
                else if(status==3) {
                    std::cout << "\nWrong input file. Exiting...\n\n";
                }
                else returnValue = true;
                break;
                
                case 'f':
                file = optarg;
                imageFile = file;
                beamFWHM /= 3600.;
                parSE.flagSearch=true;
                returnValue = true;
                break;
                
                case 'd':
                printDefaults();
                returnValue = false;
                break;
                
                case 't':
                createTemplate();
                returnValue = false;
                break;
                
                default :
                helpscreen();
                returnValue = false;
                break;
            }
        }
    }
    
    return returnValue;
}


int Param::readParams(std::string paramfile) {
    
  /// The parameters are read in from a disk file, on the assumption that each
  /// line of the file has the format "parameter value" (eg. alphafdr 0.1)
  /// 
  /// The case of the parameter name does not matter, nor does the
  /// formatting of the spaces (it can be any amount of whitespace or
  /// tabs). 
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
            std::string arg;
            ss >> arg;
            arg = makelower(arg);
            if(arg=="fitsfile")         imageFile = readFilename(ss);
            if(arg=="fitslist")         imageList = readFilename(ss);
            if(arg=="verbose")          verbose = readFlag(ss);
            if(arg=="outfolder")        outFolder = readFilename(ss);
            if(arg=="threads")          threads = readIval(ss);
            if(arg=="debug")            debug = readFlag(ss);
            if(arg=="showbar")          showbar = readFlag(ss);
            if(arg=="beamfwhm")         beamFWHM = readFval(ss);
            if(arg=="checkchannels")    checkChannels = readFlag(ss);
                
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
        
            if(arg=="2dfit")     flagRing = readFlag(ss);
            if (arg=="ellprof")  flagEllProf = readFlag(ss);
        
            // SHARED PARAMETERS BETWEEN GALMOD, GALFIT AND GALWIND
            if(arg=="galfit")    parGF.flagGALFIT  = readFlag(ss);
            if(arg=="3dfit")     parGF.flagGALFIT  = readFlag(ss);
            if(arg=="galmod")    parGM.flagGALMOD  = readFlag(ss);
            if(arg=="galwind")   parGW.flagGALWIND = readFlag(ss);
            if(arg=="spacepar")  flagSpace  = readFlag(ss);
            if(arg=="box")       readVec<int>(ss,BOX,6);
            if(arg=="nradii")    parGM.NRADII = parGF.NRADII               = readIval(ss);
            if(arg=="radii")     parGM.RADII  = parGF.RADII                = readFilename(ss);
            if(arg=="radsep")    parGM.RADSEP = parGF.RADSEP               = readDval(ss);
            if(arg=="vrot")      parGM.VROT   = parGF.VROT                 = readFilename(ss);
            if(arg=="vrad")      parGM.VRAD   = parGF.VRAD                 = readFilename(ss);
            if(arg=="z0")        parGM.Z0     = parGF.Z0                   = readFilename(ss);
            if(arg=="vvert")     parGM.VVERT                               = readFilename(ss);
            if(arg=="dvdz")      parGM.DVDZ                                = readFilename(ss);
            if(arg=="zcyl")      parGM.ZCYL                                = readFilename(ss);
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
            if(arg=="nlines")    parGM.NLINES = parGF.NLINES               = readIval(ss);
            
            // GALFIT ONLY PARAMETERS
            if(arg=="deltainc")  parGF.DELTAINC   = readFval(ss);
            if(arg=="deltapa")   parGF.DELTAPHI   = readFval(ss);
            if(arg=="ftype")     parGF.FTYPE      = readIval(ss);
            if(arg=="wfunc")     parGF.WFUNC      = readIval(ss);
            if(arg=="tol")       parGF.TOL        = readDval(ss);
            if(arg=="free")      parGF.FREE       = readFilename(ss);
            if(arg=="side")      parGF.SIDE       = readFilename(ss);
            if(arg=="mask")      parGF.MASK       = readFilename(ss);
            if(arg=="bweight")   parGF.BWEIGHT    = readIval(ss);
            if(arg=="twostage")  parGF.TWOSTAGE   = readFlag(ss);
            if(arg=="flagerrors")parGF.flagERRORS = readFlag(ss);
            if(arg=="norm")      parGF.NORM       = makeupper(readFilename(ss));
            if(arg=="polyn")     parGF.POLYN      = readFilename(ss);
            if(arg=="startrad")  parGF.STARTRAD   = readIval(ss);
            if(arg=="redshift")  parGF.REDSHIFT   = readDval(ss);
            if(arg=="restwave")  parGF.RESTWAVE   = readDval(ss);
            if(arg=="distance")  parGF.DISTANCE   = readFval(ss);

            // GALWIND ONLY PARAMETERS
            if(arg=="htot")      parGW.HTOT       = readFval(ss);
            if(arg=="openang")   parGW.OPENANG    = readFval(ss);
            if(arg=="ntot")      parGW.NTOT       = readIval(ss);
            if(arg=="vwind")     parGW.VWIND      = readFilename(ss);
            if(arg=="denstype")  parGW.DENSTYPE   = readIval(ss);
 
            if(arg=="p1")        P1 = readFilename(ss);
            if(arg=="p1par")     readVec<float>(ss,P1p,3);
            if(arg=="p2")        P2 = readFilename(ss);
            if(arg=="p2par")     readVec<float>(ss,P2p,3);
            
            if(arg=="smooth")    flagSmooth = readFlag(ss);
            if(arg=="box")       readVec<int>(ss,BOX,6);
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

            if(arg=="slitfit")    flagSlitfit = readFlag(ss);
            if(arg=="wavefile")   wavefile = readFilename(ss);
            if(arg=="ivarfile")   ivarfile = readFilename(ss);
            if(arg=="linetofit")  linetofit = readFilename(ss);

            if (arg=="flagpv")    flagPV = readFlag(ss);
            if (arg=="xpos_pv")   XPOS_PV = readFval(ss);
            if (arg=="ypos_pv")   YPOS_PV = readFval(ss);
            if (arg=="pa_pv")     PA_PV = readFval(ss);
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
    

    beamFWHM /= 3600.;

    if (!checkPars()) return 3;

    return 0;

}
 
  
bool Param::checkPars() {
    
    bool good = true;
    
    checkHome(imageFile);
    checkHome(outFolder);
    if (outFolder!="" && outFolder[outFolder.size()-1]!='/') outFolder.append("/");

    if (!verbose || threads>1) showbar=false;

    // Checking MASK parameter
    std::string maskstr = makeupper(parGF.MASK);
    if (maskstr=="NONE" || maskstr=="SMOOTH" || maskstr=="SEARCH" ||
        maskstr=="THRESHOLD" || maskstr=="NEGATIVE") {
        parGF.MASK = maskstr;
    }
    else if (maskstr.find("FILE(")!=std::string::npos) {
        size_t first = maskstr.find_first_of("(");
        std::string sub1 = maskstr.substr(0,first);
        std::string sub2 = parGF.MASK.substr(first);
        parGF.MASK = sub1+sub2;
    }
    else {
        std::cout << "\n ERROR: Unknown type of mask: " << parGF.MASK << std::endl;
        std::cout << " Setting to SMOOTH" << std::endl;
        parGF.MASK = "SMOOTH";
    }


    // Checking parameters for source finder
    if (parSE.flagSearch) {
        if(parSE.searchType != "spatial" && parSE.searchType != "spectral"){
            std::cout << "You have requested a search type of \""<<parSE.searchType<<"\".\n"
                      << "Only \"spectral\" and \"spatial\" are accepted. Setting to \"spectral\".\n";
            parSE.searchType = "spectral";
        }
    
        if(parSE.flagGrowth){
            if(parSE.UserThreshold && ((parSE.threshold<parSE.growthThreshold)||
              (parSE.snrCut<parSE.growthCut))) {
                  std::cout << "Your \"growthThreshold\" parameter" << parSE.growthThreshold
                            <<" is larger than your \"threshold\""  << parSE.threshold << std::endl;
                  good = false;
            }
            if(!parSE.UserThreshold && (parSE.snrCut<parSE.growthCut)) {
                std::cout << "Your \"growthCut\" parameter " << parSE.growthCut
                          << " is larger than your \"snrCut\"" << parSE.snrCut << std::endl;
                good = false;
            }
        
            if(!good) {
                std::cout << "The growth function is being turned off\n.";
                parSE.flagGrowth=false;
                good = true;
            }
        }
    }
    
    // Checking parameters for 3DFIT and GALMOD
    if (parGF.flagGALFIT || flagSpace || parGM.flagGALMOD || flagSlitfit) {
        std::string str =" is not an optional parameter. Please specify it in the input file or set flagSearch=true";

        if (parGF.NRADII==-1 && parGF.RADII=="-1")  {
            std::cout << "3DFIT error: NRADII or RADII" << str << std::endl;
            good = false;
        }

        if (parGM.flagGALMOD) {
            if (parGM.XPOS=="-1") {
                std::cout << "3DFIT error: XPOS" << str << std::endl;
                good = false;
            }
            if (parGM.YPOS=="-1") {
                std::cout << "3DFIT error: YPOS" << str << std::endl;
                good = false;
            }
            if (parGM.RADSEP==-1 && parGM.RADII=="-1")  {
                std::cout << "3DIT error: RADSEP" << str << std::endl;
                good = false;
            }
            if (parGM.VSYS=="-1") {
                std::cout << "3DFIT error: VSYS" << str << std::endl;
                good = false;
            }
            if (parGM.VROT=="-1") {
                std::cout << "3DFIT error: VROT" << str << std::endl;
                good = false;
            }
            if (parGM.VDISP=="-1")    {
                std::cout << "3DFIT error: VDISP" << str << std::endl;
                good = false;
            }
            if (parGM.INC=="-1")  {
                std::cout << "3DFIT error: INC" << str << std::endl;
                good = false;
            }
            if (parGM.PHI=="-1")  {
                std::cout << "3DFIT error: PHI" << str << std::endl;
                good = false;
            }
            if (parGM.Z0=="-1")   {
                std::cout << "3DFIT error: Z0" << str << std::endl;
                good = false;
            }
        }

        if (parGF.NORM!="NONE" && parGF.NORM!="AZIM" && parGF.NORM!="LOCAL") {
            if (!(parGM.flagGALMOD && parGF.NORM=="BOTH")) {
                std::cout << " ERROR: Unknown type of normalization: " << parGF.NORM << std::endl;
                std::cout << "Setting to LOCAL" << std::endl;
                parGF.NORM="LOCAL";
            }
        }

        if (parGF.flagGALFIT) {
            if (parGF.FREE=="") {
                std::cout << "3DFIT error: FREE" << str << std::endl;
                good = false;
            }

            if (parGF.FTYPE<1 || parGF.FTYPE>4) {
                std::cout << "3DFIT warning: ";
                std::cout << "Not valid argument for FTYPE parameter. ";
                std::cout << "Assuming 2 (|mod-obs|).\n";
                parGF.FTYPE = 2;
            }

            if (parGF.WFUNC<0 || parGF.WFUNC>2) {
                std::cout << "3DFIT warning: ";
                std::cout << "Not valid argument for WFUNC parameter. ";
                std::cout << "Assuming 1 (|cos(θ)| weighting function).\n";
                parGF.WFUNC = 1;
            }
            
            string SIDE = parGF.SIDE;
            if (SIDE=="A"||SIDE=="APP"||SIDE=="APPR"||SIDE=="APPROACHING") SIDE = "A";
            else if (SIDE=="R"||SIDE=="REC"||SIDE=="RECED"||SIDE=="RECEDING") SIDE = "R";
            else if (SIDE=="B"||SIDE=="BOTH") SIDE = "B";
            else if (SIDE=="BS"||SIDE=="S"||SIDE=="SINGLE"||SIDE=="BOTH SINGLE") SIDE = "S";
            else {
                std::cout << "3DFIT warning: ";
                std::cout << "Not valid argument for SIDE parameter. ";
                std::cout << "Assuming B (fitting both sides of the galaxy).\n";
                parGF.SIDE = "B";
            }

        }

        if (parGF.LTYPE<1 || parGF.LTYPE>5) {
            std::cout << "3DFIT warning: ";
            std::cout << "Not valid argument for LTYPE parameter. ";
            std::cout << "Assuming 1 (gaussian layer).\n";
            parGF.LTYPE = 1;
        }

         if ((parGF.RESTWAVE!=-1 && parGF.REDSHIFT==-1) || (parGF.RESTWAVE==-1 && parGF.REDSHIFT!=-1)) {
            std::cout<< "3DFIT warning: Restwave and Redshift must be set both. Exiting...\n";
            std::terminate();
         }
        
        if (flagSlitfit) {
            checkHome(wavefile);
            checkHome(ivarfile);
            std::cout<< "3DFIT warning: Galfit and Slitfit cannot be run at the same time. "
                     << "Switching off Slitfit. \n";
            flagSlitfit = false;
        }
            
    }
    
    // Check parameters for GALWIND task
    if (parGW.flagGALWIND) {
        std::string str =" is not an optional parameter. Please specify it in the input file.";
        if (parGW.XPOS=="-1") {
            std::cout << "GALWIND error: XPOS" << str << std::endl;
            good = false;
        }
        if (parGW.YPOS=="-1") {
            std::cout << "GALWIND error: YPOS" << str << std::endl;
            good = false;
        }
        if (parGW.VSYS=="-1") {
            std::cout << "GALWIND error: VSYS" << str << std::endl;
            good = false;
        }
        if (parGW.VDISP=="-1")    {
            std::cout << "GALWIND error: VDISP" << str << std::endl;
            good = false;
        }
        if (parGW.INC=="-1")  {
            std::cout << "GALWIND error: INC" << str << std::endl;
            good = false;
        }
        if (parGW.PHI=="-1")  {
            std::cout << "GALWIND error: PHI" << str << std::endl;
            good = false;
        }
        if (parGW.DENS=="-1")   {
            std::cout << "GALWIND error: DENS" << str << std::endl;
            good = false;
        }
        if (parGW.OPENANG==-1)   {
            std::cout << "GALWIND error: OPENANG" << str << std::endl;
            good = false;
        }
        if (parGW.HTOT==-1)   {
            std::cout << "GALWIND error: HTOT" << str << std::endl;
            good = false;
        }
        if (parGW.VWIND=="-1")   {
            std::cout << "GALWIND error: VWIND" << str << std::endl;
            good = false;
        }        
    }
    
    // Checking parameters for SMOOTH
    if (flagSmooth) {
        if (bmaj==-1 && bmin==-1 && linear==-1 && factor==-1) {
            std::cout << "SMOOTH error: "
                      << "you need to specify either the new beam (BMAJ, BMIN, BPA) "
                      << "or the parameters LINEAR and DISTANCE.\n";
            good = false;
        }
        if (linear!=-1 && parGF.DISTANCE==-1) {
            std::cout << "SMOOTH error: "
                      << "with LINEAR parameter you must specify also the DISTANCE. ";
            good = false;
        }
        if (bmaj!=-1 && bmin==-1) {
            bmin = bmaj;
        }
    }
    
    return good;
}


void Param::printDefaults (std::ostream& theStream) {
    
    Param par;
    
    theStream.setf(std::ios::left);
    theStream  <<"\n--------------------- Default parameters ---------------------\n"<<std::endl;
    theStream  << std::setfill('.');
   
    recordParam(theStream, "[FITSFILE]", "FITS file to be analysed", par.getImageFile());
    recordParam(theStream, "[FITSLIST]", "List of FITS files to be analysed", par.getImageList());
    recordParam(theStream, "[THREADS]", "Number of threads", par.getThreads());
    
    theStream  <<"--------------"<<std::endl;

    recordParam(theStream, "[checkChannels]", "Checking for bad channels in the cube", stringize(par.getCheckCh()));
    recordParam(theStream, "[beamFWHM]", "Size of the beam (arcsec)", par.getBeamFWHM());
    recordParam(theStream, "[flagRobustStats]", "Using Robust statistics?", stringize(par.getFlagRobustStats()));
    
    recordParam(theStream, "[SEARCH]", "Searching for sources in cube?", stringize(par.getParSE().flagSearch));
    recordParam(theStream, "[searchType]", "   Type of searching performed", par.getParSE().searchType);
    recordParam(theStream, "[minPix]", "   Minimum # Pixels in a detection", par.getParSE().minPix);
    recordParam(theStream, "[minChannels]", "   Minimum # Channels in a detection", par.getParSE().minChannels);
    recordParam(theStream, "[minVoxels]", "   Minimum # Voxels in a detection", par.getParSE().minVoxels);
    recordParam(theStream, "[maxChannels]", "   Maximum # Channels in a detection", par.getParSE().maxChannels);
    recordParam(theStream, "[maxAngsize]", "   Max angular size of a detection in arcmin", par.getParSE().maxAngSize);
    recordParam(theStream, "[flagAdjacent]", "   Using Adjacent-pixel criterion?", stringize(par.getParSE().flagAdjacent));
    recordParam(theStream, "[threshSpatial]", "   Max. spatial separation for merging", par.getParSE().threshSpatial);
    recordParam(theStream, "[threshVelocity]", "   Max. velocity separation for merging", par.getParSE().threshVelocity);
    recordParam(theStream, "[RejectBeforeMerge]", "   Reject objects before merging?", stringize(par.getParSE().RejectBeforeMerge));
    recordParam(theStream, "[TwoStageMerging]", "   Merge objects in two stages?", stringize(par.getParSE().TwoStageMerging));
    recordParam(theStream, "[threshold]", "   Detection Threshold", par.getParSE().threshold);
    recordParam(theStream, "[snrCut]", "   SNR Threshold (in sigma)", par.getParSE().snrCut);
    recordParam(theStream, "[flagGrowth]", "   Growing objects after detection?", stringize(par.getParSE().flagGrowth));
    recordParam(theStream, "[growthCut]", "   SNR Threshold for growth", par.getParSE().growthCut);
    recordParam(theStream, "[growthThreshold]", "   Threshold for growth", par.getParSE().growthThreshold);
    
    recordParam(theStream, "[2DFIT]", "Fitting velocity field with a ring model?", stringize(par.getFlagRing()));
    
    recordParam(theStream, "[globalProfile]", "Saving the global profile?", stringize(par.getGlobProf()));
    recordParam(theStream, "[totalMap]", "Saving 0th moment map to FITS file?", stringize(par.getTotalMap()));
    recordParam(theStream, "[massdensMap]", "Saving HI mass density map to FITS file?", stringize(par.getMassDensMap()));
    recordParam(theStream, "[velocityMap]", "Saving 1st moment map to FITS file?", stringize(par.getVelMap()));
    recordParam(theStream, "[dispersionMap]", "Saving 2th moment map to FITS file?", stringize(par.getDispMap()));
    recordParam(theStream, "[rmsMap]", "Saving RMS map to FITS file?", stringize(par.getRMSMap()));
    recordParam(theStream, "[blankCut]", "SNR clipping cut for blanked areas", par.getBlankCut());
    
    recordParam(theStream, "[SMOOTH]", "Smoothing the datacube?", stringize(par.getflagSmooth()));
    recordParam(theStream, "[FFT]", "Using FFT for convolution?", stringize(par.getflagFFT()));
    recordParam(theStream, "[REDUCE]", "Reducing datacube?", stringize(par.getflagReduce()));
    recordParam(theStream, "[BOX]", "Sub-region to be used?", "entire dataset");
    
    recordParam(theStream, "[3DFIT]", "Fitting a 3D model to the datacube?", stringize(par.getParGF().flagGALFIT));
    recordParam(theStream, "[GALMOD]", "Writing a 3D model?", stringize(par.getflagGalMod()));
    recordParam(theStream, "[DENS]", "   Global column density of gas (atoms/cm2)", par.getParGF().DENS);
    recordParam(theStream, "[LTYPE]", "   Layer type along z direction", "gaussian");
    recordParam(theStream, "[FTYPE]", "   Function to be minimized", "|m-o|");
    recordParam(theStream, "[WFUNC]", "   Weighting function", "|cos(θ)|");
    recordParam(theStream, "[TOL]", "   Minimization tolerance", par.getParGF().TOL); 
    recordParam(theStream, "[MASK]", "   Type of mask", (par.getParGF().MASK));
    recordParam(theStream, "[FREE]", "   Parameters to be minimized", par.getParGF().FREE);
    recordParam(theStream, "[SIDE]", "   What side of the galaxy to be used", (par.getParGF().SIDE)); 
    recordParam(theStream, "[TWOSTAGE]", "   Two stages minimization?", stringize(par.getParGF().TWOSTAGE));
    recordParam(theStream, "[POLYN]", "     Degree of polynomial fitting angles?", par.getParGF().POLYN);
    recordParam(theStream, "[flagErrors]", "   Estimating errors?", stringize(par.getParGF().flagERRORS));
    
    recordParam(theStream, "[GALWIND]", "Generating a 3D datacube with a wind model?", stringize(par.getParGW().flagGALWIND));
    recordParam(theStream, "[NTOT]", "   Number of layers/cylinder for each cone", par.getParGW().NTOT);
    recordParam(theStream, "[DENSTYPE]", "   How to distribute density in layers", par.getParGW().DENSTYPE);
    

    theStream  << std::endl <<"-----------------------------";
    theStream  << "------------------------------\n\n";
    theStream  << std::setfill(' ');
    theStream.unsetf(std::ios::left);
    
}


void Param::createTemplate() {
    
    using namespace std;
    
    int m = 12;
    std::ofstream parf;
    parf.open("param.par");
    
    parf << "// This is a template input file for the Galfit utility.\n";
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

    parf << "// OPTIONAL: Using a mask for the minimization (default is " << parGF.MASK << "):\n";
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
    
    std::cout << "\n A template parameter input file (\"param.par\") for BBarolo " 
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
    /// Lists the parameters, a description of them, and their value.

   
    theStream.setf(std::ios::left);
    theStream  <<"\n-------------------------- Parameters -------------------------\n"<<std::endl;
    theStream  << std::setfill('.');
   
    if (par.getImageList()=="NONE") 
        recordParam(theStream, "[FITSFILE]", "FITS file to be analysed", par.getImageFile());
    else recordParam(theStream, "[FITSLIST]", "List of FITS files to be analysed", par.getImageList());
    recordParam(theStream, "[THREADS]", "Number of threads", par.getThreads());
    
    theStream  <<"--------------"<<std::endl;
  
    recordParam(theStream, "[checkChannels]", "Checking for bad channels in the cube", stringize(par.getCheckCh()));
    recordParam(theStream, "[flagRobustStats]", "Using Robust statistics?", stringize(par.getFlagRobustStats()));
    if (par.getParGF().DISTANCE!=-1)
        recordParam(theStream, "[DISTANCE]", "Distance of the galaxy (Mpc)?", par.getParGF().DISTANCE);
  
    // PARAMETERS FOR SEARCH TASK
    recordParam(theStream, "[SEARCH]", "Searching for sources in cube?", stringize(par.getParSE().flagSearch));
    if (par.getParSE().flagSearch) {
        recordParam(theStream, "[searchType]", "   Type of searching performed", par.getParSE().searchType);
        recordParam(theStream, "[minPix]", "   Minimum # Pixels in a detection", par.getParSE().minPix);
        recordParam(theStream, "[minChannels]", "   Minimum # Channels in a detection", par.getParSE().minChannels);
        recordParam(theStream, "[minVoxels]", "   Minimum # Voxels in a detection", par.getParSE().minVoxels);
        recordParam(theStream, "[maxChannels]", "   Maximum # Channels in a detection", par.getParSE().maxChannels);
        recordParam(theStream, "[maxAngsize]", "   Max angular size of a detection in arcmin", par.getParSE().maxAngSize);
        recordParam(theStream, "[flagAdjacent]", "   Using Adjacent-pixel criterion?", stringize(par.getParSE().flagAdjacent));
        if(!par.getParSE().flagAdjacent)
            recordParam(theStream, "[threshSpatial]", "   Max. spatial separation for merging", par.getParSE().threshSpatial);
        recordParam(theStream, "[threshVelocity]", "   Max. velocity separation for merging", par.getParSE().threshVelocity);
        recordParam(theStream, "[RejectBeforeMerge]", "   Reject objects before merging?", stringize(par.getParSE().RejectBeforeMerge));
        recordParam(theStream, "[TwoStageMerging]", "   Merge objects in two stages?", stringize(par.getParSE().TwoStageMerging));
        if(par.getParSE().UserThreshold){
            recordParam(theStream, "[threshold]", "   Detection Threshold", par.getParSE().threshold);
        }
        else {
            recordParam(theStream, "[snrCut]", "   SNR Threshold (in sigma)", par.getParSE().snrCut);
        }
        recordParam(theStream, "[flagGrowth]", "   Growing objects after detection?", stringize(par.getParSE().flagGrowth));
        if(par.getParSE().flagGrowth) {                  
            if(par.getParSE().flagUserGrowthT){
                recordParam(theStream, "[growthThreshold]", "   Threshold for growth", par.getParSE().growthThreshold);
            }
            else{
                recordParam(theStream, "[growthCut]", "   SNR Threshold for growth", par.getParSE().growthCut);
            }
        }
    }
    

    // PARAMETERS FOR SMOOTH
    recordParam(theStream, "[SMOOTH]", "Smoothing the datacube?", stringize(par.getflagSmooth()));
    if (par.getflagSmooth()) {
        std::string box;
        for (int i=0;i<6;i++) if (par.getBOX(i)!=-1) box += to_string<int>(par.getBOX(i))+" ";
        if (box=="") box = "NONE";
        recordParam(theStream, "[BOX]", "   Sub-region to be used?", box);
        if (par.getLinear()!=-1) {
            recordParam(theStream, "[LINEAR]", "   New linear resolution (kpc)", par.getLinear());
        }
        else if (par.getBmaj()!=-1 && par.getBmin()!=-1){
            recordParam(theStream, "[BMAJ]", "   New major beam (arcsec)", par.getBmaj());
            recordParam(theStream, "[BMIN]", "   New minor beam (arcsec)", par.getBmin());
            recordParam(theStream, "[BPA]", "   New beam position angle (degree)", par.getBpa());
        }
        else recordParam(theStream, "[FACTOR]", "   New beam factor (times old beam)", par.getFactor());

        recordParam(theStream, "[FFT]", "   Using FFT for convolution?", stringize(par.getflagFFT()));
        recordParam(theStream, "[REDUCE]", "   Reducing datacube?", stringize(par.getflagReduce()));
    }
    
    


    recordParam(theStream, "[3DFIT]", "Fitting a 3D model to the datacube?", stringize(par.getflagGalFit()));
    recordParam(theStream, "[GALMOD]", "Writing a 3D model?", stringize(par.getflagGalMod()));

    if (par.getflagGalFit() || par.getflagGalMod()) {
        if (par.getParGF().RADII!="-1") {
            recordParam(theStream, "[RADII]", "   Radii for rings", par.getParGF().NRADII);
        }
        else {
            recordParam(theStream, "[NRADII]", "   Number of radii", par.getParGF().NRADII);
            recordParam(theStream, "[RADSEP]", "   Separation between radii (arcsec)", par.getParGF().RADSEP);
        }
        
        recordParam(theStream, "[XPOS]", "   X center of the galaxy (pixel)", par.getParGF().XPOS);
        recordParam(theStream, "[YPOS]", "   Y center of the galaxy (pixel)", par.getParGF().YPOS);
        recordParam(theStream, "[VSYS]", "   Systemic velocity of the galaxy (km/s)", par.getParGF().VSYS);
        recordParam(theStream, "[VROT]", "   Initial global rotation velocity (km/s)", par.getParGF().VROT);
        recordParam(theStream, "[VRAD]", "   Initial global radial velocity (km/s)", par.getParGF().VRAD);
        recordParam(theStream, "[VDISP]", "   Initial global velocity dispersion (km/s)", par.getParGF().VDISP);
        recordParam(theStream, "[INC]", "   Initial global inclination (degrees)", par.getParGF().INC);
        recordParam(theStream, "[PA]", "   Initial global position angle (degrees)", par.getParGF().PHI);
        recordParam(theStream, "[Z0]", "   Scale height of the disk (arcsec)", par.getParGF().Z0);
        recordParam(theStream, "[DENS]", "   Global column density of gas (atoms/cm2)", par.getParGF().DENS);
        recordParam(theStream, "[FREE]", "   Parameters to be minimized", par.getParGF().FREE);
        recordParam(theStream, "[MASK]", "   Type of mask", par.getParGF().MASK);
        recordParam(theStream, "[SIDE]", "   Side of the galaxy to be used", (par.getParGF().SIDE)); 
        recordParam(theStream, "[NORM]", "   Type of normalization", (par.getParGF().NORM));

        
        std::string ltype;
        switch (par.getParGF().LTYPE) {
            case 1:
                ltype = "gaussian";
                break;
            case 2:
                ltype = "sech2";
                break;
            case 3:
                ltype = "exponential";
                break;
            case 4:
                ltype = "Lorentzian";
                break;
            case 5:
                ltype = "box";
                break;
            default:
                ltype = "";
                break;
        }
        recordParam(theStream, "[LTYPE]", "   Layer type along z direction", ltype);
        std::string ftype;
        switch (par.getParGF().FTYPE) {
            case 1:
                ftype = "chi-squared";
                break;
            case 2:
                ftype = "|m-o|";
                break;
            case 3:
                ftype = "|m-o|/|m+o|";
                break;
            case 4:
                ftype = "(m-o)^2";
                break;
            default:
                ftype = "";
                break;
            
        }
        recordParam(theStream, "[FTYPE]", "   Function to be minimized", ftype);
        std::string wfunc;
        switch (par.getParGF().WFUNC) {
            case 0:
                wfunc = "uniform";
                break;
            case 1:
                wfunc = "|cos(θ)|";
                break;
            case 2:
                wfunc = "cos(θ)^2";
                break;
            default:
                wfunc = "";
                break;
            
        }
        recordParam(theStream, "[WFUNC]", "   Weighting function", wfunc);
        
        recordParam(theStream, "[TOL]", "   Minimization tolerance", par.getParGF().TOL); 
        recordParam(theStream, "[SIDE]", "   What side of the galaxy to be used", (par.getParGF().SIDE)); 
        recordParam(theStream, "[TWOSTAGE]", "   Two stages minimization?", stringize(par.getParGF().TWOSTAGE));
        if (par.getParGF().TWOSTAGE)
            recordParam(theStream, "[POLYN]", "     Degree of polynomial fitting angles?", par.getParGF().POLYN);
        recordParam(theStream, "[flagErrors]", "   Estimating errors?", stringize(par.getParGF().flagERRORS));

    }
    
    // PARAMETERS FOR GALWIND 
    recordParam(theStream, "[GALWIND]", "Generating a 3D datacube with a wind model?", stringize(par.getParGW().flagGALWIND));
    if (par.getParGW().flagGALWIND) {
        recordParam(theStream, "[VWIND]",   "   Radial velocity of the wind (km/s)", par.getParGW().VWIND);
        recordParam(theStream, "[OPENANG]", "   Wind opening angle (degrees)", par.getParGW().OPENANG);
        recordParam(theStream, "[HTOT]",    "   Wind maximum distance (arcsec)", par.getParGW().HTOT);
        recordParam(theStream, "[XPOS]",    "   X center of the galaxy (pixel)", par.getParGW().XPOS);
        recordParam(theStream, "[YPOS]",    "   Y center of the galaxy (pixel)", par.getParGW().YPOS);
        recordParam(theStream, "[VSYS]",    "   Systemic velocity of the galaxy (km/s)", par.getParGW().VSYS);
        recordParam(theStream, "[VDISP]",   "   Global velocity dispersion (km/s)", par.getParGW().VDISP);
        recordParam(theStream, "[INC]",     "   Global inclination (degrees)", par.getParGW().INC);
        recordParam(theStream, "[PA]",      "   Global position angle (degrees)", par.getParGW().PHI);
        recordParam(theStream, "[DENS]",    "   Global column density of gas (atoms/cm2)", par.getParGW().DENS);
        recordParam(theStream, "[NTOT]",    "   Number of layers/cylinder for each cone", par.getParGW().NTOT);
        recordParam(theStream, "[DENSTYPE]","   How to distribute density in layers", par.getParGW().DENSTYPE);
    }
    
    // PARAMETERS FOR 2DFIT
    recordParam(theStream, "[2DFIT]", "Fitting velocity field with a ring model?", stringize(par.getFlagRing()));
    if (par.getFlagRing()) {
        if (par.getParGF().RADII!="-1") {
            recordParam(theStream, "[RADII]", "   Radii for rings", par.getParGF().NRADII);
        }
        else {
            recordParam(theStream, "[NRADII]", "   Number of radii", par.getParGF().NRADII);
            recordParam(theStream, "[RADSEP]", "   Separation between radii (arcsec)", par.getParGF().RADSEP);
        }
        recordParam(theStream, "[XPOS]", "   X center of the galaxy (pixel)", par.getParGF().XPOS);
        recordParam(theStream, "[YPOS]", "   Y center of the galaxy (pixel)", par.getParGF().YPOS);
        recordParam(theStream, "[VSYS]", "   Systemic velocity of the galaxy (km/s)", par.getParGF().VSYS);
        recordParam(theStream, "[VROT]", "   Rotation velocity (km/s)", par.getParGF().VROT);
        recordParam(theStream, "[VRAD]", "   Radial velocity (km/s)", par.getParGF().VRAD);
        recordParam(theStream, "[INC]",  "   Inclination angle (degrees)", par.getParGF().INC);
        recordParam(theStream, "[PA]",   "   Position angle (degrees)", par.getParGF().PHI);
        recordParam(theStream, "[FREE]", "   Parameters to be fit", par.getParGF().FREE);
        recordParam(theStream, "[MASK]", "   Type of mask for velocity map", par.getParGF().MASK);
        recordParam(theStream, "[SIDE]", "   Side of the galaxy to be used", (par.getParGF().SIDE)); 
    }
    
    // PARAMETERS FOR ELLPROF
    recordParam(theStream, "[ELLPROF]", "Deriving radial intensity profile?", stringize(par.getFlagEllProf()));
    if (par.getFlagEllProf()) {
        if (par.getParGF().RADII!="-1") {
            recordParam(theStream, "[RADII]", "   Radii for rings", par.getParGF().NRADII);
        }
        else {
            recordParam(theStream, "[NRADII]", "   Number of radii", par.getParGF().NRADII);
            recordParam(theStream, "[RADSEP]", "   Separation between radii (arcsec)", par.getParGF().RADSEP);
        }
        recordParam(theStream, "[XPOS]", "   X center of the galaxy (pixel)", par.getParGF().XPOS);
        recordParam(theStream, "[YPOS]", "   Y center of the galaxy (pixel)", par.getParGF().YPOS);
        recordParam(theStream, "[INC]",  "   Inclination angle (degrees)", par.getParGF().INC);
        recordParam(theStream, "[PA]",   "   Position angle (degrees)", par.getParGF().PHI);
        recordParam(theStream, "[MASK]", "   Type of mask for intensity map", par.getParGF().MASK);
        recordParam(theStream, "[SIDE]", "   Side of the galaxy to be used", (par.getParGF().SIDE)); 
    }
    
    // PARAMETERS FOR MOMENT MAPS
    if (par.getGlobProf())
        recordParam(theStream, "[globalProfile]", "Saving the global profile?", stringize(par.getGlobProf()));
    if (par.getTotalMap())
        recordParam(theStream, "[totalMap]",      "Saving 0th moment map to FITS file?", stringize(par.getTotalMap()));
    if (par.getMassDensMap())
        recordParam(theStream, "[massdensMap]",   "Saving HI mass density map to FITS file?", stringize(par.getMassDensMap()));
    if (par.getVelMap())
        recordParam(theStream, "[velocityMap]",   "Saving 1st moment map to FITS file?", stringize(par.getVelMap()));
    if (par.getDispMap())
        recordParam(theStream, "[dispersionMap]", "Saving 2th moment map to FITS file?", stringize(par.getDispMap()));
    if (par.getRMSMap())
        recordParam(theStream, "[rmsMap]", "Saving RMS map to FITS file?", stringize(par.getRMSMap()));
    if (par.getMaps()) 
        recordParam(theStream, "[MASK]",          "   Mask used for maps and profile?", par.getMASK());
    
    
    
    theStream  << std::endl <<"-----------------------------";
    theStream  << "------------------------------\n\n";
    theStream  << std::setfill(' ');
    theStream.unsetf(std::ios::left);
 
    return theStream;
}


void helpscreen(std::ostream& Str) {
    
    using std::cout;
    using std::endl;
    int m=26;

    Str << endl << endl
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
        << "                                                          `\"\"\"\"\"\"`   \n\n\n\n "
        << "  Usage:                 BBarolo option [file] \n\n\n"
        << " Options: \n\n"
        << setw(m) << left << "        -f" 
        << "BBarolo runs in default and automatic mode \n"
        << setw(m) << left << " " 
        << "with all the default parameters. [file]   \n"
        << setw(m) << left << " " 
        << "parameter is mandatory and it is the name \n"
        << setw(m) << left << " " 
        << "of the fitsfile to be analysed. BBarolo   \n"
        << setw(m) << left << " " 
        << "will find sources in the cube, estimate   \n" 
        << setw(m) << left << " " 
        << "initial parameters and fit a 3D tilted-   \n"
        << setw(m) << left << " " 
        << "ring model.\n\n"
        << setw(m) << left << " "  
        << "Example:       BBarolo -f myfitsfile.fits \n"
        << endl << endl
        << setw(m) << left << "        -p" 
        << "BBarolo runs with a parameter file. [file] \n"
        << setw(m) << left << " " 
        << "is mandatory and it is the name of the    \n"
        << setw(m) << left << " " 
        << "file where all the wanted parameters have \n" 
        << setw(m) << left << " " 
        << "been listed. We recommend to run BBarolo  \n"
        << setw(m) << left << " " 
        << "this way to achieve better results.\n\n"
        << setw(m) << left << " " 
        << "Example:       BBarolo -p param.par       \n"
        << endl << endl
        << setw(m) << left << "        -d" 
        << "Prints on the screen a list all the availa-\n"
        << setw(m) << left << " " 
        << "ble parameters and their default values.  \n"
        << endl << endl
        << setw(m) << left << "        -t" 
        << "Creates a template input parameter file    \n"
        << setw(m) << left << " " 
        << "named param.par.                          \n"
        << endl << endl;
    }

