// -----------------------------------------------------------------------
// param.cpp: Member functions for the Param class.
// -----------------------------------------------------------------------

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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <unistd.h>
#include "param.hh"
#include "../Utilities/utils.hh"



Param::Param() {
    
    defaultValues();
}


void Param::defaultValues() {
	
	imageFile         	= "";
    imageList			= "NONE";
    outFolder			= "";
    verbose				= true;
    showbar             = true;
    beamFWHM			= 30.;
    checkChannels		= false;
	flagRobustStats   	= true;

    flagSearch			= false;
    searchType        	= "spatial";
    snrCut            	= 5.;
    threshold         	= 0.;
    flagUserThreshold 	= false;
    flagAdjacent      	= true;
    threshSpatial     	= -1;
    threshVelocity    	= -1;
    minChannels       	= -1;
    minPix            	= -1;
    minVoxels         	= -1;
    maxChannels         = -1;
    maxAngSize          = -1;
    RejectBeforeMerge 	= true;
    TwoStageMerging		= true;        
    flagGrowth			= true;
    growthCut			= 3.;
    flagUserGrowthT		= false;
    growthThreshold		= 0.;

    globprof			= false;
    totalmap			= false;
    velocitymap			= false;
    dispersionmap		= false;
    useBlankCube		= false;
    blankCut			= 3;

    flagRing			= false;
    interactive			= false;
    nrings				= 30;

    flagGalFit			= true;
    flagGalMod			= false;
	NRADII				= -1;
	RADII				= "-1";
	XPOS				= "-1";							
	YPOS				= "-1";						
	RADSEP				= -1;						
	VSYS				= "-1";						
	VROT				= "-1";						
	VDISP				= "-1";						
	INC					= "-1";						
    DELTAINC			= 5;
	PHI					= "-1";
    DELTAPHI			= 15;
	Z0					= "-1";							
    DENS				= "-1";
    CDENS				= 10;
	LTYPE				= 1;
	FTYPE				= 2;
    WFUNC				= 2;
	NV					= -1;							
	TOL					= 1.0E-3;						
	FREE				= "VROT VDISP INC PA";
    MASK				= "SMOOTH";
	SIDE				= "B";
	SM					= true;
    NORM                = "LOCAL";
    BWEIGHT				= 2;
    startRAD            = 0;
	TwoStage			= true;
	flagErrors			= false;
    POLYN				= "bezier";
    flagSpace			= false;
	distance			= -1;	

    flagSmooth			= false;
	flagFFT				= true;
	bmaj				= -1;						
	bmin				= -1;
    bpa					= 0;
    obmaj				= -1;
    obmin				= -1;
    obpa				= 0;
	linear				= -1;
    factor				= 2;
	flagReduce			= false;
    smo_out             = "NONE";
    for (int i=0; i<6; i++) BOX[i] = -1;

    flagSlitfit         = false;
    wavefile            = "NONE";
    ivarfile            = "NONE";
    linetofit           = "Ha";
    redshift            = 0.;

    flagPV              = false;
    XPOS_PV             = 0;
    YPOS_PV             = 0;
    PA_PV               = 0;
}

  
Param::Param (const Param& p) {
	
    operator=(p);
}

  
Param& Param::operator= (const Param& p) {
	
    if(this == &p) return *this;
    this->imageFile         = p.imageFile;
    this->imageList			= p.imageList;
    for (unsigned int i=0; i<images.size(); i++) 
		this->images[i] 	= p.images[i];
    this->outFolder			= p.outFolder;
    this->beamFWHM			= p.beamFWHM;
    this->checkChannels		= p.checkChannels;
    this->verbose           = p.verbose; 
    this->showbar           = p.showbar;
    this->flagRobustStats   = p.flagRobustStats;
    this->snrCut            = p.snrCut;
    this->threshold         = p.threshold;
	
	this->flagSearch		= p.flagSearch;
    this->searchType        = p.searchType;
    this->flagAdjacent      = p.flagAdjacent;
    this->threshSpatial     = p.threshSpatial;
    this->threshVelocity    = p.threshVelocity;
    this->minChannels       = p.minChannels;
    this->minPix            = p.minPix;
    this->minVoxels         = p.minVoxels;
    this->maxChannels       = p.maxChannels;
    this->maxAngSize        = p.maxAngSize;
    this->RejectBeforeMerge = p.RejectBeforeMerge;
    this->TwoStageMerging 	= p.TwoStageMerging;
    this->flagGrowth		= p.flagGrowth;
    this->growthCut		= p.growthCut;
    this->growthThreshold	= p.growthThreshold;
    this->flagUserGrowthT	= p.flagUserGrowthT;
    this->flagUserThreshold= p.flagUserThreshold;
    
    this->globprof			= p.globprof;
    this->velocitymap		= p.velocitymap;
    this->totalmap			= p.totalmap;
    this->dispersionmap		= p.dispersionmap;
    this->useBlankCube		= p.useBlankCube;
    this->blankCut			= p.blankCut;
         
    this->flagRing			= p.flagRing;
    this->interactive		= p.interactive;
    this->nrings			= p.nrings;
 
	this->flagGalFit		= p.flagGalFit;
	this->flagGalMod		= p.flagGalMod;
	for (int i=0; i<6; i++) this->BOX[i] = p.BOX[i];
	this->NRADII			= p.NRADII;
	this->RADII			= p.RADII;
	this->XPOS				= p.XPOS;							
	this->YPOS				= p.YPOS;						
	this->RADSEP			= p.RADSEP;						
	this->VSYS				= p.VSYS;						
	this->VROT				= p.VROT;						
	this->VDISP				= p.VDISP;						
	this->INC				= p.INC;						
	this->DELTAINC			= p.DELTAINC;						
	this->PHI				= p.PHI;	
	this->DELTAPHI			= p.DELTAPHI;					
	this->Z0				= p.Z0;							
	this->DENS				= p.DENS;						
	this->CDENS				= p.CDENS;						
	this->LTYPE				= p.LTYPE;	
	this->FTYPE				= p.FTYPE;						
	this->NV				= p.NV;							
	this->TOL				= p.TOL;
	this->TwoStage			= p.TwoStage;
    this->flagErrors		= p.flagErrors;
	this->POLYN				= p.POLYN;						
	this->FREE				= p.FREE;	
	this->MASK				= p.MASK;
	this->SIDE				= p.SIDE;
	this->SM				= p.SM;
	this->WFUNC				= p.WFUNC;
	this->BWEIGHT			= p.BWEIGHT;
    this->NORM              = p.NORM;
    this->startRAD          = p.startRAD;
	
	this->flagSpace		= p.flagSpace;
	this->P1				= p.P1;
	this->P2				= p.P2;
	for (int i=0; i<3; i++) {
		this->P1p[i]		= p.P1p[i];
		this->P2p[i]		= p.P2p[i];
	}	
	
	this->distance			= p.distance;
	
	this->flagSmooth			= p.flagSmooth;
	this->flagFFT			= p.flagFFT;
	this->bmaj				= p.bmaj;						
	this->bmin				= p.bmin;
    this->bpa				= p.bpa;
    this->obmaj				= p.obmaj;
    this->obmin				= p.obmin;
    this->obpa				= p.obpa;
	this->linear			= p.linear;	
	this->factor			= p.factor;	
	this->flagReduce		= p.flagReduce;
    this->smo_out           = p.smo_out;

    this->flagSlitfit       = p.flagSlitfit;
    this->wavefile          = p.wavefile;
    this->ivarfile          = p.ivarfile;
    this->linetofit         = p.linetofit;
    this->redshift          = p.redshift;
    
    this->flagPV            = p.flagPV;
    this->XPOS_PV           = p.XPOS_PV;
    this->YPOS_PV           = p.YPOS_PV;
    this->PA_PV             = p.PA_PV;


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
                flagSearch=true;
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
  /// \param paramfile 	A std::string containing the parameter filename.
  /// 
  /// \return 			1 if the parameter file does not exist. SUCCESS if
  /// 					it is able to read it.


    std::ifstream fin(paramfile.c_str());
    if(!fin.is_open()) return 1;
    std::string line;
    while(!std::getline(fin,line,'\n').eof()){
		if(line[0]!='#'){
			std::stringstream ss;
			ss.str(line);
			std::string arg;
			ss >> arg;
			arg = makelower(arg);
            if(arg=="fitsfile")      	imageFile = readFilename(ss);
            if(arg=="fitslist")         imageList = readFilename(ss);
			if(arg=="verbose")         	verbose = readFlag(ss);
            if(arg=="outfolder")        outFolder = readFilename(ss);
            if(arg=="showbar")          showbar = readFlag(ss);
			if(arg=="beamfwhm") 		beamFWHM = readFval(ss);
			if(arg=="checkchannels")	checkChannels = readFlag(ss);
				
			if(arg=="flagrobuststats") 	flagRobustStats = readFlag(ss); 
			if(arg=="snrcut")          	snrCut = readFval(ss); 
			if(arg=="threshold"){
				threshold = readFval(ss);
				flagUserThreshold = true;
			}
		
			if(arg=="flagsearch")		flagSearch = readFlag(ss);
            if(arg=="searchtype")      	 searchType = readSval(ss);
            if(arg=="flagadjacent")    	 flagAdjacent = readFlag(ss);
            if(arg=="threshspatial")   	 threshSpatial = readFval(ss);
            if(arg=="threshvelocity")  	 threshVelocity = readFval(ss);
            if(arg=="minchannels")     	 minChannels = readIval(ss);
            if(arg=="minvoxels")       	 minVoxels = readIval(ss);
            if(arg=="minpix")          	 minPix = readIval(ss);
            if(arg=="maxchannels")     	 maxChannels = readIval(ss);
            if(arg=="maxangsize")      	 maxAngSize = readFval(ss);
            if(arg=="rejectbeforemerge") RejectBeforeMerge = readFlag(ss);
            if(arg=="twostagemerging") 	 TwoStageMerging = readFlag(ss);
            if(arg=="flaggrowth")      	 flagGrowth = readFlag(ss);
            if(arg=="growthcut")       	 growthCut = readFval(ss);
            if(arg=="growththreshold"){
                growthThreshold = readFval(ss);
                flagUserGrowthT  = true;
            }

			if(arg=="globalprofile") 	globprof = readFlag(ss);
			if(arg=="totalmap") 	   	totalmap = readFlag(ss);
			if(arg=="velocitymap") 	   	velocitymap = readFlag(ss);
			if(arg=="dispersionmap")   	dispersionmap = readFlag(ss);
			if(arg=="blankcube")	   	useBlankCube = readFlag(ss);
			if(arg=="blankcut")			blankCut = readFval(ss);
		
			if(arg=="flagring")      	flagRing = readFlag(ss);
			if(arg=="interactive")      interactive = readFlag(ss);
			if(arg=="numrings")			nrings = readIval(ss);
		
			if(arg=="distance")			distance = readFval(ss);
            if(arg=="galfit")			flagGalFit = readFlag(ss);
            if(arg=="galmod")			flagGalMod = readFlag(ss);
			if(arg=="spacepar")			flagSpace  = readFlag(ss);
            if (arg=="box")				readVec<int>(ss,BOX,6);
            if (arg=="nradii")			NRADII = readIval(ss);
            if (arg=="radii")			RADII = readFilename(ss);
            if (arg=="xpos")			XPOS = makelower(readFilename(ss));
            if (arg=="ypos")			YPOS = makelower(readFilename(ss));
            if (arg=="radsep")          RADSEP = readDval(ss);
            if (arg=="vsys")			VSYS = readFilename(ss);
            if (arg=="vrot")			VROT = readFilename(ss);
            if (arg=="vdisp")			VDISP = readFilename(ss);
            if (arg=="inc")             INC = readFilename(ss);
            if (arg=="deltainc")		DELTAINC = readFval(ss);
            if (arg=="pa")				PHI = readFilename(ss);
            if (arg=="deltapa")         DELTAPHI = readFval(ss);
            if (arg=="z0")				Z0 = readFilename(ss);
            if (arg=="dens")			DENS = readFilename(ss);
            if (arg=="cdens")			CDENS = readFval(ss);
            if (arg=="ltype")			LTYPE = readIval(ss);
            if (arg=="ftype")			FTYPE = readIval(ss);
            if (arg=="wfunc")			WFUNC = readIval(ss);
            if (arg=="nv")				NV = readIval(ss);
            if (arg=="tol")				TOL = readDval(ss);
            if (arg=="free")			FREE = readFilename(ss);
            if (arg=="side")			SIDE = readFilename(ss);
            if (arg=="mask")			MASK = makeupper(readFilename(ss));
            if (arg=="bweight")			BWEIGHT = readIval(ss);
            if (arg=="twostage")		TwoStage = readFlag(ss);
            if (arg=="flagerrors")		flagErrors = readFlag(ss);
            if (arg=="norm")        	NORM = makeupper(readFilename(ss));
            if (arg=="polyn")           POLYN = readFilename(ss);
            if (arg=="sm")				SM = readFlag(ss);
            if (arg=="startrad")        startRAD = readIval(ss);

            if (arg=="p1")				P1 = readFilename(ss);
            if (arg=="p1par")			readVec<float>(ss,P1p,3);
            if (arg=="p2")				P2 = readFilename(ss);
            if (arg=="p2par")			readVec<float>(ss,P2p,3);
			
			if (arg=="smooth")				flagSmooth = readFlag(ss);
			if (arg=="box")					readVec<int>(ss,BOX,6);
			if (arg=="bmaj")				bmaj = readDval(ss);
			if (arg=="bmin")				bmin = readDval(ss);
			if (arg=="bpa")					bpa = readDval(ss);
            if (arg=="obmaj")				obmaj = readDval(ss);
            if (arg=="obmin")				obmin = readDval(ss);
            if (arg=="obpa")				obpa = readDval(ss);
			if (arg=="linear")				linear = readDval(ss);
			if (arg=="factor")				factor = readDval(ss);
			if (arg=="fft")					flagFFT = readFlag(ss);
			if (arg=="reduce")				flagReduce = readFlag(ss);
            if (arg=="smoothoutput")        smo_out = readFilename(ss);

            if (arg=="slitfit")             flagSlitfit = readFlag(ss);
            if (arg=="wavefile")            wavefile = readFilename(ss);
            if (arg=="ivarfile")            ivarfile = readFilename(ss);
            if (arg=="linetofit")           linetofit = readFilename(ss);
            if (arg=="redshift")            redshift = readDval(ss);

            if (arg=="flagpv")              flagPV = readFlag(ss);
            if (arg=="xpos_pv")             XPOS_PV = readFval(ss);
            if (arg=="ypos_pv")             YPOS_PV = readFval(ss);
            if (arg=="pa_pv")               PA_PV = readFval(ss);

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
	
	checkHome(imageFile);
    checkHome(outFolder);
    if (outFolder!="" && outFolder[outFolder.size()-1]!='/') outFolder.append("/");

    if (!verbose) showbar=false;

    // Can only have "spatial" or "spectral" as search types.
	if(flagSearch && searchType != "spatial" && searchType != "spectral"){
        std::cout<<"You have requested a search type of \""<<searchType<<"\".\n"
				<< "Only \"spectral\" and \"spatial\" are accepted. Setting to \"spectral\".\n";
		searchType = "spectral";
    }
	
	if(flagGrowth){
		bool good = true;
		if(flagUserThreshold && ((threshold<growthThreshold)||(snrCut<growthCut))) {
            std::cout << "Your \"growthThreshold\" parameter" << growthThreshold
					  <<" is larger than your \"threshold\""  << threshold << std::endl;
			good = false;
		}
		if(!flagUserThreshold && (snrCut<growthCut)) {
            std::cout << "Your \"growthCut\" parameter " << growthCut
					  << " is larger than your \"snrCut\"" << snrCut << std::endl;
			good = false;
		}
		
		if(!good) {
            std::cout << "The growth function is being turned off\n.";
			flagGrowth=false;
			good = true;
		}
	}

    if (flagGalFit || flagSpace || flagGalMod) {
		bool good = true;
		std::string str =" is not an optional parameter. Please specify it in the input file or set flagSearch=true";

        if (NRADII==0) {
            std::cout << "3DFIT error: NRADII cannot be 0!" << std::endl;
            good = false;
        }

        if (flagGalMod) {
			if (NRADII==-1 && RADII=="-1")	{
                std::cout << "3DFIT error: NRADII or RADII" << str << std::endl;
				good = false;
			}
			if (XPOS=="-1")	{
                std::cout << "3DFIT error: XPOS" << str << std::endl;
				good = false;
			}
			if (YPOS=="-1")	{
                std::cout << "3DFIT error: YPOS" << str << std::endl;
				good = false;
			}
			if (RADSEP==-1 && RADII=="-1")	{
                std::cout << "3DIT error: RADSEP" << str << std::endl;
				good = false;
			}
			if (VSYS=="-1")	{
                std::cout << "3DFIT error: VSYS" << str << std::endl;
				good = false;
			}
			if (VROT=="-1")	{
                std::cout << "3DFIT error: VROT" << str << std::endl;
				good = false;
            }
            if (VDISP=="-1")	{
                std::cout << "3DFIT error: VDISP" << str << std::endl;
				good = false;
            }
			if (INC=="-1")	{
                std::cout << "3DFIT error: INC" << str << std::endl;
				good = false;
			}
			if (PHI=="-1")	{
                std::cout << "3DFIT error: PHI" << str << std::endl;
				good = false;
            }
            if (Z0=="-1")	{
                std::cout << "3DFIT error: Z0" << str << std::endl;
				good = false;
            }
		}
		
        if (flagGalFit) {
            if (FREE=="") {
                std::cout << "3DFIT error: FREE" << str << std::endl;
                good = false;
            }

            if (FTYPE<1 || FTYPE>4) {
                std::cout << "3DFIT warning: ";
                std::cout << "Not valid argument for FTYPE parameter. ";
                std::cout << "Assuming 2 (|mod-obs|).\n";
                FTYPE = 2;
            }

            if (WFUNC<0 || WFUNC>2) {
                std::cout << "3DFIT warning: ";
                std::cout << "Not valid argument for WFUNC parameter. ";
                std::cout << "Assuming 1 (|cos(θ)| weighting function).\n";
                WFUNC = 1;
            }

            if (SIDE=="A"||SIDE=="APP"||SIDE=="APPR"||SIDE=="APPROACHING") SIDE = "A";
            else if (SIDE=="R"||SIDE=="REC"||SIDE=="RECED"||SIDE=="RECEDING") SIDE = "R";
            else if (SIDE=="B"||SIDE=="BOTH") SIDE = "B";
            else if (SIDE=="BS"||SIDE=="S"||SIDE=="SINGLE"||SIDE=="BOTH SINGLE") SIDE = "S";
            else {
                std::cout << "3DFIT warning: ";
                std::cout << "Not valid argument for SIDE parameter. ";
                std::cout << "Assuming B (fitting both sides of the galaxy).\n";
                SIDE = "B";
            }

        }

		if (LTYPE<1 || LTYPE>5) {
            std::cout << "3DFIT warning: ";
            std::cout << "Not valid argument for LTYPE parameter. ";
            std::cout << "Assuming 1 (gaussian layer).\n";
			LTYPE = 1;
		}
		

        if (flagSlitfit==true) {
            std::cout<< "3DFIT warning: Galfit and Splitfit cannot be run at the same time. "
                     << "Switching off Slitfit. \n";
            flagSlitfit=false;
        }
			
        return good;
	}
	
	if (flagSmooth) {
		if (bmaj==-1 && bmin==-1 && linear==-1 && factor==-1) {
            std::cout << "SMOOTH error: "
					  << "you need to specify either the new beam (BMAJ, BMIN, BPA) "
					  << "or the parameters LINEAR and DISTANCE.\n";
			return false;
		}
		if (linear!=-1 && distance==-1) {
            std::cout << "SMOOTH error: "
					  << "with LINEAR parameter you must specify also the DISTANCE. ";
			return false;
		}
		if (bmaj!=-1 && bmin==-1) {
			bmin = bmaj;
		}
	}

    if (flagSlitfit) {
        checkHome(wavefile);
        checkHome(ivarfile);
    }
	
    if (NORM!="NONE" && NORM!="AZIM" && NORM!="LOCAL") {
        std::cout << " ERROR: Unknown type of normalization: " << NORM << std::endl;
        std::cout << "Setting to LOCAL" << std::endl;
        NORM="LOCAL";
    }

    if (!(MASK=="NONE" || MASK=="SMOOTH" || MASK=="SEARCH" || MASK=="THRESHOLD" || MASK=="NEGATIVE")) {
        std::cout << " ERROR: Unknown type of mask: " << MASK << std::endl;
        std::cout << "Setting to SMOOTH" << std::endl;
        MASK="SMOOTH";
    }

	return true;
}


void Param::printDefaults (std::ostream& theStream) {
	
	Param par;
	
	theStream.setf(std::ios::left);
    theStream  <<"\n--------------------- Defaults parameters --------------------\n"<<std::endl;
    theStream  << std::setfill('.');
   
    recordParam(theStream, "[fitsFile]", "Image to be analysed", par.getImageFile());
    recordParam(theStream, "[fitsList]", "List of images to be analysed", par.getImageList());
    
    theStream  <<"--------------"<<std::endl;

	recordParam(theStream, "[checkChannels]", "Checking for bad channels in the cube", stringize(par.getCheckCh()));
	recordParam(theStream, "[beamFWHM]", "Size of the beam (arcsec)", par.getBeamFWHM());
    recordParam(theStream, "[flagRobustStats]", "Using Robust statistics?", stringize(par.getFlagRobustStats()));
    recordParam(theStream, "[flagSearch]", "Searching for sources in cube?", stringize(par.getSearch()));
    recordParam(theStream, "[searchType]", "   Type of searching performed", par.getSearchType());
    recordParam(theStream, "[minPix]", "   Minimum # Pixels in a detection", par.getMinPix());
    recordParam(theStream, "[minChannels]", "   Minimum # Channels in a detection", par.getMinChannels());
    recordParam(theStream, "[minVoxels]", "   Minimum # Voxels in a detection", par.getMinVoxels());
    recordParam(theStream, "[maxChannels]", "   Maximum # Channels in a detection", par.getMaxChannels());
    recordParam(theStream, "[maxAngsize]", "   Max angular size of a detection in arcmin", par.getMaxAngSize());
    recordParam(theStream, "[flagAdjacent]", "   Using Adjacent-pixel criterion?", stringize(par.getFlagAdjacent()));
    recordParam(theStream, "[threshSpatial]", "   Max. spatial separation for merging", par.getThreshS());
	recordParam(theStream, "[threshVelocity]", "   Max. velocity separation for merging", par.getThreshV());
    recordParam(theStream, "[RejectBeforeMerge]", "   Reject objects before merging?", stringize(par.getRejectBeforeMerge()));
    recordParam(theStream, "[TwoStageMerging]", "   Merge objects in two stages?", stringize(par.getTwoStageMerging()));
	recordParam(theStream, "[threshold]", "   Detection Threshold", par.getThreshold());
	recordParam(theStream, "[snrCut]", "   SNR Threshold (in sigma)", par.getCut());
	recordParam(theStream, "[flagGrowth]", "   Growing objects after detection?", stringize(par.getFlagGrowth()));
	recordParam(theStream, "[growthCut]", "   SNR Threshold for growth", par.getGrowthCut());
	
    recordParam(theStream, "[flagRing]", "Fitting velocity field with a ring model?", stringize(par.getFlagRing()));
	recordParam(theStream, "[interactive]", "   Using interactive mode during the fit?", stringize(par.isInteractive()));
	recordParam(theStream, "[numRings]",  "   Number of rings for the model",  par.getNrings());
    
    recordParam(theStream, "[globalProfile]", "Saving the global profile?", stringize(par.getGlobProf()));
    recordParam(theStream, "[totalMap]", "Saving 0th moment map to FITS file?", stringize(par.getTotalMap()));
    recordParam(theStream, "[velocityMap]", "Saving 1st moment map to FITS file?", stringize(par.getVelMap()));
    recordParam(theStream, "[dispersionMap]", "Saving 2th moment map to FITS file?", stringize(par.getDispMap()));
	recordParam(theStream, "[blankCube]", "Using a blanked cube for maps and profile?", stringize(par.getBlankCube()));
	recordParam(theStream, "[blankCut]", "SNR clipping cut for blanked areas", par.getBlankCut());
    
    recordParam(theStream, "[SMOOTH]", "Smoothing the datacube?", stringize(par.getflagSmooth()));
    recordParam(theStream, "[FFT]", "Using FFT for convolution?", stringize(par.getflagFFT()));
	recordParam(theStream, "[REDUCE]", "Reducing datacube?", stringize(par.getflagReduce()));
    recordParam(theStream, "[BOX]", "Sub-region to be used?", "entire dataset");
    
    recordParam(theStream, "[GALFIT]", "Fitting a 3D model to the datacube?", stringize(par.getflagGalFit()));
	recordParam(theStream, "[DENS]", "   Global column density of gas (atoms/cm2)", par.getDENS());
	recordParam(theStream, "[LTYPE]", "   Layer type along z direction", "gaussian");
	recordParam(theStream, "[FTYPE]", "   Function to be minimized", "|m-o|");
	recordParam(theStream, "[WFUNC]", "   Weighting function", "|cos(θ)|");
	recordParam(theStream, "[TOL]", "   Minimization tolerance", par.getTOL());	
    recordParam(theStream, "[MASK]", "   Type of mask", (par.getMASK()));
	recordParam(theStream, "[SIDE]", "   What side of the galaxy to be used", (par.getSIDE()));	
	recordParam(theStream, "[TWOSTAGE]", "   Two stages minimization?", stringize(par.getTwoStage()));
	recordParam(theStream, "[POLYN]", "     Degree of polynomial fitting angles?", par.getPOLYN());
    recordParam(theStream, "[flagErrors]", "   Estimating errors?", stringize(par.getflagErrors()));

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
	
	parf << "// Using the Galfit utility? Must be true!!\n";
	parf << setw(m) << left << "GALFIT";
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
	parf << setw(m) << left << "VDISP" << setw(m) << left << "10" << endl;
	parf << setw(m) << left << "INC"   << setw(m) << left << "60" << endl;
	parf << setw(m) << left << "PA"    << setw(m) << left << "100"<< endl;
	parf << setw(m) << left << "Z0"    << setw(m) << left << "20\n" << endl;
	
	parf << "// Free parameters for the minimization.\n";
	parf << setw(m) << left << "FREE" << setw(m) << left << "VROT VDISP INC PA\n" << endl;
	
    parf << "// OPTIONAL: Function to be minimezed (default is 2):\n";
	parf <<	"// = 1: chi-squared.\n";
	parf <<	"// = 2: |mod-obs|.\n";
	parf <<	"// = 3: |mod-obs|/|mod+obs|.\n";
	parf <<	"// = 4: (mod-obs)^2.\n";
	parf << setw(m) << left << "FTYPE" << setw(m) << left << "1\n"<< endl;

    parf << "// OPTIONAL: Weighting function (default is 1):\n";
	parf <<	"// = 0: uniform weight.\n";
	parf <<	"// = 1: |cos(θ)|.\n";
	parf <<	"// = 2: cos(θ)^2.\n";
	parf <<	"// θ is the azimuthal angle.\n";
	parf << setw(m) << left << "WFUNC" << setw(m) << left << "1\n"<< endl;

    parf << "// OPTIONAL: Layer type along z (default is gaussian):\n";
    parf <<	"// = 1: gaussian layer.\n";
    parf <<	"// = 2: sech2 layer.\n";
    parf <<	"// = 3: exponential layer.\n";
    parf <<	"// = 4: Lorentzian layer.\n";
    parf <<	"// = 5: box layer.;\n";
    parf << setw(m) << left << "LTYPE" << setw(m) << left << "1\n"<< endl;

	parf << "// OPTIONAL: Number of subcloud in a velocity profile.\n";
	parf << "// (default is = total number of channels).\n";
	parf << setw(m) << left << "NV" << setw(m) << left << "60\n"<< endl;
	
	parf << "// OPTIONAL: Surface density of clouds in the plane of ring (1e20)..\n";
	parf << "// (default is = 1).\n";
	parf << setw(m) << left << "CDENS" << setw(m) << left << "10\n"<< endl;

    parf << "// OPTIONAL: Tolerance for the minimization (default is 0.001).\n";
    parf << setw(m) << left << "TOL" << setw(m) << left << "1E-03\n"<< endl;

	parf << "// OPTIONAL: Using a mask for the minimization (default=true).\n";
	parf << setw(m) << left << "MASK" << setw(m) << left << "TRUE\n"<< endl;
	
	parf << "// OPTIONAL: Side of the galaxy to be fitted (default=both):.\n";
	parf <<	"// = A: Approaching.\n";
	parf <<	"// = R: Receding.\n";
	parf <<	"// = B: Both.\n";
//	parf <<	"// = S: Both but separated.\n";
	parf << setw(m) << left << "SIDE" << setw(m) << left << "B\n"<< endl;
	
	parf << "// OPTIONAL: Using a two stages minimization (default=false).\n";
	parf << setw(m) << left << "TWOSTAGE" << setw(m) << left << "false\n"<< endl;
	
	parf << "// OPTIONAL: Degree of polynomial fitting angles (default=2).\n";
	parf << setw(m) << left << "POLYN" << setw(m) << left << "3\n"<< endl;
	
    parf << "// OPTIONAL: Enabling error estimation (default=false).\n";
    parf << setw(m) << left << "flagErrors" << setw(m) << left << "true\n"<< endl;
	
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
        recordParam(theStream, "[FitsFile]", "Image to be analysed", par.getImageFile());
    else recordParam(theStream, "[FitsList]", "List of images to be analysed", par.getImageList());
    
    theStream  <<"--------------"<<std::endl;
  
	recordParam(theStream, "[checkChannels]", "Checking for bad channels in the cube", stringize(par.getCheckCh()));
    recordParam(theStream, "[flagRobustStats]", "Using Robust statistics?", stringize(par.getFlagRobustStats()));
    
    recordParam(theStream, "[flagSearch]", "Searching for sources in cube?", stringize(par.getSearch()));
    if (par.getSearch()) {
		recordParam(theStream, "[searchType]", "   Type of searching performed", par.getSearchType());
		recordParam(theStream, "[minPix]", "   Minimum # Pixels in a detection", par.getMinPix());
		recordParam(theStream, "[minChannels]", "   Minimum # Channels in a detection", par.getMinChannels());
		recordParam(theStream, "[minVoxels]", "   Minimum # Voxels in a detection", par.getMinVoxels());
        recordParam(theStream, "[maxChannels]", "   Maximum # Channels in a detection", par.getMaxChannels());
        recordParam(theStream, "[maxAngsize]", "   Max angular size of a detection in arcmin", par.getMaxAngSize());
		recordParam(theStream, "[flagAdjacent]", "   Using Adjacent-pixel criterion?", stringize(par.getFlagAdjacent()));
		if(!par.getFlagAdjacent()){
			recordParam(theStream, "[threshSpatial]", "   Max. spatial separation for merging", par.getThreshS());
		}	
		recordParam(theStream, "[threshVelocity]", "   Max. velocity separation for merging", par.getThreshV());
		recordParam(theStream, "[RejectBeforeMerge]", "   Reject objects before merging?", stringize(par.getRejectBeforeMerge()));
		recordParam(theStream, "[TwoStageMerging]", "   Merge objects in two stages?", stringize(par.getTwoStageMerging()));
		if(par.getFlagUserThreshold()){
			recordParam(theStream, "[threshold]", "   Detection Threshold", par.getThreshold());
		}
		else {
			recordParam(theStream, "[snrCut]", "   SNR Threshold (in sigma)", par.getCut());
		}
		recordParam(theStream, "[flagGrowth]", "   Growing objects after detection?", stringize(par.getFlagGrowth()));
		if(par.getFlagGrowth()) {			       
			if(par.getFlagUserGrowthThreshold()){
				recordParam(theStream, "[growthThreshold]", "   Threshold for growth", par.getGrowthThreshold());
			}
			else{
				recordParam(theStream, "[growthCut]", "   SNR Threshold for growth", par.getGrowthCut());
			}
		}
	}
            
    recordParam(theStream, "[flagRing]", "Fitting velocity field with a ring model?", stringize(par.getFlagRing()));
    if (par.getFlagRing()) {
		recordParam(theStream, "[interactive]", "   Using interactive mode during the fit?", stringize(par.isInteractive()));
		if (!par.isInteractive()) 
			recordParam(theStream, "[numRings]",  "   Number of rings for the model",  par.getNrings());
	}
    
    recordParam(theStream, "[globalProfile]", "Saving the global profile?", stringize(par.getGlobProf()));
    recordParam(theStream, "[totalMap]", "Saving 0th moment map to FITS file?", stringize(par.getTotalMap()));
    recordParam(theStream, "[velocityMap]", "Saving 1st moment map to FITS file?", stringize(par.getVelMap()));
    recordParam(theStream, "[dispersionMap]", "Saving 2th moment map to FITS file?", stringize(par.getDispMap()));
    if (par.getTotalMap() || par.getDispMap() || par.getVelMap()) {
		recordParam(theStream, "[blankCube]", "   Using a blanked cube for maps and profile?", stringize(par.getBlankCube()));
		if (par.getBlankCube())	
			if (par.getBmaj()!=-1) recordParam(theStream, "[Bmaj]", "   Beam major axis of smoothed cube for mask", par.getBmaj());
			if (par.getBmin()!=-1) recordParam(theStream, "[Bmin]", "   Beam minor axis of smoothed cube for mask", par.getBmin());
			if (par.getBpa()!=-1) recordParam(theStream,  "[Bpa]", "   Beam position angle of smoothed cube for mask", par.getBpa());
			recordParam(theStream, "[blankCut]", "   SNR clipping cut for blanked areas", par.getBlankCut());
    }
    
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
		recordParam(theStream, "[REDUCE]", "Reducing datacube?", stringize(par.getflagReduce()));
	}
    
    
    if (par.getDistance()!=-1)
		recordParam(theStream, "[Distance]", "Distance of the galaxy (Mpc)?", par.getDistance());

    recordParam(theStream, "[GALFIT]", "Fitting a 3D model to the datacube?", stringize(par.getflagGalFit()));
    recordParam(theStream, "[GALMOD]", "Writing a 3D model to the datacube?", stringize(par.getflagGalMod()));

    if (par.getflagGalFit() || par.getflagGalMod()) {
		std::string box;
		for (int i=0;i<6;i++) if (par.getBOX(i)!=-1) box += to_string<int>(par.getBOX(i))+" ";
		if (box=="") box = "NONE";
		recordParam(theStream, "[BOX]", "   Sub-region to be used?", box);
		recordParam(theStream, "[NRADII]", "   Number of radii", par.getNRADII());
		recordParam(theStream, "[RADSEP]", "   Separation between radii (arcsec)", par.getRADSEP());
		recordParam(theStream, "[XPOS]", "   X center of the galaxy (pixel)", par.getXPOS());
		recordParam(theStream, "[YPOS]", "   Y center of the galaxy (pixel)", par.getYPOS());
		recordParam(theStream, "[VSYS]", "   Systemic velocity of the galaxy (km/s)", par.getVSYS());
		recordParam(theStream, "[VROT]", "   Initial global rotation velocity (km/s)", par.getVROT());
		recordParam(theStream, "[VDISP]", "   Initial global velocity dispersion (km/s)", par.getVDISP());
		recordParam(theStream, "[INC]", "   Initial global inclination (degrees)", par.getINC());
		recordParam(theStream, "[PA]", "   Initial global position angle (degrees)", par.getPHI());
		recordParam(theStream, "[Z0]", "   Scale height of the disk (arcsec)", par.getZ0());
		recordParam(theStream, "[DENS]", "   Global column density of gas (atoms/cm2)", par.getDENS());
		recordParam(theStream, "[FREE]", "   Parameters to be minimized", par.getFREE());
        recordParam(theStream, "[MASK]", "   Type of mask?", par.getMASK());
		recordParam(theStream, "[SIDE]", "   What side of the galaxy to be used", (par.getSIDE()));	
        recordParam(theStream, "[NORM]", "   Type of normalization?", (par.getNORM()));

		
		std::string ltype;
		switch (par.getLTYPE()) {
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
		switch (par.getFTYPE()) {
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
		switch (par.getWFUNC()) {
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
		
		recordParam(theStream, "[TOL]", "   Minimization tolerance", par.getTOL());	
		
		recordParam(theStream, "[TWOSTAGE]", "   Two stages minimization?", stringize(par.getTwoStage()));
		if (par.getTwoStage())
			recordParam(theStream, "[POLYN]", "     Degree of polynomial fitting angles?", par.getPOLYN());
		recordParam(theStream, "[flagErrors]", "   Estimating errors?", stringize(par.getflagErrors()));

	}
	
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

	Str	<< endl << endl
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

