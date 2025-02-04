/////////////////////////////////////////////
// bbarolo.cpp: Main source file of BBarolo
/////////////////////////////////////////////

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
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>
#include <Arrays/param.hh>
#include <Arrays/cube.hh>
#include <Arrays/stats.hh>
#include <Map/detection.hh>
#include <Tasks/moment.hh>
#include <Tasks/ringmodel.hh>
#include <Tasks/smooth3D.hh>
#include <Tasks/galfit.hh>
#include <Tasks/galwind.hh>
#include <Tasks/ellprof.hh>
#include <Tasks/spacepar.hh>
#include <Tasks/rendering3D.hh>
#include <Utilities/utils.hh>
#include <Utilities/paramguess.hh>

#ifdef MPI
#include <mpi.h>
#include <csignal>
#endif

#ifdef DOUBLE_PRECISION
using BBreal = double;
#else
using BBreal = float;
#endif


bool BBcore (Param *par);
bool BBauto (Cube<BBreal> *c);
void sigsegv_handler(int signal);
int BBarolo_MPI (int argc, char *argv[]);

/*
#include<signal.h>
struct sigaction osa;
void action_sigint(int sig_no)
{
    printf("\nI tap SIGINT and returns back \n");
    sigaction(SIGINT,&osa,NULL);
    kill(0,SIGINT);
}
*/

int main (int argc, char *argv[]) {

//    struct sigaction act;
//    act.sa_handler = action_sigint;
//    sigemptyset(&act.sa_mask);
//    act.sa_flags = 0;
//    sigaction(SIGINT, &act, 0);

#ifdef MPI
    return BBarolo_MPI(argc,argv);

#else
    struct timeval begin, end;
    gettimeofday(&begin, NULL);

    // Reading in parameters.
    Param *par = new Param;
    if (!par->getopts(argc, argv)) return EXIT_FAILURE;

    bool verbose = par->isVerbose();
    if (verbose && !par->getLogFile()) {
        welcomeMessage();
        std::cout << *par;
    }

    for (int im=0; im<par->getListSize(); im++) {

        if (par->getListSize()>1 && verbose && !par->getLogFile()) {
            std::cout << setfill('_') << std::endl;
            std::cout << setw(70) << "" << std::endl << std::endl;
            std::string s = "Working on "+ par->getImage(im)+" ";
            std::cout << setfill(' ') << right << setw(70) << s;
            std::cout << std::endl << left;
            std::cout << std::endl << " File "<< im+1 << " of " 
                      << par->getListSize()<< std::endl << std::endl;
        }

        par->setImageFile(par->getImage(im));
        
        if (!BBcore(par)) {
            if(par->getListSize()-im>1) std::cout << "Skipping to next file...\n";
            else {std::cout << "Exiting ...\n\n"; return EXIT_FAILURE;}
        }
    }

    delete par;
    
    gettimeofday(&end, NULL);
    double time = (end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0);
    if (verbose) {
        if (time<1) std::cout << "\nExecution time: " << setprecision(2) << time << " sec.\n";
        else std::cout << "\nExecution time: " << int(time/60) << " min and " << int(time)%60 << " sec.\n";
    }

    return EXIT_SUCCESS;
#endif
}


bool BBcore (Param *par) {

    // Interface for all BBarolo's tasks

    Cube<BBreal> *c = new Cube<BBreal>;
    c->saveParam(*par);

    // Determining the output directory
    std::string outfolder = c->pars().getOutfolder();
    if (outfolder=="") {
       std::string path = get_currentpath();
       Header h;
       h.setWarning(false);
       
       if (!h.header_read(par->getImageFile())) {
           std::cout << par->getImageFile() << " is not a readable FITS file!\n";
           delete c;
           return false;
       }
       outfolder = path+"/output/"+h.Obname()+"/";
       c->pars().setOutfolder(outfolder);
    }
    mkdirp(outfolder.c_str());

    // Redirecting output to a log file if requested
    std::ofstream outf;
    if (par->getLogFile()) {
        std::string logfile = outfolder+"BBlog_0.txt";
        int i = 0;
        while (fexists(logfile)) 
            logfile = outfolder+"BBlog_"+to_string(++i)+".txt";
        outf.open(logfile);
        std::cout.rdbuf(outf.rdbuf());
        std::cerr.rdbuf(outf.rdbuf());
        if (par->isVerbose()) {
            welcomeMessage();
            std::cout << *par;
        }
        c->pars().setShowbar(false);
    }

    // Reading in FITS file
    if (!c->readCube(par->getImageFile())) {
        std::cout << par->getImageFile() << " is not a readable FITS file!\n";
        delete c;
        return false;
    }

    // Fully automated run
    if (par->getFlagAuto()) {
        //return BBauto(c);
        par->getParGF().flagGALFIT = true;
    }

    // Continuum subtraction
    if (par->getFlatContsub()) {
        c->continuumSubtract();
        c->fitswrite_3d((outfolder+c->Head().Name()+"_contsub.fits").c_str(),true);
    }


    if (par->getCheckCube()) c->CheckCube(par->getCheckCube());


    // Statistic utility ------------------------------------------
    if (par->getFlagStats()) {
        c->setCubeStats();
        c->stat().tofile(c->pars().getOutfolder()+c->Head().Name()+"_stats.txt");
    }
    // --------------------------------------------------------------
    
    
    // Mask making utility ----------------------------------------
    if (par->getMakeMask()) {
        if (par->getMASK().find("LARGEST")!=std::string::npos) c->BlankMask(NULL,true);
        else c->BlankMask(NULL,false);
    }
    // --------------------------------------------------------------


    // Spatial smoothing utility ------------------------------------------
    if (par->getflagSmooth()) {
        Smooth3D<BBreal> *sm = new Smooth3D<BBreal>;
        sm->cubesmooth(c);
        sm->fitswrite();
        delete sm;
    }
    // --------------------------------------------------------------


    // Spectral smoothing utility -----------------------------------
    if (par->getflagSmoothSpectral()) {
        SpectralSmooth3D<BBreal> *sm = new SpectralSmooth3D<BBreal>(par->getWindowType(),par->getWindowSize());
        sm->smooth(c);
        sm->fitswrite(c);
        delete sm;
    }
    // --------------------------------------------------------------


    // Source finding utility --------------------------------------
    if (par->getflagSearch()) {
        c->search();
        std::ofstream detout((outfolder+"detections.txt").c_str());
        c->printDetections(detout);
        c->writeDetections();
        if (par->getParSE().cubelets) {
            c->writeCubelets();
            c->plotDetections();
        }
    }
    // --------------------------------------------------------------


    // 3D Cube Fitting task -----------------------------------------
    if (par->getflagGalFit()) {
        Model::Galfit<BBreal> *fit = new Model::Galfit<BBreal>(c);
        fit->galfit();
        if (par->getParGF().TWOSTAGE) fit->SecondStage();
        if (par->getFlagDebug()) fit->writeModel("BOTH",par->getFlagPlots());
        else fit->writeModel(fit->Norms(2),par->getFlagPlots());
        delete fit;
    }
    // --------------------------------------------------------------


    // Cube Model task -----------------------------------------------
    if (par->getflagGalMod()) {
        Model::Galfit<BBreal> *fit = new Model::Galfit<BBreal>(c);
        if (par->getFlagDebug()) fit->writeModel("BOTH",false);
        else fit->writeModel(par->getParGF().NORM,false);
        delete fit;
    }
    //----------------------------------------------------------------


    // Full parameter space task --------------------------------------
    if (par->getflagSpace()) {
        Spacepar<BBreal> *sp = new Spacepar<BBreal>(c);
        sp->calculate();
        sp->plotAll_Python();
        delete sp;
    }
    //----------------------------------------------------------------


    // GalWind task --------------------------------------------------
    if (par->getParGW().flagGALWIND) {
        GalWind<BBreal> *w = new GalWind<BBreal>(c);
        w->compute();
        if (par->getParGF().SM) w->smooth();
        w->writeFITS();
        w->writeMomentMaps();
        w->writePV();
        if (par->getFlagPlots()) w->makePlots();
        delete w;
    }
    //----------------------------------------------------------------


    // 2D tilted-ring fitting task -----------------------------------
    if (par->getFlagRing()) {
        Ringmodel<BBreal> *trmod = new Ringmodel<BBreal>(c);
        trmod->ringfit(c->pars().getThreads(),c->pars().isVerbose(),c->pars().getShowbar());
        std::string fout = c->pars().getOutfolder()+c->Head().Name()+"_2dtrm.txt";
        std::ofstream fileo(fout.c_str());
        trmod->printfinal(fileo,c->Head());
        fileo.close();
        trmod->printfinal(std::cout,c->Head());
        trmod->writeModel(c->pars().getOutfolder()+c->Head().Name()+"_2d_mod.fits",c->Head());
        delete trmod;
    }
    //-----------------------------------------------------------------


    // Moment maps task -----------------------------------------------
    if (par->getMaps()) {
        std::string s = outfolder+c->Head().Name();
        bool masking = par->getMASK()=="NONE" ? false : true;
        MomentMap<BBreal> map;
        map.input(c);
        if (par->getParMA().totalmap) {
            map.ZeroMoment(masking,par->getParMA().maptype);
            map.fitswrite_2d((s+"map_0th.fits").c_str());
            if (par->getParMA().SNmap) map.SNMap(masking);
        }
        if (par->getParMA().velocitymap) {
            map.FirstMoment(masking,par->getParMA().maptype);
            map.fitswrite_2d((s+"map_1st.fits").c_str());
        }
        if (par->getParMA().dispersionmap) {
            map.SecondMoment(masking,par->getParMA().maptype);
            map.fitswrite_2d((s+"map_2nd.fits").c_str());
        }
        if (par->getParMA().massdensmap) {
            map.HIMassDensityMap(masking);
            map.fitswrite_2d((s+"map_massdens.fits").c_str());
        }
        if (par->getParMA().rmsmap) {
            map.RMSMap(masking);
            map.fitswrite_2d((s+"map_RMS.fits").c_str());
        }
        if (par->getParMA().globprof) {
            Image2D<BBreal> spectrum;
            spectrum.extractGlobalSpectrum(c);
            std::string specfname = c->pars().getOutfolder()+"spectrum.txt";
            std::ofstream specfile; specfile.open(specfname.c_str());
            for (int i=0; i<c->DimZ(); i++)
                specfile << c->getZphys(i) << " " << spectrum.Array(i) << std::endl;
        }
    }
    //-----------------------------------------------------------------


    // PVs extraction task --------------------------------------------
    if (par->getFlagPV()) {
        std::string s = outfolder+c->Head().Name();
        PvSlice<BBreal> *pv = new PvSlice<BBreal>(c);
        pv->slice();
        //pv->slice_old();
        pv->fitswrite_2d((s+"pv_"+to_string(par->getPA_PV(),0)+".fits").c_str(),true);
        delete pv;
    }
    
    //-----------------------------------------------------------------


    // Repixeling task ------------------------------------------------
    if (par->getflagReduce() && !(par->getflagSmooth() || par->getflagSmoothSpectral())) {
        std::string name = c->pars().getOutfolder()+c->Head().Name()+"_red.fits";
        std::string redtype = par->getReduce()=="spectral" ? "spectral" : "spatial";
        Cube<BBreal> *red = c->Reduce(floor(c->pars().getFactor()),redtype);
        red->fitswrite_3d(name.c_str(),true);
        delete red;
    }
    //-----------------------------------------------------------------


    // Kinematics fitting to slit data --------------------------------
    if (par->getFlagSlitfit()) {
        Model::Galfit<BBreal> *sfit = new Model::Galfit<BBreal>;
        sfit->slit_init(c);
        sfit->galfit();
        if (par->getParGF().TWOSTAGE) sfit->SecondStage();
        sfit->writeModel_slit();
    }
    //-----------------------------------------------------------------


    // Radial profile ------------------------------------------------
    if (par->getFlagEllProf()) {
        Tasks::Ellprof<BBreal> *ell = new Tasks::Ellprof<BBreal>(c);
        ell->RadialProfile();
        std::string fout = c->pars().getOutfolder()+c->Head().Name()+"_densprof.txt";
        std::ofstream fileo(fout.c_str());
        ell->printProfile(fileo,ell->getNseg()-1);
        ell->printProfile(std::cout,ell->getNseg()-1);
        ell->writeMap(c->pars().getOutfolder()+c->Head().Name()+"_densmap.fits");
        fileo.close();
        delete ell;
    }
    //-----------------------------------------------------------------


    // 3D Rendering ---------------------------------------------------
    if (par->getFlagRend3D()) {
        Rendering3D<BBreal> *r3d = new Rendering3D<BBreal>(c);
        r3d->compute(par->getRendAngle());
        r3d->writefits();
        delete r3d;
    }
    //-----------------------------------------------------------------

    delete c;

    outf.close();

    return true;
}

bool BBauto (Cube<BBreal> *c) {

    Param &p = c->pars();
    // First thing, searching the cube for sources
    c->search();

    int starts[3];

    for (int i=0; i<c->getNumObj(); i++) {
        Detection<BBreal> *obj = c->pObject(i);
        obj->calcWCSparams(c->Head());

        std::cout << p.getParSE().edges << std::endl;

        Cube<BBreal> *cl = c->extractCubelet(obj,p.getParSE().edges,starts);
        cl->saveStats(c->stat());
        //cl->search();

        Detection<BBreal> newdec(*obj);
        newdec.addOffsets(starts[0],starts[1],starts[2]);
//        cl->saveSources(*se);
        //c->writeDetections();
        //ParamGuess<BBreal> pg(c,&newdec);
        //pg.findAll();
       //pg.findCentre();
        std::cout << c->getZphys(0) << std::endl;
        std::cout << cl->getZphys(0) << std::endl;
        ParamGuess<BBreal> pg2(c,obj);
        pg2.findAll();
//        pg2.findCentre();
        std::cout << starts[0] << " " << starts[1] << " " << starts[2] << std::endl;
       // std::cout << pg.xcentre << " " << pg.ycentre << " " << pg.vsystem << std::endl;
        std::cout << pg2.xcentre << " " << pg2.ycentre << " " << pg2.vsystem << " " << pg2.posang << " " << pg2.inclin << std::endl;
        //pg.plotGuess("pippo1.pdf");
        //std::cin.ignore();
        pg2.plotGuess("pippo3.pdf");

        pg2.tuneWithTiltedRing();
        pg2.plotGuess("pippo5.pdf");
        c->writeCubelets();


        // Now fitting a 3dmodel
        //Model::Galfit<BBreal> *fit = new Model::Galfit<BBreal>(cl);
        //fit->galfit();
        //if (p.getParGF().TWOSTAGE) fit->SecondStage();
        //if (p.getFlagDebug()) fit->writeModel("BOTH",p.getFlagPlots());
        //else fit->writeModel(p.getParGF().NORM,p.getFlagPlots());
        //delete fit;

    }

    return true;
}


#ifdef MPI

void sigsegv_handler(int signal) {
    // this signal handler is not POSIX-compliant, but apparently signals
    // in MPI are a mess anyway, so I guess don't worry too much...
    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cerr << "Segmentation fault (SIGSEGV) in process rank " << rank << "." << std::endl;
    std::cerr << "Process terminated (remaining processes will continue)." << std::endl;
    std::cerr << "Exited with segmentation fault (SIGSEGV) error!" << std::endl;
    MPI_Finalize();
    _Exit(0);
}

int BBarolo_MPI (int argc, char *argv[]) {
  
    std::signal(SIGSEGV, sigsegv_handler);
  
    if (argc!=3 || std::string(argv[1])!="-l") {
        std::cerr << "Usage: BBarolo_MPI -l listfile \n";
        return EXIT_FAILURE;
    }
    
    int rank=0, nprocs=1;
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    std::string listfile = std::string(argv[2]);
    std::vector<std::string> parfiles;
  
    // Reading parameter list file
    std::ifstream file(listfile.c_str());
    if(!file) std::cerr << "File " << listfile << " does not exist! \n";
    std::string s;
    while(file.good()) {
        getline(file, s);
        if (s!="") parfiles.push_back(s);
    }
    file.close();
  
    // Dividing parameter files among MPI processes
    int lsize = parfiles.size()/nprocs;
    int res_lsize = parfiles.size()%nprocs;
  
    int start[nprocs+1];
    for (int i=0; i<nprocs; i++) {
        if (i<res_lsize) start[i]=(lsize+1)*i;
        else start[i]=i*lsize+res_lsize;
    }
    start[nprocs]=parfiles.size();
  
    // Redirecting std::cout to files
    std::filebuf buf;
    std::string pout = "pout."+to_string(rank);
    buf.open(pout, std::ios::out );
    auto oldbuf = std::cout.rdbuf(&buf);
  
    // Main loop over parameter files
    for (auto i=start[rank]; i<start[rank+1]; i++) {
    
        Param *par = new Param;
        par->readParamFile(parfiles[i]);
        std::cout << *par;
        if (!BBcore(par)) {
            if(par->getListSize()-i>1) std::cout << "Skipping to next file...\n";
            else {std::cout << "Exiting ...\n\n"; return EXIT_FAILURE;}
            delete par;
            continue;
        }
    
        delete par;
    }

    // Back to std::cout
    std::cout.rdbuf(oldbuf);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

#endif
