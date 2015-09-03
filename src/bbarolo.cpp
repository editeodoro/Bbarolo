/////////////////////////////////////////////
// bbarolo.cpp: Main source file of BBarolo
/////////////////////////////////////////////

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
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <sys/stat.h>
#include "./Arrays/param.hh"
#include "./Arrays/cube.hh"
#include "./Arrays/stats.hh"
#include "./Utilities/moment.hh"
#include "./Utilities/ringmodel.hh"
#include "./Utilities/smooth3D.hh"
#include "./Utilities/galfit.hh"
#include "./Utilities/spacepar.hh"
#include "./Utilities/utils.hh"

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

	clock_t start = clock();
    Param *par = new Param;

	if (!par->getopts(argc, argv)) return EXIT_FAILURE;
	if (par->getImageList()=="NONE") par->setImage(par->getImageFile());
	std::cout << *par;

	mkdir ("./output/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
	for (int im=0; im<par->getListSize(); im++) {
		
		if (par->getListSize()>1) {
			std::cout << setfill('_') << std::endl;
			std::cout << setw(70) << "" << std::endl << std::endl;
			std::string s = "Working on "+ par->getImage(im)+" ";
			std::cout << setfill(' ') << right << setw(70) << s;
			std::cout << std::endl << left;
			std::cout << std::endl << " File "<< im+1
					  << " of " << par->getListSize()<<std::endl<<std::endl;
		}
	
		par->setImageFile(par->getImage(im));
		std::string fname = par->getImage(im);
		int found = fname.find("[");
		if (found>=0) fname.erase(found, fname.size()-found); 
		if (!fexists(fname)) {
			std::cout << "\nError reading " << par->getImage(im) 
				  	  << " : the file doesn't exist!\n";
			if(par->getListSize()-im>1) std::cout << "Skipping to next file...\n";
            else {std::cout << "Exiting ...\n\n"; return EXIT_FAILURE;}
			continue;
		}



        Cube<float> *c = new Cube<float>;
		c->saveParam(*par);
		
		if (!c->readCube(par->getImageFile())) {
			std::cout << par->getImageFile() << " is not a readable FITS image!\n";
			if(par->getListSize()-im>1) std::cout << "Skipping to next file...\n";
			else std::cout << "Exiting ...\n\n";
			delete c;
			continue;
        }

		std::string outfolder = "./output/"+c->Head().Obname()+"/";
		c->pars().setOutfolder(outfolder);
		mkdir (outfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
		if (par->getCheckCh()) c->CheckChannels();
	
		if (par->getflagSmooth()) {
            Smooth3D<float> *sm = new Smooth3D<float>;
			sm->cubesmooth(c);
			sm->fitswrite();
			delete sm;
		}
	
		///<<<<< Searching stuff
		if (par->getSearch()) {			
			c->Search();
			c->plotDetections();
            std::ofstream detout((outfolder+"detections.txt").c_str());
            c->printDetections(detout);

        }
	
		///<<<<< Cube Fitting
		if (par->getflagGalFit()) {
            Model::Galfit<float> *fit = new Model::Galfit<float>(c);
            fit->galfit();
            if (par->getNORM()=="AZIM") fit->writeModel_azim();
            if (par->getTwoStage()) {
                fit->SecondStage();
                fit->writeModel_norm();
            } else if (!par->getTwoStage() && par->getNORM()!="AZIM") fit->writeModel_norm();
			delete fit;
		}
		///----------------------
	
		
		///<<<<< Cube Model
		if (par->getflagGalMod()) {
            Model::Galfit<float> *fit = new Model::Galfit<float>(c);
            if (par->getNORM()=="LOCAL") fit->writeModel_norm();
            if (par->getNORM()=="AZIM") fit->writeModel_azim() ;
			else {

                Model::Galmod<float> *mod = fit->getModel();
///////////////
//                long axis[3] = {mod->Out()->DimX(),mod->Out()->DimY(),mod->Out()->DimZ()};
//                float *ar = new float[c->NumPix()];
//                for (int i=0; i<c->NumPix();i++) ar[i] = mod->Out()->Array(i);
//                FitsWrite_3D("diocaro.fits",ar,axis);
/////////////////


				std::string mfile = c->pars().getOutfolder()+c->Head().Name()+"mod.fits";
				mod->Out()->Head().setMinMax(0.,0.);
				mod->Out()->Head().setName(c->Head().Name()+"mod");
				mod->Out()->fitswrite_3d(mfile.c_str());
                fit->plotParam(mod->Out(), mod->Out());
				delete mod;
			}
			delete fit;
		}
		///----------------------

		if (par->getflagSpace()) {
            Spacepar<float> *sp = new Spacepar<float>(c);
			sp->calculate();
			delete sp;
        }

        //<<<<<<< 2D Tilted-Ring Model
        if (par->getFlagRing()) {
            Ringmodel *trmod = new Ringmodel(c);
            trmod->ringfit();
            trmod->printfinal(std::cout);
            delete trmod;
        }
        ///----------------------

		if (par->getMaps()) {
            MomentMap<float> map;
			map.input(c);
			std::string s = outfolder+c->Head().Name();
			if (par->getTotalMap()) {
				map.ZeroMoment(par->getBlankCube());
				map.fitswrite_2d((s+"map_0th.fits").c_str());
			}
			if (par->getVelMap()) {
				map.FirstMoment(par->getBlankCube());
				map.fitswrite_2d((s+"map_1st.fits").c_str());
			}
			if (par->getDispMap()) {
				map.SecondMoment(par->getBlankCube());
				map.fitswrite_2d((s+"map_2nd.fits").c_str());
			}
		}	
		
		if (par->getGlobProf()) {
            Image2D<float> spectrum;
			spectrum.extractGlobalSpectrum(c);
			std::string specfname = c->pars().getOutfolder()+"spectrum.txt";
			std::ofstream specfile; specfile.open(specfname.c_str());
			for (int i=0; i<c->DimZ(); i++) 
				specfile << c->getZphys(i) << " " << spectrum.Array(i) << std::endl;
		}
		
		if (par->getflagReduce() && !par->getflagSmooth()) {
			std::string name = c->pars().getOutfolder()+c->Head().Name()+"_red.fits";
            Cube<float> *red = c->Reduce(floor(c->pars().getFactor()));
			red->fitswrite_3d(name.c_str(),true);
			delete red;
		}
		
		
        // Giusto per prendere i PV lungo angoli
        if (par->getFlagPV()) {
            float xpos = par->getXPOS_PV();
            float ypos = par->getYPOS_PV();
            float ang = par->getPA_PV();
            std::string s = outfolder+c->Head().Name();
            Image2D<float> *pv = PositionVelocity(c,xpos,ypos,ang);
            pv->fitswrite_2d((s+"pv_"+to_string(ang,0)+".fits").c_str());
            delete pv;
        }
		


        if (par->getFlagSlitfit()) {
            Model::Galfit<float> *sfit = new Model::Galfit<float>;
            sfit->slit_init(c);
            sfit->galfit();
            if (par->getTwoStage()) {
                sfit->SecondStage();
                sfit->writeModel_slit();
            }
            else sfit->writeModel_slit();
        }


		delete c;
	
	}
	
	delete par;
	
	clock_t end = clock();
	double time=(double(end-start))/CLOCKS_PER_SEC;
	std::cout << "\nExecution time: " << int(time/60) 
	   	  	  << " min and " << int(time)%60 << " sec.\n";
	
	
	return EXIT_SUCCESS;
	
}
