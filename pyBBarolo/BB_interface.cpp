#include<iostream>
#include<Arrays/param.hh>
#include<Arrays/cube.hh>
#include<Tasks/galmod.hh>
#include<Tasks/galfit.hh>
#include<Tasks/galwind.hh>

using namespace Model;
using namespace std;

// Wrapped C functions cannot be stopped in python with CTRL-C.
// Using handler to catch the SIGINT. Use signal(SIGINT, signalHandler)
// in function that have long execution time.
void signalHandler(int signum) {std::cerr << "Killed by the user.\n"; exit(signum);}


extern "C" {

// Interface for the Cube class ///////////////////////////////////////////////////////
Cube<float>* Cube_new(const char* fname) {return new Cube<float>(string(fname));}
void Cube_delete (Cube<float> *c) {delete c;}
int* Cube_axisdim(Cube<float> *c) {return c->AxisDim();};
float* Cube_array(Cube<float> *c) {return c->Array();}
////////////////////////////////////////////////////////////////////////////////////////

// Interface for the Rings struct //////////////////////////////////////////////////////
Rings<float>* Rings_new() {return new Rings<float>;};
void Rings_set(Rings<float>* r, int size, float* radii, float* xpos, float* ypos, float* vsys, float* vrot, float* vdisp, 
               float* vrad, float* vvert, float* dvdz, float* zcyl, float* dens, float* z0, float* inc, float* phi, int* nv)
                   {r->setRings(size,radii,xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens,z0,inc,phi,nv);}
////////////////////////////////////////////////////////////////////////////////////////



// Interface for Galfit class //////////////////////////////////////////////////////////
Galfit<float>* Galfit_new(Cube<float>* c) {return new Galfit<float>(c);}
Galfit<float>* Galfit_new_all(Cube<float> *c, Rings<float> *inrings, float DELTAINC, float DELTAPHI, int LTYPE, 
                              int FTYPE, int WFUNC, int BWEIGHT, int NV, double TOL, int CDENS, int STARTRAD, 
                              const char* MASK, const char* NORM, const char* FREE, const char* SIDE, bool TWOSTAGE, 
                              const char* POLYN, bool ERRORS, bool SMOOTH, float DISTANCE, double REDSHIFT, 
                              double RESTWAVE, const char* OUTFOLD) {
                              return new Galfit<float>(c,inrings,DELTAINC,DELTAPHI,LTYPE,FTYPE,WFUNC,BWEIGHT,NV,TOL,
                              CDENS,STARTRAD,string(MASK),string(NORM),string(FREE),string(SIDE),TWOSTAGE,string(POLYN),
                              ERRORS,SMOOTH,DISTANCE,REDSHIFT,RESTWAVE,string(OUTFOLD));}

void Galfit_delete(Galfit<float> *g) {delete g;}
bool Galfit_galfit(Galfit<float> *g) {signal(SIGINT, signalHandler); g->galfit(); return true;}
bool Galfit_secondStage(Galfit<float> *g) {signal(SIGINT, signalHandler); return g->SecondStage();}
void Galfit_writeModel(Galfit<float> *g, const char* norm) {signal(SIGINT, signalHandler); g->writeModel(string(norm));}

////////////////////////////////////////////////////////////////////////////////////////
 

// Interface for Galfit class //////////////////////////////////////////////////////////
GalWind<float>* Galwind_new (Cube<float> *c, double x0, double y0, double pa, double inc, double disp, 
            double dens, double vsys, double vw, double openang, double htot, 
            int denstype, int ntot, int cdens, int nv) { 
            return new GalWind<float>(c,x0,y0,pa,inc,disp,dens,vsys,vw,openang,htot,denstype,ntot,cdens,nv);}
    
void Galwind_delete(GalWind<float> *gw) {delete gw;} 
float* Galwind_array(GalWind<float> *gw) {return gw->getArray();}
bool Galwind_compute(GalWind<float> *gw) {signal(SIGINT, signalHandler); return gw->compute();}
bool Galwind_smooth(GalWind<float> *gw) {signal(SIGINT, signalHandler); return gw->smooth();}
bool Galwind_writeFITS(GalWind<float> *gw) {return gw->writeFITS();}
bool Galwind_writeMomentMaps(GalWind<float> *gw) {return gw->writeMomentMaps();}
////////////////////////////////////////////////////////////////////////////////////////
}