#include<iostream>
#include<csignal>
#include<Arrays/param.hh>
#include<Arrays/cube.hh>
#include<Arrays/rings.hh>
#include<Tasks/galmod.hh>
#include<Tasks/galfit.hh>
#include<Tasks/galwind.hh>
#include<Tasks/ringmodel.hh>
#include<Tasks/ellprof.hh>
#include<Tasks/smooth3D.hh>
#include<Utilities/paramguess.hh>


using namespace Model;
using namespace std;
using namespace Tasks;

// Wrapped C functions cannot be stopped in python with CTRL-C.
// Using handler to catch the SIGINT. Use signal(SIGINT, signalHandler)
// in functions that have long execution time.
void signalHandler(int signum) {std::cerr << "Killed by the user.\n"; exit(signum);}

 
extern "C" {

void delete_array(int *a) {delete [] a;}

// Interface for the Param class //////////////////////////////////////////////////////
Param* Param_new() {return new Param;}
void Param_setfromfile(Param *p, const char* pfile) {p->readParamFile(string(pfile));}
void Param_setfromstr(Param *p, const char* pstr) {p->readParamString(string(pstr));}
void Param_delete (Param *p) {delete p;}
////////////////////////////////////////////////////////////////////////////////////////


// Interface for the Cube class ///////////////////////////////////////////////////////
Cube<float>* Cube_new(const char* fname) {return new Cube<float>(string(fname),false);}
void Cube_delete (Cube<float> *c) {delete c;}
int* Cube_axisdim(Cube<float> *c) {return c->AxisDim();}
float* Cube_array(Cube<float> *c) {return c->Array();}
void Cube_setBeam(Cube<float> *c, float bmaj, float bmin, float bpa) {c->setBeam(bmaj,bmin,bpa);}
float* Cube_getBeam(Cube<float> *c) {return c->getBeam();}
bool* Cube_getMask(Cube<float> *c) {return c->Mask();}
////////////////////////////////////////////////////////////////////////////////////////


// Interface for the Rings struct //////////////////////////////////////////////////////
Rings<float>* Rings_new() {return new Rings<float>;}
void Rings_set(Rings<float>* r, int size, float* radii, float* xpos, float* ypos, float* vsys, float* vrot, float* vdisp, 
               float* vrad, float* vvert, float* dvdz, float* zcyl, float* dens, float* z0, float* inc, float* phi)
                   {r->setRings(size,radii,xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens,z0,inc,phi);}
void Rings_delete(Rings<float>* r) {delete r;}
////////////////////////////////////////////////////////////////////////////////////////


// Interface for Galmod class //////////////////////////////////////////////////////////
Galmod<float>* Galmod_new(Cube<float> *c, Rings<float> *r, int NV, int LTYPE, int CMODE, float CDENS, int ISEED) 
                          {Galmod<float> *g = new Galmod<float>; g->input(c,r,NV,LTYPE,CMODE,CDENS,ISEED); return g;}
Galmod<float>* Galmod_new_par(Cube<float> *c, Rings<float> *r, Param *p) {Galmod<float> *g = new Galmod<float>;
                              c->pars() = *p; c->Head().setRedshift(p->getRedshift()); c->Head().setWave0(p->getRestwave());
                              g->input(c,r,p->getParGM().NV,p->getParGM().LTYPE,p->getParGM().CMODE,
                              p->getParGM().CDENS,p->getParGM().ISEED); return g;}
void Galmod_delete(Galmod<float> *g) {delete g;}
float* Galmod_array(Galmod<float> *g) {return g->getArray();}
void Galmod_set_array(Galmod<float> *g, float *a) {return g->setArray(a);}
bool Galmod_compute(Galmod<float> *g) {signal(SIGINT, signalHandler); return g->calculate();}
bool Galmod_smooth(Galmod<float> *g) {signal(SIGINT, signalHandler); return g->smooth();}
////////////////////////////////////////////////////////////////////////////////////////


// Interface for Galfit class //////////////////////////////////////////////////////////
Galfit<float>* Galfit_new(Cube<float>* c) {return new Galfit<float>(c);}
Galfit<float>* Galfit_new_par(Cube<float> *c, Rings<float> *inrings, Param *p) {
                              return new Galfit<float>(c,inrings,p);}
void Galfit_delete(Galfit<float> *g) {delete g;}
float* Galfit_initialGuesses(Cube<float> *c, const char* xpos, const char* ypos, const char* inc, const char* pa) {
                             GALFIT_PAR p; p.XPOS=string(xpos); p.YPOS=string(ypos); p.INC=string(inc); p.PHI=string(pa);
                             ParamGuess<float> *ip = EstimateInitial(c,&p); float *r =  new float[8]; r[0]=ip->nrings; 
                             r[1]=ip->radsep; r[2]=ip->xcentre; r[3]=ip->ycentre; r[4]=ip->vsystem; r[5]=ip->vrot; 
                             r[6]=ip->inclin; r[7]=ip->posang; delete ip; return r;}
bool Galfit_galfit(Galfit<float> *g) {signal(SIGINT, signalHandler); g->galfit(); return true;}
bool Galfit_secondStage(Galfit<float> *g) {signal(SIGINT, signalHandler); return g->SecondStage();}
float Galfit_calcresiduals(Galfit<float> *g, Rings<float> *r) {return g->calculateResiduals(r);}
void Galfit_getModelSize(Galfit<float> *g, Rings<float> *r, int *bhi, int *blo) {g->getModelSize(r,blo,bhi);}
Galmod<float>* Galfit_getModel(Galfit<float> *g, Rings<float> *r, int *bhi, int *blo) {signal(SIGINT, signalHandler); 
                                                                        return g->getModel(r,bhi,blo,nullptr,false);}
float* Galfit_getModel_BBB(Galfit<float> *g, Rings<float> *r, int *bhi, int *blo, int iseed) {signal(SIGINT, signalHandler); 
                                                                        return g->getModel_BBB(r,bhi,blo,iseed);}
void Galfit_writeModel(Galfit<float> *g, const char* norm, bool plots) {signal(SIGINT, signalHandler); g->writeModel(string(norm),plots);}
void Galfit_writeOutputs(Galfit<float> *g, Galmod<float> *m, Ellprof<float> *e, bool plots) {signal(SIGINT, signalHandler); g->writeOutputs(m->Out(),e,plots);}
void Galfit_setOutRings(Galfit<float> *g, Rings<float> *r) {g->setOutRings(r); g->writeRingFile("rings_final1.txt",r);}
int Galfit_plotModel(Galfit<float> *g) {signal(SIGINT, signalHandler); return g->plotAll_Python();}
////////////////////////////////////////////////////////////////////////////////////////
 

// Interface for Galwind class //////////////////////////////////////////////////////////
GalWind<float>* Galwind_new (Cube<float> *c, float x0, float y0, float pa, float inc, float disp, 
                             float dens, float vsys, float vw, float openang, float htot, 
                             int denstype, int ntot, int cdens, int nv, int NTHREADS) { 
            return new GalWind<float>(c,x0,y0,pa,inc,disp,dens,vsys,vw,openang,htot,denstype,ntot,cdens,nv, NTHREADS);}
    
void Galwind_delete(GalWind<float> *gw) {delete gw;} 
float* Galwind_array(GalWind<float> *gw) {return gw->getArray();}
bool Galwind_compute(GalWind<float> *gw) {signal(SIGINT, signalHandler); return gw->compute();}
bool Galwind_smooth(GalWind<float> *gw) {signal(SIGINT, signalHandler); return gw->smooth();}
bool Galwind_writeFITS(GalWind<float> *gw) {return gw->writeFITS();}
bool Galwind_writeMomentMaps(GalWind<float> *gw) {return gw->writeMomentMaps();}
//////////////////////////////////////////////////////////////////////////////////////////

// Interface for Search class //////////////////////////////////////////////////////////
void Search_search(Cube<float> *c, const char* searchtype, float snrCut, float threshold, bool adjacent, 
                   int threshSpatial, int threshVelocity, int minPixels, int minChannels,
                   int minVoxels, int maxChannels, float maxAngSize, bool flagGrowth,
                   float growthCut, float growthThreshold, bool RejectBefore, bool TwoStage,int NTHREADS) 
                   {signal(SIGINT, signalHandler); c->search(string(searchtype),snrCut,threshold,
                    adjacent,threshSpatial,threshVelocity,minPixels,minChannels,minVoxels,maxChannels,
                    maxAngSize,flagGrowth,growthCut,growthThreshold,RejectBefore,TwoStage,NTHREADS);}
//////////////////////////////////////////////////////////////////////////////////////////
                    

// Interface for 2DFIT class //////////////////////////////////////////////////////////
Ringmodel<float>* Fit2D_new(Cube<float> *c, Rings<float> *r, const char* mask, const char* free, const char* side, int wfunc, int NTHREADS) {
                    c->pars().setMASK(string(mask)); c->pars().getParGF().FREE = string(free); 
                    c->pars().getParGF().SIDE = string(side); c->pars().getParGF().WFUNC= wfunc; c->pars().setThreads(NTHREADS);
                     Ringmodel<float> *rm = new Ringmodel<float>(); rm->setfromCube(c,r); return rm;}
void Fit2D_delete(Ringmodel<float> *rm) {delete rm;}
void Fit2D_compute(Ringmodel<float> *rm) {signal(SIGINT, signalHandler); rm->ringfit();}
void Fit2D_write(Ringmodel<float> *rm, Cube<float> *c, const char *fout) {std::ofstream fileo(fout); rm->printfinal(fileo, c->Head());}
//////////////////////////////////////////////////////////////////////////////////////////


// Interface for ELLPROF class //////////////////////////////////////////////////////////
Ellprof<float>* Ellprof_new(Cube<float> *c, Rings<float> *r, const char* mask, const char* side, int NTHREADS) {
                            c->pars().setMASK(string(mask)); c->pars().getParGF().SIDE = string(side);
                            c->pars().setThreads(NTHREADS); return new Ellprof<float>(c,r);}
Ellprof<float>* Ellprof_new_alt(Cube<float> *c, Rings<float> *r) {return new Ellprof<float>(c,r);} 
void Ellprof_delete(Ellprof<float> *e) {delete e;}
void Ellprof_compute(Ellprof<float> *e) {signal(SIGINT, signalHandler); e->RadialProfile();}
void Ellprof_write(Ellprof<float> *e, const char *fout) {std::ofstream fileo(fout); e->printProfile(fileo);}
double* Ellprof_dens_array(Ellprof<float> *e) { double *d = new double[e->getNrad()]; 
                                                for (int i=e->getNrad(); i--;) d[i] = e->getMedian(i); return d;}
void Ellprof_update_rings(Ellprof<float> *e, Rings<float> *r) {e->update_rings(r);}
//////////////////////////////////////////////////////////////////////////////////////////                  


// Interface for SpectralSmooth3D class /////////////////////////////////////////////////
SpectralSmooth3D<float>* SpectralSmooth3D_new(const char *wtype, size_t wsize) {return new SpectralSmooth3D<float>(string(wtype),wsize);}
void SpectralSmooth3D_delete(SpectralSmooth3D<float> *h) {delete h;}
void SpectralSmooth3D_compute(SpectralSmooth3D<float> *h, Cube<float> *c, int NTHREADS) {
                              signal(SIGINT, signalHandler); c->pars().setThreads(NTHREADS); h->smooth(c);}
float* SpectralSmooth3D_array(SpectralSmooth3D<float> *h) {return h->Array();}
void SpectralSmooth3D_write(SpectralSmooth3D<float> *h, Cube<float> *c, const char *fout, bool average) {
                            c->pars().setflagReduce(average); h->fitswrite(c, string(fout));}
/////////////////////////////////////////////////////////////////////////////////////////
}
