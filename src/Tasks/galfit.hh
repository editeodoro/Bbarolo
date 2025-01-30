//--------------------------------------------------------------------
// galfit.hh: Definition of the Galfit class.
//--------------------------------------------------------------------

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

#ifndef GALFIT_HH_
#define GALFIT_HH_

#include <iostream>
#include <Arrays/param.hh>
#include <Arrays/cube.hh>
#include <Arrays/image.hh>
#include <Arrays/rings.hh>
#include <Tasks/galmod.hh>
#include <Tasks/ellprof.hh>
#include <Utilities/paramguess.hh>

namespace Model {

enum ALLPARS {VROT, VDISP, DENS, Z0, INC, 
              PA, XPOS, YPOS, VSYS, VRAD, MAXPAR};

template <class T>
class Galfit
{
public:
    // Constructors:
    // 1) Empty constructor
    // 2) Give an object Cube with parameters recorded in c->pars()
    // 3) Give an object Cube + Rings + a GALWIND_PAR structure
    // 4) Give all paramaters separately
    Galfit() {}
    Galfit(Cube<T> *c);
    Galfit(Cube<T> *c, Rings<T> *inrings, GALFIT_PAR *p) {setup(c,inrings,p);}
    Galfit(Cube<T> *c, Rings<T> *inrings, Param *p);
    virtual ~Galfit();
    Galfit(const Galfit &g) {operator=(g);}     //< Copy constructor.
    Galfit& operator=(const Galfit &g);         //< Copy operator.
    
    Cube<T>  *In () {return in;}
    Rings<T> *Inrings  () {return inr;}
    Rings<T> *Outrings () {return outr;}
    void setOutRings (Rings<T>* r) {*outr = *r;}

    /// Functions defined in galfit.cpp
    void setup (Cube<T> *c, Rings<T> *inrings, GALFIT_PAR *p);
    void galfit();
    bool SecondStage();
    double calculateResiduals(Rings<T> *r) {return getFuncValue(r);}
    void writeRingFile(std::string filename, Rings<T> *r, T ***errors=nullptr);
    Galmod<T>* getModel(Rings<T> *dr, int *bhi, int* blo, Model::Galmod<T> *modsoFar=nullptr, bool finalModel=false);
    bool AsymmetricDrift(T *rad, T *densprof, T *dispprof, T *rotcur, T *inc, int nn);

    // Functions defined in galfit_min.cpp
    void   getModelSize(Rings<T> *dring, int *blo, int *bhi);

    /// Functions defined in galfit_out.cpp
    void writeModel(std::string normtype, bool makeplots=true);
    void writeOutputs (Cube<T> *mod, Tasks::Ellprof<T> *e, bool makeplots=true);
    void writePVs(Cube<T> *mod, std::string suffix="");
    void writeKinematicMaps(Cube<T> *mod, std::string suffix="");
    std::vector<std::string> writeScripts_Python();
    int  plotAll_Python ();

    /// Functions defined in slitfit.cpp
    void slit_init(Cube<T> *c);
    void writeModel_slit();

    template <class Type> friend class Galmod;

protected:
    GALFIT_PAR par;                         //< Container for GALFIT parameters
    Cube<T>  *in;                           //< A pointer to the input cube.
    Rings<T> *inr;                          //< A pointer to the initial rings.
    Rings<T> *outr;                         //< Output rings.
    bool     inDefined = false;             //< Wheter inr have been defined inside Galfit.
    bool     outDefined = false;            //< Whether the output rings have been defined.
    bool     mpar[MAXPAR];                  //< Mask for parameters to be used.
    bool     *mask;                         //< Mask for areas to be used for chi2 calculation.
    bool     *mask2D = nullptr;             //< Mask for derived maps.
    double   arcconv;                       //< Conversion factor to arcsec.
    int      nfree;                         //< Number of free parameters.
    float    distance;                      //< Distance of the galaxy in Mpc.
    T        maxs[MAXPAR];                  //< Maximum values allowed for parameters.
    T        mins[MAXPAR];                  //< Minimum values allowed for parameters.
    double   *cfield;                       //< Convolution field.
    bool     cfieldAllocated = false;
    int      NconX;                         //< Convolution field X-dimension.
    int      NconY;                         //< Convolution field Y-dimensione
    int      wpow = 1;                      //< Weighing function power.
    bool     second = false;
    Cube<T>  *line_im;                      //< Line Image;
    bool     line_imDefined = false;
    float    *chan_noise;                   //< Noise in each channel map.
    T        data_noise;                    //< Global noise in the cube. 
    bool     chan_noiseAllocated = false;
    bool     global = false;                //< Whether to fit all parameters at once.
    bool     reverse = false;               //< Using reverse cumulative fitting
    bool     verb = true;
    

    /// Pointer to the function to be minimized (3d or 2d slit)
    typedef double (Galfit<T>::*funcPtr) (Rings<T> *, T *, int*, int*);
    funcPtr func_norm = &Model::Galfit<T>::norm_local;
    
    bool setCfield ();
    void setFree();
    void fit_straight(T ***errors, bool *fitok, ostream &fout);
    void fit_reverse(T ***errors, bool *fitok, ostream &fout);
    bool regularizeParams(std::vector<T> x, std::vector<T> y, std::vector<T> &yout, int rtype);


    /// Functions defined in galfit_min.cpp
    bool   minimize(Rings<T> *dring, T &minimum, T *pmin, Galmod<T> *modsoFar=nullptr);
    double mtry(Rings<T> *dring, T **p, T *y, T *psum, const int ihi, const double fac, Galmod<T> *modsoFar=nullptr);
    double func3D(Rings<T> *dring, T *zpar, Galmod<T> *modsoFar=nullptr);
    double getFuncValue(Rings<T> *dring, Galmod<T> *modsoFar=nullptr);
    void   Convolve(T *array, int *bsize);
    void   Convolve_fft(T *array, int *bsize);
    double norm_local(Rings<T> *dring, T *array, int *bhi, int *blo);
    double norm_azim (Rings<T> *dring, T *array, int *bhi, int *blo);
    double norm_none (Rings<T> *dring, T *array, int *bhi, int *blo);

    double slitfunc (Rings<T> *dring, T *array, int *bhi, int *blo);
    bool IsIn (int x, int y, int *blo, Rings<T> *dr, double &th);
    inline bool getSide (double theta);
    inline double getResValue(T obs, T mod, double weight, double noise_weight);
    std::vector<Pixel<T> >* getRingRegion (Rings<T> *dring, int *bhi, int *blo);

    /// Functions defined in galfit_out.cpp
    void printDetails  (Rings<T> *dr, T fmin, long pix, std::ostream& str=std::cout);
    void showInitial (Rings<T> *inr, std::ostream& Stream);
    void printInitial (Rings<T> *inr, std::string outfile);
    void DensityProfile (T *surf_dens, int *count);
    int* getErrorColumns();
    int  plotParam() {plotPar_Gnuplot(); return plotAll_Python();}
    void plotPVs_Gnuplot(Image2D<T> *pva_d, Image2D<T> *pvb_d, Image2D<T> *pva_m, Image2D<T> *pvb_m);
    void plotPar_Gnuplot();

    /// Functions defined in galfit_erros.cpp
    void getErrors(Rings<T> *dring, T** err, int ir, T minimum);

};


/// Some function to conveniently write Galfit rings
void writeHeader(std::ostream &fout, bool *mpar, bool writeErrors, bool writeBadRings);
template <class T>
void writeRing(std::ostream &fout, Rings<T> *r, int i, double toKpc, int nfree, bool writeErrors, T ***errors, bool writeBadRings, bool fitOK);
template <class T>
void printRing(std::ostream &fout, Rings<T> *r, int i, double minimum, double toKpc, bool *mpar, int nthreads);

void WarningMessage(std::ostream &fout, std::string msg);
}

#endif
