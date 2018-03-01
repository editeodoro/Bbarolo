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
    Galfit() {defaults();}
    Galfit(Cube<T> *c);
    Galfit(Cube<T> *c, Rings<T> *inrings, GALFIT_PAR *p) {setup(c,inrings,p);}
    Galfit (Cube<T> *c, Rings<T> *inrings, float DELTAINC=5, float DELTAPHI=15, 
            int LTYPE=1, int FTYPE=2, int WFUNC=1, int BWEIGHT=1, int NV=-1, double TOL=1E-03, 
            int CDENS=10, int STARTRAD=0, std::string MASK="SMOOTH", std::string NORM="LOCAL", 
            std::string FREE="VROT VDISP INC PA", std::string SIDE="B", bool TWOSTAGE=true, 
            std::string POLYN="bezier", bool ERRORS=false, bool SMOOTH=true, float DISTANCE=-1, 
            double redshift=-1, double RESTWAVE=-1, std::string OUTFOLD="./", int THREADS=1);
    virtual ~Galfit();
    Galfit(const Galfit &g) {operator=(g);}     //< Copy constructor.
    Galfit& operator=(const Galfit &g);         //< Copy operator.
    
    Cube<T>  *In () {return in;}
    Rings<T> *Inrings  () {return inr;}
    Rings<T> *Outrings () {return outr;}

    /// Functions defined in galfit.cpp
    void setup (Cube<T> *c, Rings<T> *inrings, GALFIT_PAR *p);
    void galfit();
    Galmod<T>* getModel();
    bool SecondStage();

    /// Functions defined in galfit_out.cpp
    void writeModel(std::string normtype, bool makeplots=true);
    void writePVs(Cube<T> *mod, std::string suffix="");
    int  plotAll_Python ();
    

    /// Functions defined in slitfit.cpp
    void slit_init(Cube<T> *c);
    void writeModel_slit();

    template <class Type> friend class Galmod;

protected:
    Cube<T>  *in;               //< A pointer to the input cube.
    Rings<T> *inr;              //< A pointer to the initial rings.
    Rings<T> *outr;             //< Output rings.
    bool     inDefined;         //< Wheter inr have been defined inside Galfit.
    bool     outDefined;        //< Whether the output rings have been defined.
    bool     *mpar;             //< Mask for parameters to be used.
    bool     *mask;             //< Mask for areas to be used for chi2 calculation.
    double   arcconv;           //< Conversion factor to arcsec.
    int      nfree;             //< Number of free parameters.
    float    distance;          //< Distance of the galaxy in Mpc.  
    T        *maxs;             //< Maximum values allowed for parameters.
    T        *mins;             //< Minimum values allowed for parameters.
    double   *cfield;           //< Convolution field.  
    bool     cfieldAllocated;
    int      NconX;             //< Convolution field X-dimension.
    int      NconY;             //< Convolution field Y-dimensione 
    int      wpow;              //< Weighing function power.
    int      anglepar;          //< Number of parameter for polynomial fit of angles.
    bool     second;
    bool     verb;
    Cube<T>  *line_im;          //< Line Image;
    bool     line_imDefined;
    float    *chan_noise;        //< Noise in each channel map.
    bool     chan_noiseAllocated;
    bool     global;             //< Whether to fit all parameters at once.

    GALFIT_PAR par;             // Container for GALFIT parameters

    /// Pointer to the function to be minimized (3d or 2d slit)
    typedef double (Galfit<T>::*funcPtr) (Rings<T> *, T *, int*, int*);
    funcPtr func_norm;
    
    void defaults ();
    bool setCfield ();
    void setFree();

    /// Functions defined in galfit_min.cpp
    bool   minimize(Rings<T> *dring, T &minimum, T *pmin);
    T      mtry(Rings<T> *dring, T **p, T *y, T *psum, const int ihi, const double fac);
    T      func3D(Rings<T> *dring, T *zpar);
    T      model(Rings<T> *dring);
    void   Convolve(T *array, int *bsize);
    void   Convolve_fft(T *array, int *bsize);
    double norm_local(Rings<T> *dring, T *array, int *bhi, int *blo);
    double norm_azim (Rings<T> *dring, T *array, int *bhi, int *blo);
    double norm_none (Rings<T> *dring, T *array, int *bhi, int *blo);

    double slitfunc  (Rings<T> *dring, T *array, int *bhi, int *blo);
    inline bool IsIn (int x, int y, int *blo, Rings<T> *dr, double &th);
    inline bool getSide (double theta);
    inline double getFuncValue(T obs, T mod, double weight, double noise_weight);
    std::vector<Pixel<T> >* getRingRegion (Rings<T> *dring, int *bhi, int *blo);

    /// Functions defined in galfit_out.cpp
    void printDetails  (Rings<T> *dr, T fmin, long pix, std::ostream& str=std::cout);
    void showInitial (Rings<T> *inr, std::ostream& Stream);
    void printInitial (Rings<T> *inr);
    void DensityProfile (T *surf_dens, int *count);
    int* getErrorColumns();
    int  plotParam() {plotPar_Gnuplot();return plotAll_Python();}
    void plotPVs_Gnuplot(Image2D<T> *pva_d, Image2D<T> *pvb_d, Image2D<T> *pva_m, Image2D<T> *pvb_m);
    void plotPar_Gnuplot();

    /// Functions defined in galfit_erros.cpp
    void getErrors(Rings<T> *dring, T** err, int ir, T minimum);

};

}

template <class T> T polyn (T*, T*, int);
template <class T> void polynd (T*, T*, T*, int);
template <class T> void fpolyn (T , T*, int, T&, T*);

#endif
