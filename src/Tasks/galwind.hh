// -----------------------------------------------------------------------
// galwind.hh: Definition and functions for the GalWind class.
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

#ifndef GALWIND_HH
#define GALWIND_HH

#include <iostream>
#include <string>
#include <Arrays/cube.hh>

template <class T>   
class GalWind
{
// GalWind is derived from Federico Lelli algorithm.
// This class build a bi-conical outflow model with constant outflow velocity 
// and velocity dispersion. The cone is build adding N cylinders with the same 
// thickness and symmetry axis but increasing diameter. After adding the cylinders, 
// the projected bi-cone is convoluted with the observed PSF (assumed to be Gaussian).
// IMPORTANT NOTE:
// The conversion from column density to flux/beam does NOT depend on 
// frequency or any other property. For Dens = 1.25e20 atoms/cm^-2 = 1 Msun/pc^2, it 
// returns 0.0135 Jy/pixel^2 for any emission line. For HI data, to get the flux in
// Jy/pixel^2, each channel maps should be divided by 0.005*(f/c)**2/cdelt. where lightspeed
// in m/s, f = mid frequency of that channel in Hz, and cdelt = channel width in Hz.

public:
    // Constructors:
    // 1) Give an object Cube with parameters recorded in c->pars()
    // 2) Give an object Cube + a GALWIND_PAR structure
    // 3) Give all paramaters separately
    GalWind(Cube<T> *c) : GalWind(c,c->pars().getParGW()) {}
    GalWind(Cube<T> *c, GALWIND_PAR &p) {in=c, par=p; in->checkBeam();}
    GalWind(Cube<T> *c, double x0, double y0, double pa, double inc, double disp, 
            double dens, double vsys, double vw, double openang, double htot, 
            int denstype, int ntot=25, int cdens=10, int nv=10);
    // Copy constructor and destructor
    GalWind (const GalWind<T>& gw) {operator=(gw);}
    GalWind& operator= (const GalWind<T>& gw);
    ~GalWind() {if (outDefined) delete out;}
    
    // Overloaded () operator to easily access output array
    inline T& operator() (size_t npix) {return out->Array(npix);}
    inline T& operator() (size_t x, size_t y,size_t z) {return out->Array(out->nPix(x,y,z));}
    
    // Obvious inline functions
    Cube<T>*    getOut() {return out;}
    T*          getArray() {return out->Array();}
    GALWIND_PAR getPar() {return par;}
    
    // Functions to calculate the model, smooth it and write it to FITS files
    bool compute();
    bool smooth(bool scalefac=true);
    bool writeFITS (std::string fname="", bool fullHead=false);
    bool writeMomentMaps();

private:
    Cube<T>     *in;                //< A pointer to the input cube.
    Cube<T>     *out;               //< The Cube containing the model.
    bool        outDefined = false; //< Whether the out Cube is defined
    GALWIND_PAR par;                //< Parameters of the task
};

#endif
