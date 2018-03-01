//--------------------------------------------------------------------
// galmod.hh: Definition of the Galmod class.
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

// Galmod makes model observations of the HI gas in a spiral galaxy.
// The class is derived from the GALMOD task from GIPSY software.
//
// The correct way to call the class is the following:
//
// 1) call input(...) function: 
//      -Cube *c:       the Cube object from which the model is built.
//      -int *Bowup:    an array with upper coordinate limits.
//      -int *Bowlow:   an array with lower coordinate limits.
//      -Rings *ring:   a Rings object containing the input rings. 
//                      In input, distances have to be in ARCSEC, 
//                      velocities in KM/S, angles in DEGREE and
//                      column densities in HI atoms/cm^2. Required
//                      quantities are nr, xpos, ypos, radii, vrot,
//                      vsys, vdisp, dens, z0, inc and phi.
//      -int NV:        Number of subclouds in the velocity profile of 
//                      a single cloud. Can be set equal to the number 
//                      of subset. 
//      -int LTYPE[1]:  Layertype. The density profile in the direction 
//                      perpendicular to the plane of the rings. Values:
//                          = 1: gaussian layer.
//                          = 2: sech2 layer.
//                          = 3: exponential layer.
//                          = 4: Lorentzian layer.
//                          = 5: box layer.
//      -int CMODE[1]:  It determines the dependence of the number of 
//                      clouds on the surface density of the HI. 
//      -int CDENS[1]:  Surface density of clouds in the plane of the 
//                      rings per area of a pixel. Default is 1.
//      -int ISEED[-1]: Number to call the random number generator. 
//
// 2) call calculate() function.
//
//
//
//
//

#ifndef GALMOD_HH_
#define GALMOD_HH_

#include <iostream>
#include <vector>
#include <random>
#include <Arrays/cube.hh>
#include <Arrays/param.hh>
#include <Arrays/rings.hh>
#include <Utilities/utils.hh>

namespace Model {
    
template <class Type>   
class Galmod
{

public:

    Galmod();                                   //< Default constructor.
    virtual ~Galmod();                          //< Destructor.
    Galmod(const Galmod &g);                    //< Copy constructor.
    Galmod& operator=(const Galmod &g);         //< Copy operator.
    void defaults();
    
    // Obvious inline functions 
    Cube<Type>  *In () {return in;}
    Cube<Type>  *Out() {return out;}
    Rings<Type> *Ring() {return r;}
    float *getArray() {return out->Array();}
    
    void input(Cube<Type> *c, int *Boxup, int *Boxlow, Rings<Type> *rings, 
               int NV=-1, int LTYPE=1, int CMODE=1, float CDENS=1.0, int ISEED=-1);
    
    void input(Cube<Type> *c, Rings<Type> *rings, int NV=-1, int LTYPE=1, int CMODE=1, 
               float CDENS=1.0, int ISEED=-1);
               
    bool calculate();
    bool smooth(bool usescalefac=true);
    void normalize();


protected:

    Cube<Type>  *in;                        //< A pointer to the input cube.
    Cube<Type>  *out;                       //< The Cube containing the model.
    bool    outDefined;
    double  crpix3, crval3;                 //< Header keywords.
    double drval3, cdelt3;                  //<
    double  cdelt[2];                       //<
    std::string cunit3, ctype3;         
    std::string ctype[2], cunit[2]; 
    
    int     blo[2], bhi[2], bsize[2];       //< Boxes.
    double  pixsize[2], pixarea;            //< Pixels information.
    int     nsubs;                          //< Number of subsets.
    int     velsys;                         //< Type of vsys.   
    int     axtyp;                          //< Type of axes.
    float   *cd2i;                          //< Conversion from cd to intensity. 
    bool    subAllocated;                   //< 
    double  freq0;                          //< Rest frequency.
    double  chwidth;                        //< Velocity width of the channels in M/S.
    double  sig_instr;                      //< Instrumental dispersion in M/S.
    float   crota2;                         //< Rotation angle.
        
    Rings<Type> *r;                         //< Set of rings.       
    bool    ringDefined;            
    
    int     ltype;                          //< Layer type.
    int     cmode;                           
    float   cdens;                          //< Surf. dens. of clouds per area of a pixel.
    int     iseed;  
    float   arcmconv;                       //< Conversion to arcmin.
    bool    readytomod;
    bool    modCalculated;
    
    // Random number engines
    std::mt19937 generator;
    std::uniform_real_distribution<float> uniform;
    std::normal_distribution<float> gaussia;
    
    /// Private functions.
    
    void    initialize(Cube<Type> *c, int *Boxup, int *Boxlow);
    void    ringIO(Rings<Type> *rings);
    void    galmod();
    void    NHItoRAD();
    double  velgrid(double v);
    double  fdev(int &idum);
    
};

}

//#include "galmod.cpp"

#endif
