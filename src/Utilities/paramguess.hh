//--------------------------------------------------------------------
// paramguess.hh: A class to estimate initial parameters for 3D Fit
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

#ifndef PARAMGUESS_HH_
#define PARAMGUESS_HH_

#include <iostream>
#include <Arrays/cube.hh>
#include <Map/detection.hh>


template <class T>
class ParamGuess
// ParamGuess is a class to estimate geometrical and kinematical parameter
// for a galaxy from an emission-line datacube. These parameters can be used
// for example as initial guesses for 3DFIT task in BBarolo.
// 
// The class uses the total integrated spectrum, the intensity map and the  
// velocity field. 
//
// After calling the constructor, call findAll() to estimate all parameters.
// Individual parameters can be estimated as well, but they must follow this
// order: 1) findSystemicVelocity(), 2) findCentre(), 3) findPositionAngle(),
// 4) findInclination(), 5) findRings() and findRotationVelocity()
//
{
public:
    T   xcentre = 0;        //< Centre on X-axis
    T   ycentre = 0;        //< Centre on Y-axis
    T   vsystem = 0;        //< Systemic velocity
    T   inclin  = 0;        //< Inclination angle
    T   posang  = 0;        //< Position angle
    T   Rmax    = 0;        //< Maximum radius
    T   vrot    = 0;        //< Rotation velocity
    int nrings  = 0;        //< Number of rings
    T   radsep  = 0;        //< Ring width
    T*      Intmap;         //< Intensity map
    T*      Vemap;          //< Velocity field

    // Constructor and destructor
    ParamGuess(Cube<T> *c, Detection<T> *object);
    ~ParamGuess() {delete [] Vemap; delete [] Intmap;}

    void setCentre(T xcen, T ycen) {xcentre=xcen; ycentre=ycen;}
    void setPosang(T pa) {posang=pa;}
    
    void findAll();
    void findCentre();
    void findSystemicVelocity();
    void findPositionAngle(int algorithm=1);
    void findInclination(int algorithm=2);
    void findRotationVelocity();
    void findRings();
    void tuneWithTiltedRing();
    int  plotGuess(std::string outfile="initial_guesses.pdf");
    
private:
    Cube<T> *in;            //< Pointer to input datacube
    Detection<T> *obj;      //< Pointer to the detection to be analysed

    
    int     major_max[2];   //< Upper pixel coordinates on the major axis
    int     minor_max[2];   //< Upper pixel coordinates on the minor axis
    int     major_min[2];   //< Lower pixel coordinates on the major axis
    int     minor_min[2];   //< Lower pixel coordinates on the minor axis
    float   totflux_obs;    //< Total flux

    void  setAxesLine(T xcen, T ycen, T pa, float *maj, float *min);
    T     findAxisLength(float *lpar, int *coords_up, int *coords_low);
    T     funcEllipse(T *mypar);
    T     funcIncfromMap(T *mypar);
    bool  fitSimplex(int ndim, T **p);  
    
    typedef T (ParamGuess<T>::*funcPtr) (T *);
    funcPtr func = &ParamGuess<T>::funcEllipse;
};




#endif
