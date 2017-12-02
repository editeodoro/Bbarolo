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
//                      It must be negative.
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
#include <Arrays/cube.hh>


template <class T>
struct Ring {
    
    Ring (T radius, T xpos, T ypos, T vsys, T vrot, T vdisp, T vrad, T vvert, T dvdz,
          T zcyl, T dens, T z0, T inc, T phi, T pa, int nv) : radius(radius), xpos(xpos),
          ypos(ypos), vsys(vsys), vrot(vrot), vdisp(vdisp), vrad(vrad), vvert(vvert),
          dvdz(dvdz), zcyl(zcyl), dens(dens), z0(z0), inc(inc), phi(phi), pa(pa), nv(nv) {}
    
    T    radius;         //< Radius
    T    xpos;           //< X-center 
    T    ypos;           //< Y-center
    T    vsys;           //< Systemic velocity 
    T    vrot;           //< Rotational velocity.
    T    vdisp;          //< Velocity dispersion.
    T    vrad;           //< Radial velocity.
    T    vvert;          //< Vertical velocity.
    T    dvdz;           //< Vertical rotational gradient (km/s/arcsec).
    T    zcyl;           //< Height where the rotational gradient starts.
    T    dens;           //< Column densities.
    T    z0;             //< Scaleheights of the HI-layer.
    T    inc;            //< Inclination angles.
    T    phi;            //< Position angles.
    T    pa;             //< Position angles+rotation angles.
    int  nv;             //< Number of subclouds. 
};

template <class T>
struct Rings {
    int     nr;                     //< Number of rings.
    double radsep;                  //< Separation between rings.
    std::vector<T> radii;
    std::vector<T> xpos;
    std::vector<T> ypos;
    std::vector<T> vsys;
    std::vector<T> vrot;
    std::vector<T> vdisp;
    std::vector<T> vrad;
    std::vector<T> vvert;
    std::vector<T> dvdz;
    std::vector<T> zcyl;
    std::vector<T> dens;
    std::vector<T> z0;
    std::vector<T> inc;
    std::vector<T> phi;
    std::vector<T> pa;
    std::vector<int>  nv;
    
    
    void ClearAll () {
        xpos.clear(); ypos.clear(); vsys.clear(); radii.clear(); vrot.clear(); vdisp.clear();
        vrad.clear(); vvert.clear(); dvdz.clear(); zcyl.clear(); dens.clear(); z0.clear();
        inc.clear(); phi.clear(); pa.clear(); nv.clear();
    } 
    
    void setRings (int size, T* radii, T* xpos, T* ypos, T* vsys, T* vrot, T* vdisp, T* vrad, 
                   T* vvert, T* dvdz, T* zcyl, T* dens, T* z0, T* inc, T* phi, int* nv) {
        this->nr = size;
        for (int i=0; i<nr; i++) {
            this->radii.push_back(radii[i]);            
            this->xpos.push_back(xpos[i]);
            this->ypos.push_back(ypos[i]);
            this->vsys.push_back(vsys[i]);
            this->vrot.push_back(vrot[i]);
            this->vdisp.push_back(vdisp[i]);
            this->vrad.push_back(vrad[i]);
            this->vvert.push_back(vvert[i]);
            this->dvdz.push_back(dvdz[i]);
            this->zcyl.push_back(zcyl[i]);
            this->dens.push_back(dens[i]);
            this->z0.push_back(z0[i]);
            this->inc.push_back(inc[i]);
            this->phi.push_back(phi[i]);
            this->pa.push_back(0);
            this->nv.push_back(nv[i]);
        }
        this->radsep = 0;
        if (nr>1) this->radsep = radii[1]-radii[0];
    }
};


template <class T>
struct newRings {
    int     nr;                     //< Number of rings.
    double radsep;                  //< Separation between rings.
    std::vector<Ring<T> > rings;     
    inline Ring<T>& operator[] (size_t n) {return rings[n];}

};



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
        
    
    void input(Cube<Type> *c, int *Boxup, int *Boxlow, Rings<Type> *rings, 
               int NV, int LTYPE=1, int CMODE=1, float CDENS=1.0, int ISEED=-1);
    
    void input(Cube<Type> *c, Rings<Type> *rings, int LTYPE);
               
    void calculate();
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
        
    /// Private functions.
    
    void    initialize(Cube<Type> *c, int *Boxup, int *Boxlow);
    void    ringIO(Rings<Type> *rings);
    void    galmod();
    void    NHItoRAD();
    double  velgrid(double v);
    double  fdev(int& idum);
    double  gasdev(int &idum);
    int     iran(int &idum);

};

}

//#include "galmod.cpp"

#endif
