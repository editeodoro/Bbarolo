//---------------------------------------------------------------
// rings.hh: Definitions of ring containers.
//---------------------------------------------------------------

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

#ifndef RINGS_HH_
#define RINGS_HH_

#include <iostream>
#include <vector>

template <class T>
struct Ring {
    
    Ring (T radius, T xpos, T ypos, T vsys, T vrot, T vdisp, T vrad, T vvert, T dvdz,
          T zcyl, T dens, T z0, T inc, T phi, T pa) : radius(radius), xpos(xpos),
          ypos(ypos), vsys(vsys), vrot(vrot), vdisp(vdisp), vrad(vrad), vvert(vvert),
          dvdz(dvdz), zcyl(zcyl), dens(dens), z0(z0), inc(inc), phi(phi), pa(pa) {}
    
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
};


template <class T>
struct Rings {
    int     nr;                 //< Number of rings.
    double radsep;              //< Separation between rings.
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
    int id;                     // ID of the ring (for 3dFIT)
    
    
    void ClearAll () {
        xpos.clear(); ypos.clear(); vsys.clear(); radii.clear(); vrot.clear(); vdisp.clear();
        vrad.clear(); vvert.clear(); dvdz.clear(); zcyl.clear(); dens.clear(); z0.clear();
        inc.clear(); phi.clear(); pa.clear();
    } 
    
    void setRings (int size, T* radii, T* xpos, T* ypos, T* vsys, T* vrot, T* vdisp, T* vrad, 
                   T* vvert, T* dvdz, T* zcyl, T* dens, T* z0, T* inc, T* phi) {
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
        }
        this->radsep = 0;
        if (nr>1) this->radsep = radii[1]-radii[0];
    }
    
    void setRings (T minrad, T maxrad, T radsep, T xpos, T ypos, T vsys, T vrot, T vdisp, 
                   T vrad, T vvert, T dvdz, T zcyl, T dens, T z0, T inc, T phi) {
        
        this->nr = (maxrad-minrad)/radsep + 1;
        this->radsep = radsep;    
        for (int i=0; i<this->nr; i++) {
            this->radii.push_back(minrad+i*radsep);            
            this->xpos.push_back(xpos);
            this->ypos.push_back(ypos);
            this->vsys.push_back(vsys);
            this->vrot.push_back(vrot);
            this->vdisp.push_back(vdisp);
            this->vrad.push_back(vrad);
            this->vvert.push_back(vvert);
            this->dvdz.push_back(dvdz);
            this->zcyl.push_back(zcyl);
            this->dens.push_back(dens);
            this->z0.push_back(z0);
            this->inc.push_back(inc);
            this->phi.push_back(phi);
            this->pa.push_back(0);
        }
    }
};


template <class T>
struct newRings {
    int     nr;                     //< Number of rings.
    double radsep;                  //< Separation between rings.
    std::vector<Ring<T> > rings;     
    inline Ring<T>& operator[] (size_t n) {return rings[n];}

};

// A structure to define a spherical shell
template <class T>
struct Shells {
    int    ns;                  //< Number of shells.
    double sep;                 //< Separation between shells.
    std::vector<T> radii;       //< Radii
    std::vector<T> xpos;        //< Xcenter
    std::vector<T> ypos;        //< Ycenter
    std::vector<T> vsys;        //< Systemic velocity
    std::vector<T> vrot;        //< Azimuthal velocity
    std::vector<T> vsph;        //< Radial velocity (spherical)
    std::vector<T> vdisp;       //< Velocity dispersion
    std::vector<T> dens;        //< Column density
    std::vector<T> inc;         //< Inclination angle
    std::vector<T> pa;          //< Position angle
    std::vector<T> openang;     //< Opening angle
};

#endif
    