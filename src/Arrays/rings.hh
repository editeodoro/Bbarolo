//---------------------------------------------------------------
// rings.hh: Definitions of ring/rings/shell containers.
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
#include <iomanip>

using namespace std;

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
class Rings 
{
public:
    int     nr;                 //< Number of rings.
    double radsep = -1;         //< Separation between rings.
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
    
    Rings() {nr=radsep=0;}
    ~Rings() {ClearAll();}
    
    void ClearAll () {
        xpos.clear(); ypos.clear(); vsys.clear(); radii.clear(); vrot.clear(); vdisp.clear();
        vrad.clear(); vvert.clear(); dvdz.clear(); zcyl.clear(); dens.clear(); z0.clear();
        inc.clear(); phi.clear(); pa.clear();
    } 
    
    
    void addRing (T radii, T xpos, T ypos, T vsys, T vrot, T vdisp, T vrad, 
                  T vvert, T dvdz, T zcyl, T dens, T z0, T inc, T phi) {
        
        this->radii.push_back(radii);
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
        
        this->nr = this->radii.size();
    }
    
    void insertRing (unsigned where, T radii, T xpos, T ypos, T vsys, T vrot, T vdisp, 
                     T vrad, T vvert, T dvdz, T zcyl, T dens, T z0, T inc, T phi) {
        
        this->radii.insert(this->radii.begin()+where,radii);
        this->xpos.insert(this->xpos.begin()+where,xpos);
        this->ypos.insert(this->ypos.begin()+where,ypos);
        this->vsys.insert(this->vsys.begin()+where,vsys);
        this->vrot.insert(this->vrot.begin()+where,vrot);
        this->vdisp.insert(this->vdisp.begin()+where,vdisp);
        this->vrad.insert(this->vrad.begin()+where,vrad);
        this->vvert.insert(this->vvert.begin()+where,vvert);
        this->dvdz.insert(this->dvdz.begin()+where,dvdz);
        this->zcyl.insert(this->zcyl.begin()+where,zcyl);
        this->dens.insert(this->dens.begin()+where,dens);
        this->z0.insert(this->z0.begin()+where,z0);
        this->inc.insert(this->inc.begin()+where,inc);
        this->phi.insert(this->phi.begin()+where,phi);
        //this->pa.insert(this->pa.begin()+where,0);
        
        this->nr = this->radii.size();
    }

    void addRings(int size, T *radii, T *xpos, T *ypos, T *vsys, T *vrot, T *vdisp, 
                  T *vrad, T *vvert, T *dvdz, T *zcyl, T *dens, T *z0, T *inc, T *phi) {
        for (int i=0; i<size; i++) {
            this->addRing(radii[i],xpos[i],ypos[i],vsys[i],vrot[i],vdisp[i],vrad[i],
                          vvert[i],dvdz[i],zcyl[i],dens[i],z0[i],inc[i],phi[i]);
        }
    }
    
    
    void addRings(int size, T *radii, T xpos, T ypos, T vsys, T vrot, T vdisp, 
                      T vrad, T vvert, T dvdz, T zcyl, T dens, T z0, T inc, T phi) {
            for (int i=0; i<size; i++) {
                this->addRing(radii[i],xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens,z0,inc,phi);
            }
        }
    
    
    void setRings (int size, T* radii, T* xpos, T* ypos, T* vsys, T* vrot, T* vdisp, T* vrad, 
                   T* vvert, T* dvdz, T* zcyl, T* dens, T* z0, T* inc, T* phi) {
        this->addRings(size,radii,xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens,z0,inc,phi);
        this->radsep = this->nr>1 ? (this->radii[1]-this->radii[0]) : 0;
    }
    
    
    void setRings (T minrad, T maxrad, T radsep, T xpos, T ypos, T vsys, T vrot, T vdisp, 
                   T vrad, T vvert, T dvdz, T zcyl, T dens, T z0, T inc, T phi) {
        
        int size = (maxrad-minrad)/radsep + 1;
        for (int i=0; i<size; i++) {
            this->addRing(minrad+i*radsep,xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens,z0,inc,phi);
        }
        this->radsep = radsep;    
    }

    
    void deleteRing(int nring) {

        this->radii.erase(this->radii.begin()+nring);
        this->xpos.erase(this->xpos.begin()+nring);
        this->ypos.erase(this->ypos.begin()+nring);
        this->vsys.erase(this->vsys.begin()+nring);
        this->vrot.erase(this->vrot.begin()+nring);
        this->vdisp.erase(this->vdisp.begin()+nring);
        this->vrad.erase(this->vrad.begin()+nring);
        this->vvert.erase(this->vvert.begin()+nring);
        this->dvdz.erase(this->dvdz.begin()+nring);
        this->zcyl.erase(this->zcyl.begin()+nring);
        this->dens.erase(this->dens.begin()+nring);
        this->z0.erase(this->z0.begin()+nring);
        this->inc.erase(this->inc.begin()+nring);
        this->phi.erase(this->phi.begin()+nring);
        //this->pa.erase(this->pa.begin()+nring);

        this->nr = this->radii.size();
    }

    void printRing(int nring, ostream& theStream = cout) {
        
        int m=8, n=11;
        
        theStream << "\n Ring #" << nring+1 << " (idx=" << nring <<") at radius " << this->radii[nring] << " arcsec \n";
        theStream << fixed << setprecision(2);
        
        string s;
        s = "    Vrot";
        theStream << setw(n) << left << s << setw(3) << right << "= "
                  << setw(m) << this->vrot[nring] << left << setw(m) << "  km/s";

        s = "        Disp";
        theStream << setw(n+4) << left << s << setw(3) << right << "= "
                  << setw(m-1) << this->vdisp[nring]
                  << left << setw(m) << "  km/s" << endl;

        s = "    Vrad";
        theStream << setw(n) << left << s << setw(3) << right << "= "
                  << setw(m) << this->vrad[nring] << left << setw(m) << "  km/s";
        
        s = "        Vsys";
        theStream << setw(n+4) << left << s << setw(3) << right << "= "
                  << setw(m-1) << this->vsys[nring] << left << setw(m) << "  km/s" << endl;
            
        s = "    Inc";
        theStream << setw(n) << left << s << setw(3) << right << "= "
                  << setw(m) << this->inc[nring] << left << setw(m) << "  deg";

        s = "        PA";
        theStream << setw(n+4) << left << s << setw(3) << right << "= "
                  << setw(m-1) << this->phi[nring] << left << setw(m) << "  deg" << endl;

        s = "    Xpos";
        theStream << setw(n) << left << s << setw(3) << right << "= "
                  << setw(m) << this->xpos[nring] << left << setw(m) << "  pix";

        s = "        Ypos";
        theStream << setw(n+4) << left << s << setw(3) << right << "= "
                  << setw(m-1) << this->ypos[nring] << left << setw(m) << "  pix" << endl;
            
        s = "    Z0";
        theStream << setw(n) << left << s << setw(3) << right << "= "
                  << setw(m) << this->z0[nring] << left << setw(m) << "  arcs";
        
        s = "        Dens";
        theStream << setw(n+4) << left << s << setw(3) << right << "= "
                  << setw(m-1) << this->dens[nring]/1E20 << left << setw(m) << " 1E20 1/cm2" << endl;
        
        s = "    VVERT";
        theStream << setw(n) << left << s << setw(3) << right << "= "
                  << setw(m) << this->vvert[nring] << left << setw(m) << "  km/s";
        
        s = "        ZCYL";
        theStream << setw(n+4) << left << s << setw(3) << right << "= "
                  << setw(m-1) << this->zcyl[nring] << left << setw(m) << " arcs" << endl;
        
        s = "    DVDZ";
        theStream << setw(n) << left << s << setw(3) << right << "= "
                  << setw(m) << this->dvdz[nring] << left << setw(m) << "  km/s/arcs";
        
        theStream << endl;
    }
    
    void printRings(ostream& theStream = cout){for (int i=0; i<this->nr; i++) this->printRing(i,theStream);}
    
    template <class M>
    friend std::ostream& operator<< (ostream& theStream, Rings<M>& r);

};

template <class M>
ostream& operator<< (ostream& theStream, Rings<M>& r) {r.printRings(theStream); return theStream;}


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
    
