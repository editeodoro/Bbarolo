//------------------------------------------------------------
// moment.hh: Definition of the MomentMap and PvSlice classes.
//------------------------------------------------------------

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

#ifndef MOMENTMAP_HH_
#define MOMENTMAP_HH_


#include <Arrays/cube.hh>
#include <Arrays/image.hh>


//////////////////////////////////////////////////////////////////////////
/// MomentMap class to extract a moment maps from a datacube
//////////////////////////////////////////////////////////////////////////
template <class T> 
class MomentMap : public Image2D <T> 
{
public:
    MomentMap();
    ~MomentMap() {}
    MomentMap(const MomentMap &i);      
    MomentMap& operator=(const MomentMap &i);

    void input (Cube<T> *c, int *Blo, int *Bhi, bool *m=nullptr);
    void input (Cube<T> *c, bool *m=nullptr);
    void SumMap (bool msk);
    void HIMassDensityMap (bool msk);
    void ZeroMoment  (bool msk, std::string mtype="MOMENT") {storeMap(msk,0,mtype);}
    void FirstMoment (bool msk, std::string mtype="MOMENT") {storeMap(msk,1,mtype);}
    void SecondMoment(bool msk, std::string mtype="MOMENT") {storeMap(msk,2,mtype);}
    void RMSMap (float level=0.1, float sncut = 1.5);
    void SNMap(bool msk);
    bool setHead(int type); 
    
    bool fitSpectrum (size_t x, size_t y, bool msk, double *bestfitpar);
    bool calculateMoments (size_t x, size_t y, bool msk, double *moments);
    void storeMap(bool msk, int whichmap, std::string map_type);

private:
    Cube<T> *in;
    int blo[3],bhi[3];
    int nsubs;
    bool *mask = nullptr;
    
    typedef bool (MomentMap<T>::*funcPtr) (size_t, size_t, bool, double*);
    funcPtr map_Type = &MomentMap<T>::calculateMoments;
    
};


// A function to extract all kinematic maps at the same time
template <class T>
std::vector< MomentMap<T> > getAllMoments(Cube<T> *c, bool usemask=true, bool *mask=nullptr, std::string mtype="MOMENT");


/////////////////////////////////////////////////////////////////////////////////////
/// A class to extract position-velocity slices
/////////////////////////////////////////////////////////////////////////////////////
template <class T>
class PvSlice : public Image2D<T>
/// PvSlicer class to extract a position-velocity slice from a datacube
/// The slice is always a straigh line and can be defined:
///  1) Through two points (constructor #1)
///  2) Through one point and an angle (constructor #2)
///
/// In the first case, the slice can include just a part of the datacube,
/// while in the second the entire cube is sliced.
/// Usage: call one of the constructors and the slice().
{
public:
    PvSlice(Cube<T> *c);
    PvSlice(Cube<T> *c, double x1, double y1, double x2, double y2, int width=0) 
            : in(c), x1(x1), x2(x2), y1(y1), y2(y2), isAngle(false), width(width) {}
    PvSlice(Cube<T> *c, double x0, double y0, double angle, int width=0) 
            : in(c), x0(x0), y0(y0), angle(angle), isAngle(true), width(width) {}
    ~PvSlice(){if (locusAllocated) {delete [] x_locus; delete [] y_locus;}}
    bool slice ();
    bool slice_old ();
    
private:
    Cube<T> *in;                    //< A pointer to the input datacube
    int    xpix, ypix, zpix;        //< Size of the datacube to slice
    double x1, x2, y1, y2;          //< Defining slice with two points
    double x0, y0, angle;           //< Defining slice with a point and a angle (from North->West)
    bool   isAngle;                 //< Whether slice is defined with angle
    double *x_locus, *y_locus;      //< The slice
    bool  locusAllocated = false;   //< Whether x_locus and y_locus are allocated
    int   num_points;               //< Number of pixels along the slice
    int   width = 0;                //< Half width of the slice in pixels
    float nalias = 0.5;             //< Type of anti-aliasing.
    
    double weight (double x, double y, double cx, double cy) {return fabs((1-(x-cx))*(1-(y-cy)));}
     
    bool  define_slice();
    bool  check_bounds (double *blx, double *bly, double *Trx, double *Try);
    bool  pvslice ();
    void  define_header();
};

// A front-end function to extract PV-diagrams with a center and a angle.
template <class T>
PvSlice<T>* PositionVelocity (Cube<T> *c, float x0, float y0, float Phi, bool oldmethod=true);

#endif
