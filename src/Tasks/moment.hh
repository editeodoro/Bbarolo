//----------------------------------------------------------
// moment.hh: Definition of the MomentMap class.
//----------------------------------------------------------

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



template <class T> 
class MomentMap : public Image2D <T> 
{
public:
    MomentMap();
    ~MomentMap() {}
    MomentMap(const MomentMap &i);      
    MomentMap& operator=(const MomentMap &i);

    void input (Cube<T> *c, int *Blo, int *Bhi);
    void input (Cube<T> *c);
    void SumMap (bool msk);
    void HIMassDensityMap (bool msk);
    void ZeroMoment  (bool msk, std::string mtype="MOMENT") {storeMap(msk,0,mtype);}
    void FirstMoment (bool msk, std::string mtype="MOMENT") {storeMap(msk,1,mtype);}
    void SecondMoment(bool msk, std::string mtype="MOMENT") {storeMap(msk,2,mtype);}
    void RMSMap (float level=0.1, float sncut = 1.5);
    bool setHead(int type); 
    
private:
    Cube<T> *in;
    int blo[3],bhi[3];
    int nsubs;
    
    typedef bool (MomentMap<T>::*funcPtr) (size_t, size_t, bool, double*);
    funcPtr map_Type = &MomentMap<T>::calculateMoments;
    
    bool fitSpectrum (size_t x, size_t y, bool msk, double *bestfitpar);
    bool calculateMoments (size_t x, size_t y, bool msk, double *moments);
    void storeMap(bool msk, int whichmap, std::string map_type);
    
    
};


// A function to extract PV-diagrams
template <class T> Image2D<T>* PositionVelocity (Cube<T> *c, float x0, float y0, float Phi);


#endif
