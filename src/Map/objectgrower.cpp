// -----------------------------------------------------------------------
// objectgrower.cpp: Implementation of the object growing functions
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

#include <iostream>
#include <vector>
#include <algorithm>
#include <Map/objectgrower.hh>
#include <Map/detection.hh>
#include <Arrays/cube.hh>
#include <Arrays/stats.hh>
#include <Map/voxel.hh>

template <class T>
ObjectGrower<T>::ObjectGrower() {}
template ObjectGrower<short>::ObjectGrower();
template ObjectGrower<int>::ObjectGrower();
template ObjectGrower<long>::ObjectGrower();
template ObjectGrower<float>::ObjectGrower();
template ObjectGrower<double>::ObjectGrower();


template <class T>
ObjectGrower<T>::ObjectGrower(ObjectGrower<T> &o) {
    this->operator=(o);
}
template ObjectGrower<short>::ObjectGrower(ObjectGrower<short>&);
template ObjectGrower<int>::ObjectGrower(ObjectGrower<int>&);
template ObjectGrower<long>::ObjectGrower(ObjectGrower<long>&);
template ObjectGrower<float>::ObjectGrower(ObjectGrower<float>&);
template ObjectGrower<double>::ObjectGrower(ObjectGrower<double>&);


template <class T>
ObjectGrower<T>& ObjectGrower<T>::operator=(const ObjectGrower<T> &o) {
    
    if(this == &o) return *this;
    this->itsFlagArray = o.itsFlagArray;
    this->itsArrayDim = o.itsArrayDim; 
    this->itsGrowthStats = o.itsGrowthStats;
    this->itsSpatialThresh = o.itsSpatialThresh;
    this->itsVelocityThresh = o.itsVelocityThresh;
    this->itsFluxArray = o.itsFluxArray;
    return *this;
}
template ObjectGrower<short>& ObjectGrower<short>::operator=(const ObjectGrower<short>&);
template ObjectGrower<int>& ObjectGrower<int>::operator=(const ObjectGrower<int>&);
template ObjectGrower<long>& ObjectGrower<long>::operator=(const ObjectGrower<long>&);
template ObjectGrower<float>& ObjectGrower<float>::operator=(const ObjectGrower<float>&);
template ObjectGrower<double>& ObjectGrower<double>::operator=(const ObjectGrower<double>&);


template <class T>
void ObjectGrower<T>::define(Cube<T> *theCube) {
    
    /// @details This copies all necessary information from the Cube
    /// and its parameters & statistics. It also defines the array of
    /// pixel flags, which involves looking at each object to assign
    /// them as detected, all blank & "milky-way" pixels to assign
    /// them appropriately, and all others to "available". It is only
    /// the latter that will be considered in the growing function.
    /// @param theCube A pointer to a duchamp::Cube 

    itsGrowthStats = Statistics::Stats<T>(theCube->stat()); 
    if(theCube->pars().getFlagUserGrowthThreshold())
        itsGrowthStats.setThreshold(theCube->pars().getGrowthThreshold());
    else
        itsGrowthStats.setThresholdSNR(theCube->pars().getGrowthCut());    
    itsGrowthStats.setUseFDR(false);

    itsFluxArray = theCube->Array();

    itsArrayDim = std::vector<size_t>(3);
    itsArrayDim[0]=theCube->DimX();
    itsArrayDim[1]=theCube->DimY();
    itsArrayDim[2]=theCube->DimZ();
    size_t spatsize=itsArrayDim[0]*itsArrayDim[1];
    size_t fullsize=spatsize*itsArrayDim[2];

    if(theCube->pars().getFlagAdjacent()) itsSpatialThresh = 1;
    else itsSpatialThresh = int(theCube->pars().getThreshS());
    itsVelocityThresh = int(theCube->pars().getThreshV());

    itsFlagArray = std::vector<STATE>(fullsize,AVAILABLE);

    for(int o=0;o<theCube->getNumObj();o++){
        std::vector<Voxel<T> > voxlist = theCube->getObject(o).getPixelSet();
        for(size_t i=0; i<voxlist.size(); i++){
            size_t pos=voxlist[i].getX()+voxlist[i].getY()*itsArrayDim[0]+voxlist[i].getZ()*spatsize;
            itsFlagArray[pos] = DETECTED;
        }
    }

}
template void ObjectGrower<short>::define(Cube<short>*);
template void ObjectGrower<int>::define(Cube<int>*);
template void ObjectGrower<long>::define(Cube<long>*);
template void ObjectGrower<float>::define(Cube<float>*);
template void ObjectGrower<double>::define(Cube<double>*);


template <class T>
void ObjectGrower<T>::updateDetectMap(short *map) {

    int numNondegDim=0;
    for(int i=0;i<3;i++) if(itsArrayDim[i]>1) numNondegDim++;

    if(numNondegDim>1) {
        size_t spatsize=itsArrayDim[0]*itsArrayDim[1];
        for(size_t xy=0;xy<spatsize;xy++){
            short ct=0;
            for(size_t z=0;z<itsArrayDim[2];z++){
                if(itsFlagArray[xy+z*spatsize]==DETECTED) ct++;
            }
            map[xy]=ct;
        }
    }
    else{
        for(size_t z=0;z<itsArrayDim[2];z++){
            map[z] = (itsFlagArray[z] == DETECTED) ? 1 : 0;
        }
    }

}
template void ObjectGrower<short>::updateDetectMap(short*);
template void ObjectGrower<int>::updateDetectMap(short*);
template void ObjectGrower<long>::updateDetectMap(short*);
template void ObjectGrower<float>::updateDetectMap(short*);
template void ObjectGrower<double>::updateDetectMap(short*);


template <class T>
void ObjectGrower<T>::grow(Detection<T> *theObject) {
    
    /// @details This function grows the provided object out to the
    /// secondary threshold provided in itsGrowthStats. For each pixel
    /// in an object, all surrounding pixels are considered and, if
    /// their flag is AVAILABLE, their flux is examined. If it's above
    /// the threshold, that pixel is added to the list to be looked at
    /// and their flag is changed to DETECTED. 
    /// @param theObject The duchamp::Detection object to be grown. It
    /// is returned with new pixels in place. Only the basic
    /// parameters that belong to PixelInfo::Object3D are
    /// recalculated.

    size_t spatsize=itsArrayDim[0]*itsArrayDim[1];
    long zero = 0;
    std::vector<Voxel<T> > voxlist = theObject->getPixelSet();
    size_t origSize = voxlist.size();
    long xpt,ypt,zpt;
    long xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
    size_t pos;
    for(size_t i=0; i<voxlist.size(); i++){

        xpt=voxlist[i].getX();
        ypt=voxlist[i].getY();
        zpt=voxlist[i].getZ();
      
        xmin = size_t(max(xpt-itsSpatialThresh, zero));
        xmax = size_t(min(xpt+itsSpatialThresh, long(itsArrayDim[0])-1));
        ymin = size_t(max(ypt-itsSpatialThresh, zero));
        ymax = size_t(min(ypt+itsSpatialThresh, long(itsArrayDim[1])-1));
        zmin = size_t(max(zpt-itsVelocityThresh, zero));
        zmax = size_t(min(zpt+itsVelocityThresh, long(itsArrayDim[2])-1));
      
        //loop over surrounding pixels.
        for(x=xmin; x<=xmax; x++){
            for(y=ymin; y<=ymax; y++){
                for(z=zmin; z<=zmax; z++){
                    pos=x+y*itsArrayDim[0]+z*spatsize;
                    if(((x!=xpt)||(y!=ypt)||(z!=zpt))
                        && itsFlagArray[pos]==AVAILABLE ) {
                            if(itsGrowthStats.isDetection(itsFluxArray[pos])){
                                itsFlagArray[pos]=DETECTED;
                                voxlist.push_back(Voxel<T>(x,y,z));
                            }
                    }

                } 
            }
        } 
    } 

    // Add in new pixels to the Detection
    for(size_t i=origSize; i<voxlist.size(); i++){
        theObject->addPixel(voxlist[i]);
    }
   

}
template void ObjectGrower<short>::grow(Detection<short>*);
template void ObjectGrower<int>::grow(Detection<int>*);
template void ObjectGrower<long>::grow(Detection<long>*);
template void ObjectGrower<float>::grow(Detection<float>*);
template void ObjectGrower<double>::grow(Detection<double>*);


template <class T>
std::vector<Voxel<T> > ObjectGrower<T>::growFromPixel(Voxel<T> &vox) {

    std::vector<Voxel<T> > newVoxels;

    long xpt=vox.getX();
    long ypt=vox.getY();
    long zpt=vox.getZ();
    size_t spatsize=itsArrayDim[0]*itsArrayDim[1];
    long zero = 0;

    int xmin = max(xpt - itsSpatialThresh, zero);
    int xmax = min(xpt + itsSpatialThresh, long(itsArrayDim[0])-1);
    int ymin = max(ypt - itsSpatialThresh, zero);
    int ymax = min(ypt + itsSpatialThresh, long(itsArrayDim[1])-1);
    int zmin = max(zpt - itsVelocityThresh, zero);
    int zmax = min(zpt + itsVelocityThresh, long(itsArrayDim[2])-1);
      
    size_t pos;
    Voxel<T> nvox;
    std::vector<Voxel<T> > morevox;
    for(int x=xmin; x<=xmax; x++){
        for(int y=ymin; y<=ymax; y++){
            for(int z=zmin; z<=zmax; z++){
                pos=x+y*itsArrayDim[0]+z*spatsize;
                if(((x!=xpt)||(y!=ypt)||(z!=zpt))
                    && itsFlagArray[pos]==AVAILABLE ) {
                    if(itsGrowthStats.isDetection(itsFluxArray[pos])){
                        itsFlagArray[pos]=DETECTED;
                        nvox.setXYZF(x,y,z,itsFluxArray[pos]);
                        newVoxels.push_back(nvox);        
                    }
                }
            }
        }
    } 

    typename std::vector<Voxel<T> >::iterator v,v2;
    for(v=newVoxels.begin();v<newVoxels.end();v++) {
        std::vector<Voxel<T> > morevox = growFromPixel(*v);
        for(v2=morevox.begin();v2<morevox.end();v2++) 
            newVoxels.push_back(*v2);
    }

    return newVoxels;

}
template std::vector<Voxel<short> > ObjectGrower<short>::growFromPixel(Voxel<short>&);
template std::vector<Voxel<int> > ObjectGrower<int>::growFromPixel(Voxel<int>&);
template std::vector<Voxel<long> > ObjectGrower<long>::growFromPixel(Voxel<long>&);
template std::vector<Voxel<float> > ObjectGrower<float>::growFromPixel(Voxel<float>&);
template std::vector<Voxel<double> > ObjectGrower<double>::growFromPixel(Voxel<double>&);


