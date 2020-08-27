//----------------------------------------------------------
// rendering3D.hh: Definition of the Rendering3D class.
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

#ifndef RENDERING3D_HH_
#define RENDERING3D_HH_

#include <Arrays/cube.hh>
#include <Arrays/image.hh>
#include <Tasks/smooth3D.hh>
#include <Utilities/progressbar.hh>


template <class T> 
class Rendering3D
{
public:
    Rendering3D<T>(Cube<T> *c) : in(c) {}

    ~Rendering3D () {if (out!=nullptr) delete out;}

    void compute(float azangle=360.);
    void writefits(std::string fname="", float smoothfactor=1.2);

private:
    Cube<T> *in;             /// A pointer to the input datacube.
    Cube<T> *out = nullptr;  /// A cube containing the 3D rendering.

    void rot(double *coord, double *e, double phi, double *rotated);
};


template <class T>
void Rendering3D<T>::compute(float azangle) {

    if (!in->MaskAll()) in->BlankMask();
    if (out!=nullptr) delete out;

    int size[3]={in->DimX(),in->DimY(),in->DimZ()};
    int isize = int(azangle);

    int axis[3]={size[0],size[1],isize};
    out = new Cube<float>(axis);
    out->saveHead(in->Head());
    out->saveParam(in->pars());
    for (int i=0; i<out->NumPix(); i++) out->Array(i)=0;

    double e[3] = {0,1,0};
    int nthreads = in->pars().getThreads();

    ProgressBar bar(true,in->pars().isVerbose(),in->pars().getShowbar());

#pragma omp parallel num_threads(nthreads)
{
    bar.init(" Rendering 3D... ",isize);
#pragma omp	for
    for (int i=0; i<isize; i++) {
        bar.update(i+1);
        int count=0;
        for (int z=0; z<size[2]; z++) {
            for (int y=0; y<size[1]; y++) {
                for (int x=0; x<size[0]; x++) {
                    double coord[3] = {x-size[0]/2.,y-size[1]/2.,z-size[2]/2.};
                    double rotat[3] = {0.,0.,0.};
                    rot(coord,e,i,rotat);
                    int x1 = lround(rotat[0]+size[0]/2.);
                    int y1 = lround(rotat[1]+size[1]/2.);
                    int z1 = lround(rotat[2]+size[2]/2.);
                    long nPix_rot = in->nPix(x1,y1,z1);
                    long nPix3D = out->nPix(x,y,i);
                    if (x1<0 || x1>=size[0]) continue;
                    if (y1<0 || y1>=size[1]) continue;
                    if (z1<0 || z1>=size[2]) continue;

                    if (in->Mask(nPix_rot)) {
                        //if (array[nPix_rot]>array3D[nPix3D])
                        //array3D[nPix3D] = array[nPix_rot];
                        out->Array(nPix3D)+=in->Array(nPix_rot);
                        count++;
                    }
                }
            }
        }

        for (int j=0; j<size[0]*size[1]; j++)
            out->Array(j+i*size[0]*size[1]) /= count;
    }
}

    bar.fillSpace("Done.\n");

    out->Head().setCrpix(2, 0);
    out->Head().setCrval(2, 0);
    out->Head().setCdelt(2, 1);
    out->Head().setCunit(2, "DEGREES");
    out->Head().setCtype(2, "L.O.S. ANGLE");

}


template <class T>
void Rendering3D<T>::writefits(std::string fname, float smoothfactor) {

    if (fname=="") fname = in->pars().getOutfolder()+in->Head().Obname()+"_3D.fits";

    if (smoothfactor>1) {
        Smooth3D<float> *smoothed = new Smooth3D<float>;
        Beam oldbeam = {in->Head().Bmaj()*3600,
                        in->Head().Bmin()*3600,
                        in->Head().Bpa()};
        Beam newbeam = {smoothfactor*in->Head().Bmaj()*3600,
                        smoothfactor*in->Head().Bmin()*3600,
                        in->Head().Bpa()};
        smoothed->setUseBlanks(false);
        smoothed->smooth(out, oldbeam, newbeam, out->Array(), out->Array());
        delete smoothed;
    }

    out->fitswrite_3d(fname.c_str());

}


template <class T>
void Rendering3D<T>::rot(double *coord, double *e, double phi, double *rotated) {

    double x = coord[0];
    double y = coord[1];
    double z = coord[2];

    double cosphi = cos(phi*M_PI/180.);
    double sinphi = sin(phi*M_PI/180.);

    double aa = e[0]*e[0]+(1-e[0]*e[0])*cosphi;
    double ab = (1-cosphi)*e[0]*e[1]-sinphi*e[2];
    double ac = (1-cosphi)*e[0]*e[2]+sinphi*e[1];
    double ba = (1-cosphi)*e[1]*e[0]+sinphi*e[2];
    double bb = e[1]*e[1]+(1-e[1]*e[1])*cosphi;
    double bc = (1-cosphi)*e[1]*e[2]-sinphi*e[0];
    double ca = (1-cosphi)*e[2]*e[0]-sinphi*e[1];
    double cb = (1-cosphi)*e[2]*e[1]+sinphi*e[0];
    double cc = e[2]*e[2]+(1-e[2]*e[2])*cosphi;

    rotated[0] = x*aa+y*ab+z*ac;
    rotated[1] = x*ba+y*bb+z*bc;
    rotated[2] = x*ca+y*cb+z*cc;
}


#endif
