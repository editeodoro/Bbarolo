//------------------------------------------------------------------------
// galfit_min.cpp: Members functions for minimization of the Galfit class.
//------------------------------------------------------------------------

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
#include <cfloat>
#include <cmath>
#include <functional>
#include <Arrays/cube.hh>
#include <Tasks/galmod.hh>
#include <Tasks/galfit.hh>
#include <Utilities/utils.hh>
#include <Utilities/conv2D.hh>
#include <Utilities/allocator.hpp>

namespace Model {

template <class T>
bool Galfit<T>::minimize(Rings<T> *dring, T &minimum, T *pmin, Galmod<T> *modsoFar) {

    /// This function uses the Downhill Simplex Method
    /// in multidimensions due to Nelder and Mead.

    const double TINY = 1.0e-10;
    const double tol  = par.TOL;
    
    int ndim=nfree;
    if (global) ndim=nfree*dring->nr;

    const int NMAX=200*ndim;

    int mpts=ndim+1;
    T **p = allocate_2D<T>(mpts,ndim);

    T *point = new T[ndim];
    T *dels  = new T[ndim];

    /*
    int k=0;
    if (mpar[VROT])  point[k++]=dring->vrot.front();
    if (mpar[VDISP]) point[k++]=dring->vdisp.front();
    if (mpar[DENS])  point[k++]=dring->dens.front();
    if (mpar[Z0])    point[k++]=dring->z0.front();
    if (mpar[INC])   point[k++]=dring->inc.front();
    if (mpar[PA])    point[k++]=dring->phi.front();
    if (mpar[XPOS])  point[k++]=dring->xpos.front();
    if (mpar[YPOS])  point[k++]=dring->ypos.front();
    if (mpar[VSYS])  point[k++]=dring->vsys.front();
    if (mpar[VRAD])  point[k++]=dring->vrad.front();

    /// Determine the initial simplex.
    for (int i=0; i<ndim; i++) {
        dels[i]  = 0.2*point[i];
        point[i] = point[i]-0.1*point[i];					///-<<<<<<<<<<<< Totalmente arbitrario. pensaci su.
    }
    */

    /// Determine the initial simplex (new way).
    //    if (global) {
    //        int n=dring->nr, k=0;
    //        if (mpar[VROT])  for (int j=0;j<n;j++) {dels[k]=30.; point[k++]=dring->vrot[j];}
    //        if (mpar[VDISP]) for (int j=0;j<n;j++) {dels[k]=10.; point[k++]=dring->vdisp[j];}
    //        if (mpar[DENS])  for (int j=0;j<n;j++) {dels[k]=0.2*dring->dens.front(); point[k++]=dring->dens[j];}
    //        if (mpar[Z0])    for (int j=0;j<n;j++) {dels[k]=0.2*dring->z0.front(); point[k++]=dring->z0[j];}
    //        if (mpar[INC])   for (int j=0;j<n;j++) {dels[k]=5.; point[k++]=dring->vexp[j];}
    //        if (mpar[PA])    for (int j=0;j<n;j++) {dels[k]=10.; point[k++]=dring->phi[j];}
    //        if (mpar[XPOS])  for (int j=0;j<n;j++) {dels[k]=(maxs[XPOS]-mins[XPOS])/4.; point[k++]=dring->xpos[j];}
    //        if (mpar[YPOS])  for (int j=0;j<n;j++) {dels[k]=(maxs[YPOS]-mins[YPOS])/4.; point[k++]=dring->ypos[j];}
    //        if (mpar[VSYS])  for (int j=0;j<n;j++) {dels[k]=(maxs[VSYS]-mins[VSYS])/4.; point[k++]=dring->vsys[j];}
    //        if (mpar[VRAD])  for (int j=0;j<n;j++) {dels[k]=(maxs[VRAD]-mins[VRAD])/4.; point[k++]=dring->vrad[j];}
    //    }
    //    else {
    //        int k=0;
    //        if (mpar[VROT])  {dels[k]=30.; point[k++]=dring->vrot.front();}
    //        if (mpar[VDISP]) {dels[k]=10.; point[k++]= dring->vdisp.front();}
    //        if (mpar[DENS])  {dels[k]=0.2*dring->dens.front(); point[k++]= dring->dens.front();}
    //        if (mpar[Z0])    {dels[k]=0.2*dring->z0.front(); point[k++]= dring->z0.front();}
    //        if (mpar[INC])   {dels[k]=5.; point[k++]= dring->inc.front();}
    //        if (mpar[PA])    {dels[k]=10.;point[k++]= dring->phi.front();}
    //        if (mpar[XPOS])  {dels[k]=(maxs[XPOS]-mins[XPOS])/4.;point[k++]= dring->xpos.front();}
    //        if (mpar[YPOS])  {dels[k]=(maxs[YPOS]-mins[YPOS])/4.;point[k++]= dring->ypos.front();}
    //        if (mpar[VSYS])  {dels[k]=(maxs[VSYS]-mins[VSYS])/4.;point[k++]= dring->vsys.front();}
    //        if (mpar[VRAD])  {dels[k]=(maxs[VRAD]-mins[VRAD])/4.;point[k++]= dring->vrad.front();}
    //    }
    int n=1, k=0;
    if (global) n=dring->nr;
    if (mpar[VROT]) {
        for (int j=0;j<n;j++) {
            if (dring->vrot[j]>40) dels[k]=30.;
            else dels[k]=5;
            point[k++]=dring->vrot[j];
        }
    }
    if (mpar[VDISP]) {
        for (int j=0;j<n;j++) {
            //dels[k] = 0.5*dring->vdisp.front();
            if (dring->vdisp[j]>80) dels[k]=20.;
            else if (dring->vdisp[j]>20. && dring->vdisp[j]<=80.) dels[k]=10.;
            else if (dring->vdisp[j]<=2) dels[k]=0.2;
            else dels[k]=5;
            point[k++]=dring->vdisp[j];
        }
    }
    if (mpar[DENS])  for (int j=0;j<n;j++) {dels[k]=0.2*dring->dens.front(); point[k++]=dring->dens[j];}
    if (mpar[Z0])    for (int j=0;j<n;j++) {dels[k]=0.2*dring->z0.front(); point[k++]=dring->z0[j];}
    if (mpar[INC])   for (int j=0;j<n;j++) {dels[k]=5.; point[k++]=dring->inc[j];}
    if (mpar[PA])    for (int j=0;j<n;j++) {dels[k]=10.; point[k++]=dring->phi[j];}
    if (mpar[XPOS])  for (int j=0;j<n;j++) {dels[k]=(maxs[XPOS]-mins[XPOS])/4.; point[k++]=dring->xpos[j];}
    if (mpar[YPOS])  for (int j=0;j<n;j++) {dels[k]=(maxs[YPOS]-mins[YPOS])/4.; point[k++]=dring->ypos[j];}
    if (mpar[VSYS])  for (int j=0;j<n;j++) {dels[k]=(maxs[VSYS]-mins[VSYS])/4.; point[k++]=dring->vsys[j];}
    if (mpar[VRAD])  for (int j=0;j<n;j++) {dels[k]=15; point[k++]=dring->vrad[j];}

    // Following line offsets the initial simplex such that the initial point given is at the centre
    // of the simplex instead of just one of the verteces.
    //for (int j=0; j<ndim; j++) point[j] -= dels[j]/2.;

    // Build the initial matrix.
    for (int i=0; i<mpts; i++) {
        for (int j=0; j<ndim; j++) p[i][j]=point[j];
        if (i!=0) p[i][i-1] += dels[i-1];
        // Following to print initial simplex
        //for (int j=0; j<ndim; j++) std::cout << p[i][j] << " ";
        //std::cout << std::endl;
    }
    
    delete [] point;
    delete [] dels;

    T psum[ndim], x[ndim];
    T *y = new T[mpts];

    for (int i=0; i<mpts; i++) {
        for (int j=0; j<ndim; j++) x[j]=p[i][j];
        y[i]=func3D(dring,x,modsoFar);
    }

    int nfunc=0;
    for (int j=0; j<ndim; j++) {
        T sum=0.0;
        for (int i=0; i<mpts; i++) sum += p[i][j];
        psum[j]=sum;
    }

    // Main cycle begins here.
    for (;;) {
        int ihi, inhi;
        int ilo=0;
        // First determine the highest point (worst) and the
        // lowest (best).
        ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
        for (int i=0; i<mpts; i++) {
            if (y[i]<=y[ilo]) ilo=i;
            if (y[i]>y[ihi]) {
                inhi=ihi;
                ihi=i;
            }
            else if (y[i]>y[inhi] && i!=ihi) inhi=i;
        }

        // Compute the fractional range from highest to lowest.
        double rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);

        // If it is satisfactory, put best point and value in slot 0;
        if (rtol<tol) {
            std::swap(y[0],y[ilo]);
            for (int i=0; i<ndim; i++) {
                std::swap(p[0][i],p[ilo][i]);
                pmin[i]=p[0][i];
            }
            minimum=y[0];
            deallocate_2D(p,ndim+1);
            delete [] y;
            return true;
        }

        if (nfunc>=NMAX) {
            deallocate_2D(p,ndim+1);
            delete [] y;
            return false;
        }
        nfunc += 2;

        // Try a new iteration. Try extrapolate by a factor -1 through
        // the simplex face across from the high point, i.e. reflect the
        // simplex from the high point.
        T ytry=mtry(dring, p,y,psum,ihi,-1.0,modsoFar);
        if (ytry<=y[ilo]) {
            // If it gives a value better than the best point, try an
            // additional extrapolation by a factor 2.
            ytry=mtry(dring, p,y,psum,ihi,2.0,modsoFar);
        }
        else if (ytry>=y[inhi]) {
            // Otherwise, if the reflected point is worse than the second
            // highest, look for an intermediate lower point, i.e. do a
            // on dimensional contraction.
            T ysave=y[ihi];
            ytry=mtry(dring, p,y,psum,ihi,0.5,modsoFar);
            if (ytry>=ysave) {
                // Can't seem to get rid of that high point. Better contract
                // around the lowest (best) point.
                for (int i=0; i<mpts; i++) {
                    if (i!=ilo) {
                        for (int j=0; j<ndim; j++)
                            p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                        y[i]=func3D(dring,psum,modsoFar);
                    }
                }

                nfunc += ndim;

                for (int j=0;j<ndim;j++) {
                    T sum=0.0;
                    for (int i=0; i<mpts; i++) sum += p[i][j];
                    psum[j]=sum;
                }
            }
        }
        else --nfunc;
    }

}	
template bool Galfit<float>::minimize(Rings<float>*,float&,float*,Galmod<float> *);
template bool Galfit<double>::minimize(Rings<double>*,double&,double*,Galmod<double> *);
	
	
template <class T>
T Galfit<T>::mtry(Rings<T> *dring, T **p, T *y, T *psum, const int ihi, const double fac, Galmod<T> *modsoFar) {

    // Extrapolates by a factor fac through the face of the simplex across from
    // the high point, tries it, and replaces the high point if the new point is better.

    int ndim = nfree;
    T ptry[ndim];
    double fac1=(1.0-fac)/ndim;
    double fac2=fac1-fac;

    for (int j=0; j<ndim; j++)
        ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;

    // Evaluate the function at the trial point.
    T ytry=func3D(dring, ptry, modsoFar);
    if (ytry<y[ihi]) {
        // If it is better than the highest, then replace the highest.
        y[ihi]=ytry;
        for (int j=0; j<ndim; j++) {
            psum[j] += ptry[j]-p[ihi][j];
            p[ihi][j]=ptry[j];
        }
    }
    return ytry;
}
template float Galfit<float>::mtry(Rings<float>*,float**,float*,float*,const int,const double,Galmod<float> *);
template double Galfit<double>::mtry(Rings<double>*,double**,double*,double*,const int,const double,Galmod<double> *);


template <class T>
T Galfit<T>::func3D(Rings<T> *dring, T *zpar, Galmod<T> *modsoFar) {

    static std::uniform_real_distribution<float> uniform = std::uniform_real_distribution<float>(0.,1.);
    static std::mt19937 generator(-1);
    static auto fran = std::bind(uniform, generator);
    
//    T vrot  = dring->vrot.front();
//    T vdisp = dring->vdisp.front();
//    T dens  = dring->dens.front();
//    T z0    = dring->z0.front();
//    T inc   = dring->inc.front();
//    T phi   = dring->phi.front();
//    T xpos  = dring->xpos.front();
//    T ypos  = dring->ypos.front();
//    T vsys  = dring->vsys.front();
//    T vrad  = dring->vrad.front();
    

    int n=1, np=0, w_r=dring->id;

    if (global) {n=dring->nr; w_r=0;}
    T vrot[n],vdisp[n],dens[n],z0[n],inc[n],phi[n],xpos[n],ypos[n],vsys[n],vrad[n];

    for (int j=0; j<n; j++) {
      vrot[j]  = dring->vrot[j];
      vdisp[j] = dring->vdisp[j];
      dens[j]  = dring->dens[j];
      z0[j]    = dring->z0[j];
      inc[j]   = dring->inc[j];
      phi[j]   = dring->phi[j];
      xpos[j]  = dring->xpos[j];
      ypos[j]  = dring->ypos[j];
      vsys[j]  = dring->vsys[j];
      vrad[j]  = dring->vrad[j];
    }

    // Following to print current paraemters
    // for (int i=0;i<nfree;i++) std::cout << zpar[i] << " ";
    // std::cout << std::endl;
        
    for (int i=0; i<MAXPAR; i++) {
        if (mpar[i]) {
            switch(i) {
            case VROT:
                for (int j=0; j<n; j++) {
                    double minv = (inr->vrot[w_r]-par.DELTAVROT)>mins[VROT] ? inr->vrot[w_r]-par.DELTAVROT : mins[VROT];
                    double maxv = (inr->vrot[w_r]+par.DELTAVROT)<maxs[VROT] ? inr->vrot[w_r]+par.DELTAVROT : maxs[VROT];
                    if (zpar[np]>maxv) vrot[j]= maxv - fran()*par.DELTAVROT;
                    else if (zpar[np]<minv) vrot[j]= minv + fran()*par.DELTAVROT;
                    else vrot[j] = zpar[np];
                    // Maximum rotation velocity inside the datacube
                    //maxv = fabs(AlltoVel<T>(in->getZphys(in->DimZ()-1),in->Head())-vsys[j]);
                    //minv = fabs(AlltoVel<T>(in->getZphys(0),in->Head())-vsys[j]);
                    //double maxvrot = std::max(maxv/sin(inc[j]*M_PI/180.),minv/sin(inc[j]*M_PI/180.));
                    //if (zpar[np]>maxvrot) vrot[j] = maxvrot-fran()*maxvrot;
                    zpar[np++] = vrot[j];
                }
                break;

            case VDISP:
                for (int j=0; j<n; j++) {
                    if (zpar[np]>maxs[VDISP]) vdisp[j]= maxs[VDISP]-fran()*inr->vdisp[w_r];
                    else if (zpar[np]<mins[VDISP]) vdisp[j]= mins[VDISP]+fran()*inr->vdisp[w_r];
                    else vdisp[j] = zpar[np];
                    //vdisp[j] = (zpar[np]<mins[VDISP] || zpar[np]>maxs[VDISP]) ? inr->vdisp[w_r]+0.1*mins[VDISP] : zpar[np];
                    //vdisp = mins[VDISP]+0.1*mins[VDISP];
                    //vdisp = outr->vdisp[w_r-1];
                    zpar[np++] = vdisp[j];
                }
                break;

            case DENS:
                for (int j=0; j<n; j++) dens[j] = zpar[np++];
                break;

            case Z0:
                for (int j=0; j<n; j++) {
                    z0[j] = (zpar[np]<mins[Z0] || zpar[np]>maxs[Z0]) ? inr->z0[w_r] : zpar[np];
                    zpar[np++] = z0[j];
                }
                break;

            case INC:
                for (int j=0; j<n; j++) {
                    double mini = (inr->inc[w_r]-par.DELTAINC)>mins[INC] ? inr->inc[w_r]-par.DELTAINC : mins[INC];
                    double maxi = (inr->inc[w_r]+par.DELTAINC)<maxs[INC] ? inr->inc[w_r]+par.DELTAINC : maxs[INC];
                    if (zpar[np]>maxi) inc[j]= maxi;//- fran()*par.DELTAINC;
                    else if (zpar[np]<mini) inc[j]= mini;// + fran()*par.DELTAINC;
                    else inc[j] = zpar[np];
                    zpar[np++] = inc[j];
                }
                break;

            case PA:
                for (int j=0; j<n; j++) {
                    double minp = (inr->phi[w_r]-par.DELTAPHI)>mins[PA] ? inr->phi[w_r]-par.DELTAPHI : mins[PA];
                    double maxp = (inr->phi[w_r]+par.DELTAPHI)<maxs[PA] ? inr->phi[w_r]+par.DELTAPHI : maxs[PA];
                    if (zpar[np]>maxp) phi[j]= maxp;//- fran()*par.DELTAPHI;
                    else if (zpar[np]<minp) phi[j]= minp;//+ fran()*par.DELTAPHI;
                    else phi[j] = zpar[np];
                    zpar[np++] = phi[j];
                }
                break;

            case XPOS:
                for (int j=0; j<n; j++) {
                    xpos[j] = (zpar[np]<mins[XPOS] || zpar[np]>maxs[XPOS]) ? inr->xpos[w_r] : zpar[np];
                    zpar[np++] = xpos[j];
                }
                break;

            case YPOS:
                for (int j=0; j<n; j++) {
                    ypos[j] = (zpar[np]<mins[YPOS] || zpar[np]>maxs[YPOS]) ? inr->ypos[w_r] : zpar[np];
                    zpar[np++] = ypos[j];
                }
                break;

            case VSYS:
                for (int j=0; j<n; j++) {
                    vsys[j] = (zpar[np]<mins[VSYS] || zpar[np]>maxs[VSYS]) ? inr->vsys[w_r] : zpar[np];
                    zpar[np++] = vsys[j];
                }
                break;
            
            case VRAD:
                for (int j=0; j<n; j++) {
                    vrad[j] = (zpar[np]<mins[VRAD] || zpar[np]>maxs[VRAD]) ? inr->vrad[w_r] : zpar[np];
                    zpar[np] = vrad[j];
                }
                    break;

            default:
                break;
            }
        }
    }

    /*
    bool out_of_boundaries=false;
    for (int i=0; i<MAXPAR; i++) {
        if (mpar[i]) {
            switch(i) {
                case VROT:
                    vrot = zpar[npar];
                    if (vrot<mins[VROT] || vrot>maxs[VROT]) out_of_boundaries = true;
                    npar++;
                    break;

                case VDISP:
                    vdisp = zpar[npar];
                    if (vdisp<mins[VDISP] || vdisp>maxs[VDISP]) out_of_boundaries = true;
                    npar++;
                    break;

                case DENS:
                    dens = zpar[npar];
                    npar++;
                    break;

                case Z0:
                    z0 = zpar[npar];
                    if (z0<mins[Z0] || z0>maxs[Z0]) out_of_boundaries = true;
                    npar++;
                    break;

                case INC:
                    inc = zpar[npar];
                    if (inc<mins[INC] || inc>maxs[INC]) out_of_boundaries = true;
                    npar++;
                    break;

                case PA:
                    phi = zpar[npar];
                    if (phi<mins[PA] || phi>maxs[PA]) out_of_boundaries = true;
                    npar++;
                    break;

                case XPOS:
                    xpos = zpar[npar];
                    if (xpos<mins[XPOS] || xpos>maxs[XPOS]) out_of_boundaries = true;
                    npar++;
                    break;

                case YPOS:
                    ypos = zpar[npar];
                    if (ypos<mins[YPOS] || ypos>maxs[YPOS]) out_of_boundaries = true;
                    npar++;
                    break;

                case VSYS:
                    vsys = zpar[npar];
                    if (vsys<mins[VSYS] || vsys>maxs[VSYS]) out_of_boundaries = true;
                    npar++;
                    break;
                
                case VRAD:
                    vrad = zpar[npar];
                    if (vrad<mins[VRAD] || vrad>maxs[VRAD]) out_of_boundaries = true;
                    npar++;
                    break;

                    default:
                    break;
            }
        }
    }

    if (out_of_boundaries) return FLT_MAX;
*/


    for (int i=0,j=0; i<dring->nr; i++) {
        j=i*global;
        dring->vrot[i] 	= vrot[j];
        dring->vdisp[i]	= vdisp[j];
        dring->z0[i]	= z0[j];
        dring->dens[i]	= dens[j];
        dring->inc[i]	= inc[j];
        dring->phi[i]	= phi[j];
        dring->xpos[i]	= xpos[j];
        dring->ypos[i]	= ypos[j];
        dring->vsys[i]	= vsys[j];
        dring->vrad[i]	= vrad[j];
    }

    return getFuncValue(dring,modsoFar);

}
template float Galfit<float>::func3D(Rings<float>*,float*,Galmod<float> *);
template double Galfit<double>::func3D(Rings<double>*,double*,Galmod<double> *);


template <class T>
T Galfit<T>::getFuncValue(Rings<T> *dring, Galmod<T> *modsoFar) {

    // Getting the sizes of model cube based on last ring
    int blo[2], bhi[2], bsize[2];
    if (reverse) getModelSize(outr,blo,bhi,bsize);
    else getModelSize(dring,blo,bhi,bsize);

    int nv = par.NV<0 ? in->DimZ() : par.NV;

    // Calculating the model
    Model::Galmod<T> *mod = new Model::Galmod<T>;
    mod->input(in,bhi,blo,dring,nv,par.LTYPE,1,par.CDENS);
    mod->calculate();

    // Adding up the "sofar" model, if requested
    T *modp = mod->Out()->Array();
    if (modsoFar!=nullptr) {
        for (auto i=mod->Out()->NumPix(); i--;) modp[i] += modsoFar->Out()->Array()[i];
    }
    
    //<<<<< Convolution....
    if (par.SM) {
        if (in->pars().getflagFFT()) Convolve_fft(modp, bsize);
        else Convolve(modp, bsize);
    }
    
    //<<<<< Normalizing & calculating the residuals....
    double minfunc = (this->*func_norm)(dring,modp,bhi,blo);
    
    delete mod;
	return minfunc; 
}
template float Galfit<float>::getFuncValue(Rings<float>*,Galmod<float> *);
template double Galfit<double>::getFuncValue(Rings<double>*, Galmod<double> *);


template <class T>
void Galfit<T>::Convolve(T *array, int *bsize) {

    if (cfieldAllocated) {
        int ndx = (bsize[0]+NconX-1);
        int ndy = (bsize[1]+NconY-1);
        T *beforeCON = new T[ndx*ndy];
        T *afterCON  = new T[ndx*ndy];
        for (int z=0; z<in->DimZ(); z++) {
            for (int x=0; x<ndx; x++) {
                for (int y=0; y<ndy; y++) {
                    long nPix = x+y*ndx;
                    int mXpos = x-(NconX-1)/2;
                    int mYpos = y-(NconY-1)/2;
                    long mPix = mXpos+mYpos*bsize[0]+z*bsize[0]*bsize[1];
                    afterCON[nPix] = beforeCON[nPix] = 0;
                    if (x>=(NconX-1)/2 && x<=(bsize[0]+(NconX-1)/2) &&
                            y>=(NconY-1)/2 && y<=(bsize[1]+(NconY-1)/2)) {
                        beforeCON[nPix] = array[mPix];
                    }
                }
            }

            for (int yc=0; yc<NconY; yc++) {
                for (int xc=0; xc<NconX; xc++) {
                    T cf = cfield[NconX-xc-1+(NconY-yc-1)*NconX];
                    if (cf!=0.0) {
                        for (int y=0; y<bsize[1]; y++) {
                            T *v1 = &beforeCON[(yc+y)*ndx+xc];
                            T *v2 = &afterCON[(NconY/2+y)*ndx+NconX/2];
                            for (int x=0; x<bsize[0]; x++) v2[x] = v2[x]+cf*v1[x];
                        }
                    }
                }
            }

            for (int x=(NconX-1)/2; x<(bsize[0]+(NconX-1)/2); x++) {
                for (int y=(NconY-1)/2; y<(bsize[1]+(NconY-1)/2); y++) {
                    long nPix = x+y*(bsize[0]+NconX-1);
                    long mPix = (x-(NconX-1)/2)+(y-(NconY-1)/2)*bsize[0]+z*bsize[0]*bsize[1];
                    //if (IsIn(x-(NconX-1)/2,y-(NconY-1)/2,blo,dring))
                    array[mPix] = afterCON[nPix];
                    //else modp[mPix] = 0;
                }
            }

        }

        delete [] beforeCON;
        delete [] afterCON;
    }

}
template void Galfit<float>::Convolve(float*,int*);
template void Galfit<double>::Convolve(double*,int*);


template <class T>
void Galfit<T>::Convolve_fft(T *array, int *bsize) {

    if (cfieldAllocated) {
        long size = bsize[0]*bsize[1];
        double *beforeCON = new double[size];
        Conv2D cfft;
        init_Conv2D(cfft,LINEAR_SAME, bsize[0], bsize[1], NconX, NconY);
        
        for (uint z=in->DimZ(); z--;) {
            T *ptr = &array[z*size];
            for (uint i=size; i--;) beforeCON[i] = isNaN(ptr[i]) ? 0 : ptr[i];
            convolve (cfft,beforeCON, cfield);
            for (uint i=size; i--;)
                ptr[i] = (cfft.dst[i]<1.E-12) ? 0. : cfft.dst[i];	//<<<< Un po' arbitrario, non mi piace.
        }

        clear_Conv2D(cfft);
        delete [] beforeCON;
    }

}
template void Galfit<float>::Convolve_fft(float*,int*);
template void Galfit<double>::Convolve_fft(double*,int*);

///*
// ANELLO 2D
template <class T>
double Galfit<T>::norm_local (Rings<T> *dring, T *array, int *bhi, int *blo) {

    int bsize[2] = {bhi[0]-blo[0], bhi[1]-blo[1]};

    int numPix_ring=0, numBlanks=0, numPix_tot=0;
    double minfunc = 0;
    int bweight = par.BWEIGHT;
    
    // Weighting function can be either a cos(theta)^n or a sin(theta)^n
    double (*wfunc)(double) = cos;
    if (wpow<0) wfunc = sin;
    
    for (uint y=bsize[1]; y--;) {
        for (uint x=bsize[0]; x--;) {
            double theta;
            if (!IsIn(x,y,blo,dring,theta)) continue;
            if (!getSide(theta)) continue;
            numPix_ring++;

            //< Factor for normalization.
            T modSum=0, obsSum = 0, factor=0;
            for (uint z=in->DimZ(); z--;) {
                long modPix = x+y*bsize[0]+z*bsize[0]*bsize[1];
                long obsPix = in->nPix(x+blo[0],y+blo[1],z);
                modSum += array[modPix];
                if (in->Array(obsPix)>0) obsSum += in->Array(obsPix)*mask[obsPix];
            }
            if (modSum!=0) factor = obsSum/modSum;
            else factor=0;

            //double costh = fabs(cos(theta*M_PI/180.));
            double wf = fabs(wfunc(theta*M_PI/180.));
            double wi = std::pow(wf, fabs(wpow));

            // Normalizing and residuals.
            for (uint z=in->DimZ(); z--;) {
                long modPix = x+y*bsize[0]+z*bsize[0]*bsize[1];
                long obsPix = in->nPix(x+blo[0],y+blo[1],z);
                array[modPix] *= factor;
                T obs = in->Array(obsPix)>0 ? in->Array(obsPix) : data_noise;
                T mod = array[modPix];

                if (mask[obsPix]==0) {
                    if (mod==0) continue;
                    else numBlanks++;
                }

                numPix_tot++;
                minfunc += getResValue(obs,mod,wi,chan_noise[z]);
            }
        }
    }
    //numPix_ring=1;
    //return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/numPix_ring;
    return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/((numPix_tot-numBlanks));

}
template double Galfit<float>::norm_local(Rings<float>*,float*,int*,int*);
template double Galfit<double>::norm_local(Rings<double>*,double*,int*,int*);
//*/

/*
//ANELLO 3D
template <class T>
double Galfit<T>::norm_local (Rings<T> *dring, T *array, int *bhi, int *blo) {

    std::cout << "RING3D_LOCAL" << std::endl;
    int bsize[2] = {bhi[0]-blo[0], bhi[1]-blo[1]};

    std::vector<Pixel<T> > *anulus = getRingRegion(dring, bhi, blo);

    int numPix_ring=0, numBlanks=0, numPix_tot=0;
    double minfunc = 0;
    int bweight = par.BWEIGHT;

    typename std::vector<Pixel<T> >::iterator pix;
    for(pix=anulus->begin();pix<anulus->end();pix++) {
        numPix_ring++;
        long x = pix->getX();
        long y = pix->getY();
        T theta = pix->getF();
        if (!getSide(theta)) continue;

        //< Factor for normalization.
        T modSum=0, obsSum = 0, factor=0;
        for (uint z=in->DimZ(); z--;) {
            long modPix = x+y*bsize[0]+z*bsize[0]*bsize[1];
            long obsPix = in->nPix(x+blo[0],y+blo[1],z);
            modSum += array[modPix];
            if (in->Array(obsPix)>0) obsSum += in->Array(obsPix)*mask[obsPix];
        }
        if (modSum!=0) factor = obsSum/modSum;
        else factor=0;

        double costh = fabs(cos(theta*M_PI/180.));
        double wi = std::pow(costh, double(wpow));

        // Normalizing and residuals.
        for (uint z=in->DimZ(); z--;) {
            long modPix = x+y*bsize[0]+z*bsize[0]*bsize[1];
            long obsPix = in->nPix(x+blo[0],y+blo[1],z);
            array[modPix] *= factor;
            T obs = in->Array(obsPix)>0 ? in->Array(obsPix) : data_noise;
            T mod = array[modPix];

            if (mask[obsPix]==0 && mod==0) continue;
            else if (mask[obsPix]==0 && mod!=0) {
                numBlanks++;
                obs = data_noise;
            }

            numPix_tot++;
            minfunc += getResValue(obs,mod,wi,chan_noise[z]);

        }
    }

    //numPix_ring=1;
    //return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/numPix_ring;
    return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/((numPix_tot-numBlanks));

}
template double Galfit<float>::norm_local(Rings<float>*,float*,int*,int*);
template double Galfit<double>::norm_local(Rings<double>*,double*,int*,int*);
*/

///*
//ANELLO 2D
template <class T>
double Galfit<T>::norm_azim (Rings<T> *dring, T *array, int *bhi, int *blo) {

    int bsize[2] = {bhi[0]-blo[0], bhi[1]-blo[1]};

    int numPix_ring=0, numBlanks=0, numPix_tot=0;
    double minfunc = 0;
    int bweight = par.BWEIGHT;
    
    // Weighting function can be either a cos(theta)^n or a sin(theta)^n
    double (*wfunc)(double) = cos;
    if (wpow<0) wfunc = sin;

    //< Factor for normalization.
    T obsSum=0, modSum=0, factor=1;
    for (uint y=bsize[1]; y--;) {
        for (uint x=bsize[0]; x--;) {
            double theta;
            if (!IsIn(x,y,blo,dring,theta)) continue;
            if (!getSide(theta)) continue;

            for (uint z=in->DimZ(); z--;) {
                long modPix = x+y*bsize[0]+z*bsize[0]*bsize[1];
                long obsPix = in->nPix(x+blo[0],y+blo[1],z);
                modSum += array[modPix];
                obsSum += in->Array(obsPix)*mask[obsPix];
            }
        }
    }

    if (modSum!=0) factor = obsSum/modSum;
    else factor=0;
    
    for (uint y=bsize[1]; y--;) {
        for (uint x=bsize[0]; x--;) {
            double theta;
            if (!IsIn(x,y,blo,dring,theta)) continue;
            if (!getSide(theta)) continue;
            numPix_ring++;
            //double costh = fabs(cos(theta*M_PI/180.));
            double wf = fabs(wfunc(theta*M_PI/180.));
            double wi = std::pow(wf, fabs(wpow));

            // Normalizing and residuals.
            for (uint z=in->DimZ(); z--;) {
                long modPix = x+y*bsize[0]+z*bsize[0]*bsize[1];
                long obsPix = in->nPix(x+blo[0],y+blo[1],z);
                array[modPix] *= factor;
                T obs = in->Array(obsPix)>0 ? in->Array(obsPix) : data_noise;
                T mod = array[modPix];

                if (mask[obsPix]==0) {
                    if (mod==0) continue;
                    else numBlanks++;
                }

                numPix_tot++;
                minfunc += getResValue(obs,mod,wi,chan_noise[z]);
            }
        }
    }
    
    //numPix_ring=1;
    //return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/numPix_ring;
    return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/((numPix_tot-numBlanks));

}
template double Galfit<float>::norm_azim(Rings<float>*,float*,int*,int*);
template double Galfit<double>::norm_azim(Rings<double>*,double*,int*,int*);
//*/

/*
//ANELLO 3D
template <class T>
double Galfit<T>::norm_azim (Rings<T> *dring, T *array, int *bhi, int *blo) {

    std::cout << "RING3D_AZIM" << std::endl;
    int bsize[2] = {bhi[0]-blo[0], bhi[1]-blo[1]};

    // Weighting function can be either a cos(theta)^n or a sin(theta)^n
    double (*wfunc)(double) = cos;
    if (wpow<0) wfunc = sin;

    std::vector<Pixel<T> > *anulus = getRingRegion(dring, bhi, blo);

    int numPix_ring=0, numBlanks=0, numPix_tot=0;
    double minfunc = 0;
    int bweight = par.BWEIGHT;

    //< Factor for normalization.
    T modSum=0, obsSum = 0, factor=1;
    typename std::vector<Pixel<T> >::iterator pix;
    for(pix=anulus->begin();pix<anulus->end();pix++) {
        numPix_ring++;
        long x = pix->getX();
        long y = pix->getY();
        T theta = pix->getF();
        if (!getSide(theta)) continue;

        for (uint z=in->DimZ(); z--;) {
            long modPix = x+y*bsize[0]+z*bsize[0]*bsize[1];
            long obsPix = in->nPix(x+blo[0],y+blo[1],z);
            modSum += array[modPix];
            obsSum += in->Array(obsPix)*mask[obsPix];
        }
        if (modSum!=0) factor = obsSum/modSum;
        else factor=0;
    }
    
    if (modSum!=0) factor = obsSum/modSum;
    else factor=0;
    
    for(pix=anulus->begin();pix<anulus->end();pix++) {
        long x = pix->getX();
        long y = pix->getY();
        T theta = pix->getF();
        if (!getSide(theta)) continue;
        
        double wf = fabs(wfunc(theta*M_PI/180.));
        double wi = std::pow(wf, fabs(wpow));
        
        // Normalizing and residuals.
        for (uint z=in->DimZ(); z--;) {
            long modPix = x+y*bsize[0]+z*bsize[0]*bsize[1];
            long obsPix = in->nPix(x+blo[0],y+blo[1],z);
            array[modPix] *= factor;
            T obs = in->Array(obsPix)>0 ? in->Array(obsPix) : data_noise;
            T mod = array[modPix];

            if (mask[obsPix]==0) {
                if (mod==0) continue;
                else numBlanks++;
            }

            numPix_tot++;
            minfunc += getResValue(obs,mod,wi,chan_noise[z]);
        }
    }

    //numPix_ring=1;
    //return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/numPix_ring;
    return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/((numPix_tot-numBlanks));

}
template double Galfit<float>::norm_azim(Rings<float>*,float*,int*,int*);
template double Galfit<double>::norm_azim(Rings<double>*,double*,int*,int*);
*/


template <class T>
double Galfit<T>::norm_none (Rings<T> *dring, T *array, int *bhi, int *blo) {

    int bsize[2] = {bhi[0]-blo[0], bhi[1]-blo[1]};

    int numPix_ring=0, numBlanks=0, numPix_tot=0;
    double minfunc = 0;
    int bweight = par.BWEIGHT;

    // Weighting function can be either a cos(theta)^n or a sin(theta)^n
    double (*wfunc)(double) = cos;
    if (wpow<0) wfunc = sin;
    
    for (uint y=bsize[1]; y--;) {
        for (uint x=bsize[0]; x--;) {
            double theta;
            if (!IsIn(x,y,blo,dring,theta)) continue;
            if (!getSide(theta)) continue;

            numPix_ring++;

            //double costh = fabs(cos(theta*M_PI/180.));
            double wf = fabs(wfunc(theta*M_PI/180.));
            double wi = std::pow(wf, fabs(wpow));
            
            // Normalizing and residuals.
            for (uint z=in->DimZ(); z--;) {
                long modPix = x+y*bsize[0]+z*bsize[0]*bsize[1];
                long obsPix = in->nPix(x+blo[0],y+blo[1],z);
                T obs = in->Array(obsPix)>0 ? in->Array(obsPix) : data_noise;
                T mod = array[modPix];

                if (mask[obsPix]==0) {
                    if (mod==0) continue;
                    else numBlanks++;
                }

                numPix_tot++;
                minfunc += getResValue(obs,mod,wi,chan_noise[z]);
            }
        }
    }
    //numPix_ring=1;
    //return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/numPix_ring;
    return std::pow((1+numBlanks/T(numPix_tot)),bweight)*minfunc/((numPix_tot-numBlanks));

}
template double Galfit<float>::norm_none(Rings<float>*,float*,int*,int*);
template double Galfit<double>::norm_none(Rings<double>*,double*,int*,int*);


template <class T>
void Galfit<T>::getModelSize(Rings<T> *dring, int *blo, int *bhi,int* bsize) {
    
    // Calculating sizes of model cube (only extends to the last ring)
    int xdis = ceil((dring->radii.back()+3*dring->z0.back())/(fabs(in->Head().Cdelt(0))*arcconv));
    int ydis = ceil((dring->radii.back()+3*dring->z0.back())/(fabs(in->Head().Cdelt(1))*arcconv));
    int xpos = ceil(dring->xpos.back());
    int ypos = ceil(dring->ypos.back());
    blo[0] = xpos-xdis;
    blo[1] = ypos-ydis;
    bhi[0] = xpos+xdis+1;
    bhi[1] = ypos+ydis+1;
    if (blo[0]<0) blo[0] = 0;
    if (blo[1]<0) blo[1] = 0;
    if (bhi[0]>in->DimX()) bhi[0] = in->DimX();	
    if (bhi[1]>in->DimY()) bhi[1] = in->DimY();
    bsize[0] = bhi[0]-blo[0];
    bsize[1] = bhi[1]-blo[1];	

}
template void Galfit<float>::getModelSize(Rings<float>*,int*,int*,int*);
template void Galfit<double>::getModelSize(Rings<double>*,int*,int*,int*);


template <class T> 
bool Galfit<T>::IsIn (int x, int y, int *blo, Rings<T> *dr, double &th) {

    // This function verifies that we are inside the rings.
    // Return also the value of azimutal angle th of (x,y) coordinates.

    double F = M_PI/180.;
    double pixScale = ((fabs(in->Head().Cdelt(0))*arcconv)+(fabs(in->Head().Cdelt(1))*arcconv))/2.;
    T inc = dr->inc.back();
    T phi = dr->phi.back();
    double r1 = dr->radii.front()/pixScale;
    double r2 = dr->radii.back()/pixScale;
    double x0 = dr->xpos.back()-blo[0];
    double y0 = dr->ypos.back()-blo[1];
    double xr =  -(x-x0)*sin(F*phi)+(y-y0)*cos(F*phi);
    double yr = (-(x-x0)*cos(F*phi)-(y-y0)*sin(F*phi))/cos(F*inc);
    double r = sqrt(xr*xr+yr*yr);
    if (r<0.1) th = 0.0;
    else th = atan2(yr, xr)/F;
    return r>=r1 && r<=r2;

}
template bool Galfit<float>::IsIn(int,int,int*,Rings<float>*,double&);
template bool Galfit<double>::IsIn(int,int,int*,Rings<double>*,double&);


template <class T>
inline bool Galfit<T>::getSide (double theta) {

    if (par.SIDE=="R") return (fabs(theta)<=90.0);
    else if (par.SIDE=="A") return (fabs(theta)>=90.0);
    else return true;

}
template bool Galfit<float>::getSide(double);
template bool Galfit<double>::getSide(double);


template <class T>
inline double Galfit<T>::getResValue(T obs, T mod, double weight, double noise_weight) {

    double value = 0;
    switch(par.FTYPE) {
        case 1:
            value = weight*std::pow(mod-obs,2)/std::sqrt(obs)/noise_weight;
            break;
        case 2:
            value = weight*fabs(mod-obs)/noise_weight;
            break;
        case 3:
            value = weight*fabs(mod-obs)/(mod+obs)/noise_weight;
            break;
        case 4:
            value = weight*std::pow(mod-obs,2)/noise_weight;
            break;
        default:
            break;
    }
    return value;

}
template double Galfit<float>::getResValue(float,float,double,double);
template double Galfit<double>::getResValue(double,double,double,double);


template <class T>
std::vector<Pixel<T> >* Galfit<T>::getRingRegion (Rings<T> *dring, int *bhi, int *blo) {

    std::vector<Pixel<T> > *anulus= new std::vector<Pixel<T> >;
    long bsize[2] = {bhi[0]-blo[0], bhi[1]-blo[1]};
    bool *an = new bool[bsize[0]*bsize[1]];
    for (int i=0;i<bsize[0]*bsize[1];i++) an[i]=false;

    T R1  = dring->radii.front()/(in->Head().PixScale()*arcconv);
    T R2  = dring->radii.back()/(in->Head().PixScale()*arcconv);
    T phi = dring->phi.back();
    T inc = dring->inc.back();
    T psi = 0.;
    T z0  = dring->z0.back()/(in->Head().PixScale()*arcconv); ///prima prendevo 3*dring->....
    T x0  = dring->xpos.back()-blo[0]-1;
    T y0  = dring->ypos.back()-blo[1]-1;

    double **matrices = RotMatrices(inc,psi,-phi-90);
    int size[2] = {3,3};
    double *rotmatrix = MatrixProduct(&matrices[2][0], size, &matrices[0][0],size);

    int xyrange = lround(R2);
    int zrange = lround(z0);
    int sizecoord[2] = {3,1};
    for (int z=-zrange; z<=zrange; z++) {
        for (int y=-xyrange; y<=xyrange; y++) {
            for(int x=-xyrange; x<=xyrange; x++) {
                double r = sqrt(x*x+y*y);
                if (r<=R2 && r>=R1) {
                    double coord[3]={double(x),double(y),double(z)};
                    double *coordrot = MatrixProduct(rotmatrix,size,coord,sizecoord);
                    int xrot = lround(coordrot[0]+x0);
                    int yrot = lround(coordrot[1]+y0);
                    if (xrot>=0 && xrot<bsize[0] &&
                            yrot>=0 && yrot<bsize[1]) {
                        double theta;
                        if (r<0.1) theta = 0.0;
                        else theta = atan2(y, x)/M_PI*180.;
                        if(!an[xrot+yrot*bsize[0]]) {
                            an[xrot+yrot*bsize[0]] = true;
                            Pixel<T> pix(xrot,yrot,theta);
                            anulus->push_back(pix);
                        }
                    }
                }
            }
        }
    }

    delete [] an;
    deallocate_2D<double>(matrices,3);
    delete [] rotmatrix;
    return anulus;

}
template std::vector<Pixel<float> >* Galfit<float>::getRingRegion (Rings<float>*,int*,int*);
template std::vector<Pixel<double> >* Galfit<double>::getRingRegion (Rings<double>*,int*,int*);


}


