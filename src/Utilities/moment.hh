//----------------------------------------------------------
// mmaps.hh: Definition of the MomentMap class.
//----------------------------------------------------------

/*-----------------------------------------------------------------------
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 Bbarp;p is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with Bbarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning Bbarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#ifndef MOMENTMAP_HH_
#define MOMENTMAP_HH_

#include <iostream>
#include <fitsio.h>

#include "../Arrays/cube.hh"
#include "../Arrays/image.hh"
#include "../Utilities/utils.hh"
#include "../Utilities/progressbar.hh"


template <class T> 
class MomentMap : public Image2D <T> 
{
public:
	MomentMap();
	~MomentMap() {};
    MomentMap(const MomentMap &i);		
    MomentMap& operator=(const MomentMap &i);

	void input (Cube<T> *c, int *Blo, int *Bhi);
	void input (Cube<T> *c);
	void ZeroMoment (bool msk);
	void FirstMoment (bool msk);
	void SecondMoment (bool msk);
	
private:
	Cube<T> *in;
	int	blo[3],bhi[3];
	int	nsubs;

	bool setHead(int type); 
	
};


template <class T>
MomentMap<T>::MomentMap() {
	
	this->numPix = 0;
	this->numAxes = 2;
	this->arrayAllocated = false;
	this->headDefined	   = false;
	this->statsDefined   = false;
	
}


template <class T> 
void MomentMap<T>::input (Cube<T> *c, int *Blo, int *Bhi) {
	
	in = c;
	
	for (int i=0; i<3; i++) {
		blo[i] = Blo[i];
		bhi[i] = Bhi[i];
	}
	
	this->axisDim[0] = bhi[0]-blo[0];
	this->axisDim[1] = bhi[1]-blo[1];
	nsubs = bhi[2]-blo[2]; 
	
	this->numPix = this->axisDim[0]*this->axisDim[1];
	this->array = new T[this->numPix];
	this->arrayAllocated = true;
	
}


template <class T> 
void MomentMap<T>::input (Cube<T> *c) {
	
	in = c;
	
	for (int i=0; i<3; i++) {
		blo[i] = 0;
		bhi[i] = in->AxesDim(i);
	}
	
	this->axisDim[0] = bhi[0]-blo[0];
	this->axisDim[1] = bhi[1]-blo[1];
	nsubs = bhi[2]-blo[2]; 
	
	this->numPix = this->axisDim[0]*this->axisDim[1];
	this->array = new T[this->numPix];
	this->arrayAllocated = true;
	
}


template <class T>
void MomentMap<T>::ZeroMoment (bool msk) {
	
	if (!this->arrayAllocated) {
		std::cout << "MOMENT MAPS error: ";
		std::cout << "Array not allocated. Call 'input' first!!\n";
		abort();
	}
		
	if (!(this->headDefined=setHead(0))) {
		std::cout<< "MOMENT MAPS warning: cannot create new header.\n";
	}
	
	bool isVerbose = in->pars().isVerbose();
	
	if(msk && !in->MaskAll()) in->BlankMask();
	
	float deltaV = 1;
	if (in->HeadDef()) deltaV = fabs(DeltaVel<T>(in->Head()));
	
	ProgressBar bar(" Extracting 0th moment map... ", true);
    bar.setShowbar(in->pars().getShowbar());
	if (isVerbose) bar.init(this->axisDim[0]);

    for (int x=0; x<this->axisDim[0]; x++) {
		if (isVerbose) bar.update(x+1);
		for (int y=0; y<this->axisDim[1]; y++) {
			float fluxsum = 0;
			for (int z=0; z<nsubs; z++) {
				long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
				if (msk)fluxsum += in->Array(npix)*in->Mask(npix);
				else fluxsum += in->Array(npix);
			}
            if (in->HeadDef())
                this->array[x+y*this->axisDim[0]] = FluxtoJy(fluxsum, in->Head())*deltaV;
            else this->array[x+y*this->axisDim[0]] = fluxsum;
		}
	}

	if (isVerbose) bar.fillSpace(" Done.\n");
	
}


template <class T>
void MomentMap<T>::FirstMoment (bool msk) {
	
	
	if (!this->arrayAllocated) {
		std::cout << "MOMENT MAPS error: ";
		std::cout << "Array not allocated. Call 'input' first!!\n";
		abort();
	}
		
	if (!(this->headDefined=setHead(1))) {
		std::cout<< "MOMENT MAPS warning: cannot create new header.\n";
	}
	
	bool isVerbose = in->pars().isVerbose();
	
	if(msk && !in->MaskAll()) in->BlankMask();
			
	ProgressBar bar(" Extracting 1st moment map... ", true);
    bar.setShowbar(in->pars().getShowbar());
	if (isVerbose) bar.init(this->axisDim[0]);
	
	
	for (int x=0; x<this->axisDim[0]; x++) {		
		if (isVerbose) bar.update(x+1);
		for (int y=0; y<this->axisDim[1]; y++) {
			T num = 0;
			T denom = 0;
			for (int z=0; z<nsubs; z++) {
				long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
				T VEL;
				if (in->HeadDef()) {
					double crpix2 = in->Head().Crpix(2);
					double cdelt2 = in->Head().Cdelt(2);
					double crval2 = in->Head().Crval(2);
					T zval = (z+1+blo[2]-crpix2)*cdelt2+crval2;
                    VEL = AlltoVel<T>(zval, in->Head());
				}
				else VEL = z+blo[2];
				if (msk) {
					num += in->Array(npix)*VEL*in->Mask(npix);
					denom += in->Array(npix)*in->Mask(npix);
				}
				else {
					num += in->Array(npix)*VEL;
					denom += in->Array(npix);	
				}
			}
			this->array[x+y*this->axisDim[0]]=num/denom;	
		}
	}
	
	if (isVerbose) bar.fillSpace(" Done.\n");

}


template <class T>
void MomentMap<T>::SecondMoment (bool msk) {
	
	if (!this->arrayAllocated) {
		std::cout << "MOMENT MAPS error: ";
		std::cout << "Array not allocated. Call 'input' first!!\n";
		abort();
	}
		
	if (!(this->headDefined=setHead(2))) {
		std::cout<< "MOMENT MAPS warning: cannot create new header.\n";
	}
	
	bool isVerbose = in->pars().isVerbose();
	
	if(msk && !in->MaskAll()) in->BlankMask();

	ProgressBar bar(" Extracting 2nd moment map... ", true);
    bar.setShowbar(in->pars().getShowbar());
	if (isVerbose) bar.init(this->axisDim[0]);
	
	for (int x=0; x<this->axisDim[0]; x++) {		
		if (isVerbose) bar.update(x+1);
		for (int y=0; y<this->axisDim[1]; y++) {
			T num=0, denom=0;
			T numfrst=0, denomfrst=0;
			T firstmoment=0;
			for (int z=0; z<nsubs; z++) {
				long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
				T VEL;
				if (in->HeadDef()) {
					double crpix2 = in->Head().Crpix(2);
					double cdelt2 = in->Head().Cdelt(2);
					double crval2 = in->Head().Crval(2);
					T zval = (z+1+blo[2]-crpix2)*cdelt2+crval2;
                    VEL = AlltoVel<T>(zval, in->Head());
				}
				else VEL = z+blo[2];
				if (msk) {	
					if (in->Mask(npix)) {
						numfrst += in->Array(npix)*VEL;
						denomfrst += in->Array(npix);				
					}
				}
				else {
					numfrst += in->Array(npix)*VEL;
					denomfrst += in->Array(npix);	
				}
			}
			if(denomfrst!=0) firstmoment = numfrst/denomfrst;
			else firstmoment = 0;
			for (int z=0; z<nsubs; z++) {
				long npix = in->nPix(x+blo[0],y+blo[1],z+blo[2]);
				T VEL;
				if (in->HeadDef()) {
					double crpix2 = in->Head().Crpix(2);
					double cdelt2 = in->Head().Cdelt(2);
					double crval2 = in->Head().Crval(2);
					T zval = (z+1+blo[2]-crpix2)*cdelt2+crval2;
                    VEL = AlltoVel<T>(zval, in->Head());
				}
				else VEL = z+blo[2];
				if (msk) {
					if (in->Mask(npix)) { 
						num += in->Array(npix)*(VEL-firstmoment)*(VEL-firstmoment);
						denom += in->Array(npix);
					}
				}
				else {
					num += in->Array(npix)*(VEL-firstmoment)*(VEL-firstmoment);
					denom += in->Array(npix);
				}
			}
			this->array[x+y*this->axisDim[0]]=sqrt(num/denom);			
		}
	}
	if (isVerbose) bar.fillSpace(" Done.\n");
	
}


template <class T> 
bool MomentMap<T>::setHead(int type) {
	
	if (!in->HeadDef()) return false; 
	else { 
		this->copyHeader(in->Head());
		this->head.setCrpix(0, in->Head().Crpix(0)-blo[0]);
		this->head.setCrpix(1, in->Head().Crpix(1)-blo[1]);
		if (type==0) {				
			this->head.setBtype("intensity");
			if (in->Head().BeamArea()!=0) 
				this->head.setBunit("JY * KM/S");
			else this->head.setBunit("JY/BEAM * KM/S");
		}
		else if (type==1 || type==2) {
			this->head.setBtype("velocity");
			this->head.setBunit("KM/S");
		}
	}	
	
	return true;
	
}


template <class T>
Image2D<T>* PositionVelocity (Cube<T> *c, int x0, int y0, T phi) {
	
	while(phi>=360) phi -= 360;
	while(phi<0) phi += 360;
	
	double P = phi*M_PI/180.;
	int dim[2] = {0, c->DimZ()}; 
	Image2D<T> *pv;
	int xdim=c->DimX(), ydim=c->DimY();
	int xmax=xdim, ymax=ydim;
	int xmin=0, ymin=0; 
	
	if (phi==90 || phi==270) {
		dim[0] = xdim;
		pv  = new Image2D<T>(dim);
		for (int x=0; x<dim[0]; x++) 
			for (int z=0; z<dim[1]; z++) 
				pv->Array()[x+z*dim[0]] = c->Array(c->nPix(x,y0,z));				
	}
	else {	
		std::vector<int> xx, yy; 
		double mx=0, my=0;
		
		if (phi<90) my = tan(P+M_PI_2);
		else if (phi>90 && phi<270) my = tan(P-M_PI_2);
		else if (phi>270) my = tan(P-3*M_PI_2);

		mx = 1./my;
		
		int x_1 = lround(x0+(c->DimY()-y0)/my); 
		int x_0 = lround(x0-y0/my);
		xmax = my>0 ? x_1 : x_0;
		xmin = my>0 ? x_0 : x_1;
		if (xmax<xmin) std::swap(xmax,xmin);
		if (xmax>c->DimX()) xmax=c->DimX();
		if (xmin<0) xmin=0;
		xdim = fabs(xmax-xmin);
				
		int y_1 = lround(y0+(c->DimX()-x0)/mx); 
		int y_0 = lround(y0-x0/mx);
		ymax = mx>0 ? y_1 : y_0;
		ymin = mx>0 ? y_0 : y_1;
		if (ymax<ymin) std::swap(ymax,ymin);
		if (ymax>c->DimY()) ymax=c->DimY();
		if (ymin<0) ymin=0;
		ydim = fabs(ymax-ymin);
		
		if (xdim>=ydim) {
			int nxdim=0;
			for (int x=0; x<c->DimX(); x++) { 
				int y1 = lround(my*(x-x0)+y0);
				bool isin = y1>=0 && y1<c->DimY();	
				if (isin) {
					xx.push_back(x);
					yy.push_back(y1);
					nxdim++;
				}
			}
			dim[0]=nxdim;
		}
		else {
			int nxdim=0;
			for (int y=0; y<c->DimY(); y++) { 
				int x1 = lround(mx*(y-y0)+x0);
				bool isin = x1>=0 && x1<c->DimX();	
				if (isin) {
					xx.push_back(x1);
					yy.push_back(y);
					nxdim++;
				}
			}
			dim[0]=nxdim;
		}

		pv  = new Image2D<T>(dim);	
		for (int i=0; i<dim[0]; i++) { 
			for (int z=0; z<c->DimZ(); z++) {
				pv->Array()[i+z*dim[0]] = c->Array(c->nPix(xx[i],yy[i],z));		
			} 
		}
	}
	
	Header &h = c->Head();
	pv->copyHeader(h);
	float crpix0 = xdim>ydim ? x0-xmin : y0-ymin;
	float xdom = xdim*h.Cdelt(0);
	float ydom = ydim*h.Cdelt(1);
	float crdelt0 = sqrt(xdom*xdom+ydom*ydom)/pv->DimX();
    pv->Head().setCrpix(0, crpix0+1);
	pv->Head().setCrval(0, 0);
	pv->Head().setCdelt(0, crdelt0);
	pv->Head().setCtype(0, "Offset");
	pv->Head().setCrpix(1, h.Crpix(2));
	pv->Head().setCdelt(1, h.Cdelt(2));
	pv->Head().setCrval(1, h.Crval(2));
	pv->Head().setCunit(1, h.Cunit(2));
	pv->Head().setCtype(1, h.Ctype(2));
	pv->Head().setDataMax(0.);
	pv->Head().setDataMin(0.);
	
	std::string name = c->Head().Name()+"_pv"+to_string(phi);
	pv->Head().setName(name);

	
	return pv;
	
}



#endif
