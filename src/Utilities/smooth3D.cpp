//---------------------------------------------------------------
// smooth3D.cpp: Member functions of the Smooth3D class.
//---------------------------------------------------------------

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

#include <iostream>
#include <cstring>
#include <iomanip>
#include "../Arrays/cube.hh"
#include "smooth3D.hh"
#include "utils.hh"
#include "progressbar.hh"
#include "conv2D.hh"

#define BLANK 0xff800000	

template <class T>
void Smooth3D<T>::defaults() {
	
	cutoffratio=1.0/10000.0;
	maxconv=512*512;
	arrayAllocated 	= false;
	confieAllocated	= false;
	beamDefined 	= false;
	usescalefac		= true;
	fft				= true;
	
}
template void Smooth3D<short>::defaults();
template void Smooth3D<int>::defaults();
template void Smooth3D<long>::defaults();
template void Smooth3D<float>::defaults();
template void Smooth3D<double>::defaults();


template <class T>
Smooth3D<T>::Smooth3D() {
	
	defaults();
}
template Smooth3D<short>::Smooth3D();
template Smooth3D<int>::Smooth3D();
template Smooth3D<long>::Smooth3D();
template Smooth3D<float>::Smooth3D();
template Smooth3D<double>::Smooth3D();


template <class T>
Smooth3D<T>::~Smooth3D() {
	
	if (arrayAllocated) delete [] array;	
	if (confieAllocated) delete [] confie;		
}
template Smooth3D<short>::~Smooth3D();
template Smooth3D<int>::~Smooth3D();
template Smooth3D<long>::~Smooth3D();
template Smooth3D<float>::~Smooth3D();
template Smooth3D<double>::~Smooth3D();


template <class T>
Smooth3D<T>::Smooth3D(const Smooth3D<T> &s) {

	operator=(s);
}
template Smooth3D<short>::Smooth3D(const Smooth3D<short> &);
template Smooth3D<int>::Smooth3D(const Smooth3D<int> &);
template Smooth3D<long>::Smooth3D(const Smooth3D<long> &);
template Smooth3D<float>::Smooth3D(const Smooth3D<float> &);
template Smooth3D<double>::Smooth3D(const Smooth3D<double> &);


template <class T>
Smooth3D<T>& Smooth3D<T>::operator=(const Smooth3D<T> &s) {
	
	if(this==&s) return *this;
	
	this->in	 = s.in;
	this->NdatX = s.NdatX;
	this->NdatY = s.NdatY;
	this->NdatZ = s.NdatZ;
	this->NconX = s.NconX;
	this->NconY = s.NconY;
	this->cutoffratio = s.cutoffratio;
	this->maxconv = s.maxconv;
	this->crota	  = s.crota;
	this->oldbeam = s.oldbeam;
	this->newbeam = s.newbeam;
	this->fft	  = s.fft;
	
	this->beamDefined = s.beamDefined;
	if (beamDefined) {
		this->conbeam = s.conbeam;
		this->scalefac= s.scalefac;
		this->usescalefac= s.usescalefac;
	}
	
	for (int i=0; i<2; i++) {
		this->bhi[i] = s.bhi[i];
		this->blo[i] = s.blo[i];
		this->fhi[i] = s.fhi[i];
		this->flo[i] = s.fhi[i];
		this->dimAxes[i] = s.dimAxes[i];
		this->gridspace[i]=s.gridspace[i];
	}
	
	if (this->arrayAllocated) delete [] array;
	this->arrayAllocated = s.arrayAllocated;
	if (this->arrayAllocated) {
		this->array = new T[NdatX*NdatY*NdatZ];
		for (int i=0; i<NdatX*NdatY*NdatZ; i++)
			this->array[i] = s.array[i];
	}

	if (this->confieAllocated) delete [] confie;
	this->confieAllocated = s.confieAllocated;
	if (this->confieAllocated) {
		this->confie = new double[maxconv];
		for (int i=0; i<maxconv; i++)
			this->confie[i] = s.confie[i];
	}
	
    return *this;
}
template Smooth3D<short>& Smooth3D<short>::operator=(const Smooth3D<short> &);
template Smooth3D<int>& Smooth3D<int>::operator=(const Smooth3D<int> &);
template Smooth3D<long>& Smooth3D<long>::operator=(const Smooth3D<long> &);
template Smooth3D<float>& Smooth3D<float>::operator=(const Smooth3D<float> &);
template Smooth3D<double>& Smooth3D<double>::operator=(const Smooth3D<double> &);


template <class T>
void Smooth3D<T>::cubesmooth(Cube<T> *c) {
		
	std::string cunit = c->Head().Cunit(0);
	double unittoarc=0;
	if (cunit=="DEGREE" || cunit=="DEGREES" || cunit=="DEG" ||
		cunit=="degree" || cunit=="degrees" || cunit=="deg") 
			unittoarc = 3600;
	else if (cunit=="ARCSEC" || cunit=="ARCS" ||
			 cunit=="arcsec" || cunit=="arcs") 
			unittoarc = 1;
	else if (cunit=="ARCMIN" || cunit=="ARCM" ||
			 cunit=="arcmin" || cunit=="arcm") 
			unittoarc = 60;
	else {
		std::cout << "SMOOTH error (unknown CUNIT for RA-DEC): ";
		std::cout << "cannot convert to ARCSEC.\n";
		std::cout << cunit;
		abort();
	}
	
    Beam OB = {c->pars().getOBmaj(), c->pars().getOBmin(), c->pars().getOBpa()};


    if (OB.bmaj==-1) OB.bmaj = c->Head().Bmaj()*unittoarc;
    if (OB.bmin==-1) OB.bmin = c->Head().Bmin()*unittoarc;
    if (OB.bpa==0) OB.bpa = c->Head().Bpa();


    if (OB.bmaj==0 || OB.bmin==0) {
        std::cout << "SMOOTH warning: BMAJ, BMIN keywords are not in the FITS header."
                  << "Assuming a " << c->pars().getBeamFWHM()*3600 << " arcsec beam."
                  << "Change this value with beamFWHM parameter.\n";
        OB.bmaj=OB.bmin=c->pars().getBeamFWHM()*3600;
	}
	
	Beam NB;
	double npa = c->pars().getBpa();
	if(npa==-1) npa = OB.bpa;
	
	float linear = c->pars().getLinear();
	float factor = c->pars().getFactor();
	if (linear!=-1) {
		NB.bmaj = linear/KpcPerArc(c->pars().getDistance());
		NB.bmin = NB.bmaj;
		NB.bpa  = OB.bpa;
	}
	else {
        NB.bmaj = c->pars().getBmaj()==-1 ? factor*OB.bmaj : c->pars().getBmaj();
        NB.bmin = c->pars().getBmin()==-1 ? factor*OB.bmin : c->pars().getBmin();
		NB.bpa  = npa;
	}
	
	int Bhi[3], Blo[3];
	
	for (int i=0; i<6; i+=2) {
		int j = i/2;
		Blo[j] = c->pars().getBOX(i);
		Bhi[j] = c->pars().getBOX(i+1);
		if (Blo[j]<0 || Blo[j]>c->AxesDim(j)) Blo[j]=0;
		if (Bhi[j]<0 || Bhi[j]>c->AxesDim(j)) Bhi[j]=c->AxesDim(j);
	}
	
	fft = c->pars().getflagFFT();
	
	smooth(c, Bhi, Blo, OB, NB);
}
template void Smooth3D<short>::cubesmooth(Cube<short>*);
template void Smooth3D<int>::cubesmooth(Cube<int>*);
template void Smooth3D<long>::cubesmooth(Cube<long>*);
template void Smooth3D<float>::cubesmooth(Cube<float>*);
template void Smooth3D<double>::cubesmooth(Cube<double>*);


template <class T>
void Smooth3D<T>::smooth(Cube<T> *c, Beam Oldbeam, Beam Newbeam) {
	
	int Bhi[3] = {c->DimX(), c->DimY(), c->DimZ()};
	int Blo[3] = {0,0,0};
	
	smooth(c, Bhi, Blo, Oldbeam, Newbeam);
}
template void Smooth3D<short>::smooth(Cube<short>*,Beam, Beam);
template void Smooth3D<int>::smooth(Cube<int>*,Beam, Beam);
template void Smooth3D<long>::smooth(Cube<long>*,Beam, Beam);
template void Smooth3D<float>::smooth(Cube<float>*,Beam, Beam);
template void Smooth3D<double>::smooth(Cube<double>*,Beam, Beam);


template <class T>
void Smooth3D<T>::smooth(Cube<T> *c, int *Bhi, int *Blo, Beam Oldbeam, Beam Newbeam) {
	
	/// Front end function for Smooth3D. It initializes all needed 
	/// variable. 
		
	using namespace std;
	in = c;
	
	for (int i=0; i<3; i++) {
		dimAxes[i]	= c->AxesDim(i);
		bhi[i] 	 	= Bhi[i];
		blo[i]		= Blo[i];
		fhi[i]		= dimAxes[i];
		flo[i]		= 0;
	}

	NdatX = bhi[0]-blo[0];
	NdatY = bhi[1]-blo[1];
	NdatZ = bhi[2]-blo[2];
	
	crota = c->Head().Crota();
	
	beamDefined = defineBeam(Oldbeam, Newbeam);
	if (!beamDefined) abort();
	
	if (c->pars().isVerbose()) {

		int m=20;
		int n=7;	
		cout << showpoint << fixed << setprecision(2);
	
		cout << setfill('=') << setw(30) << right << " SMOOTH " << setw(24) << " ";
		cout << setfill(' ') << std::endl;
		cout << setw(m) << right << "Convert old beam: " << setw(n) 
				  << oldbeam.bmaj << "," << setw(n) << oldbeam.bmin 
				  << setw(n) << " arcsec" << setw(n) << oldbeam.bpa 
				  << " deg" << endl;
	
	
		cout << setw(m) << right << "to the new beam: " << setw(n) 
				  << newbeam.bmaj << "," << setw(n) << newbeam.bmin 
				  << setw(n) << " arcsec" << setw(n) << newbeam.bpa 
				  << " deg" << endl;
	
		cout << setw(m) << right << "with conv. beam: " << setw(n)
			      << conbeam.bmaj << "," << setw(n) << conbeam.bmin 
				  << setw(n) << " arcsec" << setw(n) << conbeam.bpa+90 
				  << " deg" << endl << endl;
	
		cout << setw(m) << "Scale factor: " << setw(n) 
				  << setprecision(4) << scalefac << endl;
				  
		cout << setfill('=') << setw(54) << " " << endl << endl;
		cout << setfill(' '); 
	}
	
	array = new T [NdatX*NdatY*NdatZ];
	arrayAllocated =true;
		
	bool allOK;

	if (fft) allOK = calculatefft(c->Array(), array);
	else allOK = calculate(c->Array(), array);
	if (!allOK) {
		std::cout << "SMOOTH error: cannot smooth data\n";
		abort();
	}

}
template void Smooth3D<short>::smooth(Cube<short>*,int*,int*,Beam, Beam);
template void Smooth3D<int>::smooth(Cube<int>*,int*,int*,Beam, Beam);
template void Smooth3D<long>::smooth(Cube<long>*,int*,int*,Beam, Beam);
template void Smooth3D<float>::smooth(Cube<float>*,int*,int*,Beam, Beam);
template void Smooth3D<double>::smooth(Cube<double>*,int*,int*,Beam, Beam);


template <class T>
bool Smooth3D<T>::smooth(Cube<T> *c, Beam Oldbeam, Beam Newbeam, T *OldArray, T *NewArray) {
	
	/// Front end function for Smooth3D without boxes. It gets
	/// information from Cube *c (dimAxes, cunit etc...) and
	/// smooths the OldArray in the NewArray. All arrays must 
	/// have the same X-Y-Z dimensions of Cube object. 
	/// This function does not write the smoothed array on the 
	/// Smooth3D::array variable!!!!
	

	using namespace std;
	in = c;
	
	for (int i=0; i<3; i++) {
		dimAxes[i]	= c->AxesDim(i);
		fhi[i]		= dimAxes[i];
		flo[i]		= 0;
		bhi[i] 	 	= fhi[i];
		blo[i]		= flo[i];
	}

	NdatX = bhi[0]-blo[0];
	NdatY = bhi[1]-blo[1];
	NdatZ = bhi[2]-blo[2];
	
	crota = c->Head().Crota();
	
	beamDefined = defineBeam(Oldbeam, Newbeam);	
	if (!beamDefined) return false;
		
	bool allOK;
	if (fft) allOK = calculatefft(OldArray, NewArray);
	else allOK = calculate(OldArray, NewArray);
	if (!allOK) {
		std::cout << "SMOOTH error: cannot smooth data\n";
		return false;
	}
	
	return true;
}
template bool Smooth3D<short>::smooth(Cube<short>*,Beam,Beam,short*,short*);
template bool Smooth3D<int>::smooth(Cube<int>*,Beam,Beam,int*,int*);
template bool Smooth3D<long>::smooth(Cube<long>*,Beam,Beam,long*,long*);
template bool Smooth3D<float>::smooth(Cube<float>*,Beam,Beam,float*,float*);
template bool Smooth3D<double>::smooth(Cube<double>*,Beam,Beam,double*,double*);


template <class T>
bool Smooth3D<T>::defineBeam(Beam Oldbeam, Beam Newbeam) {
	
	/// A function that defines the convolution beam, the 
	/// convolution field and the scale factor for smoothing.
		
	oldbeam = Oldbeam;
	newbeam = Newbeam;
	
	if (oldbeam.bmaj<oldbeam.bmin) {
		std::cout << "Major axis < minor axis. Inverting...";
		double dum = oldbeam.bmaj;
		oldbeam.bmaj = oldbeam.bmin;
		oldbeam.bmin = dum; 
	}
	
	bool agreed = ((newbeam.bmaj>=oldbeam.bmaj) && (newbeam.bmin>=oldbeam.bmin)); 
 
	if (!agreed) {
		std::cout << "SMOOTH error: new beam smaller than old beam\n";
		return false;
	}
	if (newbeam.bmaj<newbeam.bmin) {
		std::cout << "Major axis < minor axis. Inverting...";
		double dum = newbeam.bmaj;
		newbeam.bmaj = newbeam.bmin;
		newbeam.bmin = dum; 
	}

	if (!Convpars()) {
		std::cout << "SMOOTHING error: cannot calculate convolution parameters!\n";
		return false;
	}
	
	
	std::string cunit = in->Head().Cunit(0);
	double unittoarc=0;
	if (cunit=="DEGREE" || cunit=="DEGREES" || cunit=="DEG" ||
		cunit=="degree" || cunit=="degrees" || cunit=="deg") 
			unittoarc = 3600;
	else if (cunit=="ARCSEC" || cunit=="ARCS" ||
			 cunit=="arcsec" || cunit=="arcs") 
			unittoarc = 1;
	else if (cunit=="ARCMIN" || cunit=="ARCM" ||
			 cunit=="arcmin" || cunit=="arcm") 
			unittoarc = 60;
	else {
		std::cout << "SMOOTH error (unknown CUNIT for RA-DEC): ";
		std::cout << "cannot convert to ARCSEC.\n";
		std::cout << cunit;
		return false; 
	}
	
	double *datadummy = new double[NdatX*NdatY];
	T *dataI 	= new T[NdatX*NdatY];
    T *dataO	 	= new T[NdatX*NdatY];
	
	gridspace[0]=in->Head().Cdelt(0)*unittoarc; 
	gridspace[1]=in->Head().Cdelt(1)*unittoarc;

	confie = new double[maxconv];
	confieAllocated=true;
		
	conbeam.bpa -= crota;
	if (!Fillgauss2d(conbeam, 1.0, true, NconX, NconY, confie)) {
		std::cout << "SMOOTH error: Convolution region to big!\n";
		return false;
	}
	conbeam.bpa += crota;

	int maxbufsize=dimAxes[0]*NconY;
	if (maxbufsize>NdatX*NdatY) {
		cout << maxbufsize << "  " << NdatX << "  " << NdatY << " " << newbeam.bmaj << " " << NdatX*NdatY << endl;
		std::cout<<"SMOOTH: Convolution field too big or buffer too small.\n";
		return false;
	}


	int dumNconX, dumNconY;
	oldbeam.bpa=oldbeam.bpa+90-crota;
	if (!Fillgauss2d(oldbeam, 1.0, false, dumNconX, dumNconY, datadummy)) {
		std::cout << "SMOOTH error: Cannot calculate a default scale!\n";
		return false ;
	}
    else {
		int mx = max(NconX+10, dumNconX)/2;
		int my = max(NconY+10, dumNconY)/2;
		if ((mx*my)>NdatX*NdatY) {
			std::cout << "SMOOTH error: Gaussians too big for buffer!\n";
			return false;
		}
		int jp = 0;
		int ip = 0;
		int lx = dumNconX/2;
		int ly = dumNconY/2;

		for (int iy=-my; iy<=my; iy++) {
			for (int ix=-mx; ix<=mx; ix++) {
				ip++;
				if((ix<-lx) || (ix>lx) || (iy<-ly) || (iy>ly)) 
					dataI[ip] = 0.0;
				else {
					jp++;
					dataI[ip] = datadummy[jp-1];
				}
			}
		}
		lx = 2*mx+1;
		ly = 2*my+1;
		int Iresult = Convolve(confie, NconX, NconY, dataI, dataO, lx, ly);

        if (Iresult!=0) {
			std::cout << "Give scale factor for output: ";
			std::cin >> scalefac;
			std::cout << std::endl;
		}
        else scalefac = 1.0/dataO[(lx*ly)/2+1];
	}
	oldbeam.bpa=oldbeam.bpa-90+crota;

    delete [] datadummy;
	delete [] dataI;
	delete [] dataO;
	
	return true;
}
template bool Smooth3D<short>::defineBeam(Beam,Beam);
template bool Smooth3D<int>::defineBeam(Beam,Beam);
template bool Smooth3D<long>::defineBeam(Beam,Beam);
template bool Smooth3D<float>::defineBeam(Beam,Beam);
template bool Smooth3D<double>::defineBeam(Beam,Beam);


template <class T>
bool Smooth3D<T>::calculate(T *OldArray, T *NewArray) {
	
	if (!beamDefined) {
		std::cout << "SMOOTH error: Convolution beam is not set.\n";
		return false;
	}
	
	long size = (NdatX+NconX-1)*(NdatY+NconY-1);
	T *beforeCON = new T[size];
	T *afterCON  = new T[size];
		
	bool verb = in->pars().isVerbose();
    ProgressBar bar("Smoothing... ",true);
    bar.setShowbar(in->pars().getShowbar());
    if (verb) bar.init(NdatZ);

	if (!usescalefac) scalefac=1.0;

	for (int z=0; z<NdatZ; z++) {
		if (verb) bar.update(z+1);	
		for (int x=0; x<(NdatX+NconX-1); x++) {
			for (int y=0; y<(NdatY+NconY-1); y++) {
				long nPix = x+y*(NdatX+NconX-1);
				int oXpos = (x+blo[0]-(NconX-1)/2);
				int oYpos = (y+blo[1]-(NconY-1)/2);
				long oPix = oXpos+oYpos*dimAxes[0]+(z+blo[2])*dimAxes[0]*dimAxes[1];	
				if (x>=(NconX-1)/2 && x<=(NdatX+(NconX-1)/2) && 
					y>=(NconY-1)/2 && y<=(NdatY+(NconY-1)/2)) {
					beforeCON[nPix] = OldArray[oPix];
				}
				else beforeCON[nPix] = 0;	
			}
		}
		
		int Iresult = Convolve(confie, NconX, NconY, beforeCON, 
								afterCON, NdatX+NconX-1, NdatY+NconY-1);
		
		 
		if (Iresult!=0) {
			std::cout<<"SMOOTH error: cannot convolve requested functions (code="
					 << Iresult << ")" << std::endl;
			return false;
		}
       

		for (int x=(NconX-1)/2; x<(NdatX+(NconX-1)/2); x++) {
			for (int y=(NconY-1)/2; y<(NdatY+(NconY-1)/2); y++) {
				long nPix = x+y*(NdatX+NconX-1);
				long oPix = (x-(NconX-1)/2)+(y-(NconY-1)/2)*NdatX+z*NdatX*NdatY;	
				NewArray[oPix] = afterCON[nPix]*scalefac;
			}
		}
		
		
	}	

	if (verb) bar.fillSpace("OK.\n");

	delete [] beforeCON;
	delete [] afterCON;
	
	return true;
}
template bool Smooth3D<short>::calculate(short*,short*);
template bool Smooth3D<int>::calculate(int*,int*);
template bool Smooth3D<long>::calculate(long*,long*);
template bool Smooth3D<float>::calculate(float*,float*);
template bool Smooth3D<double>::calculate(double*,double*);


template <class T>
bool Smooth3D<T>::calculatefft(T *OldArray, T *NewArray) {
	
	if (!beamDefined) {
		std::cout << "SMOOTH error: Convolution beam is not set.\n";
		return false;
	}
	
	double *beforeCON = new double[NdatX*NdatY];
		
	bool verb = in->pars().isVerbose();
    ProgressBar bar(" Smoothing... ",false);
    bar.setShowbar(in->pars().getShowbar());
    if (verb) bar.init(NdatZ);
	
	Conv2D conv_fft;
	init_Conv2D (conv_fft, LINEAR_SAME, NdatX, NdatY, NconX, NconY);
	
	if (!usescalefac) scalefac=1.0;

	for (int z=0; z<NdatZ; z++) {
		if (verb) bar.update(z+1);	
		for (int x=0; x<NdatX; x++) {
			for (int y=0; y<NdatY; y++) {
				long nPix = x+y*NdatX;
				long oPix = (x+blo[0])+(y+blo[1])*dimAxes[0]+(z+blo[2])*dimAxes[0]*dimAxes[1];	
				beforeCON[nPix] = OldArray[oPix];	
			}
		}
		
		convolve (conv_fft, beforeCON, confie);	

		for (int x=0; x<NdatX; x++) {
			for (int y=0; y<NdatY; y++) {
				long nPix = x+y*NdatX;
				long oPix = x+y*NdatX+z*NdatX*NdatY;	
				NewArray[oPix] = conv_fft.dst[nPix]*scalefac;
			}
		}
		
	}

	if (verb) bar.fillSpace(" OK.\n");

	clear_Conv2D(conv_fft);
	
	delete [] beforeCON;

	return true;
}
template bool Smooth3D<short>::calculatefft(short*,short*);
template bool Smooth3D<int>::calculatefft(int*,int*);
template bool Smooth3D<long>::calculatefft(long*,long*);
template bool Smooth3D<float>::calculatefft(float*,float*);
template bool Smooth3D<double>::calculatefft(double*,double*);


template <class T>	
int Smooth3D<T>::Convolve(double *cfie, int ncx, int ncy, T *dat1, T *dat2, int ndx, int ndy) {
	
	/// Convolves data with two dimensional convolution function.
    /// Takes care of BLANKs.
    ///
    /// INPUTS:	the convolution field array cfie with dimensions 
    /// ncx and ncy and the array to be convolved dat1 with 
    /// dimensions ndx and ndy.
    ///
    /// OUTPUT: The convoluted array dat2.
    ///
    /// RETURN:         0 : succesfull convolution.
    ///                 1 : NDATX < NCONX.
    ///                 2 : NDATY < NCONY.
    ///                 3 : NCONX not an odd number.
    ///                 4 : NCONY not an odd number.
    
	int bcx = ncx/2;
	int bcy = ncy/2;
	int icx = ncx/2+1;
	int icy = ncy/2+1;
	int nblank = 0;
	int ntx = ndx-ncx+1;
	int nty = ndy-ncy+1;
	int ier = 0;
	int x,y,k;

	if (ntx<1) ier=1;                	// x size too small 
	if (nty<1) ier=2;                	// y size too small 
	if ((icx+bcx)!=ncx) ier=4;    		// odd x size of convolution function 
	if ((icy+bcy)!=ncy) ier=8;		   	// odd y size of convolution function 
	if (ier!=0) return ier;		         
	for (k=0, y=0; y<ndy; y++) {
		for (x=0; x<ndx; x++) {
			//if (dat1[k]==BLANK) nblank++;	
			//if (x<bcx || x>(ndx-icx)) dat2[k++] = BLANK;
			//else if (y<bcy || y>(ndy-icy)) dat2[k++] = BLANK;
			//else 
			dat2[k++] = 0.0;
		}
    }
   
	for (y=0; y<ncy; y++) {
		for (x=0; x<ncx; x++) {
			T cf = cfie[ncx-x-1+(ncy-y-1)*ncx];
			if (cf!=0.0) {
				for (int y2=0; y2<nty; y2++) {
					int y1 = y+y2;
					T *v1 = &dat1[y1*ndx+x];;
					T *v2 = &dat2[(bcy+y2)*ndx+bcx];
					if (!nblank) 
						for (int i=0; i<ntx; i++) v2[i] = v2[i]+cf*v1[i];
					else  {     
						for (int i=0; i<ntx; i++) { 
							//if (v1[i]!=BLANK && v2[i]!=BLANK) 
								v2[i] = v2[i]+cf*v1[i];
							//else v2[i] = BLANK;
						}
					}
				}
			}	
		}
   }
  
   return ier;
}
template int Smooth3D<short>::Convolve(double*,int,int,short*,short*,int,int);
template int Smooth3D<int>::Convolve(double*,int,int,int*,int*,int,int);
template int Smooth3D<long>::Convolve(double*,int,int,long*,long*,int,int);
template int Smooth3D<float>::Convolve(double*,int,int,float*,float*,int,int);
template int Smooth3D<double>::Convolve(double*,int,int,double*,double*,int,int);


template <class T>
bool Smooth3D<T>::Convpars() {
	      
	/// Determine parameters conbeam={bmaj, bmin, bpa} of 2d-gaussian 
	/// needed to convolve 'oldbeam'={oldbmaj, oldbmin, oldpa} to 
	/// 'newbeam'={newbmaj, newbmin, newpa}.
	/// Note that all angles are in DEGREES and wrt. pos. X-axis!
	/// Convoluted beam values are calculated with formulas given by 
	/// J. P. Wild in Aust. J. Phys.,1970,23,113-115. 
	///
	/// A distribution F2 is convolved with a gaussian F1 and the 
	/// result is a new distribution F0. For the distributions we can 
	/// write: F0 = F1*F2 (where * denotes a 2-dim convolution).
	/// This function finds the beam and angle of F1. In what follows a 
	/// and b denote the major and minor SEMI-axes and 'th' is the angle 
	/// of the major axis with the positive X-axis. 

	const double degtorad = atan(1.)/45.; 
	const double radtodeg = 45./atan(1.);

	double a2  = oldbeam.bmaj/2;
	double b2  = oldbeam.bmin/2;
	double a0  = newbeam.bmaj/2;
	double b0  = newbeam.bmin/2;
	double D0  = a0*a0-b0*b0;
	double D2  = a2*a2-b2*b2;    
	double th2 = oldbeam.bpa*degtorad;
	double th0 = newbeam.bpa*degtorad;
	double D1  = sqrt(D0*D0+D2*D2-2*D0*D2*cos(2*(th0-th2)));	
	
	double a1, b1, th1;
	
	double arg = 0.5*(a0*a0+b0*b0-a2*a2-b2*b2+D1); 
	if (arg<0) {
		std::cout << "SMOOTHING error: unsuitable new beam parameters!\n";
	return false;
	}
	else a1 = sqrt(arg);
      
	arg = 0.5*(a0*a0+b0*b0-a2*a2-b2*b2-D1); 
	if (arg<0) {
		std::cout << "SMOOTHING error: unsuitable new beam parameters!\n";
	return false;
	}
	else b1 = sqrt(arg);
	
	double nom   = D0*sin(2*th0)-D2*sin(2*th2);
    double denom = D0*cos(2*th0)-D2*cos(2*th2); 
    if (denom==0 && nom==0) th1=0;
	else {
		T twoth1 = atan2(nom,denom);
		th1 = twoth1/2;          
    }
      
    conbeam.bmaj = 2*a1;
    conbeam.bmin = 2*b1;
    conbeam.bpa  = th1*radtodeg;

	return true;
 
}
template bool Smooth3D<short>::Convpars();
template bool Smooth3D<int>::Convpars();
template bool Smooth3D<long>::Convpars();
template bool Smooth3D<float>::Convpars();
template bool Smooth3D<double>::Convpars();


template <class T>
bool Smooth3D<T>::Fillgauss2d(Beam varbeam, float ampl, bool norm, int &nconx, int &ncony, double *cfie) {

	/// Fill cfie with values of (rotated) 2-d Gaussian funtion.
	///
	/// -Beam varbeam(INPUT): 	the Beam object of the 2d Gaussian.
	/// -float ampl(INPUT):		the amplitude of 2d Gaussian. Only used if 
	///							norm=false.
    /// -bool norm(INPUT):		the 2d has to be normalized?
    /// -int nconx(OUTPUT):		X-dimension of the 2d Gaussian.
    /// -int ncony(OUTPUT):		Y-dimension of the 2d Gaussian.                              
	/// -float *cfie(OUTPUT):	An array containing the 2d Gaussian.
	///
	/// 
	/// Returns true if all ok, otherwise false.
	 
	
	const double degtorad = atan(1.)/45.; 
	double amplitude;
   
	if (norm) amplitude = 1.0; 
	else amplitude = ampl; 

	double phi = varbeam.bpa*degtorad;        
	double cs = cos(phi);
	double sn = sin(phi);  
	double beam[2] = {fabs(varbeam.bmaj), fabs(varbeam.bmin)};

	double xr = 0.5*beam[0];
	double yr = 0.5*beam[1];

	double extend = sqrt(-1.0*log(cutoffratio)/log(2.0));
	xr *= extend;
	yr *= extend;
   
	double x1 = fabs(xr*cs-0.0*sn);
	double y1 = fabs(xr*sn+0.0*cs);
	double x2 = fabs(0.0*cs-yr*sn);
	double y2 = fabs(0.0*sn+yr*cs);
	double x  = (x2>x1 ? x2:x1);
	double y  = (y2>y1 ? y2:y1);

	double gridspac[2] = {fabs(gridspace[0]), fabs(gridspace[1])}; 
   
	int Xmax = round(x/gridspac[0]); 
	int Ymax = round(y/gridspac[1]); 

	int nconX = 2*Xmax+1;    
	int nconY = 2*Ymax+1; 

	if ((nconX*nconY)>maxconv) {
		std::cout << "SMOOTH error: Convolution function too big for buffer.\n";
		return false;
	} else {
		nconx = nconX;
		ncony = nconY;         
	}
      
	double argfac = -4.0 * log(2.0);

	double totalarea = 0;
	for (int j=-Ymax; j<=Ymax; j++) {
		for (int i=-Xmax; i<=Xmax; i++) {
			int pos = (j+Ymax)*nconX+(i+Xmax);    
			x = i*gridspac[0];
			y = j*gridspac[1];
			xr = x*cs + y*sn;
			yr = -1.0*x*sn + y*cs;
			double argX=0;
			double argY=0;
			if (beam[0]!=0) argX = xr/beam[0];
			if (beam[1]!=0) argY = yr/beam[1];
			double arg = argfac*(argX*argX+argY*argY);
			double c = exp(arg);
			if (c>=cutoffratio) {
				c *= amplitude;
				cfie[pos] = c;
				totalarea += c;
			} 
			else cfie[pos] = 0;
		}
	}
   

	if (norm) {
		for (int i=0; i<(nconX*nconY); i++) {
			cfie[i] = cfie[i]/totalarea;
		}
	}
	
	return true;   
}
template bool Smooth3D<short>::Fillgauss2d(Beam,float,bool,int&,int&,double*);
template bool Smooth3D<int>::Fillgauss2d(Beam,float,bool,int&,int&,double*);
template bool Smooth3D<long>::Fillgauss2d(Beam,float,bool,int&,int&,double*);
template bool Smooth3D<float>::Fillgauss2d(Beam,float,bool,int&,int&,double*);
template bool Smooth3D<double>::Fillgauss2d(Beam,float,bool,int&,int&,double*);


template <class T> 
void Smooth3D<T>::fitswrite() {
	
	int ax[3] = {NdatX, NdatY, NdatZ};
	Cube<T> *out = new Cube<T>(ax);
	for (int i=0; i<out->NumPix(); i++) 
		out->Array()[i] = array[i];
 	out->saveHead(in->Head());
 	out->saveParam(in->pars());
 	out->Head().setDimAx(0, NdatX);
 	out->Head().setDimAx(1, NdatY);
 	out->Head().setDimAx(2, NdatZ);
 	out->Head().setCrpix(0, in->Head().Crpix(0)-blo[0]);
 	out->Head().setCrpix(1, in->Head().Crpix(1)-blo[1]);
 	out->Head().setCrpix(2, in->Head().Crpix(2)-blo[2]);
	out->Head().setBmaj(newbeam.bmaj/3600.);
	out->Head().setBmin(newbeam.bmin/3600.);
	out->Head().setBpa(newbeam.bpa);
	out->Head().calcArea();
    /*
	std::string name = in->pars().getImageFile();	
	name = in->pars().getImageFile();
	int found = name.find(".fits");
	name.insert(found, "_s"+to_string<int>(lround(newbeam.bmaj)));
    */

    std::string name = in->pars().getSmoothOut();
    if (name=="NONE") {
        in->pars().getOutfolder()+in->Head().Name();
        name += "_s"+to_string<int>(lround(newbeam.bmaj));
        if (in->pars().getflagReduce()) name+="red";
        name+=".fits";
    }

	std::string nbeam = to_string(newbeam.bmaj)+" x "+to_string(newbeam.bmin)+" arcsec, ";
	nbeam += "pa = "+to_string(newbeam.bpa)+" deg";
	out->Head().Keys().push_back("HISTORY BBAROLO SMOOTHING: New beam "+nbeam);
	std::string obeam = to_string(oldbeam.bmaj)+" x "+to_string(oldbeam.bmin)+" arcsec, ";
	obeam += "pa = "+to_string(oldbeam.bpa)+" deg";
	out->Head().Keys().push_back("HISTORY BBAROLO SMOOTHING: Old beam "+obeam);
	
	int factor = floor(in->pars().getFactor());
	if (in->pars().getflagReduce() && factor>1) {
		Cube<T> *red = out->Reduce(factor);	
        red->fitswrite_3d(name.c_str(),true);
		delete red;
	}
	else {
		T minn,maxx;
		findMinMax<T>(out->Array(), out->NumPix(), minn, maxx);
		out->Head().setDataMax(double(maxx));
		out->Head().setDataMin(double(minn));
		out->fitswrite_3d(name.c_str(),true);
	}
	delete out;
}
template void Smooth3D<short>::fitswrite();
template void Smooth3D<int>::fitswrite();
template void Smooth3D<long>::fitswrite();
template void Smooth3D<float>::fitswrite();
template void Smooth3D<double>::fitswrite();	
