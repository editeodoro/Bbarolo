//---------------------------------------------------------------
// utils.cpp: Generic utility functions.
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
 along with BBarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning BBarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "utils.hh"
#include "../Arrays/header.hh"
#include "../Map/voxel.hh"
#include <fitsio.h>
#include <wcslib/wcs.h>


double KpcPerArc(double d) {return 2*d*tan(M_PI/2/180)/3.6;}
double VeltoDist(double vsys) {return vsys/70.;}
double RedtoDist(double redshift) {return redshift*299792.458/70.;};

template <class Type> 
bool isNaN (Type n) {
	volatile Type d = n; 
	return d != d;
}
template bool isNaN (short);
template bool isNaN (int);
template bool isNaN (long);
template bool isNaN (float);
template bool isNaN (double);


bool fexists(std::string filename) {std::ifstream ifile(filename.c_str()); return ifile;}


void FitsWrite_2D (const char *filename, float *image, long xsize, long ysize) {
	
	fitsfile *fptr;
	int status, numAxes = 2;
	long firstPix = 1, numPix;
	long dimAxes[2] = {xsize, ysize};
	numPix = xsize*ysize;
	
	remove (filename);
	
	//Create the new file
	status = 0;
	if (fits_create_file (&fptr, filename, &status)) {
		fits_report_error (stderr, status);
		}
	
	//Create the primary array image	
    if (fits_create_img (fptr, FLOAT_IMG, numAxes, dimAxes, &status)){
		fits_report_error (stderr, status);
	}
	
	
    if (fits_write_img (fptr, TFLOAT, firstPix, numPix, image, &status)){
		fits_report_error (stderr, status);
	}
	
	// Close the FITS File
    if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
	}
	
}



void FitsWrite_2D (const char *filename, double *image, long xsize, long ysize) {

    fitsfile *fptr;
    int status, numAxes = 2;
    long firstPix = 1, numPix;
    long dimAxes[2] = {xsize, ysize};
    numPix = xsize*ysize;

    remove (filename);

    //Create the new file
    status = 0;
    if (fits_create_file (&fptr, filename, &status)) {
        fits_report_error (stderr, status);
        }

    //Create the primary array image
    if (fits_create_img (fptr, DOUBLE_IMG, numAxes, dimAxes, &status)){
        fits_report_error (stderr, status);
    }


    if (fits_write_img (fptr, TDOUBLE, firstPix, numPix, image, &status)){
        fits_report_error (stderr, status);
    }

    // Close the FITS File
    if (fits_close_file(fptr, &status)){
        fits_report_error(stderr, status);
    }

}


void FitsWrite_3D (const char *outfile, float *outcube, long *dimAxes) {
	
	fitsfile *fptr;
	int status = 0;      
	long  fpixel = 1, naxis = 3, nelements;
	long dnaxes[3] = {dimAxes[0], dimAxes[1], dimAxes[2]};
	long naxes[3] = {dimAxes[2], dimAxes[1], dimAxes[0]};
 
	remove(outfile);             
	status = 0;   
	fits_create_file(&fptr, outfile, &status);  
	fits_create_img(fptr, FLOAT_IMG, naxis, dnaxes, &status);
    	
	nelements = naxes[0]*naxes[1]*naxes[2]; 
	
	fits_write_img(fptr, TFLOAT, fpixel, nelements, outcube, &status);

	fits_close_file(fptr, &status);      
  
	fits_report_error(stderr, status);  
}


void FitsWrite_3D (const char *outfile, short *outcube, long *dimAxes) {
	
	fitsfile *fptr;
	int status = 0;      
	long  fpixel = 1, naxis = 3, nelements;
	long dnaxes[3] = {dimAxes[0], dimAxes[1], dimAxes[2]};
	long naxes[3] = {dimAxes[2], dimAxes[1], dimAxes[0]};
 
	remove(outfile);             
	status = 0;   
	fits_create_file(&fptr, outfile, &status);  
	fits_create_img(fptr, SHORT_IMG, naxis, dnaxes, &status);
    	
	nelements = naxes[0]*naxes[1]*naxes[2]; 
	
	fits_write_img(fptr, TSHORT, fpixel, nelements, outcube, &status);

	fits_close_file(fptr, &status);      
  
	fits_report_error(stderr, status);  
}


template <class T> 
T AlltoVel (T in, Header &h) {
	
  /// This function convert a spectral input value "in" in 
  /// a output value in units of KM/S.
	
    std::string cunit2 = h.Cunit(h.NumAx()-1);
    if (h.NumAx()>3) cunit2 = h.Cunit(2);

    double freq0 = h.Freq0();
    const double c = 299792458;
    T vel_km_s;


	if (cunit2=="KM/S" || cunit2=="Km/s" || cunit2=="km/s" || cunit2=="kms") return in;
    else if (cunit2=="M/S" || cunit2=="m/s" || cunit2=="M/s" || cunit2=="ms") {
		vel_km_s = in/1000.;
	}
	else if (cunit2=="HZ" || cunit2=="Hz" || cunit2=="hz") {
		const double HIrest = freq0;
		T vel_m_s = c*(HIrest*HIrest-in*in)/(HIrest*HIrest+in*in);
		vel_km_s = vel_m_s/1000.;
	}
	else if (cunit2=="MHZ" || cunit2=="MHz" || cunit2=="Mhz") {
        const double HIrest = freq0;
		T vel_m_s = c*(HIrest*HIrest-in*in)/(HIrest*HIrest+in*in);
		vel_km_s = vel_m_s/1000.;
	}
    else if (cunit2=="MUM" || cunit2=="mum" || cunit2=="Mum" || cunit2=="um" ||
             cunit2=="A" || cunit2=="a" || cunit2=="Ang" || cunit2=="ang" ||
             cunit2=="Angstrom" || cunit2=="angstrom" || cunit2=="ANGSTROM") {          // Micron or Angstrom for Hi-Z
        //int z_cent = floor(h.DimAx(2)/2.)-1;
        int z_cent = h.Crpix(2)-1;
        double line_wave = (z_cent+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
        T vel_m_s = c*(in*in-line_wave*line_wave)/(in*in+line_wave*line_wave);
		vel_km_s = vel_m_s/1000.;
    }
	else return in;
	
	return vel_km_s;
	
}
template short AlltoVel (short, Header &);
template int AlltoVel (int, Header &);
template long AlltoVel (long, Header &);
template float AlltoVel (float, Header &);
template double AlltoVel (double, Header &);


template <class T> 
T DeltaVel (Header &h) {

 /// This function convert the spectral cdelt value 
 /// in a output cdelt value in units of KM/S.	
	
	T deltaV=0;
	long zdim = h.DimAx(2);
	
	if (h.Cunit(2)=="KM/S" || h.Cunit(2)=="Km/s" || h.Cunit(2)=="km/s" || h.Cunit(2)=="kms") return h.Cdelt(2);
	else if (h.Cunit(2)=="M/S" || h.Cunit(2)=="m/s" || h.Cunit(2)=="M/s") {
		deltaV =h.Cdelt(2)/1000.;
	}
	else if (h.Cunit(2)=="HZ" || h.Cunit(2)=="Hz" || h.Cunit(2)=="hz") {
		
		const double HIrest = h.Freq0();
		const double c = 299792458;

		T sum=0;
		for (int i=0; i<zdim-1; i++) {
			T ipix = (i+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
			T ivel = c*(HIrest*HIrest-ipix*ipix)/(HIrest*HIrest+ipix*ipix);
			T spix = (i+2-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
			T svel = c*(HIrest*HIrest-spix*spix)/(HIrest*HIrest+spix*spix);
			T diff = svel - ivel;
			sum += diff;
		}
	
		deltaV = sum/(zdim-1);
		deltaV /= 1000;
	}	
	else if (h.Cunit(2)=="MHZ" || h.Cunit(2)=="MHz" || h.Cunit(2)=="Mhz") {
		
        const double HIrest = h.Freq0();
		const double c = 299792458;
		
		T sum=0;
		for (int i=0; i<zdim-1; i++) {
			T ipix = (i+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
			T ivel = c*(HIrest*HIrest-ipix*ipix)/(HIrest*HIrest+ipix*ipix);
			T spix = (i+2-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
			T svel = c*(HIrest*HIrest-spix*spix)/(HIrest*HIrest+spix*spix);
			T diff = svel - ivel;
			sum += diff;
		}
	
		deltaV = sum/(zdim-1);
		deltaV /= 1000;
	}
    else if (h.Cunit(2)=="MUM" || h.Cunit(2)=="mum" || h.Cunit(2)=="Mum" || h.Cunit(2)=="um" ||
             h.Cunit(2)=="A" || h.Cunit(2)=="a" || h.Cunit(2)=="Ang" || h.Cunit(2)=="ang" ||
             h.Cunit(2)=="Angstrom" || h.Cunit(2)=="angstrom" || h.Cunit(2)=="ANGSTROM") {            // Micron or Angstrom for Hi-Z
        const double c = 299792458;
        int z_cent = floor(h.DimAx(2)/2.)-1;
        double line_wave = (z_cent+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
		
		T sum=0;
		for (int i=0; i<zdim-1; i++) {
			T ipix = (i+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
            T ivel = c*(ipix*ipix-line_wave*line_wave)/(ipix*ipix+line_wave*line_wave);
			T spix = (i+2-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
            T svel = c*(spix*spix-line_wave*line_wave)/(spix*spix+line_wave*line_wave);
			T diff = svel - ivel;
			sum += diff;
		}

        deltaV = sum/(zdim-1);
		deltaV /= 1000;
	}
	else return h.Cdelt(2);
	
	return deltaV;
}
template short DeltaVel (Header&);
template int DeltaVel (Header&);
template long DeltaVel (Header&);
template float DeltaVel (Header&);
template double DeltaVel (Header&);


template <class T> 
T FluxtoJy (T in, Header &h) {
	
 /// This function convert the input flux value  
 /// in a output flux value in units of Jy.	
 
	T fluxvalue, fluxJY;
	
	fluxvalue = in;
	
	if (h.Bunit()=="W.U." || h.Bunit()=="w.u." || h.Bunit()=="W.u.") {
		fluxJY = 5E-3*fluxvalue;
		if (h.BeamArea()!=0) fluxJY /= h.BeamArea();
	}
	else if (h.Bunit()=="JY/BEAM" || h.Bunit()=="Jy/beam" ||
		     h.Bunit()=="jy/beam" || h.Bunit()=="Jy/Beam") {
		fluxJY = fluxvalue;
		if (h.BeamArea()!=0) fluxJY /= h.BeamArea();
	}
	else if (h.Bunit()=="JY" || h.Bunit()=="Jy" || h.Bunit()=="jy") {
		return fluxvalue; 
	}
	else return fluxvalue;
	
	return fluxJY;
	
}
template short FluxtoJy (short, Header&);
template int FluxtoJy (int, Header&);
template long FluxtoJy (long, Header&);
template float FluxtoJy (float, Header&);
template double FluxtoJy (double, Header&);



std::string decToDMS(const double dec, std::string type, int decPrecision) {
  /** 
   *Converts a decimal angle (in degrees) to a format reflecting the axis type:
   *  RA   (right ascension):     hh:mm:ss.ss, with dec modulo 360. (24hrs)
   *  DEC  (declination):        sdd:mm:ss.ss  (with sign, either + or -)
   *  GLON (galactic longitude): ddd:mm:ss.ss, with dec made modulo 360.
   *  GLAT (galactic latitude):  sdd:mm:ss.ss  (with sign, either + or -)
   *    Any other type defaults to RA, and prints warning.
   *
   * \param dec Decimal value of the angle, in degrees.
   * \param type String indicating desired type of output. Options RA, DEC, 
   *              GLON, GLAT
   * \return String with angle in desired format.
   */

	double dec_abs,degD,minD,minDint,sec;
	int deg,min;
	const double minPerHour=60.;
	const double degPerHour=15.;
	double thisDec = dec;
	std::string sign="";
	int degSize = 2; // number of figures in the degrees part of the output.
	
	int precision=std::max(0,decPrecision);
	std::string Type = type;	
    int found = Type.find("RA");
	if(found>=0) Type = "RA";
	else {
		found = Type.find("DEC");
		if(found>=0) Type = "DEC";
		else {
			found = Type.find("GLON");
			if (found>=0) Type = "GLON";
			else {
				found = Type.find("GLAT");
				if (found>=0) Type = "GLAT";
				else {	// UNKNOWN TYPE -- DEFAULT TO RA.
					std::cout << "WARNING <decToDMS> : Unknown axis type ("
							  << type << "). Defaulting to using RA.\n";
					while (thisDec<0.) thisDec += 360.;
					while (thisDec>=360.) thisDec -= 360.;
					thisDec /= degPerHour;
				}
			}
		}
	}
	
	if (Type=="RA") precision++;
	if ((Type=="RA")||(Type=="GLON")) {
		if(Type=="GLON")  degSize = 3; // longitude has three figures in degrees.
		while (thisDec<0.) thisDec += 360.;
		while (thisDec>=360.) thisDec -= 360.;
		if(Type=="RA") thisDec /= degPerHour;  // Convert to hours.
	}
	else if((Type=="DEC")||(Type=="GLAT")) {
		if(thisDec<0.) sign = "-";
		else sign = "+";
	}
 
	dec_abs = fabs(thisDec);
	minD = modf(dec_abs, &degD) * minPerHour;
	sec = modf(minD, &minDint) * minPerHour;
	deg = int(degD);
	min = int(minDint);

	if(fabs(sec-minPerHour)<std::pow(double(10),-double(precision))){ // to prevent rounding errors stuffing things up
		sec=0.;
		min++;
		if(min==60){
			min=0;
			deg++;
		}
	}

	std::stringstream ss(std::stringstream::out);
	ss.setf(std::ios::showpoint);
	ss.setf(std::ios::fixed);
	ss << sign;
	ss << std::setw(degSize)<<std::setfill('0')<<deg<<":";
	ss<<std::setw(2)<<std::setfill('0')<<min<<":";
	if(precision>0)
		ss<<std::setw(precision+3)<<std::setprecision(precision)<<sec;
	else {
		ss << std::setw(2) << int(sec);
	}

	return ss.str();
}


double dmsToDec(std::string dms) {
  /** 
   *  double dmsToDec(string)
   *   Converts a std::string in the format +12:23:34.45 to a decimal angle in degrees.
   *   Assumes the angle given is in degrees, so if passing RA as the argument,
   *   need to multiply by 15 to get the result in degrees rather than hours.
   *   The sign of the angle is preserved, if present.
   */

	bool isNeg = false;
	if(dms[0]=='-') isNeg = true;

	std::stringstream ss;
	ss.str(dms);
	std::string deg,min,sec;
	getline(ss,deg,':');
	getline(ss,min,':');
	getline(ss,sec);
	char *end;
	double d = strtod(deg.c_str(),&end);
	double m = strtod(min.c_str(),&end);
	double s = strtod(sec.c_str(),&end);  

	double dec = fabs(d) + m/60. + s/3600.;
	if(isNeg) dec = dec * -1.;

	return dec;

}


template <class T> 
T angularSeparation(T &ra1, T &dec1, T &ra2, T &dec2) {
  	/*
   *  Enter ra & dec for two positions. (all positions in degrees)
   *  Returns the angular separation in degrees.
   */
	
	const long double degToRadian=M_PI/180.;
	long double dra = (ra1-ra2)*degToRadian;
	long double d1 = dec1*degToRadian;
	long double d2 = dec2*degToRadian;
	long double angsep;
    if((fabs(ra1-ra2) < 1./3600.)&&(fabs(dec1-dec2)<1./3600.))
		return sqrt(dra*dra + (d1-d2)*(d1-d2)) / degToRadian;
	else {
		if(fabs(ra1-ra2) < 1./3600.)
		angsep = cos(d1)*cos(d2) - dra*dra*cos(d1)*cos(d2)/2. + sin(d1)*sin(d2);
		else angsep = cos(dra)*cos(d1)*cos(d2) + sin(d1)*sin(d2);
		double dangsep = acos(angsep) / degToRadian;
		return dangsep;
  	}

}
template short angularSeparation(short&,short&,short&,short&);
template int angularSeparation(int&,int&,int&,int&);
template long angularSeparation(long&,long&,long&,long&);
template float angularSeparation(float&,float&,float&,float&);
template double angularSeparation(double&,double&,double&,double&);


double arcsconv(std::string cunit) {
	
	if (cunit=="DEGREE" || cunit=="DEGREES" || cunit=="DEG" ||
		cunit=="degree" || cunit=="degrees" || cunit=="deg") 
		return 3600.;
	else if (cunit=="ARCSEC" || cunit=="ARCS" ||
			 cunit=="arcsec" || cunit=="arcs") 
		return 1.;
	else if (cunit=="ARCMIN" || cunit=="ARCM" ||
			 cunit=="arcmin" || cunit=="arcm") 
		return 60.;
	else {
		std::cout << "Conversion error (unknown CUNIT for RA-DEC): ";
		std::cout << "cannot convert to ARCSEC.\n";
		std::cout << cunit;
		std::terminate(); 
	}
	
}


double degconv(std::string cunit) {

    if (cunit=="DEGREE" || cunit=="DEGREES" || cunit=="DEG" ||
        cunit=="degree" || cunit=="degrees" || cunit=="deg")
        return 1.;
    else if (cunit=="ARCSEC" || cunit=="ARCS" ||
        cunit=="arcsec" || cunit=="arcs")
        return 1/3600;
    else if (cunit=="ARCMIN" || cunit=="ARCM" ||
        cunit=="arcmin" || cunit=="arcm")
        return 1/60.;
    else {
        std::cout << "Conversion error (unknown CUNIT for RA-DEC): ";
        std::cout << "cannot convert to DEGREE.\n";
        std::cout << cunit;
        std::terminate();
    }
}



template <class T>
T Pbcor (PixelInfo::Voxel<T> &v, Header &h) {
	
	bool haveCorr = h.Telesc()=="WSRT" || h.Telesc()=="VLA" || h.Telesc()=="ATCA"
					|| h.Telesc()=="FST" || h.Telesc()=="GMRT";
	if (!haveCorr) return v.getF();
	
	const double degtorad = M_PI/180.;
	T fluxcorr=0;
	
	float pcRA = h.Crval(0);				/// Pointing center R.A. (in degrees).
	float pcDEC = h.Crval(1);				/// Pointing center DEC  (in degrees).
	
	if (pcDEC>=0) pcDEC = 90 - pcDEC;			/// Conversion of DEC to polar angle.
	else pcDEC = 90 + fabs(pcDEC);
	
	pcRA  *= degtorad;							/// Conversions to radians.
	pcDEC *= degtorad;
	
	float posRA = ((v.getX()+1-h.Crpix(0))*h.Cdelt(0)+h.Crval(0));   
	float posDEC = ((v.getY()+1-h.Crpix(1))*h.Cdelt(1)+h.Crval(1));
	
	if (posDEC>=0) posDEC = 90 - posDEC;
	else posDEC = 90 + fabs(posDEC);
	if (posRA<0) posRA += 360;
	
	posRA  *= degtorad;
	posDEC *= degtorad;
	
	float acarg = cos(pcDEC)*cos(posDEC) + sin(pcDEC)*sin(posDEC)*cos(pcRA-posRA);
	acarg = std::max<float>(-1., std::min<float>(1.,acarg));
	
	float angdist = acos(acarg)/degtorad;
	
	float freq=0; 								/// Frequency in GHz
	
	if (h.Cunit(2)=="KM/S" || h.Cunit(2)=="Km/s" || h.Cunit(2)=="km/s") {
		const float HIrest = 1.420405751;
		const float c = 299792.458;
		float vel = ((v.getZ()+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2));
		freq = HIrest*sqrt((1-vel/c)/(1+vel/c));
	}
	else if (h.Cunit(2)=="M/S" || h.Cunit(2)=="m/s" || h.Cunit(2)=="M/s") {
		const float HIrest = 1.420405751;
		const float c = 299792458;
		float vel = ((v.getZ()+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2));
		freq = HIrest*sqrt((1-vel/c)/(1+vel/c));
	}
	else if (h.Cunit(2)=="HZ" || h.Cunit(2)=="Hz" || h.Cunit(2)=="hz") {
		freq=((v.getZ()+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2))/(1.e09);
	}	
	else if (h.Cunit(2)=="MHZ" || h.Cunit(2)=="MHz" || h.Cunit(2)=="Mhz") {
		freq=((v.getZ()+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2))/(1000);
	}
	else return -2;
	
	if (h.Telesc()=="WSRT") {
		const float calib = 61.18;
		float cutoff = 0.023;
		float pbc = pow(cos(calib*freq*angdist*degtorad), 6);
		
		if (pbc<cutoff) return 0;
		else fluxcorr = v.getF()/pbc;
	}
	else if (h.Telesc()=="VLA") {
		float RF = angdist*60*freq;					/// The distance is now in arcmin.
		float cutoff = 0.023;
		float a0, a1, a2, a3, a4;
		
		if (freq<1.43) {
			a0 = 1.;
			a1 = -1.329E-03;
			a2 = +6.445E-07; 
			a3 = -1.146E-10;
			a4 = 0.;
		}
		else if (freq>=1.43 && freq<=1.73) {
			a0 = 1.;
			a1 = -1.343E-03;
			a2 = +6.579E-07; 
			a3 = -1.186E-10;
			a4 = 0.;
		}
		else {
			a0 = +0.9920378;
			a1 = +0.9956885E-03;
			a2 = +0.3814573E-05; 
			a3 = -0.5311695E-08;
			a4 = +0.3980963E-11;
		}
		
		float polyn = a0 + a1*pow(RF,2) + a2*pow(RF,4) + a3*pow(RF,6) + a4*pow(RF,8);
		float pbc = 1/polyn;
		pbc = 1.0/std::max<float>(1, pbc);
		
		
		if (pbc<cutoff) return 0;
		else fluxcorr = v.getF()/pbc;	
	}
	else if (h.Telesc()=="ATCA") {
		float RF = angdist*60*freq;					/// The distance is now in arcmin.
		if (RF>50) return 0;						/// Cut-off at 50 arcmin*GHz
		const float a0 = 1.;
		const float a1 = 8.99E-04;
		const float a2 = 2.15E-06; 
		const float a3 = -2.23E-09;
		const float a4 = 1.56E-12;
		float polyn = a0 + a1*pow(RF,2) + a2*pow(RF,4) + a3*pow(RF,6) +a4*pow(RF,8);
		float pbc = 1/polyn;
		pbc = 1.0/std::max<float>(1, pbc);
		
		fluxcorr = v.getF()/pbc;	
	}
	else if (h.Telesc()=="GMRT") {
		float RF = angdist*60*freq;					/// The distance is now in arcmin.
		const float a0 = -2.27961E-03;	
		const float a1 = 21.4611E-07;
		const float a2 = -9.7929E-10; 
		const float a3 = 1.80153E-13;
		float polyn = 1 + a0*pow(RF,2) + a1*pow(RF,4) + a2*pow(RF,6) +a3*pow(RF,8);
		float pbc = 1/polyn;
		pbc = 1.0/std::max<float>(1, pbc);
		
		fluxcorr = v.getF()/pbc;	
	}
	else if (h.Telesc()=="FST") {
		float RF = angdist*freq;
		if (angdist>2.8) return 0;				/// Cut-off at 2.8 degrees.
		const float a0 = 0.8031;
		
		float pbc = exp(-a0*RF*RF);
		
		fluxcorr = v.getF()/pbc;
	}
	
	
	return fluxcorr;
}
template short Pbcor (PixelInfo::Voxel<short>&, Header&);
template int Pbcor (PixelInfo::Voxel<int>&, Header&);
template long Pbcor (PixelInfo::Voxel<long>&, Header&);
template float Pbcor (PixelInfo::Voxel<float>&, Header&);
template double Pbcor (PixelInfo::Voxel<double>&, Header&);




template <class T> 
void Pbcor (long x, long y, long z, T &flux, Header &h) {
	
	PixelInfo::Voxel<T> vox (x, y, z, flux);
	flux = Pbcor<T>(vox, h);
	
} 
template void Pbcor (long, long, long, short&, Header&);
template void Pbcor (long, long, long, int&, Header&);
template void Pbcor (long, long, long, long&, Header&);
template void Pbcor (long, long, long, float&, Header&);
template void Pbcor (long, long, long, double&, Header&);


template <class T>
bool getData (std::vector<std::vector<T> > &allData, std::string file, bool err_verbose) {
	
	std::ifstream filein(file.c_str());
	if (!filein) {
        if (err_verbose) std::cout << "\n ERROR: " << file << " doesn't exist!" << std::endl;
		return false;
	}
	std::string line;
	while (getline(filein,line)) {
		std::vector<T> lineData;
		float val;
		std::istringstream lineStream(line);
		while (lineStream >> val) {
			lineData.push_back(val);
		}
		allData.push_back(lineData);
	}
	filein.close();
	
	return true;
	
}
template bool getData (std::vector<std::vector<short> > &,std::string,bool);
template bool getData (std::vector<std::vector<int> > &,std::string,bool);
template bool getData (std::vector<std::vector<long> > &,std::string,bool);
template bool getData (std::vector<std::vector<float> > &,std::string,bool);
template bool getData (std::vector<std::vector<double> > &,std::string,bool);


template <class T>
bool getDataColumn (std::vector<T> &data, std::string filestring) {
	
	// This function read a column into the array "data". 
	// The input file string should as follow:
	//
	// "file(file.dat,n,m:k)"
	//
	// where file.dat is the file to open, n is the number of the
	// column to be read and m:k is the range of rows to be read.
	
	
	int found = filestring.find("ile");
	if (found==-1) return false;	
	
	filestring.erase(0,5);
	
	std::string filename = filestring;
	std::string column = filestring;
	std::string rows = filestring;	
	int col=0, row_low=0, row_high=INT_MAX;
	
	found = filename.find(",");
	if (found==-1) filename.erase(filename.size()-1);
	else {
		filename.erase(found, filename.size()-1);
		column.erase(0,found+1);
		rows.erase(0,found+1);
		found = column.find(",");
		if (found==-1) column.erase(column.size()-1);
		else {
			column.erase(found, column.size()-1);
			rows.erase(0,found+1);
			rows.erase(rows.size()-1);
			found = rows.find(":");
			if (found==-1) row_low = atoi(rows.c_str())-1;
			else {
				std::string rows_low  = rows.substr(0,found);
				std::string rows_high = rows.substr(found+1, rows.size()-1);
				row_low = atoi(rows_low.c_str())-1;
				if (rows_high.size()!=0) row_high = atoi(rows_high.c_str())-1; 
			}
		}
		col = atoi(column.c_str())-1;
	}
	
	filename = deblank(filename);
	checkHome(filename);
	std::ifstream filein(filename.c_str());
	if (!filein) {
		std::cout << "\n ERROR: " << filename << " doesn't exist!";
		return false;
	}
	
	int row_n=0, col_n=0;
	std::string line;
	while (std::getline(filein,line)) {
		if (row_n>=row_low && row_n<=row_high) {
			T val;
			std::istringstream lineStream(line);
			while (lineStream>>val) {
				if (col_n++==col) data.push_back(val);
			}
			col_n=0;
		}
		row_n++;
	}
	filein.close();
	
	return true;
}
template bool getDataColumn (std::vector<short> &,std::string);
template bool getDataColumn (std::vector<int> &,std::string);
template bool getDataColumn (std::vector<long> &,std::string);
template bool getDataColumn (std::vector<float> &,std::string);
template bool getDataColumn (std::vector<double> &,std::string);	



template <> int selectBitpix<short>() {return SHORT_IMG;};
template <> int selectBitpix<int>() {return SHORT_IMG;};
template <> int selectBitpix<long>() {return LONG_IMG;};
template <> int selectBitpix<float>() {return FLOAT_IMG;};
template <> int selectBitpix<double>() {return DOUBLE_IMG;};

template <> int selectDatatype<short>() {return TSHORT;};
template <> int selectDatatype<int>() {return TINT;};
template <> int selectDatatype<long>() {return TLONG;};
template <> int selectDatatype<float>() {return TFLOAT;};
template <> int selectDatatype<double>() {return TDOUBLE;};
