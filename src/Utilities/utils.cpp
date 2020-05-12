//---------------------------------------------------------------
// utils.cpp: Generic utility functions.
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <Utilities/utils.hh>
#include <Arrays/header.hh>
#include <Arrays/param.hh>
#include <Arrays/rings.hh>
#include <Map/voxel.hh>
#include <fitsio.h>
#include <wcslib/wcs.h>
#include <sys/stat.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>


std::string get_currentpath() {
    std::string path = "";
    char pathbuf[2048];
    if (getcwd(pathbuf, sizeof(pathbuf))==NULL) std::cerr << "Could not find the current PID \n";
    else path = pathbuf;
    return path;
}

double KpcPerArc(double d) {return 2*d*tan(M_PI/2/180)/3.6;}
double VeltoDist(double vsys) {return vsys/70.;}
double RedtoDist(double redshift) {return redshift*299792.458/70.;}


bool fexists(std::string filename) {std::ifstream ifile(filename.c_str()); return ifile.good();}


bool mkdirp(const char* path, mode_t mode) {
    /** Utility function to create directory tree */

    if(path[0] == '\0') return false;

    char* p = const_cast<char*>(path);
    while (*p != '\0') {
        p++;
        while(*p != '\0' && *p != '/') p++;
        char v = *p;
        *p = '\0';

        if(mkdir(path, mode) != 0 && errno != EEXIST) {
            *p = v;
            return false;
        }
        *p = v;
    }

    return true;
}


template <class T> 
T AlltoVel (T in, Header &h) {
    
  /// This function convert a spectral input value "in" in 
  /// a output value in units of KM/S.
  /// No errors if can not convert, just return input value
    
    std::string cunit2 = makelower(h.Cunit(h.NumAx()-1));
    if (h.NumAx()>3) cunit2 = makelower(h.Cunit(2));

    const double c = 299792458;
    T vel_km_s = in;

    if (cunit2=="km/s" || cunit2=="kms") vel_km_s = in;
    else if (cunit2=="m/s" || cunit2=="ms") vel_km_s = in/1000.;
    else if (cunit2=="hz" || cunit2=="mhz") {
        double frest = h.Freq0()/(1+h.Redshift());
        T vel_m_s = c*(frest*frest-in*in)/(frest*frest+in*in);
        vel_km_s = vel_m_s/1000.;
    }
    else if (cunit2=="mum" || cunit2=="um" || cunit2=="micron" ||
             cunit2=="a" || cunit2=="ang"  || cunit2=="angstrom") {          // Micron or Angstrom for Hi-Z
        //int z_cent = floor(h.DimAx(2)/2.)-1;
        double z_cent = h.Crpix(2)-1;
        double restw = h.Wave0(), reds = h.Redshift();
        if (restw!=-1) {
            z_cent = (restw*(1+reds)-h.Crval(2))/h.Cdelt(2)+h.Crpix(2)-1;
        }
        double line_wave = (z_cent+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
        
        T vel_m_s = c*(in*in-line_wave*line_wave)/(in*in+line_wave*line_wave);
        vel_km_s = vel_m_s/1000.;
    }
    
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
    
    T deltaV = h.Cdelt(2);
    long zdim = h.DimAx(2);
    std::string cunit2 = makelower(h.Cunit(2));

    const double c = 299792458;

    if (cunit2=="km/s" || cunit2=="kms") deltaV = h.Cdelt(2);
    else if (cunit2=="m/s" || cunit2=="ms") {
        deltaV =h.Cdelt(2)/1000.;
    }
    else if (cunit2=="hz" || cunit2=="mhz") {
        double frest = h.Freq0()/(1+h.Redshift());
            
        T sum=0;
        for (int i=0; i<zdim-1; i++) {
            T ipix = (i+1-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
            T ivel = c*(frest*frest-ipix*ipix)/(frest*frest+ipix*ipix);
            T spix = (i+2-h.Crpix(2))*h.Cdelt(2)+h.Crval(2);
            T svel = c*(frest*frest-spix*spix)/(frest*frest+spix*spix);
            T diff = svel - ivel;
            sum += diff;
        }
        deltaV = sum/(zdim-1);
        deltaV /= 1000; 
    }
    else if (cunit2=="mum" || cunit2=="um" || cunit2=="micron" ||
             cunit2=="a" || cunit2=="ang"  || cunit2=="angstrom") {            // Micron or Angstrom for Hi-Z

        double z_cent = h.DimAx(2)/2.-1;
        double restw = h.Wave0(), reds = h.Redshift();
        if (restw!=-1) {
            z_cent = (restw*(1+reds)-h.Crval(2))/h.Cdelt(2)+h.Crpix(2)-1;
        }

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
 
    T fluxJY = in;
    std::string b = deblankAll(makelower(h.Bunit()));
    size_t f = std::string::npos;
    
    if (b.find("w.u.")!=f || b.find("wu")!=f) {
        fluxJY *= 5E-3;
        if (h.BeamArea()!=0) fluxJY /= h.BeamArea();
    }
    else if (b.find("jy/b")!=f || b.find("j/b")!=f) {
        if (h.BeamArea()!=0) fluxJY /= h.BeamArea();
    }
    else if (b.find("jy")!=f) {
        return fluxJY; 
    }
    else if (b=="k") {
        // Converting from Kelvin -> Jy/Beam -> Jy
        fluxJY = in*(h.Bmaj()*3600.*h.Bmin()*3600.)/(1360.*21.106*21.106);
        if (h.BeamArea()!=0) fluxJY /= h.BeamArea();
    }
    
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
                else {  // UNKNOWN TYPE -- DEFAULT TO RA.
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
template float angularSeparation(float&,float&,float&,float&);
template double angularSeparation(double&,double&,double&,double&);


double arcsconv(std::string cunit) {
    
    std::string Cunit = makelower(cunit);
    if (Cunit=="degree" || Cunit=="degrees" || Cunit=="deg") return 3600.;
    else if (Cunit=="arcmin" || Cunit=="arcm") return 60.;
    else if (Cunit=="arcsec" || Cunit=="arcs") return 1.;
    else {
        std::cout << "Conversion error (unknown CUNIT for RA-DEC): ";
        std::cout << "cannot convert to ARCSEC.\n";
        std::cout << cunit;
        std::terminate(); 
    }
}


double degconv(std::string cunit) {

    std::string Cunit = makelower(cunit);
    if (Cunit=="degree" || Cunit=="degrees" || Cunit=="deg") return 1.;
    else if (Cunit=="arcmin" || Cunit=="arcm") return 1./60.;
    else if (Cunit=="arcsec" || Cunit=="arcs") return 1./3600.;
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
    
    float pcRA = h.Crval(0);                /// Pointing center R.A. (in degrees).
    float pcDEC = h.Crval(1);               /// Pointing center DEC  (in degrees).
    
    if (pcDEC>=0) pcDEC = 90 - pcDEC;           /// Conversion of DEC to polar angle.
    else pcDEC = 90 + fabs(pcDEC);
    
    pcRA  *= degtorad;                          /// Conversions to radians.
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
    
    float freq=0;                               /// Frequency in GHz
    
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
        float RF = angdist*60*freq;                 /// The distance is now in arcmin.
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
        float RF = angdist*60*freq;                 /// The distance is now in arcmin.
        if (RF>50) return 0;                        /// Cut-off at 50 arcmin*GHz
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
        float RF = angdist*60*freq;                 /// The distance is now in arcmin.
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
        if (angdist>2.8) return 0;              /// Cut-off at 2.8 degrees.
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


template <class T> 
Rings<T>* readRings(GALFIT_PAR &par, Header &h) {
    
    const short NPAR = 14;
    
    // Try to read ring information from an input file
    Rings<T> fr;
    bool radii_b = getDataColumn(fr.radii,par.RADII);
    bool xpos_b  = getDataColumn(fr.xpos,par.XPOS);
    bool ypos_b  = getDataColumn(fr.ypos,par.YPOS);
    bool vsys_b  = getDataColumn(fr.vsys,par.VSYS);
    bool vrot_b  = getDataColumn(fr.vrot,par.VROT);
    bool vrad_b  = getDataColumn(fr.vrad,par.VRAD);
    bool vvert_b = getDataColumn(fr.vvert,par.VVERT);
    bool dvdz_b  = getDataColumn(fr.dvdz,par.DVDZ);
    bool zcyl_b  = getDataColumn(fr.zcyl,par.ZCYL);
    bool vdisp_b = getDataColumn(fr.vdisp,par.VDISP);
    bool z0_b    = getDataColumn(fr.z0,par.Z0);
    bool dens_b  = getDataColumn(fr.dens,par.DENS); 
    bool inc_b   = getDataColumn(fr.inc,par.INC);
    bool pa_b    = getDataColumn(fr.phi,par.PHI);

    // Determine maximum number of rings (= minimum number of entries)
    size_t s[NPAR] = {fr.radii.size(),fr.xpos.size(), fr.ypos.size(),fr.vsys.size(),
                      fr.vrot.size(),fr.vrad.size(),fr.vvert.size(),fr.dvdz.size(),
                      fr.zcyl.size(),fr.vdisp.size(),fr.z0.size(),fr.dens.size(),
                      fr.inc.size(),fr.phi.size()};
    size_t max_size=UINT_MAX;
    for (short i=0; i<NPAR; i++) if (s[i]!=0 && s[i]<max_size) max_size=s[i];

    // Reading parameter rings without caring if they are file or not
    T vsys        = par.VSYS!="-1" ? atof(par.VSYS.c_str()) : 0;
    T vrot        = par.VROT!="-1" ? atof(par.VROT.c_str()) : 0;
    T vrad        = par.VRAD!="-1" ? atof(par.VRAD.c_str()) : 0;
    T vvert       = par.VVERT!="-1"? atof(par.VVERT.c_str()): 0;
    T dvdz        = par.DVDZ!="-1" ? atof(par.DVDZ.c_str()) : 0;
    T zcyl        = par.ZCYL!="-1" ? atof(par.ZCYL.c_str()) : 0;
    T vdisp       = par.VDISP!="-1"? atof(par.VDISP.c_str()): 0; 
    T z0          = par.Z0!="-1"   ? atof(par.Z0.c_str())   : 0;
    T dens        = par.DENS!="-1" ? atof(par.DENS.c_str()) : 1;
    T inc         = par.INC!="-1"  ? atof(par.INC.c_str())  : 0;
    T pa          = par.PHI!="-1"  ? atof(par.PHI.c_str())  : 0;

    T xpos = 0, ypos = 0;
    std::string pos[2] = {par.XPOS, par.YPOS};
    if (pos[0]!="-1" && pos[1]!="-1") {
        double *pixs  = getCenterCoordinates(pos, h);
        xpos = pixs[0];
        ypos = pixs[1];
    }
    
    // Setting number of rings and radsep
    size_t nr     = par.NRADII;
    double radsep = par.RADSEP;    
    nr = nr>0 && nr<max_size ? nr : max_size;
    if (radii_b) {
        radsep = 0;
        for (unsigned i=1; i<fr.radii.size()-1; i++)
            radsep += fr.radii[i+1]-fr.radii[i];
        radsep/=(fr.radii.size()-2);
    }
    
    // Filling rings with values
    Rings<T> *inR = new Rings<T>;
    inR->nr     = nr;
    inR->radsep = radsep;
    for (int i=0; i<inR->nr; i++) {
        if (radii_b) inR->radii.push_back(fr.radii[i]);
        else inR->radii.push_back(i*radsep+radsep/2.);
        if (vrot_b) inR->vrot.push_back(fr.vrot[i]);
        else inR->vrot.push_back(vrot);
        if (vrad_b) inR->vrad.push_back(fr.vrad[i]);
        else inR->vrad.push_back(vrad);
        if (vvert_b) inR->vvert.push_back(fr.vvert[i]);
        else inR->vvert.push_back(vvert);
        if (dvdz_b) inR->dvdz.push_back(fr.dvdz[i]);
        else inR->dvdz.push_back(dvdz);
        if (zcyl_b) inR->zcyl.push_back(fr.zcyl[i]);
        else inR->zcyl.push_back(zcyl);
        if (vdisp_b) inR->vdisp.push_back(fr.vdisp[i]);
        else inR->vdisp.push_back(vdisp);
        if (z0_b) inR->z0.push_back(fr.z0[i]);
        else inR->z0.push_back(z0);
        if (dens_b) inR->dens.push_back(fr.dens[i]*1.E20);
        else inR->dens.push_back(dens*1.E20);
        if (inc_b) inR->inc.push_back(fr.inc[i]);
        else inR->inc.push_back(inc);
        if (pa_b) inR->phi.push_back(fr.phi[i]);
        else inR->phi.push_back(pa);
        if (xpos_b) inR->xpos.push_back(fr.xpos[i]);
        else inR->xpos.push_back(xpos);
        if (ypos_b) inR->ypos.push_back(fr.ypos[i]);
        else inR->ypos.push_back(ypos);
        if (vsys_b) inR->vsys.push_back(fr.vsys[i]);
        else inR->vsys.push_back(vsys);
    }
    
    return inR;
}
template Rings<float>* readRings(GALFIT_PAR &, Header &);
template Rings<double>* readRings(GALFIT_PAR &, Header &);


double* getCenterCoordinates(std::string *pos, Header &h) {

    double *pixels = new double[h.NumAx()];
    double world[h.NumAx()];
    bool isPOS[h.NumAx()];
    for (int i=0; i<h.NumAx(); i++) {
        world[i] = 0;
        isPOS[i] = false;
    }
    
    for (int i=0; i<2; i++) {
        std::string coord_str = pos[i];
        std::string coord_typ = makelower(h.Ctype(i));
        if (coord_str.find('d')!=std::string::npos) {
            std::string substr = coord_str.erase(coord_str.find('d'),coord_str.size()-1);
            world[i] = atof(substr.c_str());
        }
        else if (coord_str.find(':')!=std::string::npos) {      // Found sexagesimal
            double pos_deg = dmsToDec(coord_str);
            if (coord_typ.find("ra")!=std::string::npos) pos_deg*=15;
            world[i] = pos_deg;
        }
        else isPOS[i]=true;
    }

    if (isPOS[0] && isPOS[1]) {
        pixels[0]=atof(pos[0].c_str());
        pixels[1]=atof(pos[1].c_str());
    }
    else if (!isPOS[0] && !isPOS[1]) {
        wcsToPixSingle(h.WCS(),world,pixels);
    }
    else {
        std::cerr << "WCS ERROR: please provide both center coordinates. \n";
        std::terminate();
    }
    return pixels;
}


template <class T>
T* RingRegion (Rings<T> *r, Header &h) {
	
    // Given a set of rings, return the 2D region covered by these rings
    
	long bsize[2] = {h.DimAx(0),h.DimAx(1)};
    float pscale = h.PixScale()*arcsconv(h.Cunit(0));
    T *ringregion = new T[bsize[0]*bsize[1]];
	for (int i=0;i<bsize[0]*bsize[1];i++) ringregion[i]=log(-1);	
	
    T R1  = std::max((r->radii.front()-r->radsep/2.)/pscale,0.); //#+sqrt(in->Head().BeamArea()/M_PI);
    T R2  = (r->radii.back()+r->radsep/2.)/pscale;    
	T phi = r->phi.back();
	T inc = r->inc.back();
	T psi = 0.;
    T z0  = 3*r->z0.back()/(pscale); //prima prendevo 3*dring->....
	T x0  = r->xpos.back()-1;
	T y0  = r->ypos.back()-1;
	
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
						if(isNaN(ringregion[xrot+yrot*bsize[0]])) {
							ringregion[xrot+yrot*bsize[0]] = theta;
						}
					}
				}
			}
		}
	}

	deallocate_2D<double>(matrices,3);
	delete [] rotmatrix;
	return ringregion;
	
}
template float* RingRegion (Rings<float>*,Header&);
template double* RingRegion (Rings<double>*,Header&);


template <class T>
T* SimulateNoise(double stddev, size_t size) {
    
    T *noise = new T[size];
    // Random engine generator
    auto const seed = std::random_device()();
    std::mt19937 gen(seed);
    // Gaussian distribution draw
    std::normal_distribution<T> gaussian(0.,stddev);
    // Filling the vector with random draws
    for (size_t i=size; i--;) noise[i] = gaussian(gen);
    return noise;
}
template float* SimulateNoise(double,size_t);
template double* SimulateNoise(double,size_t);


template <class T>
T* Smooth1D(T *inarray,size_t npts,std::string windowType,size_t windowSize) {
    
    // Performs smoothing on a single 1D array. Accepted windows are below 
    
    bool known_window = windowType=="HANNING" || windowType=="HANNING2" ||
                        windowType=="BOXCAR"  || windowType=="TOPHAT"   || 
                        windowType=="FLATTOP" || windowType=="BARTLETT" ||
                        windowType=="WELCH"   || windowType=="BLACKMAN";
    if (!known_window) {
         std::cerr << "Smoothing 1D: window type unknown "
                   << "Changing "<< windowType << " to \"HANNING\" .\n";
         windowType = "HANNING";
    }
    // Check window size and type
    if(windowSize%2==0){ 
      std::cerr << "Smoothing 1D: need an odd number for the window size. "
  	            << "Changing "<< windowSize << " to " << windowSize+1<<".\n";
      windowSize++;
    }
   
    // Defining coefficients for smoothing
    double *coeff = new double[windowSize];
    float scale = (windowSize+1.)/2.;
    float N = windowSize-1;
    float sum = 0;
    for(size_t j=0; j<windowSize; j++) {
        float x = j-(windowSize-1)/2.;
        if (windowType=="HANNING")            // Hanning used in radio-astronomy
            coeff[j] = 0.5+0.5*cos(x*M_PI/scale);
        else if (windowType=="HANNING2")      // Classical Hanning window
            coeff[j] = 0.5-0.5*cos(j*2*M_PI/N);
        else if (windowType=="BARTLETT")      // Triangular window
            coeff[j] = 1-fabs((j-N/2.)/(N/2.));
        else if (windowType=="WELCH")
            coeff[j] = 1-std::pow((j-N/2.)/(N/2.),2);
        else if (windowType=="BLACKMAN") 
            coeff[j] = 0.42659-0.49656*cos(j*2*M_PI/N)+0.076849*cos(j*4*M_PI/N);
        else if (windowType=="FLATTOP")
            coeff[j] =  0.21557895-0.416631580*cos(j*2*M_PI/N)+0.277263158*cos(j*4*M_PI/N)
                                  -0.083578957*cos(j*6*M_PI/N)+0.006947368*cos(j*8*M_PI/N);
        else coeff[j] = 1.;
        sum += coeff[j]; 
    }
    
    // Smooth
    T *newarray = new T[npts];
    for(size_t i=0; i<npts; i++){
        newarray[i] = 0.;
        for(size_t j=0; j<windowSize; j++){
            float x = j-(windowSize-1)/2.;
            if((i+x>0)&&(i+x<npts)) newarray[i] += coeff[j]/sum*inarray[i+int(x)];
        }
    }
    delete [] coeff;
    return newarray;
}
template float* Smooth1D(float*,size_t,std::string,size_t);
template double* Smooth1D(double*,size_t,std::string,size_t);


