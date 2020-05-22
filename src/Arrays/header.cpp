//-----------------------------------------------------------------------
// header.cpp: Member functions for the Header class.
//-----------------------------------------------------------------------

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
#include <cstring>
#include <sstream>
#include <cmath>
#include <cctype>
#include <vector>
#include <algorithm>
#include <fitsio.h>
#include <wcslib/wcs.h>
#include <wcslib/wcsunits.h>
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>
#include <Arrays/header.hh>

Header::Header () {
    
    bitpix = FLOAT_IMG;
    numAxes = bmaj = bmin = bpa = beamArea = freq0 = 0.;
    wave0 = -1;
    datamin = datamax = redshift = crota = 0.;
    dunit3 = "";
    object = "NONE";
    pointAllocated = false;
    warning = true;
    wcs = new struct wcsprm;
    wcs->flag=-1;
    wcsIsGood=false;
    nwcs=0;
}


Header::~Header () {
    
    if (pointAllocated) {
        delete [] dimAxes;              
        delete [] crpix;                    
        delete [] crval;                    
        delete [] cdelt;
        delete [] ctype;
        delete [] cunit;
    }
    
    keys.clear();
    wcsvfree(&nwcs,&wcs);
}


Header::Header(const Header& h) {
    
    operator=(h);
}


Header& Header::operator=(const Header& h) {
    
    if(this == &h) return *this;
    this->numAxes = h.numAxes;
    this->bitpix  = h.bitpix;

    if (pointAllocated) {
        delete [] this->dimAxes;
        delete [] this->crpix;
        delete [] this->crval;
        delete [] this->cdelt;
        delete [] this->ctype;
        delete [] this->cunit;
    }
    
    this->pointAllocated = h.pointAllocated;
    if(this->pointAllocated) {
        this->dimAxes = new long[numAxes];
        this->crpix = new double[numAxes];
        this->crval = new double[numAxes];
        this->cdelt = new double[numAxes];
        this->cunit = new std::string[numAxes];
        this->ctype = new std::string[numAxes];
        for (int i=0; i<numAxes; i++) {
            this->dimAxes[i] = h.dimAxes[i];
            this->crpix[i] = h.crpix[i];
            this->crval[i] = h.crval[i];
            this->cdelt[i] = h.cdelt[i];
            this->ctype[i] = h.ctype[i];
            this->cunit[i] = h.cunit[i];
        }
    }   

    this->bmaj      = h.bmaj;
    this->bmin      = h.bmin;
    this->bpa       = h.bpa;
    this->bzero     = h.bzero;
    this->bscale    = h.bscale;
    this->blank     = h.blank;
    this->beamArea  = h.beamArea;
    this->epoch     = h.epoch;
    this->freq0     = h.freq0;
    this->wave0     = h.wave0;
    this->redshift  = h.redshift;
    this->fitsname  = h.fitsname;
    this->bunit     = h.bunit;  
    this->btype     = h.btype;
    this->object    = h.object;
    this->telescope = h.telescope;
    this->dunit3    = h.dunit3;
    this->drval3    = h.drval3;
    this->datamin   = h.datamin;
    this->datamax   = h.datamax;
    this->warning   = h.warning;
    
//    this->wcs = h.wcs;
#pragma omp critical (wcs_header)
{
    this->wcs = new struct wcsprm;
    this->wcs->flag = -1;
    wcsini(true, h.wcs->naxis, this->wcs);
    wcscopy(true, h.wcs, this->wcs);
    wcsset(this->wcs);
    this->nwcs      = h.nwcs;
    this->wcsIsGood = h.wcsIsGood;
}

    for (unsigned int i=0; i<h.keys.size(); i++)
        this->keys.push_back(h.keys[i]); 
            
    calcArea();
    
    return *this;

}


void Header::setNumAx (int size){

    if (pointAllocated) {
        delete [] dimAxes;
        delete [] crpix;
        delete [] crval;
        delete [] cdelt;
        delete [] cunit;
        delete [] ctype;
    }
    numAxes = size;
    dimAxes = new long[numAxes];
    crpix   = new double[numAxes];
    crval   = new double[numAxes];
    cdelt   = new double[numAxes];
    cunit   = new std::string[numAxes];
    ctype   = new std::string[numAxes];
    pointAllocated = true;

    wcsini(true, numAxes, wcs);
}


void Header::calcArea () {
    
    float AvPixScale = (fabs(cdelt[0])+fabs(cdelt[1]))/2.;
    float scalbmaj = bmaj/AvPixScale;
    float scalbmin = bmin/AvPixScale;
    beamArea =  M_PI * (scalbmaj/2.) * (scalbmin/2.) / M_LN2;

}


bool Header::header_read (std::string fname) {
        
    fitsfile *fptr;
    int status=0, nfound;
    char comment[72];
    

    fitsname = fname;
    char filename[100];
    strcpy(filename, fname.c_str());
 
    if (fits_open_file(&fptr, filename, READONLY, &status)) {
        fits_report_error(stderr, status);
        return false;
    }
    
    status=0;
    if (fits_get_img_type(fptr, &bitpix, &status)) {
        fits_report_error(stderr, status);
        bitpix = FLOAT_IMG;
    }
        
    status=0;
    if (fits_get_img_dim(fptr, &numAxes, &status))
        fits_report_error(stderr, status);

    if (!pointAllocated) {
        dimAxes = new long [numAxes];
        crpix = new double [numAxes];
        crval = new double [numAxes];
        cdelt = new double [numAxes];
        ctype = new std::string[numAxes];
        cunit = new std::string[numAxes];
    }
    pointAllocated=true;
    
    status = 0;
    if(fits_get_img_size(fptr, numAxes, dimAxes, &status)){
      fits_report_error(stderr, status);
    }

    char Bunit[20], Btype[20], name[20], Tel[20], Dunit3[20], Keys[100];
    for (int i=0; i<20; i++) {
            Bunit[i]=Btype[i]=name[i]=Tel[i]=Dunit3[i]=Keys[i]=' ';
    }
    for (int i=20; i<100; i++) Keys[i]=' ';

    int nkeys;
    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    for (int i=1; i<=nkeys; i++) {
        fits_read_record(fptr,i,Keys,&status);
        keys.push_back(Keys);
    }
    
    status=0;
    fits_read_keys_dbl (fptr, "CDELT", 1, numAxes, cdelt, &nfound, &status);
    if (nfound==0) {
        Warning("HEADER WARNING: CDELTs keywords not found.");
        for (int i=0; i<numAxes; i++) {
            std::string keyw = "CD"+std::to_string(i+1)+"_"+std::to_string(i+1);
            if (fits_read_key_dbl (fptr, keyw.c_str(), &cdelt[i], comment, &status)) {
                fits_report_error(stderr, status);
                Warning("HEADER ERROR: No pixel spacing keywords found.");
                return false;
            }
        }
        Warning("HEADER WARNING: Found CDn_n keywords instead!");
    }
    
    status=0;
    fits_read_keys_dbl (fptr, "CRPIX", 1, numAxes, crpix, &nfound, &status);
    if (nfound==0) {
        Warning("HEADER ERROR: CRPIXs keywords not found.");
        return false;
    }
        
    status=0;
    fits_read_keys_dbl (fptr, "CRVAL", 1, numAxes, crval, &nfound, &status);
    if (nfound==0) {
        Warning("HEADER ERROR: CRVALSs keywords not found.");
        return false;
    }
        
    status=0;
    if (fits_read_key_dbl (fptr, "DRVAL3", &drval3, comment, &status)) {
        if (status!=202) fits_report_error(stderr, status);
        drval3=0;
    }

    char **Ctype = new char*[numAxes];
    char **Cunit = new char*[numAxes];
    for (int i=0; i<numAxes; i++) {
        Ctype[i] = new char[25];
        Cunit[i] = new char[25];
    }

    status=0;
    fits_read_keys_str (fptr, "CTYPE", 1, numAxes, Ctype, &nfound, &status);
    if (nfound==0) {
        Warning("HEADER WARNING: CTYPEs keywords not found. Assuming [RA,DEC,VELO].");
        if (numAxes>0) ctype[0] = "RA---NCP";
        if (numAxes>1) ctype[1] = "DEC--NCP";
        if (numAxes>2) ctype[2] = "VELO-HELO";
    }
    else {
        for (int i=0; i<numAxes; i++) {
            ctype[i] = Ctype[i];
        }
    }

    if (ctype[0].find("RA")!=std::string::npos && crval[0]<0) crval[0]+=360.;

    status=0;
    fits_read_keys_str (fptr, "CUNIT", 1, numAxes, Cunit, &nfound, &status);
    if (nfound==0) {
        if (numAxes>0) cunit[0] = "DEGREE";
        if (numAxes>1) cunit[1] = "DEGREE";
        if (numAxes>2) {
            if (ctype[2]=="FREQ" || ctype[2]=="freq" || ctype[2]=="Freq") {
                cunit[2] = "HZ";
                Warning("HEADER WARNING: CUNITs keywords not found. Assuming [DEGREE,DEGREE,HZ]");
            }
            else {
                cunit[2] = "M/S";
                Warning("HEADER WARNING: CUNITs keywords not found. Assuming [DEGREE,DEGREE,M/S]");
            }
        }
    }
    else {
        for (int i=0; i<numAxes; i++) {
            cunit[i] = Cunit[i];
            // Fixing "s-1" problem in CUNIT3
            if (cunit[i].find("km s-1")!=std::string::npos) cunit[i] = "km/s";
            if (cunit[i].find("m s-1")!=std::string::npos) cunit[i] = "m/s";
            if (cunit[i]=="") {
                cunit[i] = i==2 ? "km/s" : "DEGREE";
                std::stringstream toprint;
                toprint << "HEADER WARNING: CUNIT" << i+1 << " keywords not found. Assuming " << cunit[i]; 
                Warning(toprint.str());
            }
            
        }
    }

    for (int i=0; i<numAxes; i++) {
        delete [] Ctype[i];
        delete [] Cunit[i];
    }
    delete [] Cunit;
    delete [] Ctype;
        
    status=0;
    if (fits_read_key_str (fptr, "DUNIT3", Dunit3, comment, &status)) {
        if (status!=202) fits_report_error(stderr, status);
        dunit3 = "NONE";
    }
    else dunit3 = Dunit3;

    status=0;
    if (fits_read_key_str (fptr, "BUNIT", Bunit, comment, &status)) {
        if (status==202)
            Warning("HEADER WARNING: BUNIT keyword not found.");
        else fits_report_error(stderr, status);
        bunit = "NONE";
    }
    else bunit = Bunit;
    if (bunit.find("beam-1 Jy")!=std::string::npos ||
        bunit.find("Jy beam-1")!=std::string::npos) 
            bunit = "JY/BEAM";
    
    status=0;
    if (fits_read_key_str (fptr, "BTYPE", Btype, comment, &status)) {
        btype = "NONE";
    }
    else btype = Btype;
    
    status=0;
    if (fits_read_key_flt (fptr, "BZERO", &bzero, comment, &status)) {
        bzero = 0;
    }

    status=0;
    if (fits_read_key_flt (fptr, "BSCALE", &bscale, comment, &status)) {
        bscale = 1;
    }
    
    status=0;
    if (fits_read_key_flt (fptr, "BLANK", &blank, comment, &status)) {
        blank = 0;
    }
    
    status=0;
    if (fits_read_key_flt (fptr, "EPOCH", &epoch, comment, &status)) {
        if (fits_read_key_flt (fptr, "EQUINOX", &epoch, comment, &status))
        epoch = 0;
    }
    
    status=0;
    if (fits_read_key_dbl (fptr, "DATAMIN", &datamin, comment, &status)) {
        datamin = 0;
    }

    status=0;
    if (fits_read_key_dbl (fptr, "DATAMAX", &datamax, comment, &status)) {
        datamax = 0;
    }
    
    status=0;
    if (fits_read_key_dbl (fptr, "FREQ0", &freq0, comment, &status)) {
        status=0;
        if (fits_read_key_dbl (fptr, "RESTFREQ", &freq0, comment, &status)) {
            status=0;
            if (fits_read_key_dbl (fptr, "RESTFRQ", &freq0, comment, &status)) {
                if (dunit3=="NONE" || drval3==0 || (cunit[2]!="HZ" && cunit[2]!="hz" &&
                    cunit[2]!="MHZ" && cunit[2]!="Mhz" && cunit[2]!="mhz" &&
                    cunit[2]!="GHZ" && cunit[2]!="Ghz" && cunit[2]!="ghz" &&
                    dunit3!="KM/S" && dunit3!="km/s" && dunit3!="M/S"  && dunit3!="m/s")) {
                    Warning("HEADER WARNING: FREQ0-RESTFREQ keyword not found. Assuming 1.4204057 GHz.");
                    freq0 = 0.1420405751786E10;
                }
                else {
                    double drval3ms=0.;
                    double crval3hz=0.;
                    if (dunit3=="KM/S" || dunit3=="km/s") drval3ms=drval3*1000;
                    else if (dunit3=="M/S" || dunit3=="m/s") drval3ms=drval3;
                    if (cunit[2]=="HZ" || cunit[2]=="hz") crval3hz = crval[2];
                    else if (cunit[2]=="MHZ" || cunit[2]=="Mhz" || cunit[2]=="mhz") crval3hz=crval[2]*1.E06;
                    else if (cunit[2]=="GHZ" || cunit[2]=="Ghz" || cunit[2]=="ghz") crval3hz=crval[2]*1.E09;
                    freq0 = crval3hz*sqrt((299792458.+drval3ms)/(299792458.-drval3ms));
                }
            }
        }
    }

    status=0;
    if (fits_read_key_dbl (fptr, "CROTA", &crota, comment, &status)) {
        crota = 0;
    }

    status=0;
    if (fits_read_key_str (fptr, "OBJECT", name, comment, &status)) {
        if (status==202) Warning("HEADER WARNING: OBJECT keyword not found.");
        else fits_report_error(stderr, status);
        object = "NONE";
    }
    else object = name;

    for (uint i=0; i<object.size();i++) {
        if (object[i]=='/' || object[i]=='\\') object.replace(i,1, "-");
        if (object[i]=='(' || object[i]==')')  object.replace(i,1, "-");
        if (isspace(object[i])) object.replace(i,1, "_");
    }    

    status = 0;
    if (fits_read_key_str (fptr, "TELESCOP", Tel, comment, &status)) {
        status = 0;
        if (fits_read_key_str (fptr, "INSTRUME", Tel, comment, &status)) {
            if (status==202) Warning("HEADER WARNING: TELESCOP-INSTRUME keywords not found.");
            else fits_report_error(stderr, status);
            telescope = "NONE";
        }
        else telescope = Tel;
    }
    else telescope = Tel;
    
    
    double clbmaj=0, clbmin=0;
    status=0;
    if (fits_read_key_dbl (fptr, "BMAJ", &bmaj, comment, &status)) {
        status = 0;
        if (fits_read_key_dbl (fptr, "BMMAJ", &bmaj, comment, &status)) {
            status = 0;
            if (fits_read_key_dbl (fptr, "CLBMMAJ", &clbmaj, comment, &status)) {
                if (status==202) Warning("HEADER WARNING: BMAJ-BMMAJ-CLBMMAJ keywords not found.");
                else fits_report_error(stderr, status);
                bmaj = 0;
            }
        }
    }
    
    status=0;
    if (fits_read_key_dbl (fptr, "BMIN", &bmin, comment, &status)) {
        status = 0;
        if (fits_read_key_dbl (fptr, "BMMIN", &bmin, comment, &status)) {
            status = 0;
            if (fits_read_key_dbl (fptr, "CLBMMIN", &clbmin, comment, &status)) {
                if (status==202) Warning("HEADER WARNING: BMIN-BMMIN-CLBMMIN keywords not found.");
                else fits_report_error(stderr, status);
                bmin = 0;
            }
        }
    }
    
    status=0;
    if (fits_read_key_dbl (fptr, "BPA", &bpa, comment, &status)) {
        status = 0;
        if (fits_read_key_dbl (fptr, "BMPA", &bpa, comment, &status)) {
            if (status==202) Warning("HEADER WARNING: BPA-BMPA keywords not found.");
            else fits_report_error(stderr, status);
            bpa = 0;
        }
    }
    
    if (bmaj==0 && bmin==0 && bpa==0 && clbmaj==0 && clbmin==0) {
        for (unsigned int i=0; i<keys.size(); i++) {
            int found = keys[i].find("BMAJ=");
            char *pEnd;
            if (found>=0) {
                pEnd = &keys[i].at(found+5);
                bmaj = strtod(pEnd,NULL);
            }
            found = keys[i].find("BMIN=");
            if (found>=0) {
                pEnd = &keys[i].at(found+5);
                bmin = strtod(pEnd,NULL);

            }
            found = keys[i].find("BPA=");
            if (found>=0) {
                pEnd = &keys[i].at(found+4);
                bpa = strtod(pEnd,NULL);
            }
        }
        if (bmaj!=0 && bmin!=0) {
            std::cout << "\n--------> WARNING: beam information found in HISTORY keywords: <--------\n"
                      << " BMAJ = " << bmaj << " " << cunit[0]
                      << "  BMIN = " << bmin << " " << cunit[0]
                      << "  BPA = "  << bpa  << " DEGREE\n";
            std::cout << " It is heartly recommended to check these values before going on!!\n\n";
        }
    }    
    
    char bm[30];
    status=0;
    fits_read_key_str (fptr, "BMMAJ", bm, comment, &status);
    std::string bmstr = bm;
    bool arcsecbeam = false;
    int found = bmstr.find("D");
    if (found>=0) arcsecbeam = true;
    if (arcsecbeam) bmaj /= 3600.;
    
    status=0;
    fits_read_key_str (fptr, "BMMIN", bm, comment, &status);
    bmstr = bm;
    arcsecbeam = false;
    found = bmstr.find("D");
    if (found>=0) arcsecbeam = true;
    if (arcsecbeam) bmin /= 3600.;
    
    if (clbmaj!=0) bmaj = clbmaj/3600.;
    if (clbmin!=0) bmin = clbmin/3600.;


    int noComments = 1;     // fits_hdr2str will ignore COMMENT, HISTORY etc
    int nExc = 0;
    char *hdr=0;

    // Read in the entire PHU of the FITS file to a std::string.
    // This will be read by the wcslib functions to extract the WCS.
    status = 0;
    fits_hdr2str(fptr, noComments, NULL, nExc, &hdr, &nkeys, &status);

    status = wcsini(true, numAxes, wcs);

    int relax=1; // for wcspih -- admit all recognised informal WCS extensions
    int ctrl=2;  // for wcspih -- report each rejected card and its reason for rejection
    int nreject;
    // Parse the FITS header to fill in the wcsprm structure
    status=wcspih(hdr, nkeys, relax, ctrl, &nreject, &nwcs, &wcs);


    int stat[NWCSFIX];
    // Applies all necessary corrections to the wcsprm structure
    //  (missing cards, non-standard units or spectral types, ...)
    status = wcsfix(1, (const int*)dimAxes, wcs, stat);

    // Set up the wcsprm struct. Report if something goes wrong.
    status = wcsset(wcs);
    // Re-do the corrections to account for things like NCP projections
    status = wcsfix(1, (const int*)dimAxes, wcs, stat);

    char stype[5],scode[5],sname[22],units[8],ptype,xtype;
    int restreq;
    status = spctyp(wcs->ctype[wcs->spec],stype,scode,sname,units,&ptype,&xtype,&restreq);

    if((wcs->lng!=-1) && (wcs->lat!=-1)) wcsIsGood = true;

    // Close the FITS File
    status=0;
    if (fits_close_file(fptr, &status))
        fits_report_error(stderr, status);
    
    calcArea();

    return true;

}   


void Header::headwrite_3d (fitsfile *fptr, bool fullHead) {

    int status=0;
    char com[]= "  ";
    
    fits_update_key_dbl(fptr, "CRPIX1", crpix[0], 10, com, &status);
    fits_update_key_dbl(fptr, "CRVAL1", crval[0], 10, com, &status);
    fits_update_key_dbl(fptr, "CDELT1", cdelt[0], 10, com, &status);
    fits_update_key_str(fptr, "CTYPE1", ctype[0].c_str(), com, &status);
    fits_update_key_str(fptr, "CUNIT1", cunit[0].c_str(), com, &status);
    
    fits_update_key_dbl(fptr, "CRPIX2", crpix[1], 10, com, &status);
    fits_update_key_dbl(fptr, "CRVAL2", crval[1], 10, com, &status);
    fits_update_key_dbl(fptr, "CDELT2", cdelt[1], 10, com, &status);
    fits_update_key_str(fptr, "CTYPE2", ctype[1].c_str(), com, &status);
    fits_update_key_str(fptr, "CUNIT2", cunit[1].c_str(), com, &status);    
    
    fits_update_key_dbl(fptr, "CRPIX3", crpix[2], 10, com, &status);
    fits_update_key_dbl(fptr, "CRVAL3", crval[2], 10, com, &status);
    fits_update_key_dbl(fptr, "CDELT3", cdelt[2], 10, com, &status);
    fits_update_key_str(fptr, "CTYPE3", ctype[2].c_str(), com, &status);
    fits_update_key_str(fptr, "CUNIT3", cunit[2].c_str(), com, &status);
    if (drval3!=0) fits_update_key_dbl(fptr, "DRVAL3", drval3, 10, com, &status);
    if (dunit3!="NONE") fits_update_key_str(fptr, "DUNIT3", dunit3.c_str(), com, &status);
    fits_update_key_str(fptr, "BUNIT", bunit.c_str(), com, &status);
    
    if (bmaj!=0) fits_update_key_dbl(fptr, "BMAJ", bmaj, 10, com, &status);
    if (bmin!=0) fits_update_key_dbl(fptr, "BMIN", bmin, 10, com, &status);
    fits_update_key_dbl(fptr, "BPA", bpa, 10, com, &status);
    if (btype!="NONE") fits_update_key_str(fptr, "BTYPE", btype.c_str(), com, &status);
    //fits_update_key_flt(fptr, "BLANK", blank, 10, com, &status);
    
    if (object!="NONE") fits_update_key_str(fptr, "OBJECT", object.c_str(), com, &status);
    if (epoch!=0) fits_update_key_flt(fptr, "EPOCH", epoch, 10, com, &status);
    if (telescope!="NONE") fits_update_key_str(fptr, "TELESCOP", telescope.c_str(), com, &status);
    if (freq0!=0) fits_update_key_dbl(fptr, "FREQ0", freq0, 10, com, &status);
    if (datamax!=0) fits_update_key_dbl(fptr, "DATAMAX", datamax, 10, com, &status);
    if (datamin!=0) fits_update_key_dbl(fptr, "DATAMIN", datamin, 10, com, &status);

    if (fullHead) {
        for (uint i=0; i<keys.size(); i++) {
            status=0;
            bool towrite = false;
            int hist = keys[i].find("HISTORY");
            
            if(hist>=0) towrite=true;
            else {
                int found = keys[i].find("="); 
                if (found>=0) {
                    char keyname [] = "                                                                                  ";
                    strncpy(keyname, keys[i].c_str(),found);
                    char card[100];
                    if (fits_read_card(fptr, keyname, card, &status)) towrite=true;
                }
            }
            if (towrite) {
                status=0;
                fits_write_record(fptr, keys[i].c_str(), &status);
            }
        }
        
    }
    
    
    fits_report_error(stderr, status); 
}


void Header::headwrite_2d (fitsfile *fptr, bool fullHead) {

    int status=0;
    char com[]= "  ";
    
    fits_update_key_dbl(fptr, "CRPIX1", crpix[0], 10, com, &status);
    fits_update_key_dbl(fptr, "CRVAL1", crval[0], 10, com, &status);
    fits_update_key_dbl(fptr, "CDELT1", cdelt[0], 10, com, &status);
    fits_update_key_str(fptr, "CTYPE1", ctype[0].c_str(), com, &status);
    fits_update_key_str(fptr, "CUNIT1", cunit[0].c_str(), com, &status);
    
    fits_update_key_dbl(fptr, "CRPIX2", crpix[1], 10, com, &status);
    fits_update_key_dbl(fptr, "CRVAL2", crval[1], 10, com, &status);
    fits_update_key_dbl(fptr, "CDELT2", cdelt[1], 10, com, &status);
    fits_update_key_str(fptr, "CTYPE2", ctype[1].c_str(), com, &status);
    fits_update_key_str(fptr, "CUNIT2", cunit[1].c_str(), com, &status);
    
    fits_update_key_str(fptr, "BUNIT", bunit.c_str(), com, &status);
    
    if (bmaj!=0) fits_update_key_dbl(fptr, "BMAJ", bmaj, 10, com, &status);
    if (bmin!=0) fits_update_key_dbl(fptr, "BMIN", bmin, 10, com, &status);
    if (bpa!=0)  fits_update_key_dbl(fptr, "BPA", bpa, 10, com, &status);
    if (btype!="NONE") fits_update_key_str(fptr, "BTYPE", btype.c_str(), com, &status);
    if (epoch!=0) fits_update_key_flt(fptr, "EPOCH", epoch, 10, com, &status);
    fits_update_key_str(fptr, "OBJECT", object.c_str(), com, &status);
    //fits_update_key_flt(fptr, "BZERO", bzero, 12, com, &status);
    //fits_update_key_flt(fptr, "BSCALE", bscale, 12, com, &status);
    //fits_update_key_flt(fptr, "BLANK", blank, 12, com, &status);
    
    if (datamax!=0) fits_update_key_dbl(fptr, "DATAMAX", datamax, 10, com, &status);
    if (datamin!=0) fits_update_key_dbl(fptr, "DATAMIN", datamin, 10, com, &status);
    
    if (fullHead) {
        for (uint i=0; i<keys.size(); i++) {
            status=0;
            bool towrite = false;
            int hist = keys[i].find("HISTORY");

            if(hist>=0) towrite=true;
            else {
                int found = keys[i].find("=");
                if (found>=0) {
                    char keyname [] = "                                                                                  ";
                    strncpy(keyname, keys[i].c_str(),found);
                    char card[100];
                    if (fits_read_card(fptr, keyname, card, &status)) towrite=true;
                }
            }
            if (towrite) {
                status=0;
                fits_write_record(fptr, keys[i].c_str(), &status);
            }
        }

    }
    
    fits_report_error(stderr, status);  
}


template <class T>
bool Header::read_keyword(std::string keyword, T &key, bool err) {
    
    std::cout << "HEADER ERROR: unknown type of keyword "<< keyword << std::endl;
    return false;
}

template <>
bool Header::read_keyword<int>(std::string keyword, int &key, bool err) {
    
    fitsfile *fptr;                                 
    int status=0;
    char comment[72];   
 
    if (fits_open_file(&fptr, fitsname.c_str(), READONLY, &status)) {
        fits_report_error(stderr, status);
        return false;
    }   
    status = 0;
    long dum = long(key);
    if (fits_read_key_lng (fptr, keyword.c_str(), &dum, comment, &status)) {
        if (err) fits_report_error(stderr, status);
        return false;
    }
    key=dum;
    
    status=0;
    if (fits_close_file(fptr, &status))
        fits_report_error(stderr, status);
    
    return true;
    
}

template <>
bool Header::read_keyword<long>(std::string keyword, long &key, bool err) {
    
    fitsfile *fptr;                                 
    int status=0;
    char comment[72];   
 
    if (fits_open_file(&fptr, fitsname.c_str(), READONLY, &status)) {
        fits_report_error(stderr, status);
        return false;
    }   
    status = 0;
    if (fits_read_key_lng (fptr, keyword.c_str(), &key, comment, &status)) {
        if (err) fits_report_error(stderr, status);
        return false;
    }
    status = 0;
    if (fits_close_file(fptr, &status))
        fits_report_error(stderr, status);
    
    return true;
    
}

template <>
bool Header::read_keyword<float>(std::string keyword, float &key, bool err) {
    
    fitsfile *fptr;                                 
    int status=0;
    char comment[72];   
 
    if (fits_open_file(&fptr, fitsname.c_str(), READONLY, &status)) {
        fits_report_error(stderr, status);
        return false;
    }   
    status = 0;
    if (fits_read_key_flt (fptr, keyword.c_str(), &key, comment, &status)) {
        if (err) fits_report_error(stderr, status);
        return false;
    }
    status = 0;
    if (fits_close_file(fptr, &status))
        fits_report_error(stderr, status);
    
    return true;
    
}

template <>
bool Header::read_keyword<double>(std::string keyword, double &key, bool err) {
    
    fitsfile *fptr;                                 
    int status=0;
    char comment[72];   
 
    if (fits_open_file(&fptr, fitsname.c_str(), READONLY, &status)) {
        fits_report_error(stderr, status);
        return false;
    }   
    status = 0;
    if (fits_read_key_dbl (fptr, keyword.c_str(), &key, comment, &status)) {
        if (err) fits_report_error(stderr, status);
        return false;
    }
    status = 0;
    if (fits_close_file(fptr, &status))
        fits_report_error(stderr, status);
    
    return true;
    
}

template <>
bool Header::read_keyword<std::string>(std::string keyword, std::string &key, bool err) {
    
    fitsfile *fptr;                                 
    int status=0;
    char comment[72];
        
    if (fits_open_file(&fptr, fitsname.c_str(), READONLY, &status)) {
        fits_report_error(stderr, status);
        return false;
    }
    
    char *Key = new char[200];
    status = 0;
    if (fits_read_key_str (fptr, keyword.c_str(), Key, comment, &status)) {
        if (err) fits_report_error(stderr, status);
        return false;
    }
    
    key = Key;
    
    status = 0;
    if (fits_close_file(fptr, &status))
        fits_report_error(stderr, status);
    
    delete [] Key;
    return true;
}

template <>
bool Header::read_keyword<char*>(std::string keyword, char* &key, bool err) {
    
    fitsfile *fptr;                                 
    int status=0;
    char comment[72];
        
    if (fits_open_file(&fptr, fitsname.c_str(), READONLY, &status)) {
        fits_report_error(stderr, status);
        return false;
    }
    
    status = 0;
    if (fits_read_key_str (fptr, keyword.c_str(), key, comment, &status)) {
        if (err) fits_report_error(stderr, status);
        return false;
    }
    
    status = 0;
    if (fits_close_file(fptr, &status))
        fits_report_error(stderr, status);
    
    return true;
}
