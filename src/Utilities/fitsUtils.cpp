//---------------------------------------------------------------
// fitsUtils.cpp: Utility functions for FITS Files
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
#include <cstring>
#include <fitsio.h>
#include <Utilities/utils.hh>

template <> int selectBitpix<short>() {return SHORT_IMG;}
template <> int selectBitpix<int>() {return SHORT_IMG;}
template <> int selectBitpix<long>() {return LONG_IMG;}
template <> int selectBitpix<float>() {return FLOAT_IMG;}
template <> int selectBitpix<double>() {return DOUBLE_IMG;}

template <> int selectDatatype<short>() {return TSHORT;}
template <> int selectDatatype<int>() {return TINT;}
template <> int selectDatatype<long>() {return TLONG;}
template <> int selectDatatype<float>() {return TFLOAT;}
template <> int selectDatatype<double>() {return TDOUBLE;}


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


int modhead(int argc, char *argv[]) {
    
    // This function is modified from "modhead" FITS utility from NASA
    // https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html#modhead
    
    fitsfile *fptr;        
    char card[FLEN_CARD], newcard[FLEN_CARD];
    char oldvalue[FLEN_VALUE], comment[FLEN_COMMENT];
    int status = 0;
    int iomode, keytype;

    if (argc==4) iomode = READONLY;
    else if (argc==5) iomode = READWRITE;
    else {
        std::cout << "\n BBarolo's MODHEAD FITS utility: \n\n"
                  << " Write or modify the value of a header keyword.\n"
                  << " If a newvalue is not specified, just print the current value.\n\n"
                  << " Usage:\n   BBarolo --modhead filename[ext] keyword [newvalue]\n\n"
                  << " Examples: \n"
                  << "   BBarolo --modhead in.fits cunit1      (list the CUNIT1 keyword)\n"
                  << "   BBarolo --modhead in.fits cunit1 deg  (set CUNIT1 = 'deg')\n\n"
                  << " NOTE: it may be necessary to enclose the input file name in single \n"
                  << " quote characters on some Unix shells.\n\n";
        return 0;
    }

    if (!fits_open_file(&fptr, argv[2], iomode, &status)) {
        if (fits_read_card(fptr,argv[3],card, &status)) {
            std::cerr << "Keyword does not exist\n";
            card[0] = comment[0] = '\0';
            status = 0; 
        }
        else std::cout << card << std::endl;

        if (argc==5) {      // Write or overwrite the keyword 
            // Check if this is a protected keyword that must not be changed 
            if (*card && fits_get_keyclass(card) == TYP_STRUC_KEY) 
                std::cerr << "Protected keyword cannot be modified.\n";
            else {
                
                // Retrieve the comment string 
                if (*card) fits_parse_value(card, oldvalue, comment, &status);

                // Construct template for new keyword */
                std::string newcard = std::string(argv[3]) + " = ";
                if (isdigit(argv[4][1])) newcard += std::string(argv[4]);
                else newcard += "'" + std::string(argv[4]) + "'" ;
                
                if (*comment) newcard += "/ " + std::string(comment);
                
                // Reformat the keyword string to conform to FITS rules
                fits_parse_template(&newcard[0], card, &keytype, &status);

                // Overwrite the keyword with the new value */
                fits_update_card(fptr, argv[3], card, &status);

                std::cout << "Keyword has been changed to:\n";
                std::cout << card << std::endl;
            } 
        }
        fits_close_file(fptr, &status);
    } 
    
    // If error occured, print out error message 
    if (status) fits_report_error(stderr, status);
    return status ;
}


int remhead(int argc, char *argv[]) {
    
    // This function remove a keyword from a FITS file
    
    fitsfile *fptr;
    char card[FLEN_CARD];
    int status = 0;

    if (argc!=4) {
        std::cout << "\n BBarolo's REMHEAD FITS utility: \n\n"
                  << " Delete a header keyword from a FITS file.\n\n"
                  << " Usage:\n   BBarolo --remhead filename[ext] keyword \n\n"
                  << " Examples: \n"
                  << "   BBarolo --remhead in.fits object      (remove the OBJECT keyword)\n\n"
                  << " NOTE: it may be necessary to enclose the input file name in single \n"
                  << " quote characters on some Unix shells.\n\n";
        return 0;
    }

    if (!fits_open_file(&fptr, argv[2], READWRITE, &status)) {
        if (fits_read_card(fptr,argv[3],card, &status)) {
            std::cerr << "Keyword " << makeupper(std::string(argv[3])) << " does not exist.\n";
            status = 0;
        }
        else {
            // Check if this is a protected keyword that can not be deleted 
            if (*card && fits_get_keyclass(card) == TYP_STRUC_KEY) 
                std::cerr << "Protected keyword " << makeupper(std::string(argv[3])) << " cannot be deleted.\n";
            else {
                // Delete keyword
                if (!fits_delete_key(fptr, argv[3], &status))
                    std::cout << "Keyword " << makeupper(std::string(argv[3])) << " has been deleted.\n";
            }
        }
        fits_close_file(fptr, &status);
    } 
    
    // If error occured, print out error message 
    if (status) fits_report_error(stderr, status);
    return status ;
}


int listhead(int argc, char *argv[]) {
    
    // This function is modified from "listhead" FITS utility from NASA
    // https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html#listhead
    
    fitsfile *fptr;
    char card[FLEN_CARD];
    int status = 0;
    int single = 0, hdupos, nkeys;

    if (argc != 3) {
        std::cout << "\n BBarolo's LISTHEAD FITS utility: \n\n"
                  << " List the FITS header keywords in a single extension, or, if [ext] is\n"
                  << " not given, list the keywords in all the extensions.\n\n"
                  << " Usage:\n   BBarolo --listhead filename[ext] \n\n"
                  << " Examples: \n"
                  << "   BBarolo --listhead file.fits      (list every header in the file) \n"
                  << "   BBarolo --listhead file.fits[0]   (list primary array header) \n"
                  << "   BBarolo --listhead file.fits[2]   (list header of 2nd extension) \n"
                  << "   BBarolo --listhead file.fits[GTI] (list header of GTI extension) \n\n"
                  << " NOTE: it may be necessary to enclose the input file name in single \n"
                  << " quote characters on some Unix shells.\n\n";
        return 0;
    }

    if (!fits_open_file(&fptr, argv[2], READONLY, &status)) {
        
        // Get the current HDU position
        fits_get_hdu_num(fptr, &hdupos);

        // List only a single header if a specific extension was given 
        if (hdupos != 1 || strchr(argv[2], '[')) single = 1;

        // Main loop through each extension
        for (; !status; hdupos++) { 
 
            std::cout << "Header listing for HDU #" << hdupos << "\n";
        
            // Get # of keywords
            fits_get_hdrspace(fptr, &nkeys, NULL, &status);

            // Read and print each keywords 
            for (int ii = 1; ii<=nkeys; ii++) {  
                if (fits_read_record(fptr, ii, card, &status)) break;
                std::cout << card << std::endl;
            }
            std::cout << "END\n\n";

            if (single) break;

            // Try to move to next HDU
            fits_movrel_hdu(fptr, 1, NULL, &status);  
        }

        // Reset after normal error
        if (status == END_OF_FILE) status = 0; 

        fits_close_file(fptr, &status);
    }

    if (status) fits_report_error(stderr, status);
    return status;
}


int fitscopy(int argc, char *argv[]) {
    
    // This function is modified from "fitscopy" FITS utility from NASA
    // https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html#fitscopy
    
    fitsfile *infptr, *outfptr;
    int status = 0;

    if (argc != 4) {
        std::cout << "\n BBarolo's FITSCOPY FITS utility: \n\n"
                  << " Copy an input file to an output file, optionally filtering the file in the \n"
                  << " process. Filters may be used to extract a subimage from a larger image, \n"
                  << " select rows from a table, filter a table with a GTI time extension or a SAO\n"
                  << " region file, create or delete columns in a table, create an image by binning\n"
                  << " 2 table columns, and convert IRAF format .imh or raw binary data files into \n"
                  << " FITS images (see CFITSIO User's Guide for filtering syntax).\n\n"
                  << " Usage:\n   BBarolo --fitscopy inputfile[filter] outputfile \n\n"
                  << " Examples: \n"
                  << "   BBarolo --fitscopy in.fits[11:50,21:60] out.fits      (copy a subimage)\n"
                  << "   BBarolo --fitscopy in.fits[-*,*] out.fits             (mirror reverse axis 0)\n"
                  << "   BBarolo --fitscopy iniraf.imh out.fits                (IRAF image to FITS)\n"
                  << "   BBarolo --fitscopy in.dat[i512,512] out.fit           (binary file to FITS)\n"
                  << "   BBarolo --fitscopy in.fits[events][pi>35] out.fits    (copy rows with pi>35)\n"
                  << "   BBarolo --fitscopy in.fits[events][bin X,Y] out.fits  (bin an image)\n\n"
                  << " NOTE: it may be necessary to enclose the input file name in single quote\n"
                  << " characters on some Unix shells.\n\n";
        return 0;
    }
    
    // Open the input file 
    if (!fits_open_file(&infptr, argv[2], READONLY, &status)) {
        // Create the output file 
        remove(argv[3]);
        if (!fits_create_file(&outfptr, argv[3], &status)) {
            // Copy the previous, current, and following HDUs
            fits_copy_file(infptr, outfptr, 1, 1, 1, &status);
            fits_close_file(outfptr,  &status);
        }
        fits_close_file(infptr, &status);
    }

    if (status) fits_report_error(stderr, status);
    return status;
}


int fitsarith(int argc, char *argv[]) {
    
    // This function is modified from "imarith" FITS utility from NASA
    // https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html#imarith
 
    fitsfile *afptr, *bfptr, *outfptr;
    int status = 0;
    int atype, btype, anaxis, bnaxis, check=1, op;
    long ntodo;
    long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1};
    double value = 0.0;
    int image2=1;
    
    if (argc != 6) {
        std::cout << "\n BBarolo's FITSARITH FITS utility: \n\n"
                  << " Perform an operation between two images or between an image and a number. \n"
                  << " Supported operators are add, sub, mul, div (first character required). \n\n"
                  << " Usage:\n   BBarolo --fitsarith image1 { image2 | value } oper outimage \n\n"
                  << " Examples: \n"
                  << "   BBarolo --fitsarith in1.fits in2.fits add out.fits  (add the 2 files)\n"
                  << "   BBarolo --fitsarith in1.fits 1000.0 mul out.fits    (mult in1 by 1000)\n\n"
                  << " NOTE: it may be necessary to enclose the file names in single quote\n"
                  << " characters on some Unix shells.\n\n";

        return 0;
    }

    // Open input images
    if (fits_open_file(&afptr, argv[2], READONLY, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    if (fits_open_file(&bfptr, argv[3], READONLY, &status)) {
        value = atof(argv[3]);
        if (value == 0.0) {
            std::cerr << "FITSARITH ERROR: second argument is neither an image name"
                      << "nor a valid numerical value.\n";
            return status;
        }
        image2 = 0;
        status = 0;
    }

    // Read dimensions
    fits_get_img_dim(afptr, &anaxis, &status); 
    if (image2) fits_get_img_dim(bfptr, &bnaxis, &status);
    fits_get_img_size(afptr, 3, anaxes, &status);
    if (image2) fits_get_img_size(bfptr, 3, bnaxes, &status);

    if (status) {
        fits_report_error(stderr, status);
        return status;
    }
    
    // Check that images have same sizes and no more than 3 axes
    if (anaxis > 3) {
        std::cerr << "FITSARITH ERROR: images with > 3 dimensions are not supported.\n";
        check = 0;
    }
    else if (image2 && (anaxes[0]!=bnaxes[0] || anaxes[1]!=bnaxes[1] || anaxes[2]!=bnaxes[2])) {
        std::cerr << "FITSARITH ERROR: input images don't have same size.\n";
        check = 0;
    }
    
    // Check that requested operation is acceptable
    if      (*argv[4]=='a' || *argv[4]=='A') op = 1;
    else if (*argv[4]=='s' || *argv[4]=='S') op = 2;
    else if (*argv[4]=='m' || *argv[4]=='M') op = 3;
    else if (*argv[4]=='d' || *argv[4]=='D') op = 4;
    else {
        std::cerr << "FITSARITH ERROR: unknown arithmetic operator " << argv[4] << "\n";
        check = 0;
    }

    // Create the new empty output file if the above checks are OK
    remove(argv[5]);
    if (check && !fits_create_file(&outfptr, argv[5], &status)) {
      // Copy all the header keywords from first image to new output file
      fits_copy_header(afptr, outfptr, &status);
      // Number of pixels to read in each row */
      long npixels = anaxes[0];  
      
      double *apix = new double[npixels];
      double *bpix = nullptr;
      if (image2) bpix = new double[npixels]; 

      if (apix==nullptr || (image2 && bpix==nullptr)) {
        std::cerr << "FITSARITH ERROR: Memory allocation error\n";
        return 1;
      }
      
      long firstpix[3] = {1,1,1};
      
      // Loop over all planes of the cube
      for (firstpix[2] = 1; firstpix[2]<=anaxes[2]; firstpix[2]++) {
          // Loop over all rows of the plane 
          for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
              // Read both images as doubles, regardless of actual datatype.
              // Give starting pixel coordinate and number of pixels to read.
              if (fits_read_pix(afptr,TDOUBLE,firstpix,npixels,nullptr,apix,nullptr,&status)) break;
              if (image2 && fits_read_pix(bfptr,TDOUBLE,firstpix,npixels,nullptr,bpix,nullptr,&status)) break; 

              switch (op) {
              
                  case 1:                               // Addition
                  for(int i=0; i<npixels; i++) {
                      if (image2) apix[i] += bpix[i];
                      else apix[i] += value;
                  }
                  break;
                  
                  case 2:                               // Subtraction
                  for(int i=0; i<npixels; i++) {
                      if (image2) apix[i] -= bpix[i];
                      else apix[i] -= value;
                  }
                  break;
                  
                  case 3:                               // Multiplication
                  for(int i=0; i<npixels; i++) {
                      if (image2) apix[i] *= bpix[i];
                      else apix[i] *= value;
                  }
                  break;
                  
                  case 4:                               // Division
                  for(int i=0; i<npixels; i++) {
                      if (image2) {
                          if (bpix[i]!=0.) apix[i] /= bpix[i];
                          else apix[i]=0;
                      }
                      else apix[i] /= value;
                  }
                  break;

              }
              
              // Write new values to output image */
              fits_write_pix(outfptr, TDOUBLE, firstpix, npixels, apix, &status); 
          }
      }

      fits_close_file(outfptr, &status);
      delete [] apix;
      if (image2) delete [] bpix;
    }

    fits_close_file(afptr, &status);
    if (image2) fits_close_file(bfptr, &status);
 
    if (status) fits_report_error(stderr, status);
    return status;
}
