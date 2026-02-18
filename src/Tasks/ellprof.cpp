//--------------------------------------------------------------------
// ellprof.cpp: Members functions of the Ellprof class.
//--------------------------------------------------------------------

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
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <Arrays/cube.hh>
#include <Arrays/param.hh>
#include <Tasks/ellprof.hh>
#include <Tasks/galmod.hh>
#include <Utilities/utils.hh>
#include <Utilities/progressbar.hh>


namespace Tasks {

template <class T>
void Ellprof<T>::defaults() {
    Overlap = true;
    Nseg = 1;
    maprotation = 0;
    Range[0] = -1.0*FLT_MAX;
    Range[1] = FLT_MAX;
    Rmax=0;
    subpix[0]=subpix[1]=2;
}


template <class T>
void Ellprof<T>::allocateArrays (size_t nrad, size_t nseg) {

    deallocateArrays();
    
    Radius = new T [nrad];
    Width = new T [nrad];
    Phi = new T [nrad];
    Inc = new T [nrad];
    Cosphi = new double [nrad];
    Sinphi = new double [nrad];
    Cosinc = new double [nrad];
    Annuli = new T*[nrad];
    Segments = new float [2*nseg];

    Sum = new double*[nrad];
    Sumsqr = new double*[nrad];
    Num = new long*[nrad];
    Numblanks = new long*[nrad];

    Datamin = new double*[nrad];
    Datamax = new double*[nrad];
    Var = new double*[nrad];
    MAD = new double*[nrad];
    Contrib = new long*[nrad];

    Mean = new double*[nrad];
    Median = new double*[nrad];
    Area = new double*[nrad];
    Blankarea = new double*[nrad];
    Surfdens = new double*[nrad];
    Surfdens_Bl = new double*[nrad];

    medianArray = new std::vector<double>*[nrad];

    for (size_t i=0; i<nrad; i++) {
        Annuli[i] = new T[2];
        Sum[i] = new double[nseg];
        Sumsqr[i] = new double[nseg];
        Num[i] = new long[nseg];
        Numblanks[i] = new long[nseg];
        Datamin[i] = new double[nseg];
        Datamax[i] = new double[nseg];
        Var[i] = new double[nseg];
        MAD[i] = new double[nseg];
        Contrib[i] = new long[nseg];
        Mean[i] = new double[nseg];
        Median[i] = new double[nseg];
        Area[i] = new double[nseg];
        Blankarea[i] = new double[nseg];
        Surfdens[i] = new double[nseg];
        Surfdens_Bl[i] = new double[nseg];
        medianArray[i] = new std::vector<double>[nseg];
    }

    for (size_t i=0; i<nrad; i++) {
        for (size_t j=0; j<nseg; j++) {
            Annuli[i][j]=Sum[i][j]=Sumsqr[i][j]=Num[i][j]=Numblanks[i][j]=0;
            Datamin[i][j]=Datamax[i][j]=Var[i][j]=MAD[i][j]=Contrib[i][j]=Mean[i][j]=0;
            Median[i][j]=Area[i][j]=Blankarea[i][j]=Surfdens[i][j]=Surfdens_Bl[i][j]=0;
        }
    }

    arrayAllocated = true;
}


template <class T>
void Ellprof<T>::deallocateArrays () {

    if (arrayAllocated) {
        for (size_t i=0; i<Nrad; i++) {
            delete [] Annuli[i];
            delete [] Sum[i];
            delete [] Sumsqr[i];
            delete [] Num[i];
            delete [] Numblanks[i];
            delete [] Datamin[i];
            delete [] Datamax[i];
            delete [] Var[i];
            delete [] MAD[i];
            delete [] Contrib[i];
            delete [] Mean[i];
            delete [] Median[i];
            delete [] Area[i];
            delete [] Blankarea[i];
            delete [] Surfdens[i];
            delete [] Surfdens_Bl[i];
            delete [] medianArray[i];
        }

        delete [] Radius;
        delete [] Width;
        delete [] Phi;
        delete [] Inc;
        delete [] Cosphi;
        delete [] Sinphi;
        delete [] Cosinc;
        delete [] Annuli;
        delete [] Segments;

        delete [] Sum;
        delete [] Sumsqr;
        delete [] Num;
        delete [] Numblanks;

        delete [] Datamin;
        delete [] Datamax;
        delete [] Var;
        delete [] MAD;
        delete [] Contrib;

        delete [] Mean;
        delete [] Median;
        delete [] Area;
        delete [] Blankarea;
        delete [] Surfdens;
        delete [] Surfdens_Bl;

        delete [] medianArray;
    }
}


template <class T>
Ellprof<T>::Ellprof(Cube<T> *c, Rings<T> *inR, bool mask) {

    // Reading input rings from parameter file
    bool deallocate = false;

    if (inR==nullptr) {
        inR = readRings<T>(c->pars().getParGF(), c->Head());
        deallocate = true;
    }
    setFromCube(c,inR,mask);
    if (deallocate) delete inR;
}


template <class T>
Ellprof<T>::Ellprof(MomentMap<T> *image, size_t nrad, float width, float phi, float inc, float *pos, size_t nseg, float* segments) {

    Rings<T> *r = new Rings<T>;

    r->nr = nrad;
    r->radsep = width;
    for (size_t i=0; i<nrad; i++) {
        r->radii.push_back(i*width + width/2.);
        r->phi.push_back(phi);
        r->inc.push_back(inc);
        r->xpos.push_back(pos[0]);
        r->ypos.push_back(pos[1]);
    }

    init(image, r, nseg, segments);

    delete r;
}


template <class T>
Ellprof<T>::Ellprof(MomentMap<T> *image, Rings<T> *rings, size_t nseg, float* segments) {
    init(image,rings,nseg,segments);
}


template <class T>
void Ellprof<T>::setFromCube(Cube<T> *c, Rings<T> *inR, bool mask) {
    
    // Setting other options
    T meanPA = findMean(&inR->phi[0], inR->nr);
    int nseg = 1;
    
    float segments[4] = {0, 360., 0., 0};
    if (c->pars().getParGF().SIDE=="A") {
        nseg = 2;
        segments[2]=-90;
        segments[3]=90;
    }
    else if(c->pars().getParGF().SIDE=="R") {
        nseg = 2;
        segments[2]=90;
        segments[3]=-90;
    }
    if (meanPA>180) std::swap(segments[2], segments[3]);
    
    // Extracting moment map
    if (imAllocated) delete im;
    im = new MomentMap<T>;
    im->input(c);
    
    if (c->Head().NumAx()>2 && c->DimZ()>1) im->ZeroMoment(mask);
    else {
        for (int i=0; i<c->NumPix(); i++) im->Array(i) = c->Array(i);
        im->setHead(0);
        im->Head().setBunit(c->Head().Bunit());
    }
    
    imAllocated = true;

    init(im, inR, nseg, segments);
    
}


template <class T>
void Ellprof<T>::init(MomentMap<T> *image, Rings<T> *rings, size_t nseg, float* segments) {

    defaults();

    if (!image->HeadDef()) {
        std::cerr << "\n ELLPROF ERROR: Moment map has no proper header. Exiting ...\n";
        std::terminate();
    }

    im = image;
    
    update_rings(rings, nseg, segments);

}


template <class T>
void Ellprof<T>::update_rings(Rings<T> *rings, size_t nseg, float* segments) {

    Nseg = nseg;
    Nrad = rings->nr;

    allocateArrays(Nrad, Nseg);

    if (Nseg == 1 && segments==nullptr) 
        float segments[4] = {0., 360., 0., 0.};
    

    // Initialize the rings
    for (size_t r=0; r<Nrad; r++) {
        Radius[r] = rings->radii[r];
        Phi[r] = toangle(rings->phi[r]);
        Inc[r] = rings->inc[r];
        
        // Handling the case of a single ring
        if (rings->nr==1) Width[r] = rings->radsep>0 ? rings->radsep : rings->radii[0]/2.;
        else {
            if (r==rings->nr-1) Width[r] = (rings->radii[r]-rings->radii[r-1]);
            else Width[r] = (rings->radii[r+1]-rings->radii[r]);
        }
        
        /* Setup angle arrays */
        /*---------------------------------------------------------------*/
        /* The values are needed to rotate back so invert sign of angle  */
        /* The angle is Phi, corrected for the map PA and wrt. hor. axis */
        /*---------------------------------------------------------------*/
        Cosphi[r] = cos((-1.0*(Phi[r]+maprotation+90.0))*M_PI/180.);
        Sinphi[r] = sin((-1.0*(Phi[r]+maprotation+90.0))*M_PI/180.);

        Cosinc[r] = cos(Inc[r]*M_PI/180.);
        if (fabs(Cosinc[r]) <= 0.0001) {
            //std::cerr << "ELLPROF ERROR: " << r+1 << "th inclination (" << Inc[r] << ") is illegal! (COS(inc) = 0). Setting to 89...\n";
            Inc[r] = 89.;
            Cosinc[r] = cos(Inc[r]*M_PI/180.);
         }

        Annuli[r][0] = max(Radius[r] - 0.5*Width[r], 0.0);
        Annuli[r][1] = max(Radius[r] + 0.5*Width[r], 0.0);
        
        Rmax = max(Rmax, Annuli[r][1]);
    }

    // INITIALIZING THE SEGMENTS
    /* First segment is always the complete ring */
    Segments[0] =   0.0;
    Segments[1] = 360.0;

    if (segments!=NULL) {
        size_t istart = 0;
        // If the first is not [0,360], correct it!
        if (!(segments[0] == 0.0) && (segments[1] == 360.0)) {
            Nseg += 1;
            delete [] Segments;
            Segments = new float[Nseg];
            Segments[0] =   0.0;
            Segments[1] = 360.0;
            istart = 1;
        }
        for (size_t i=istart; i<Nseg; i++) {
            /* Convert to proper range */
            Segments[2*i]   = toangle(segments[2*i]);
            Segments[2*i+1] = toangle(segments[2*i+1]);
        }
    }

    // Checking the position of the center

    Position[0] = findMean(&rings->xpos[0], rings->nr);
    Position[1] = findMean(&rings->ypos[0], rings->nr);
    if (Position[0]>=im->DimX() || Position[0]<0 ||
        Position[1]>=im->DimY() || Position[1]<0) {
            std::cerr << "ELLPROF ERROR: The center of galaxy must be inside the map. Exiting...\n ";
            abort();
    }

    //Setting spacings in arcsec
    Dx = im->Head().Cdelt(0)*arcsconv(im->Head().Cunit(0));
    Dy = im->Head().Cdelt(1)*arcsconv(im->Head().Cunit(1));
    if (Dx==0 || Dy==0) {
        std::cerr << "ELLPROF ERROR: Null spacing. Check CUNIT and CDELT keywords in the header. Exiting...\n ";
        abort();
    }

    stepxy[0] = fabs(Dx) / (float) subpix[0];
    stepxy[1] = fabs(Dy) / (float) subpix[1];
}


template <class T>
void Ellprof<T>::setOptions (bool overlap, float *range, float *subp) {
    Overlap = overlap;
    Range[0] = range[0];
    Range[1] = range[1];
    subpix[0] = subp[0];
    subpix[1] = subp[1];
}


template <class T>
void Ellprof<T>::RadialProfile () {
   
    for (size_t i=0; i<Nrad; i++) {
        for (size_t s=0; s< Nseg; s++) {
            // Reset the statistics variables 
            Sum[i][s]=Mean[i][s]=Median[i][s]=0.0;
            Var[i][s]=MAD[i][s]=Area[i][s]= 0.0;
            Num[i][s]=Numblanks[i][s]=Contrib[i][s]= 0;
            Datamin[i][s]=FLT_MAX;
            Datamax[i][s]=-FLT_MAX;
            Surfdens[i][s]=Surfdens_Bl[i][s]=0.;
        }
    }

    float Bgridlo[2], Bgridhi[2];
    Bgridlo[0] = std::max(0, int(Position[0]-Rmax/fabs(Dx) - 1));
    Bgridhi[0] = std::min(im->DimX()-1, int(Position[0]+Rmax/fabs(Dx)+1));
    Bgridlo[1] = std::max(0, int(Position[1]-Rmax/fabs(Dy)-1));
    Bgridhi[1] = std::min(im->DimY()-1, int(Position[1]+Rmax/fabs(Dy)+1));

    // Start looping over all pixels 
    
    ProgressBar bar(true,im->pars().isVerbose(),im->pars().getShowbar());
    bar.init(" Computing radial profile... ",Bgridhi[0]-Bgridlo[0]);
    
    for (size_t x = Bgridlo[0]; x <= Bgridhi[0]; x++) {
        bar.update(x-Bgridlo[0]+1);
        for (size_t y = Bgridlo[1]; y <= Bgridhi[1]; y++) {
            processpixel(x, y, im->Array(x,y));
        }
    }
    
    bar.fillSpace("Done.\n");

    float Sumtotgeo    = 0.0;        // Sum of all complete rings 
    float Sumtotgeo_bl = 0.0;
    for (size_t n=0; n<Nrad; n++) {
        // We need the sum of of the contributions for each ring. Each contribution is 
        // multiplied by the area between two radii. The contributions are the corrected means.      
        if (Num[n][0]>0) { 
            float mean   = Sum[n][0] / Num[n][0];
            float meanbl = Sum[n][0] / (Num[n][0]+Numblanks[n][0]);
            float face_on_av_surfdens = mean * Cosinc[n]/(fabs(Dx*Dy));
            float face_on_av_surfdens_bl = meanbl * Cosinc[n]/(fabs(Dx*Dy));
            float geometricalarea = M_PI * (Annuli[n][1]*Annuli[n][1] - Annuli[n][0]*Annuli[n][0]);
            Sumtotgeo    += face_on_av_surfdens * geometricalarea;
            Sumtotgeo_bl += face_on_av_surfdens_bl * geometricalarea;
        }
    }

    // For all rings and segments the sum is calculated. Do some simple statistics 
    // using this sum and the number of pixels involved. 

    float subpixtot = subpix[0]*subpix[1];
    for (size_t i=0; i<Nrad; i++) {
        float   face_on_surfdens_bl_tot = 0.0;
        float   face_on_surfdens_tot = 0.0;
        for (size_t m=0; m<Nseg; m++) {
            float   surfdens = 0.0;
            float   face_on_surfdens = 0.0;
            float   surfdens_bl = 0.0;
            float   face_on_surfdens_bl = 0.0;
            float   area, area_bl;
            Sum[i][m] /= subpixtot;

            if (Num[i][m]==0) {
                Mean[i][m] = 0;
                Median[i][m] = 0;
                Var[i][m]  = 0;
                MAD[i][m] = 0;
                Datamin[i][m] = 0;
                Datamax[i][m] = 0;
                Area[i][m] = 0;
            }
            else {
                // The 'processpixel' function calculated too much FLUX. Each intensity has to    
                // be divided by the number of 'subpixels' in a pixel. In order to get the AREA   
                // expressed in pixels, it will be divided by the same number.                    
                Area[i][m] = Num[i][m] / subpixtot;
                Mean[i][m] = Sum[i][m] / Area[i][m];
                if (Num[i][m]>1) Var[i][m] = (Sumsqr[i][m] - Num[i][m]*Mean[i][m]*Mean[i][m])/float(Num[i][m]-1);
            }

            Blankarea[i][m] = Numblanks[i][m] / subpixtot;
            Median[i][m]  = findMedian(&medianArray[i][m][0], Num[i][m]);
            MAD[i][m]  = findMADFM(&medianArray[i][m][0], Num[i][m],Median[i][m],false);
            
            area = Area[i][m] * fabs(Dx*Dy);
            if (area == 0.0)  surfdens = 0.0;
            else surfdens = Sum[i][m] / area;

            face_on_surfdens = Cosinc[i] * surfdens;
            area_bl = Blankarea[i][m] * fabs(Dx*Dy);

            if (area + area_bl == 0.0) surfdens_bl = 0.0;
            else surfdens_bl = Sum[i][m] / (area+area_bl);
            face_on_surfdens_bl = Cosinc[i] * surfdens_bl;

            face_on_surfdens_tot += face_on_surfdens;
            face_on_surfdens_bl_tot += face_on_surfdens_bl;

            Surfdens[i][m] = surfdens;
            Surfdens_Bl[i][m] = surfdens_bl;
        }

    }
    
    
}


template <class T>
void Ellprof<T>::processpixel(int x, int y,float imval) {
/*------------------------------------------------------------*/
/* PURPOSE: Given the central position of a pixel, generate   */
/*          positions in that pixel and check whether the are */
/*          inside or outside a ring/segment.                 */
/*                                                            */
/* Example of subdivision of a pixel in y direction.          */
/*                                                            */
/*       -                                                    */
/*  |    |  +   <- Last subpixel in y                         */
/*  A    -                                                    */
/*  b    |  +                                                 */
/*  s    -         <== position of Yr                         */
/*  d    |  +                                                 */
/*  y    -                                                    */
/*  |    |  +   <- First subpixel in y                        */
/*       -                                                    */
/*                                                            */
/* Note that if the subdivision is one pixel, the position    */
/* that will be examined is the central position (Xr, Yr) of  */
/* that pixel.                                                */
/*------------------------------------------------------------*/

    bool validpixel = IsInRange(imval, Range);

    bool inside[Nrad];

    bool  Contribflag[Nrad][Nseg];  /* Did a pixel already contribute to ring? */
    for (size_t rad=0; rad<Nrad; rad++)
        for (size_t seg=0; seg<Nseg; seg++)
            Contribflag[rad][seg] = false;


    /* The pixel position converted to arcsec wrt central position */
    double absdx = fabs(Dx);
    double absdy = fabs(Dy);
    float Xr = absdx*(x-Position[0]);
    float Yr = absdy*(y-Position[1]);

    for (float posX = Xr + 0.5*(stepxy[0]-absdx); posX < Xr + 0.5*absdx; posX += stepxy[0]) {
        for (float posY = Yr + 0.5*(stepxy[1]-absdy); posY < Yr + 0.5*absdy; posY += stepxy[1]) {

            int     overcount = 0;
            float   imvaloverlap;

            for (size_t rad=0; rad<Nrad; rad++) {
                if (IsInRing(posX, posY, rad)) {
                    inside[rad] = true;
                    overcount++;
                }
                else inside[rad] = false;
            }

            imvaloverlap = imval;
            if (Overlap)
                if (overcount != 0 && imval != imval) imvaloverlap = imval/float(overcount);

            for (size_t rad=0; rad<Nrad; rad++){
                if (inside[rad]) {
                    float  theta = gettheta(posX, posY, Phi[rad], maprotation);
                    for (size_t seg=0; seg < Nseg; seg++) {
                        if (IsInSegment(theta, Segments[2*seg], Segments[2*seg+1]) ) {
                            /* Now we have a pixel (in a ring, in a segment) that is either blank or not   */
                            /* blank. If it is not a blank, but its image value is not within the wanted   */
                            /* range of values, it will be treated as a blank.                             */

                            if (validpixel) {
                                Sum[rad][seg] += imvaloverlap;
                                Sumsqr[rad][seg] += imvaloverlap*imvaloverlap;
                                Num[rad][seg] += 1;
                                Contribflag[rad][seg] = true;

                                /* Overlapping or not, the data min, max are the real data min and max.  */
                                /* not of  the weighted values.                                          */
                                if (imval>Datamax[rad][seg]) Datamax[rad][seg]=imval;
                                if (imval<Datamin[rad][seg]) Datamin[rad][seg]=imval;


                                medianArray[rad][seg].push_back(imvaloverlap);

                            }
                            else Numblanks[rad][seg]+=1;
                        }
                    }
                }
            }
        }
    }

    for (size_t rad=0; rad<Nrad; rad++)
        for (size_t seg=0; seg<Nseg; seg++)
            if (Contribflag[rad][seg]) Contrib[rad][seg] += 1;

}


template <class T>
bool Ellprof<T>::IsInRange(float value, float *Range) {
    // Check whether value of pixel is within user given range (in flux)

    if (value!=value) return false;
    if (Range[0] < Range[1]) {
        if (value >= Range[0] && value <= Range[1]) return true;
    }
    else {
        if (value < Range[1] || value > Range[0]) return true;
    }
    return false;
}


template <class T>
bool Ellprof<T>::IsInRing(float Xr, float Yr, int radnr) {
    // Return true if position Xr Yr is inside a ring

    float Xx =  Xr * Cosphi[radnr] - Yr * Sinphi[radnr];
    float Yy = (Xr * Sinphi[radnr] + Yr * Cosphi[radnr]) / Cosinc[radnr];
    float R_sqr = Xx*Xx + Yy*Yy;
    return( (Annuli[radnr][0]*Annuli[radnr][0] <= R_sqr && R_sqr < Annuli[radnr][1]*Annuli[radnr][1]) );
}


template <class T>
float Ellprof<T>::gettheta(float X,float Y,float Phi,float Crota) {
    /*------------------------------------------------------------*/
    /* Convert angle in XY plane to angle wrt major axis of map.  */
    /* 'atan2': Returns in radians the arc tangent of two real    */
    /* numbers. The arguments must not both be 0.0. If number-2   */
    /* is 0.0, the absolute value of the result is pi/2. If       */
    /* number-1 is 0.0, the result is 0.0 if number-2 is positive */
    /* and pi if number-2 is negative. Otherwise, the result is   */
    /* in the range -pi, exclusive, through +pi, inclusive, and   */
    /* is calculated as follows: arctan (argument_1/argument_2)   */
    /* If number-1 is positive, the result is positive; if        */
    /* number-1 is negative, the result is negative.              */
    /*------------------------------------------------------------*/

    float theta;

    if ( X == 0.0 && Y == 0.0 ) theta = 0.0;
    else
        theta = atan2( Y, X )*180./M_PI;         /* Angle in deg. wrt pos x axis */

    /* Now correct for angle offsets: */
    /* 360 = (90-theta) + Crota + phi */
    theta = toangle(theta + 270.0  - Phi - Crota);
    return theta;
}


template <class T>
float Ellprof<T>::toangle(float Angle) {
    /*------------------------------------------------------------*/
    /* PURPOSE: Return angle between 0 and < 360.0                */
    /*------------------------------------------------------------*/
    while (Angle < 0.0)
        Angle +=360.0;
    while (Angle > 360.0)
        Angle -=360.0;
    return Angle;
}


template <class T>
bool Ellprof<T>::IsInSegment(float Angle, float Segm1, float Segm2) {

    /* Check whether Angle of pixel is within given segment.*/
    /* Return true if range >= 360 degrees.                          */

    if ( fabs(Segm2 - Segm1) >= 360.0 ) return true;
    if (Segm1 < Segm2) {
        if (Angle >= Segm1 && Angle < Segm2) return true;
        else return false;
    }
    if (Segm1 > Segm2) {
        if (Angle >= Segm1 || Angle < Segm2) return true;
        else return false;
    }

    if (Segm1 == Segm2) {
        if (Angle == Segm1) return true;
        else return false;
    }

    return false;
}


template <class T>
void Ellprof<T>::printProfile (ostream& theStream, int seg) {
    
    std::string unit = deblankAll(im->Head().Bunit());
    std::string unit_l = makelower(unit);
    
    // Checking if units are in JY
    bool isJy = false;
    if (unit_l.find("jy")!=std::string::npos) isJy = true;
    
    theStream << "# ELLPROF results for " << im->Head().Name() << std::endl;
    theStream << fixed << setprecision(5);
    theStream << "#\n# Map units: u = " << unit << std::endl;
    theStream << "# Pixel area: " << fabs(Dx*Dy) << " arcs2" << std::endl << "#\n";
    theStream << "# Columns 2-6  : ring statistics (sum, mean, median, standard deviation and median absolute deviation from the median) in map units (see above).\n";
    theStream << "# Column  7    : effective number of valid pixels within the ring.\n";
    theStream << "# Columns 8-10 : surface density, its error and face-on surface density (inclination corrected).\n";
    theStream << "# Columns 11   : effective number of blank pixels within the ring.\n";
    theStream << "# Columns 12-13: surface density and face-on surface density including area covered by blank pixels (total area: NPIX + NPIXBL).\n";
    
    if (isJy)
//        theStream << "# Columns 14-16: (ONLY FOR HI DATA!) face-on (inclination corrected) number surface density and mass surface densities with two different techniques.\n";
        theStream << "# Columns 14-15: (ONLY FOR HI DATA!) face-on (inclination corrected) number surface density and mass surface density.\n";
    
    int m=9;
    theStream << "#\n#" << setw(m-1) << "RADIUS" << " "
              << setw(m) << "SUM" << " "
              << setw(m) << "MEAN" << " "
              << setw(m) << "MEDIAN" << " "
              << setw(m) << "STDDEV" << " "
              << setw(m) << "MAD" << " "
              << setw(m) << "NPIX" << " "
              << setw(m) << "SDENS" << " "
              << setw(m) << "ERR_SD" << " "
              << setw(m) << "SDENS_FO" << " "
              << setw(m) << "NPIXBL" << " "
              << setw(m) << "SDENSBL" << " "
              << setw(m+1) << "SDENSBL_FO" << " ";
    if (isJy)
        theStream << setw(m+1) << "HIDENS_FO" << " "
                  << setw(m+1) << "HIMASSDENS" << " ";
//                  << setw(m) << "HIMASSDENS2" << " ";

    theStream << "\n#" << setw(m-1) << "arcsec" << " "
              << setw(m) << "u" << " "
              << setw(m) << "u" << " "
              << setw(m) << "u" << " "
              << setw(m) << "u" << " "
              << setw(m) << "u" << " "
              << setw(m) << "pix" << " "
              << setw(m) << "u/arcs2" << " "
              << setw(m) << "u/arcs2" << " "
              << setw(m) << "u/arcs2" << " "
              << setw(m) << "pix" << " "
              << setw(m) << "u/arcs2" << " "
              << setw(m+1) << "u/arcs2" << " ";
    if (isJy)
        theStream << setw(m+1) << "1E20/cm2" << " "
                  << setw(m+1) << "Msun/pc2" << " ";
//                  << setw(m) << "Msun/pc2" << " ";
              
    theStream << std::endl << "#" << std::endl;

    // This works only if units of map are JY/B * KM/S or JY * KM/S
    double barea = im->Head().Bmaj()*im->Head().Bmin();
    const double arctorad = 1/3600.*M_PI/180.;
    barea *= abs(arcsconv(im->Head().Cunit(0))*arcsconv(im->Head().Cunit(1))); 
    
    for (size_t i=0; i<Nrad; i++) {
        
        // Writing 
        theStream << fixed << setprecision(3) << setw(m) << Radius[i] << " "
                  << scientific << setw(m) << Sum[i][seg] << " "
                  << setw(m) << Mean[i][seg] << " "
                  << setw(m) << Median[i][seg] << " "
                  << setw(m) << sqrt(fabs(Var[i][seg])) << " "
                  << setw(m) << MAD[i][seg] << " "
                  << setw(m) << Area[i][seg] << " "
                  << setw(m) << Surfdens[i][seg] << " "
                  << setw(m) << sqrt(fabs(Var[i][seg]))/fabs(Dx*Dy) << " "
                  << setw(m) << getSurfDensFaceOn(i,seg) << " "
                  << setw(m) << Blankarea[i][seg] << " "
                  << setw(m) << Surfdens_Bl[i][seg] << " "
                  << setw(m+1) << getSurfDensFaceOn_Bl(i,seg) << " ";
        
        if (isJy) {
            // This simply calculates mass surface density for HI data
            double massdsurfdens  = 2.36E-07*getSurfDensFaceOn(i,seg)/(arctorad*arctorad);
            // Calculating Mass surface density with Roberts75 formula.
            double massdsurfdens2 = 8794*Mean[i][seg]/barea*Cosinc[i]*im->Head().BeamArea();
            if ((unit_l.find("/b")!=std::string::npos)) {
                massdsurfdens  /= im->Head().BeamArea();
                massdsurfdens2 /= im->Head().BeamArea();
            }
            // Calculating number density in 1E20 cm^{-2}
            double numberdens = massdsurfdens/(8.41185687E-02*3.0857*3.0857);
            //totmass_req  = 3.0856*3.0856*8.41185687e-22*totmass_req;
            
            theStream << setw(m+1) << numberdens << " "
                      << setw(m+1) << massdsurfdens << " ";
//                      << setw(m) << massdsurfdens2 << " ";
        }
        
        theStream << fixed << std::endl;
    }
}


// Explicit instantiation of the class
template class Ellprof<float>;
template class Ellprof<double>;
}
