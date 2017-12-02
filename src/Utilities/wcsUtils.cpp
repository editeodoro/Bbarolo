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
#include <sstream>
#include <math.h>
#include <wcslib/wcs.h>
#include <wcslib/wcsunits.h>
#include <Utilities/utils.hh>

int pixToWCSSingle(struct wcsprm *wcs, const double *pix, double *world) {

    ///  @details
    ///   Uses wcs to convert the three-dimensional pixel position referenced
    ///   by pix to world coordinates, which are placed in the array world[].
    ///   Assumes that pix only has one point with an x,y,and z pixel positions.
    ///   If there are any other axes present (eg. Stokes), these are set to the
    ///   first pixel in that axis.
    ///   Offsets these pixel values by 1 to account for the C arrays being
    ///   indexed to 0.
    ///
    /// \param wcs The wcsprm struct containing the WCS information.
    /// \param pix The array of pixel coordinates.
    /// \param world The returned array of world coordinates.

    int naxis=wcs->naxis,npts=1,status;
    int specAxis = wcs->spec;
    if(specAxis<0) specAxis=2;
    if(specAxis>=naxis) specAxis = naxis-1;

    double *newpix = new double[naxis*npts];
    // correct from 0-indexed to 1-indexed pixel array
    for(int i=0;i<naxis;i++) newpix[i] = 1.;
    newpix[wcs->lng] += pix[0];
    newpix[wcs->lat] += pix[1];
    newpix[specAxis]+= pix[2];

    int    *stat      = new int[npts];
    double *imgcrd    = new double[naxis*npts];
    double *tempworld = new double[naxis*npts];
    double *phi       = new double[npts];
    double *theta     = new double[npts];
    status=wcsp2s(wcs, npts, naxis, newpix, imgcrd,
                  phi, theta, tempworld, stat);
    if(status>0){
        std::cout << "\nCannot convert to wcs. WCS error code = " << status
                 <<": stat="<<stat[0] << " : " << wcs_errmsg[status] << std::endl;
    }

    //return just the spatial/velocity information
    world[0] = tempworld[wcs->lng];
    world[1] = tempworld[wcs->lat];
    world[2] = tempworld[specAxis];

    delete [] stat;
    delete [] imgcrd;
    delete [] tempworld;
    delete [] phi;
    delete [] theta;
    delete [] newpix;
    return status;
}


int wcsToPixSingle(struct wcsprm *wcs, const double *world, double *pix) {

  ///  @details
  ///   Uses wcs to convert the three-dimensional world coordinate position 
  ///    referenced by world to pixel coordinates, which are placed in the 
  ///    array pix[].
  ///   Assumes that world only has one point (three values eg RA Dec Velocity)
  ///   Offsets the pixel values by 1 to account for the C arrays being 
  ///    indexed to 0.
  /// 
  /// \param wcs The wcsprm struct containing the WCS information.
  /// \param world The array of world coordinates.
  /// \param pix The returned array of pixel coordinates.

  int naxis=wcs->naxis,npts=1,status;
  int specAxis = wcs->spec;
  if(specAxis<0) specAxis=2;
  if(specAxis>=naxis) specAxis = naxis-1;

  double *tempworld = new double[naxis*npts];
  for(int i=0;i<naxis;i++) tempworld[i] = wcs->crval[i];
  tempworld[wcs->lng]  = world[0];
  tempworld[wcs->lat]  = world[1];
  tempworld[specAxis] = world[2];

  int    *stat    = new int[npts];
  double *temppix = new double[naxis*npts];
  double *imgcrd  = new double[naxis*npts];
  double *phi     = new double[npts];
  double *theta   = new double[npts];

  status=wcss2p(wcs, npts, naxis, tempworld, 
        phi, theta, imgcrd, temppix, stat);

  if( status > 0 ){
      std::cout << "\nCannot convert to wcs. WCS error code = " << status
               <<": stat="<<stat[0] << " : " << wcs_errmsg[status] << std::endl;
  }

  pix[0] = temppix[wcs->lng]-1.;
  pix[1] = temppix[wcs->lat]-1.;
  pix[2] = temppix[specAxis]-1.;
  // correct from 1-indexed to 0-indexed pixel array
  //  and only return the spatial & frequency information

  delete [] stat;
  delete [] imgcrd;
  delete [] temppix;
  delete [] phi;
  delete [] theta;
  delete [] tempworld;
  return status;
}

/*
int pixToWCSMulti(struct wcsprm *wcs, const double *pix, 
          double *world, const int npts)
{
  ///  @details
  ///   Uses wcs to convert the three-dimensional pixel positions referenced 
  ///    by pix to world coordinates, which are placed in the array world[].
  ///   pix is assumed to hold the positions of npts points.
  ///   Offsets these pixel values by 1 to account for the C arrays being 
  ///    indexed to 0.
  /// 
  /// \param wcs The wcsprm struct containing the WCS information.
  /// \param pix The array of pixel coordinates.
  /// \param world The returned array of world coordinates.
  /// \param npts The number of distinct pixels in the arrays.

  int naxis=wcs->naxis,status;
  int specAxis = wcs->spec;
  if(specAxis<0) specAxis=2;
  if(specAxis>=naxis) specAxis = naxis-1;

  // correct from 0-indexed to 1-indexed pixel array
  // Add entries for any other axes that are present, keeping the 
  //   order of pixel positions the same
  double *newpix = new double[naxis*npts];
  for(int pt=0;pt<npts;pt++){
    for(int i=0;i<naxis;i++) newpix[pt*naxis+i] = 1.;
    newpix[pt*naxis+wcs->lng]  += pix[pt*3+0];
    newpix[pt*naxis+wcs->lat]  += pix[pt*3+1];
    newpix[pt*naxis+specAxis] += pix[pt*3+2];
  }

  int    *stat      = new int[npts];
  double *imgcrd    = new double[naxis*npts];
  double *tempworld = new double[naxis*npts];
  double *phi       = new double[npts];
  double *theta     = new double[npts];
  status=wcsp2s(wcs, npts, naxis, newpix, imgcrd, 
        phi, theta, tempworld, stat);
  if(status>0){
    std::stringstream statstr;
    statstr<<"|";
    for(int i=0;i<npts;i++) statstr<<stat[i]<<"|";
    DUCHAMPTHROW("pixToWCSMulti","WCS Error Code = " << status <<": stat="<<statstr.str() << " : " << wcs_errmsg[status]);
  }
  else{
    //return just the spatial/velocity information, keeping the
    //  order of the pixel positions the same.
    for(int pt=0;pt<npts;pt++){
      world[pt*3+0] = tempworld[pt*naxis + wcs->lng];
      world[pt*3+1] = tempworld[pt*naxis + wcs->lat];
      world[pt*3+2] = tempworld[pt*naxis + specAxis];
    }
    
  }

  delete [] stat;
  delete [] imgcrd;
  delete [] tempworld;
  delete [] phi;
  delete [] theta;
  delete [] newpix;
  return status;
}

int wcsToPixMulti(struct wcsprm *wcs, const double *world, 
          double *pix, const int npts)
{
  ///  @details
  ///   Uses wcs to convert the three-dimensional world coordinate position 
  ///    referenced by world to pixel coordinates, which are placed in the
  ///    array pix[].
  ///   world is assumed to hold the positions of npts points.
  ///   Offsets the pixel values by 1 to account for the C arrays being 
  ///    indexed to 0.
  /// 
  /// \param wcs The wcsprm struct containing the WCS information.
  /// \param world The array of world coordinates.
  /// \param pix The returned array of pixel coordinates.
  /// \param npts The number of distinct pixels in the arrays.

  int naxis=wcs->naxis,status=0;
  int specAxis = wcs->spec;
  if(specAxis<0) specAxis=2;
  if(specAxis>=naxis) specAxis = naxis-1;

  // Test to see if there are other axes present, eg. stokes
  if(wcs->naxis>naxis) naxis = wcs->naxis;

  // Add entries for any other axes that are present, keeping the 
  //   order of pixel positions the same
  double *tempworld = new double[naxis*npts];
  for(int pt=0;pt<npts;pt++){
    for(int axis=0;axis<naxis;axis++) 
      tempworld[pt*naxis+axis] = wcs->crval[axis];
    tempworld[pt*naxis + wcs->lng]  = world[pt*3+0];
    tempworld[pt*naxis + wcs->lat]  = world[pt*3+1];
    tempworld[pt*naxis + specAxis] = world[pt*3+2];
  }

  int    *stat   = new int[npts];
  double *temppix = new double[naxis*npts];
  double *imgcrd = new double[naxis*npts];
  double *phi    = new double[npts];
  double *theta  = new double[npts];
  status=wcss2p(wcs,npts,naxis,tempworld,phi,theta,
        imgcrd,temppix,stat);
  if(status>0){
    std::stringstream statstr;
    statstr<<"|";
    for(int i=0;i<npts;i++) statstr<<stat[i]<<"|";
    DUCHAMPTHROW("wcsToPixMulti","WCS Error Code = " << status <<": stat="<<statstr.str() << " : " << wcs_errmsg[status]);
  }
  else{
    // correct from 1-indexed to 0-indexed pixel array 
    //  and return just the spatial/velocity information, 
    //  keeping the order of the pixel positions the same.
    for(int pt=0;pt<npts;pt++){
      pix[pt*naxis + 0] = temppix[pt*naxis + wcs->lng] - 1.;
      pix[pt*naxis + 1] = temppix[pt*naxis + wcs->lat] - 1.;
      pix[pt*naxis + 2] = temppix[pt*naxis + specAxis] - 1.;
    }
  }

  delete [] stat;
  delete [] imgcrd;
  delete [] temppix;
  delete [] phi;
  delete [] theta;
  delete [] tempworld;
  return status;
}

double pixelToVelocity(struct wcsprm *wcs, double &x, double &y, double &z, 
               std::string velUnits)
{
  ///  @details
  ///   Uses wcs to convert the three-dimensional pixel position (x,y,z)
  ///   to world coordinates.
  ///   Returns the velocity in the units given by velUnits.
  /// 
  /// \param wcs  The wcsprm struct containing the WCS information.
  /// \param x The x pixel coordinate.
  /// \param y The y pixel coordinate.
  /// \param z The z pixel coordinate.
  /// \param velUnits The string containing the units of the output velocity.
  /// \return The velocity of the pixel.

  int naxis=wcs->naxis,status=0;
  int    *stat   = new int[1];
  double *pixcrd = new double[naxis];
  double *imgcrd = new double[naxis];
  double *world  = new double[naxis];
  double *phi    = new double[1];
  double *theta  = new double[1];

  int specAxis = wcs->spec;
  if(specAxis<0) specAxis=2;
  if(specAxis>=naxis) specAxis = naxis-1;

  // correct from 0-indexed to 1-indexed pixel array by adding 1 to
  // pixel values.
  for(int i=0;i<naxis;i++) pixcrd[i] = 1.;
  pixcrd[wcs->lng] += x;
  pixcrd[wcs->lat] += y;
  pixcrd[specAxis]+= z;
  status=wcsp2s(wcs, 1, naxis, pixcrd, imgcrd, 
        phi, theta, world, stat);
  if(status>0){
    std::stringstream statstr;
    statstr<<"|";
    for(int i=0;i<naxis;i++) statstr<<stat[i]<<"|";
    DUCHAMPTHROW("pixelToVelocity","WCS Error Code = " << status <<": stat="<<stat[0] << " : " << wcs_errmsg[status]);
  }

  double vel = coordToVel(wcs, world[specAxis], velUnits);

  delete [] stat;
  delete [] pixcrd;
  delete [] imgcrd;
  delete [] world;
  delete [] phi;
  delete [] theta;
  
  return vel;
}
 
double coordToVel(struct wcsprm *wcs, const double coord, 
          std::string outputUnits)
{
  ///  @details
  ///   Convert the wcs coordinate given by coord to a velocity in km/s.
  ///   Does this by checking the ztype parameter in wcs to see if it is 
  ///    FREQ or otherwise, and converts accordingly.
  /// 
  /// \param wcs  The wcsprm struct containing the WCS information.
  /// \param coord The input WCS coordinate.
  /// \param outputUnits The string containing the units of the output velocity.
  /// \return The velocity of the input coord value.

  static int errflag = 0;
  double scale, offset, power;
  int specIndex = wcs->spec;
  if(specIndex<0) specIndex = 2;
  if(specIndex>=wcs->naxis) specIndex = wcs->naxis-1;
  int status = wcsunits( wcs->cunit[specIndex], outputUnits.c_str(), 
             &scale, &offset, &power);

  if(status > 0){
    if(errflag==0){
      std::stringstream errmsg;
      errmsg << "WCSUNITS Error Code = " << status << ":"
         << wcsunits_errmsg[status];
      if(status == 10) errmsg << "\nTried to get conversion from " 
                  << wcs->cunit[specIndex]
                  << " to " << outputUnits.c_str();
      //      errmsg << "\nUsing coordinate value instead.\n";
      DUCHAMPERROR("coordToVel", errmsg);
      errflag = 1;
    }
    return coord;
  }
  else return pow( (coord*scale + offset), power);

}

double velToCoord(struct wcsprm *wcs, const float velocity, 
          std::string outputUnits)
{
  ///  @details
  ///   Convert the velocity given to the appropriate world coordinate for wcs.
  ///   Does this by checking the ztype parameter in wcs to see if it is 
  ///    FREQ or otherwise, and converts accordingly.
  /// 
  /// \param wcs  The wcsprm struct containing the WCS information.
  /// \param velocity The input velocity.
  /// \param outputUnits The string containing the units of the input velocity.
  /// \return The WCS coordinate corresponding to the input velocity.

  double scale, offset, power;
  int specIndex = wcs->spec;
  if(specIndex<0) specIndex = 2;
  if(specIndex>=wcs->naxis) specIndex = wcs->naxis-1;
  int status = wcsunits( outputUnits.c_str(), wcs->cunit[specIndex], 
             &scale, &offset, &power);

  if(status > 0){
    std::stringstream errmsg;
    errmsg << "WCSUNITS Error Code = " << status << ":" << wcsunits_errmsg[status];
// << "\nUsing coordinate value instead.\n";
    DUCHAMPERROR("velToCoord",errmsg.str());
    return velocity;
  }
  else return (pow(velocity, 1./power) - offset) / scale;
  
}
*/
