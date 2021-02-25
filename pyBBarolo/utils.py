"""
This module defines some utility functions  
"""

########################################################################
# Copyright (C) 2021 Enrico Di Teodoro
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################

import numpy as np
from astropy.io import fits
from .wrapper import BBaroloWrapper


def emptyfits(axisDim,cdelts,crpixs=None,crvals=None,ctypes=None,cunits=None,\
              data=None,beam=None,obj=None,bunit=None,dtype=np.float32):
    
    """ Defines an empty fits file.
    
        Args:
          axisDim (list):    axis sizes
          cdelts (list):     pixel sizes
          crpixs (list):     central pixels (default: axisDim/2.+1)
          crvals (list):     central values (default: 0)
          ctypes (list):     axis types (default: ra, dec, velo-helo)
          cunits (list):     axis units (default: deg, deg, km/s)
          beam (float,list): beam size in degrees
          obj (str):         name of the object
          bunit (str):       flux unit
          dtype (type):      type of data (default: float32)
    
        Returns:
          A fits.PrimaryHDU() object 
    """
    # Default cunits and ctypes
    ctypes = ['RA---SIN', 'DEC--SIN', 'VELO-HEL'] if ctypes is None else ctypes
    cunits = ['deg', 'deg', 'km/s'] if cunits is None else cunits
    axisDim = np.array(axisDim)
    
    # Define new header
    h = fits.Header()
    # Adding WCS informations
    for i in range (len(axisDim)):
        h['CRPIX%d'%(i+1)] = axisDim[i]/2.+1 if crpixs is None else crpixs[i]
        h['CRVAL%d'%(i+1)] = 0. if crpixs is None else crvals[i]
        h['CDELT%d'%(i+1)] = cdelts[i]
        h['CUNIT%d'%(i+1)] = cunits[i]
        h['CTYPE%d'%(i+1)] = ctypes[i]
    
    # Adding beam information
    if beam is not None:
        bmaj = bmin = bpa = 0
        if isinstance(beam,(int,float)): 
            bmaj, bmin, bpa = beam, beam, 0
        elif isinstance (beam,(list,tuple,np.ndarray)) and len(beam)<4:
            if   len(beam)==3: bmaj, bmin, bpa = beam
            elif len(beam)==2: bmaj, bmin, bpa = beam[0], beam[1], 0
            elif len(beam)==1: bmaj, bmin, bpa = beam[0], 0, 0
        h['BMAJ'],h['BMIN'],h['BPA'] = bmaj, bmin, bpa
    
    if obj is not None: h['OBJECT'] = obj
    if bunit is not None: h['BUNIT'] = bunit
    
    # Defining empty data
    if data is None:
        d = np.zeros(shape=np.flip(axisDim),dtype=dtype)
    else:
        ok = True
        for i in range (len(axisDim)):
            ok *= np.flip(axisDim)[i]==data.shape[i]
        if ok:  d = data
        else: raise ValueError("Array dimensions are different than axisDim")
    
    
    return fits.PrimaryHDU(data=d.astype(dtype),header=h)
    

class SimulatedGalaxy(object):
    def __init__(self):
        ...
    
    
    
def getRandomGalaxyParameters(radmin,radmax,pixsize):
    radii = np.arange(radmin,radmax,pixsize)
    
    # Generic exponential density profile with hole
    dens = 50*np.exp(-radii/400-200/(0.5*radii+100))

    # Disc scaleheight
    z0 = np.full(len(radii),10)
    #z0[radii>500] = (radii[radii>500]-500)/50.+z0[0]
    inc = int(np.random.uniform(30,80))
    inc = np.full(len(radii),inc)
    
    # Paramaeters for PA
    pa = np.full(len(radii),0.)
    iswarppa = np.random.choice([True, False])
    if iswarppa:
        rstart  = np.random.randint(0,len(radii))
        pastart = 0
        pastop  = np.random.uniform(pastart-30,pastart+30)
        mline   = (pastop-pastart)/(radii[-1]-radii[rstart])
        pa[radii>radii[rstart]] = mline*(radii[radii>radii[rstart]]-radii[-1])+pastop

    # Parameters for vrot (arctan)
    isdwarf = False#np.random.choice([True, False])
    if isdwarf:
        vflat = np.random.uniform(10,130)
        rs = np.random.uniform(pixsize,radmax/2)
    else:
        vflat = np.random.uniform(100,300)
        rs = np.random.uniform(pixsize,radmax/10)
    vrot  = 2./np.pi*vflat*np.arctan(radii/rs)

    # Parameters for vrad
    vrad = np.full(len(radii),0.)
    rstart = np.random.randint(0,2*len(radii)/3)
    vradstart = 0
    vradstop  = np.random.uniform(vradstart-30,vradstart+30)
    mline   = (vradstop-vradstart)/(radii[-1]-radii[rstart])
    vrad[radii>radii[rstart]] = mline*(radii[radii>radii[rstart]]-radii[-1])+vradstop

    return (radii,vrot,vrad,inc,pa,z0,dens)

