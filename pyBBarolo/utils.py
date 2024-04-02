"""
This module defines some utility functions and classes  
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

import os
import numpy as np
from astropy.io import fits
from .wrapper import BBaroloWrapper, BBexe
from .pyBBarolo import isIterable, isNumber


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
    crvals = [60, 0, 0] if crvals is None else crvals
    
    axisDim = np.array(axisDim)
    
    # Define new header
    h = fits.Header()
    # Adding WCS informations
    for i in range (len(axisDim)):
        h['CRPIX%d'%(i+1)] = (axisDim[i]-1)/2.+1 if crpixs is None else crpixs[i]
        h['CRVAL%d'%(i+1)] = crvals[i]
        h['CDELT%d'%(i+1)] = cdelts[i]
        h['CUNIT%d'%(i+1)] = cunits[i]
        h['CTYPE%d'%(i+1)] = ctypes[i]
    
    if h['CDELT1']>0:
        print ("WARNING, CDELT1 should be negative!")
        
    # Adding beam information
    if beam is not None:
        bmaj = bmin = bpa = 0
        if isNumber(beam): 
            bmaj, bmin, bpa = beam, beam, 0
        elif isIterable(beam) and len(beam)<4:
            if   len(beam)==3: bmaj, bmin, bpa = beam
            elif len(beam)==2: bmaj, bmin, bpa = beam[0], beam[1], 0
            elif len(beam)==1: bmaj, bmin, bpa = beam[0], 0, 0
        h['BMAJ'],h['BMIN'],h['BPA'] = bmaj, bmin, bpa
    
    h['OBJECT'] = 'NONE'
    if obj is not None: h['OBJECT'] = obj
    if bunit is not None: h['BUNIT'] = bunit
    h['EQUINOX'] = 2000.
    
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
    


class SimulatedGalaxyCube(object):
    """ 
    A class to easily simulate an emission line datacube for a  
    galaxy with some parameters
    
    Attributes
    ----------
    f: (astropy.PrimaryHDU)
      A FITS HDU representing the emission line datacube 
    
    radii: (list, np.array)
      Radii of the galaxy model
    
    dens: (list, np.array)
      Surface-density profile of the galaxy
    
    inc: (list, np.array)
      Inclination angles of the galaxy
    
    pa: (list, np.array)
      Position angles of the galaxy
    
    z0: (list, np.array)
      Scale heights of the galaxy
    
    vrot: (list, np.array)
      Rotation curve of the galaxy
    
    disp: (list, np.array)
      Velocity dispersions of the galaxy
    
    vrad: (list, np.array)
      Radial velocities of the galaxy

    xpos: (list, np.array)
      X center of the galaxy
    
    ypos: (list, np.array)
      Y center of the galaxy
    
    vsys: (list, np.array)
      Systemic velocity of the galaxy
    
    
    Methods
    -------
    define_galaxy(...):
      Define the parameters of the galaxy model
    
    savefile(fileout):
      Write galaxy parameters to a text file
    
    run(exe=None,stdout=None,noise,outfolder,smooth,**kwargs):
      Run GALMOD and build the emission-line datacue model.
    
    """
    def __init__(self,axisDim=None,cdelts=None,**kwargs):
        """ Initialize a SimulatedGalaxyCube instance.
            Here we define the template cube that will be used to create
            the galaxy model.
            If no inputs are given, a default template cube will be used.
            
        Args:
          axisDim (list):   axis dimensions of the template cube
          cdelts (list):    pixel dimensions of the template cube
          kwargs (dict):    Additional parameters to be passed to emptyfits()

        """
        # A default cube size if not given
        if not axisDim: axisDim = [256,256,64]
        if not cdelts: cdelts = [-20./3600.,20./3600.,10]
        
        self.f     = emptyfits(axisDim,cdelts,**kwargs)
        self.radii = None
        self.dens  = None
        self.inc   = None
        self.pa    = None
        self.z0    = None
        self.vrot  = None
        self.disp  = None
        self.vrad  = None
        self.xpos  = None
        self.ypos  = None
        self.vsys  = None
        
    
    
    def define_galaxy(self,radmin=None,radmax=None,radii=None,dens=None,\
                      ishole=None,xpos=None,ypos=None,vsys=None,\
                      inc=None,warpinc=None,mininc=30,maxinc=89,maxiwarp=15,\
                      pa=None,warppa=None,minpa=0,maxpa=359,maxpawarp=30,\
                      vrot=None,vrotdwarf=None,vdisp=None,z0=None,flare=False,\
                      vrad=None,radmotion=False,maxvrad=20):
        """ Define the parameters of the model galaxy and sets all the variables
            needed by the class. This function needs to be called before run().
                      
            All inputs are optional. The user can decide to set any of them to 
            control the parameters of the galaxy model. If none is set, the galaxy
            will be completely random. 
            
        Args:
          radmin (float):      Minimum radius of the galaxy (arcsec).
          radmax (float):      Maximum radius of the galaxy (arcsec).
          radii (list):        List of radii (arcsec). If None, assume cubesize/3.
          dens (float,list):   Surface density (1E20 cm-2). If None, assume a profile.
          ishole (bool):       Whether to have a central hole in dens. If None, randomly decide.
          xpos (float,list):   X center (pixel). If None, assume the center of Xsize/2.
          ypos (float,list):   Y center (pixel). If None, assume the center of Ysize/2.
          vsys (float,list):   Systemic velocity (pixel). If None, assume the Zsize/2.
          inc (float,list):    Inclination angle (deg). If None, assume a random inc.
          warpinc (bool):      Whether to set a warp in inclination. If None, randomly decide.
          mininc (float):      Minimum inclination for random draw when inc=None.
          maxinc (float):      Maximum inclination for random draw when inc=None.
          maxiwarp (float):    Maximum deviation in the inc warp if warpinc=True.
          pa (float,list):     Position angle (deg). If None, assume a random PA.
          warppa (bool):       Whether to set a warp in PA. If None, randomly decide.
          minpa (float):       Minimum PA for random draw when pa=None.
          maxpa (float):       Maximum PA for random draw when pa=None.
          maxpawarp (float):   Maximum deviation in the PA warp if warppa=True.
          vrot (float,list):   Rotation velocity (km/s). If None, assume a random curve.
          vrotdwarf (bool):    Whether to assume a dwarf-like rotation curve if vrot=None.
          vdisp (float,list):  Velocity dispersion (km/s). If None, assume a random vdisp.
          z0 (float,list):     Disk scale height (arcs). If None, assume a random z0.
          flare (float,list):  Whether to set a flare in z0. If None, randomly decide.
          vrad (float,list):   Radial velocity (km/s). If None, assume 0.
          radmotion (float,list): Whether to set some radial motion. If None, randomly decide.
          maxvrad (float):     Maximum deviation in the vrad if radmotion=True.

        """
        
        h = self.f.header
        pixsize = np.mean([np.abs(h['CDELT1']),np.abs(h['CDELT2'])])*3600.
        
        if radii is None and radmax is None:
            # Setting a maximum radius at cube dimension/3
            radmax = np.round(h['NAXIS2']*pixsize/3.)

        # Defining radii
        if isIterable(radii): self.radii = radii
        else:
            radmin = radmin if radmin else pixsize
            self.radii = np.arange(radmin,radmax+pixsize,pixsize)
        
        radii = self.radii
        radmax = radii[-1]
        
        
        # Defining X centre 
        if isIterable(xpos):
            if len(xpos)==len(radii): self.xpos = xpos
            else: 
                raise AttributeError("radii and xpos must have the same size") 
        elif isNumber(xpos):
            self.xpos = np.full(len(radii),xpos)
        else:
            # Centering the galaxy
            self.xpos = np.full(len(radii),h['NAXIS1']/2.)
            
        
        # Defining Y centre 
        if isIterable(ypos):
            if len(ypos)==len(radii): self.ypos = ypos
            else: 
                raise AttributeError("radii and ypos must have the same size") 
        elif isNumber(xpos):
            self.ypos = np.full(len(radii),ypos)
        else:
            # Centering the galaxy
            self.ypos = np.full(len(radii),h['NAXIS2']/2.)
            
        
        # Defining systemic velocity 
        if isIterable(vsys):
            if len(vsys)==len(radii): self.vsys = vsys
            else: 
                raise AttributeError("radii and ypos must have the same size") 
        elif isNumber(vsys):
            self.vsys = np.full(len(radii),vsys)
        else:
            # Centering the galaxy in velocity
            vs = (h['NAXIS3']/2.+1-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']
            self.vsys = np.full(len(radii),vs)
        
            
        # Defining density profile
        if isIterable(dens):
            if len(dens)==len(radii): self.dens = dens
            else: 
                raise AttributeError("radii and dens must have the same size") 
        elif isNumber(dens):
            self.dens = np.full(len(radii),dens)
        else:
            # Using a random densifity profile
            d0 = 50.
            R0 = radmax/5.
            if ishole is None: ishole =  np.random.choice([True, False])
            if ishole:
                self.dens = d0*np.exp(-radii/R0
                            -R0/2./(0.5*radii+R0/4.))
            else: self.dens = d0*np.exp(-radii/R0)
        
        
        # Defining inclination
        if isIterable(inc):
            if len(inc)==len(radii): self.inc = inc
            else: 
                raise AttributeError("radii and inc must have the same size")
        else:
            if isNumber(inc):
                self.inc = np.full(len(radii),inc)
            else:
                # Using a random inclination
                i = np.random.uniform(mininc,maxinc)
                self.inc = np.full(len(radii),np.round(i))
            
            # Creating a warp in inclination if requested
            if warpinc is None: warpinc = np.random.choice([True, False])
            if warpinc:
                rstart  = np.random.randint(0,len(radii)-1)
                start = self.inc[0]
                stop  = np.random.uniform(start-maxiwarp,start+maxiwarp)
                mline = (stop-start)/(radmax-radii[rstart])
                self.inc[radii>radii[rstart]] = mline*(radii[radii>radii[rstart]]-radmax)+stop
        
        
        # Defining position angle
        if isIterable(pa):
            if len(pa)==len(radii): self.pa = pa
            else: 
                raise AttributeError("radii and pa must have the same size")
        else:
            if isNumber(pa):
                self.pa = np.full(len(radii),pa)
            else:
                # Using a random position angle
                p = np.random.uniform(minpa,maxpa)
                self.pa = np.full(len(radii),np.round(p))
        
            # Creating a warp in position angle if requested
            if warppa is None: warppa = np.random.choice([True, False])
            if warppa:
                rstart  = np.random.randint(0,len(radii)-1)
                start = self.pa[0]
                stop  = np.random.uniform(start-maxpawarp,start+maxpawarp)
                mline = (stop-start)/(radmax-radii[rstart])
                self.pa[radii>radii[rstart]] = mline*(radii[radii>radii[rstart]]-radmax)+stop
                
        
        # Defining scale height
        if isIterable(z0):
            if len(z0)==len(radii): self.z0 = z0
            else: 
                raise AttributeError("radii and z0 must have the same size")
        else:
            if isNumber(z0):
                self.z0 = np.full(len(radii),z0)
            else:
                # Random scale height
                z = np.random.uniform(0,3*pixsize)
                self.z0 = np.full(len(radii),z)
            
            # Creating a flare if requested
            if flare is None: flare = np.random.choice([True, False])
            if flare:
                rstart  = np.random.randint(0,len(radii)-1)
                start = self.z0[0]
                stop  = np.random.uniform(start,start+start)
                mline = (stop-start)/(radmax-radii[rstart])
                self.z0[radii>radii[rstart]] = mline*(radii[radii>radii[rstart]]-radmax)+stop
            

        # Defining rotation curve 
        if isIterable(vrot):
            if len(vrot)==len(radii): self.vrot = vrot
            else: 
                raise AttributeError("radii and vrot must have the same size")
        elif isNumber(vrot):
            self.vrot = np.full(len(radii),vrot)
        else:
            # Random rotation curve
            if vrotdwarf is None: vrotdwarf = np.random.choice([True, False])
            if vrotdwarf:
                vflat = np.random.uniform(10,130)
                rs = np.random.uniform(pixsize,radmax/2)
            else:
                vflat = np.random.uniform(100,350)
                rs = np.random.uniform(pixsize,radmax/10)
            self.vrot  = 2./np.pi*vflat*np.arctan(radii/rs)
            
            
        # Defining velocity dispersion 
        if isIterable(vdisp):
            if len(vdisp)==len(radii): self.disp = vdisp
            else: 
                raise AttributeError("radii and vdisp must have the same size")
        elif isNumber(vdisp):
            self.disp = np.full(len(radii),np.round(vdisp))
        else:
            # Random constant dispersion
            vd = np.random.uniform(6,15)
            self.disp = np.full(len(radii),np.round(vd)) 
        
        
        # Defining radial motions
        if isIterable(vrad):
            if len(vrad)==len(radii): self.vrad = vrad
            else: 
                raise AttributeError("radii and vrad must have the same size")
        else:
            if isNumber(vrad):
                self.vrad = np.full(len(radii),vrad)
            else:
                # No radial motions 
                self.vrad = np.full(len(radii),0.)
            
            if radmotion is None: radmotion = np.random.choice([True, False])
            if radmotion:
                rstart  = np.random.randint(0,len(radii)-1)
                start = self.vrad[0]
                stop  = np.random.uniform(start-maxvrad,start+maxvrad)
                mline = (stop-start)/(radmax-radii[rstart])
                self.vrad[radii>radii[rstart]] = mline*(radii[radii>radii[rstart]]-radmax)+stop
            
            

    def savefile(self,fileout="galaxy_params.txt"):
        """ Save galaxy parameter to a text file.
            
        Args:
          fileout (str):  Name of text file
        """
        if self.radii is not None:
            np.savetxt(fileout,np.transpose([self.radii,self.vrot,self.vrad,self.disp,self.inc,
            self.pa,self.z0,self.dens,self.xpos,self.ypos,self.vsys]),
            header="RADII         VROT    VRAD   VDISP      INC       PA      Z0    DENS    XPOS    YPOS    VSYS",
            fmt='%-12.5f %7.1f %7.1f %7.1f %8.2f %8.2f %7.2f %7.2f %7.1f %7.1f %7.1f')
        
        
        
    def run(self,exe=BBexe,stdout=None,noise=None,outfolder=None,smooth=True,**kwargs):
        """ Execute BBarolo to create a galaxy model
            
        Args:
          exe (str):           BBArolo's executable
          stdout (str):        Where to print BBarolo's output (see BBaroloWrapper.run())
          noise (float, bool): Noise of the output datacube. If True, random select it.
          outfolder (str):     Output folder where to write the model datacube.
          smooth (bool):       Whether to smooth the model to the resolution in the template cube
          kwargs (dict):       Other parameters to be passed to BB's GALMOD task.
        """
        if self.radii is None:
            raise RuntimeError("Please, define a model with define_galaxy() before calling run()")
        
        # If noise is not given, put a random noise
        noiserms = 0.
        if noise==True:
            noiserms = np.random.uniform(0,0.05)
        elif isNumber(noise): noiserms = noise
            
        try: self.f.header['BMAJ']>0
        except: smooth = False
        
        # Defining outfolder
        obj = self.f.header['OBJECT'] 
        if not outfolder: outfolder = obj
        
        # Saving a text file with parameters
        p = "galaxy_params.txt"
        self.savefile(fileout=p)
        
        # Saving the empty cube
        self.f.writeto("emptycube.fits",overwrite=True)
        opts = dict(fitsfile='emptycube.fits',
                    galmod=True,
                    radii=f'file({p},1)',
                    vrot=f'file({p},2)',
                    vrad=f'file({p},3)',
                    vdisp=f'file({p},4)',
                    inc=f'file({p},5)',
                    pa=f'file({p},6)',
                    z0=f'file({p},7)',
                    dens=f'file({p},8)',
                    xpos=f'file({p},9)',
                    ypos=f'file({p},10)',
                    vsys=f'file({p},11)',
                    norm=None,
                    mask=None,
                    linear=0,
                    cdens=10,
                    noiserms=noiserms,
                    outfolder=outfolder,
                    sm=smooth)
        
        # Deleting any entry that is overriden by kwargs
        for x in kwargs: opts.pop(x,None)
        
        bb = BBaroloWrapper(**opts,**kwargs)
        bb.run(exe=exe,stdout=stdout)
        
        # Copy galaxy parameter file to outfolder
        #bb.write_parameterfile(f'{outfolder}/BBparams.par')
        os.system(f'cp {p} {outfolder}/{obj}_params.txt')
