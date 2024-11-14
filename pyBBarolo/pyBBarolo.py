"""
This module defines several classes to interface pyBBarolo with
the underlying BBarolo C++ code, with some extra functionality.
"""

########################################################################
# Copyright (C) 2017 Enrico Di Teodoro
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

from __future__ import print_function, division
import os,sys
import numpy as np
from .BB_interface import libBB
from astropy.io import fits


def reshapePointer (p, shape):
    """Take a POINTER to c_float and reshape it.
    
    Args:
      p (POINTER(c_float)): The pointer to be reshaped
      shape (tuple, list): The new shape of the array
    
    Returns:
      ndarray: The reshaped array
    
    """
    return np.ctypeslib.as_array(p, shape=tuple(shape))


def isIterable (p):
    """Check if p is an iteratable (list,tuple or numpy array). """
    return isinstance(p,(list,tuple,np.ndarray))


def isNumber (p):
    """Check if p is an iteratable (list,tuple or numpy array). """
    return isinstance(p,(float,int,np.int32,np.int64,np.float32,np.float64))


class FitsCube(object):
    def __init__(self,fitsname):
        # File name of the FITS file to read
        self.fname = fitsname
        # Pointer to the C++ Cube object
        self._cube = None
        if not os.path.isfile(fitsname):
            raise ValueError("FitsCube ERROR: file %s does not exist."%self.fname)
        self._cube = libBB.Cube_new(self.fname.encode('utf-8'))
        # Axis dimensions
        self.dim   = libBB.Cube_axisdim(self._cube)
        # Getting array data and reshaping in (z,y,x)
        self.data  = reshapePointer(libBB.Cube_array(self._cube),self.dim[::-1])
        # Beam information
        self.beam  = libBB.Cube_getBeam(self._cube)
        # Also having a pointer to an Astropy object
        self.fapy  = fits.open(fitsname)[0]
        
        
    def __del__(self):
        if self._cube: libBB.Cube_delete(self._cube)
        
    
    def setBeam(self, bmaj, bmin, bpa=0):
        """ Change the Beam parameters
        
        Args:
          bmaj (float): Beam major axis in degrees
          bmin (float): Beam minor axis in degrees
          bpa  (float): Beam position angle in degrees (default 0.0)
        
        """
        libBB.Cube_setBeam(self._cube,np.float32(bmaj),np.float32(bmin),np.float32(bpa))



class Param(object):
    """Wrapper class for C++ Param class (param.hh)"""
    
    def __init__(self,**kwargs):
        # A dictionary with parameter names and values
        self.opts = {}
        # Pointer to the C++ Param object
        self._params = libBB.Param_new();
        self.paramDefined = False
        if kwargs:
            self.add_params(**kwargs)
        
    def add_params(self,**kwargs):
        """ Add new parameter to the parameter list """
        if kwargs: self.opts.update(kwargs)

    def add_params_opts(self,**kwargs):
        """ Add new parameter to the parameter list from _opts variable"""
        for key, value in kwargs.items():
            self.opts.update({key : value[0]})
    
    def remove_param(self,toremove):
        """ Remove an option from the parameter list """
        self.opts.pop(toremove, None)
        
    def write_parameterfile(self,fileout='param.par'):
        """ Write a BBarolo's parameter file in fileout """
        with open(fileout,'w') as f:
            f.write(self.__str__())

    def make_object(self):
        s = self.__str__()
        if self.paramDefined: libBB.Param_delete(self._params)
        libBB.Param_setfromstr(self._params,s.encode('utf-8'))
        self.paramDefined = True

    def __str__(self):
        s = "##### Input parameters for BBarolo #####\n"
        for k in self.opts:
            s += f'{k.upper():18s} {self.opts[k]} \n' 
        return s
    
    def __del__(self):
        if self.paramDefined: libBB.Param_delete(self._params)



class Rings(object):
    """Wrapper class for C++ Rings structure (rings.hh)"""
    
    def __init__(self,nrings):
        # Number of rings
        if not isNumber(nrings): raise ValueError("nrings must be a number!")
        self.nr     = nrings
        # A dictionary that stores the rings values
        self.r = dict(radii=None,xpos=None,ypos=None,vsys=None,vrot=None,vdisp=None,vrad=None,\
                      vvert=None,dvdz=None,zcyl=None,dens=None,z0=None,inc=None,phi=None)
        # Pointer to the C++ Rings object
        self._rings = None

    def set_rings (self,radii,xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens,z0,inc,phi):
        """ Define rings given the input parameters 
            
           radii: List or array
           other: float or array
        """
        
        if isIterable(radii): self.modify_parameter('radii',radii)
        else: raise ValueError("radii must be an array")
        
        self.modify_parameter('xpos',xpos)
        self.modify_parameter('ypos',ypos)
        self.modify_parameter('vsys',vsys)
        self.modify_parameter('vrot',vrot)
        self.modify_parameter('vdisp',vdisp)
        self.modify_parameter('vrad',vrad)
        self.modify_parameter('vvert',vvert)
        self.modify_parameter('dvdz',dvdz)
        self.modify_parameter('zcyl',zcyl)
        self.modify_parameter('dens',dens)
        self.modify_parameter('z0',z0)
        self.modify_parameter('inc',inc)
        self.modify_parameter('phi',phi)

        self.make_object()

    def modify_parameter(self,pname,pvalue,makeobj=False):
        """ Modifies one of the ring parameters (e.g., vrot, vdisp, etc.)
            
           pname (str): Name of the parameter
           other: float or array
        """
        if pname not in self.r:
            raise ValueError("ERROR: Unknown ring parameter %s"%pname)
        self.r[pname] = np.array(pvalue,dtype=np.float32) if isIterable(pvalue) \
                                                          else np.full(self.nr,pvalue,dtype=np.float32)
        if len(self.r[pname])!=self.nr: raise ValueError("All parameters must have size = %i"%self.nr)
        if makeobj: self.make_object()

    def make_object(self):
        """ Creates and stores the C++ Rings object """
        self.__del__()
        self._rings = libBB.Rings_new();
        libBB.Rings_set(self._rings,self.nr,self.r['radii'],self.r['xpos'],self.r['ypos'],self.r['vsys'],\
                        self.r['vrot'],self.r['vdisp'],self.r['vrad'],self.r['vvert'],self.r['dvdz'],\
                        self.r['zcyl'],self.r['dens'],self.r['z0'],self.r['inc'],self.r['phi'])
    
    def __del__(self):
        if self._rings is not None: libBB.Rings_delete(self._rings)



class Task(object):
    """Superclass for all BBarolo tasks
    
    Args:
      fitsname (str): Input FITS file
    
    """
    def __init__(self,fitsname):
        # Task name
        self.taskname = None
        # Input datacube 
        self.inp = FitsCube(fitsname)
        # Mandatory arguments of the task
        self._args = {}
        # Options for the task
        self._opts = {}
        # A Param object 
        self._par = Param(fitsfile=fitsname)
        
    
    
    def init(self,*args,**kwargs):
        """Initialize a task. Refer to derived class for *args and **kwargs."""
        self._input(*args,**kwargs)
        
        
    def set_options(self,**kwargs):
        """Set options for the task. Keywords **kwargs are given in self._opts."""
        for key in kwargs:
            if key in self._opts:
                self._opts[key][0] = self._opts[key][1](kwargs[key]) 
            else:
                raise ValueError('Option %s unknown. Try show_options() for a list of keywords'%key)
        
        self._check_options()
    
        
    def show_options(self):
        """ Show current options for the task """
        print ("\nCurrent options for %s task: %s "%(self.taskname,"-"*(43-len(self.taskname))))
        for key, value in list(self._opts.items()):
            print (" %-12s = %-10s # %s  "%(key,value[0],value[2]))
        
        print ("\n BEAM = %.4f x %.4f arcs (BPA = %.4f deg)"\
               %(self.inp.beam[0]*3600,self.inp.beam[1]*3600,self.inp.beam[2]))
        print ("%s\n"%("-"*70))
        
    
    def set_arguments(self,**kwargs):
        """Set options for the task. Keywords **kwargs are given in self._opts."""
        for key in kwargs:
            if key in self._args:
                self._args[key][0] = kwargs[key] 
            else:
                raise ValueError('Argument %s unknown. Try show_arguments() for a list of keywords'%key)
        
    
    def show_arguments(self):
        """ Show needed arguments for the task """
        print ("\nArguments for %s task: %s "%(self.taskname,"-"*(49-len(self.taskname))))
        if self.taskname=='3DFIT':
            print ("If not set, BBarolo estimates the initial guesses (see documentation)\n")
        for key, value in list(self._args.items()):
            print (" %-8s # %s  "%(key,value[1]))
            
        print ("%s\n"%("-"*70))


    def compute(self,threads=1,**kwargs):
        """ Compute the model 
        
        This function needs to be called after :func:`init`.
        """
        if not isinstance(threads,int):
            raise ValueError("%s ERROR: threads must and integer."%self.taskname)
        return self._compute(threads,**kwargs)



class Model3D(Task):
    """Superclass for all 3D models, derived from Task
    
    Args:
      fitsname (str): FITS file of the galaxy to model
    
    """
    def __init__(self,fitsname):
        # A pointer to the C++ model
        self._mod = None
        self._modCalculated = False
        # A Rings object with rings
        self._inri = None
        # The output model cube (astropy PrimaryHDU)
        self.outmodel = None
        super(Model3D,self).__init__(fitsname=fitsname)
        self._opts.update ({'cdens' : [10, np.int32, "Surface density of clouds in a ring (1E20)"],
                            'nv'    : [-1, np.int32, "Number of subclouds per profile"]})
        

    def smooth(self,beam=None):
        """Smooth the model and return the smoothed array 
        
        This function needs to be called after :func:`compute`. The smoothed array is written
        in outmodel instance and returned.
        
        Args:
          beam (tuple,list): The beam in the form (bmaj,bmin,bpa), bmaj and bmin in arcs.
                                     If None, use the beam from the FITS header.
        
        Returns:
           astropy PrimaryHDU: a datacube with the output smoothed model
        
        """
        if self.outmodel is None:
            raise ValueError("%s ERROR: you need to compute the model before calling smooth()."%self.taskname)
        if beam is None:
            # Check that Beam info is in the header
            if self.inp.beam[0]==-1 or self.inp.beam[1]==-1:
                raise RuntimeError("Beam info is not available in the header.\
                                    Please pass beam=(bmaj,bmin,bpa) paramater to smooth()")
        else:
            if isIterable(beam) and len(beam)==3:
                set_beam(bmaj=beam[0],bmin=beam[1],bpa=beam[2])
            else: 
                raise ValueError("%s ERROR: beam is a list of [bmaj,bmin,bpa]."%self.taskname)            
        return self._smooth()
        

    def set_beam(self, bmaj, bmin, bpa=0):
        """ Change the Beam parameters
        
        Args:
          bmaj (float): Beam major axis in arcsec
          bmin (float): Beam minor axis in arcsec
          bpa  (float): Beam position angle in degrees (default 0.0)
        
        """
        self.inp.setBeam(bmaj/3600.,bmin/3600.,bpa)



class GalMod(Model3D):
    """Produce a simulated 3D datacube of disk galaxy
    
    Args:
      fitsname (str): FITS file of the galaxy to model
    
    """
    def __init__(self,fitsname):
        super(GalMod,self).__init__(fitsname=fitsname)
        self.taskname = "GALMOD"
        self._opts.update({'ltype' : [2, np.int32, "Layer type along z"],
                           'cmode' : [1, np.int32, "Mode of clouds-surface density"],
                           'iseed' : [-1, np.int32, "Seed for random number generator"]})
        self._args = {'radii': [None,'Radii of the model in arcsec (must be an array)'],
                      'xpos' : [None,'X-center of the galaxy in pixels'],
                      'ypos' : [None,'Y center of the galaxy in pixels'],
                      'vsys' : [None,'Systemic velocity of the galaxy in km/s'],
                      'vrot' : [None,'Rotation velocity in km/s'],
                      'vdisp': [None,'Velocity dispersion in km/s'],
                      'vrad' : [None,'Radial velocity in km/s'],
                      'vvert': [None,'Vertical velocity in km/s'],
                      'dvdz' : [None,'Vertical rotational gradient (km/s/arcsec)'],
                      'zcyl' : [None,'Height where the rotational gradient starts (arcsec)'],
                      'dens' : [None,'Surface density of gas in 1E20 atoms/cm2'],
                      'z0'   : [None,'Disk scaleheight in arcsec'],
                      'inc'  : [None,'Inclination angle in degrees'],
                      'phi'  : [None,'Position angle of the receding part of the major axis (N->W)']}

    def __del__(self):
        if self._mod: libBB.Galmod_delete(self._mod)
        
        
    def _input(self,radii,xpos,ypos,vsys,z0,inc,phi,vrot,vdisp,dens=1,vrad=0,vvert=0,dvdz=0,zcyl=0):
        """ Initialize rings for the model
        
        Args:
          radii (list):        Radii of the model in arcsec
          xpos (float, list):  X center of the galaxy in pixels
          ypos (float, list):  Y center of the galaxy in pixels
          vsys (float, list):  Systemic velocity of the galaxy in km/s
          vrot (float, list):  Rotation velocity in km/s
          vdisp (float, list): Velocity dispersion in km/s
          vrad (float, list):  Radial velocity in km/s
          vvert (float, list): Vertical velocity in km/s
          dvdz (float, list):  Vertical rotational gradient (km/s/arcsec).
          zcyl (float, list):  Height where the rotational gradient starts (arcsec).
          dens (float, list):  Surface density of gas in 1E20 atoms/cm2
          z0 (float, list):    Disk scaleheight in arcsec
          inc (float, list):   Inclination angle in degrees
          phi (float, list):   Position angle of the receding part of the major axis (N->W)
    
        """
        
        if not isIterable(radii): raise ValueError("radii must be an array")
        self._inri = Rings(len(radii))
        self._inri.set_rings(radii,xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens*1E20,z0,inc,phi)
    
        
    def _compute(self,threads=1,**kwargs):
        """ Compute the model 
        
        This function needs to be called after :func:`input`.
        
        Returns:
          astropy PrimaryHDU: a datacube with the output model
        """
        if self._inri is None: 
            raise ValueError("GALMOD ERROR: you need to set the model with init(...) before calling compute().")
        
        self._check_options()
        self._par.add_params_opts(**self._opts)
        self._par.add_params(threads=threads,**kwargs)
        self._par.make_object()
        self._mod = libBB.Galmod_new_par(self.inp._cube,self._inri._rings,self._par._params)
        self._modCalculated = libBB.Galmod_compute(self._mod)
        data_mod  = reshapePointer(libBB.Galmod_array(self._mod),self.inp.dim[::-1])
        
        self.outmodel = fits.PrimaryHDU(data_mod)
        self.outmodel.header = self.inp.fapy.header
        return self.outmodel
    
    
    def _smooth(self):
        """ Smooth the model """
        libBB.Galmod_smooth(self._mod)
        self.outmodel.data  = reshapePointer(libBB.Galmod_array(self._mod),self.inp.dim[::-1])
        return self.outmodel
        
        
    def _check_options(self):
        """ Check if current options are ok for the C++ code """
        op = self._opts
        for key in op:
            op[key][0] = op[key][1](op[key][0])
            if key=='ltype' and op[key][0] not in range(1,6):
                raise ValueError("%s ERROR: %s can only be 1, 2, 3, 4 or 5."%(self.taskname,key))
            if key=='cmode' and op[key][0] not in range(1,3):
                raise ValueError("%s ERROR: %s can only be 1, 2."%(self.taskname,key))
            if key=='cdens' and op[key][0]<=0:
                raise ValueError("%s ERROR: %s can only be > 0."%(self.taskname,key))
            if key=='iseed' and op[key][0]>=0:
                raise ValueError("%s ERROR: %s can only be < 0."%(self.taskname,key))



class GalWind(Model3D):
    """Produce a simulated 3D datacube of a biconical outflow  
    
    Args:
      fitsname (str): FITS file of the galaxy to fit
    
    """
    def __init__(self,fitsname):
       super(GalWind,self).__init__(fitsname=fitsname)
       # Task name
       self.taskname = "GALWIND"
       self._opts.update({'ntot'  : [25, np.int32, "Total number of cylinders"],
                          'dtype' : [1, np.int32, "Vertical density type"]})
       self._opts["nv"][0] = 10
       self._args = {'xpos'   : [None, 'X-center of the galaxy in pixels'],
                     'ypos'   : [None, 'Y center of the galaxy in pixels'],
                     'vsys'   : [None, 'Systemic velocity of the galaxy in km/s'],
                     'inc'    : [None, 'Inclination angle in degrees'],
                     'phi'    : [None, 'Position angle of the receding part of the major axis (N->W)'],
                     'vdisp'  : [None, 'Velocity dispersion in km/s'],
                     'dens'   : [None, 'Surface density of gas in 1E20 atoms/cm2'],
                     'vwind'  : [None, 'Radial wind velocity of cylinders in km/s'],
                     'openang': [None, 'Wind opening angle in degrees'],
                     'htot'   : [None, 'Wind truncation height in arcsec']}
       self.ready = False
        
        
    def __del__(self):
        if self._mod: libBB.Galwind_delete(self._mod)   
            
            
    def _input(self,xpos,ypos,vsys,inc,phi,vdisp,dens,vwind,openang,htot):
        """ Set the parameters of the model.
        
        The model is built by using ntot cylinders at increasing heights from the galaxy plane. 
        Velocity, velocity dispersion and density can vary with cylinder, other parameters are 
        fixed.
        
        Args:
          xpos (float):       X center of the galaxy in pixels.
          ypos (float):       Y center of the galaxy in pixels.
          phi (float):        Position angle of the galaxy in degrees (N->W)
          inc (float):        Inclination angle of the galaxy in degrees
          vdisp (float, list):Velocity dispersion of cylinders in km/s
          dens (float, list): Surface density in units of 1E20 atoms/cm2
          vwind (float,list): Radial wind velocity of cylinders
          openang(float,list): Wind opening angle
          htot (float):       Wind truncation height in arcsec 
        """        
        super(GalWind,self).set_arguments(xpos=xpos,ypos=ypos,vsys=vsys,inc=inc,phi=phi,vdisp=vdisp,\
                                          dens=dens,vwind=vwind,openang=openang,htot=htot)
        self.ready = True
    

    def _compute(self,threads=1):
        """ Compute the model 
        
        This function needs to be called after :func:`input`.
        
        Returns:
          astropy PrimaryHDU: a datacube with the output model
        """
                
        if self.ready:
            self._check_options()
            op = self._opts
            ar = self._args
            if self._mod: self.__del__()
            self._mod = libBB.Galwind_new(self.inp._cube,ar['xpos'][0],ar['ypos'][0],ar['phi'][0],\
                                          ar['inc'][0],ar['vdisp'][0],ar['dens'][0],ar['vsys'][0],\
                                          ar['vwind'][0],ar['openang'][0],ar['htot'][0],op["dtype"][0],\
                                          op["ntot"][0],op["cdens"][0],op["nv"][0],int(threads))
            
            self._modCalculated = libBB.Galwind_compute(self._mod)
            data_mod  = reshapePointer(libBB.Galwind_array(self._mod),self.inp.dim[::-1])
        else:
            raise ValueError("GALWIND ERROR: you need to set the model with init(...) before calling compute().")
        
        self.outmodel = fits.PrimaryHDU(data_mod)
        self.outmodel.header = self.inp.fapy.header
        return self.outmodel
        
    
    def _smooth(self):
        """ Smooth the model """
        libBB.Galwind_smooth(self._mod)
        self.outmodel.data  = reshapePointer(libBB.Galwind_array(self._mod),self.inp.dim[::-1])
        return self.outmodel
    
    
    def _check_options(self):
        """ Check if current options are ok for the C++ code """
        op = self._opts
        for key in op:
            op[key][0] = op[key][1](op[key][0])
            if key=='cdens' and op[key][0]<=0:
                raise ValueError("%s ERROR: %s can only be > 0."%(self.taskname,key))
            if key=='ntot' and op[key][0]<=0:
                raise ValueError("%s ERROR: %s can only be > 0."%(self.taskname,key))
            if key=='dtype' and op[key][0] not in range(0,3):
                raise ValueError("%s ERROR: %s can only be 0, 1, 2"%(self.taskname,key))



class FitMod3D(Model3D):
    """Fit a galaxy model to a 3D datacube 
    
    Args:
      fitsname (str): FITS file of the galaxy to fit
    
    """
    
    def __init__(self,fitsname):
        super(FitMod3D,self).__init__(fitsname=fitsname)
        # Task name
        self.taskname = "3DFIT"
        # The output final rings
        self.bfit = None
        # A dictionary with options and defaults values
        self._opts.update({'ltype'   : [2, np.int32, "Layer type along z"],
                           'smooth'  : [True, bool, "If false, disable smoothing"],
                           'deltainc': [5., np.float32,"Inclination angle variation (degrees)"],
                           'deltaphi': [15., np.float32, "Position angle variation (degrees)"],
                           'ftype'   : [2, np.int32, "Residual function to minimize" ],
                           'wfunc'   : [2, np.int32, "Weighting function for major axis"],
                           'bweight' : [1, np.int32, "Weighting function for Blank pixels"],
                           'tol'     : [1E-03, np.float64, "Tolerance for minimization."],
                           'mask'    : ['SMOOTH', str, "Mask type"],
                           'norm'    : ['LOCAL', str, "Normalization type"],
                           'free'    : ['VROT VDISP', str, "Free parameters"],
                           'side'    : ['B', str, "Which side of the galaxy to fit"],
                           'twostage': [True, bool, "Regularize and fit a second model"],
                           'polyn'   : ['bezier', str, "Type of regularization"],
                           'startrad': [0, np.int32, "Starting radius"],
                           'errors'  : [False, bool, "Whether estimating errors"],
                           'distance': [-1., np.float32, "Distance of the galaxy in Mpc"],
                           'redshift': [-1., np.float64, "Redshift of the galaxy"],
                           'restwave': [-1., np.float64, "Rest wavelength of observed line"],
                           'outfolder' : ['./', str, "Directory for outputs" ]})
        self._args = {'radii': [None, 'Radii of the model in arcsec (must be an array)'],
                      'xpos' : [None, 'X-center of the galaxy in pixels'],
                      'ypos' : [None, 'Y center of the galaxy in pixels'],
                      'vsys' : [None, 'Systemic velocity of the galaxy in km/s'],
                      'vrot' : [None, 'Rotation velocity in km/s'],
                      'vdisp': [None, 'Velocity dispersion in km/s'],
                      'vrad' : [None, 'Radial velocity in km/s'],
                      'z0'   : [None, 'Disk scaleheight in arcsec'],
                      'inc'  : [None, 'Inclination angle in degrees'],
                      'phi'  : [None, 'Position angle of the receding part of the major axis (N->W)']}

    def __del__(self):
        if self._mod: libBB.Galfit_delete(self._mod)
        
    
    def _input(self,radii=None,xpos=None,ypos=None,vsys=None,vrot=None,vdisp=8,vrad=0,z0=0,inc=None,phi=None):
        """ Initialize rings for the fit 
        
        Set initial guesses for the fit and fixed parameters. Parameters not set are estimated.
        
        Args:
          radii (list):        Radii of the model in arcsec
          xpos (float, list):  X center of the galaxy in pixels
          ypos (float, list):  Y center of the galaxy in pixels
          vsys (float, list):  Systemic velocity of the galaxy in km/s
          vrot (float, list):  Rotation velocity in km/s
          vdisp (float, list): Velocity dispersion in km/s
          vrad (float, list):  Radial velocity in km/s
          z0 (float, list):    Disk scaleheight in arcsec
          inc (float, list):   Inclination angle in degrees
          phi (float, list):   Position angle of the receding part of the major axis (N->W)
    
        """
        # Check if any parameter needs to be estimated
        allpars = [radii,xpos,ypos,vsys,vrot,inc,phi]
        toEstimate = [True if (not isIterable(a) and a is None) else False for a in allpars]
                
        if True in toEstimate:
            print ("\n%s: estimating initial parameters for fit %s\n"%(self.taskname,"-"*25))
            sys.stdout.flush()
            xpos_s = '%s'%np.mean(xpos) if xpos else '-1'
            ypos_s = '%s'%np.mean(ypos) if ypos else '-1'
            inc_s  = '%s'%np.mean(inc)  if inc  else '-1'
            phi_s  = '%s'%np.mean(phi)  if phi  else '-1'
            guess  = libBB.Galfit_initialGuesses(self.inp._cube,xpos_s.encode('utf-8'),ypos_s.encode('utf-8'),inc_s.encode('utf-8'),phi_s.encode('utf-8'))
            guess  = reshapePointer(guess,[8])

            print (" Estimated parameters:")
            keys = ['radii','xpos','ypos','vsys','vrot','inc','phi']
            unit = ['arcsec','pix','pix','km/s','km/s','deg','deg']
            
            for i in range (len(toEstimate)):
                if toEstimate[i]:
                    if i==0: allpars[i] = np.linspace(guess[i+1]/2.,guess[i+1]*(guess[i]+0.5),guess[i])
                    else:    allpars[i] = guess[i+1]
                    print ("\n  %6s ="%keys[i], allpars[i], "%s"%unit[i])
                    
            radii,xpos,ypos,vsys,vrot,inc,phi = allpars
            print ("%s\n"%("-"*70))


        # Now we can initialize rings
        if not isIterable(radii): raise ValueError("radii must be an array")
        self._inri = Rings(len(radii))
        self._inri.set_rings(radii,xpos,ypos,vsys,vrot,vdisp,vrad,0.,0.,0.,1.E20,z0,inc,phi)


    def _compute(self,threads=1):
        """ Fit the model.

        Run this function after having set initial parameters with :func:`init` and options with
        :func:`set_options`. If initial guesses are not defined, the function will be run in 
        automated mode, e.g., it will estimate initial parameters and then run the fit.

        Returns:
          ndarray: An array fo size 2. Index 0 is a n x m matrix where n = number of rings, 
          m = number of parameters. Index 1 is an astropy PrimaryHDU object datacube 
          with the model.

        """
        if self._inri is None: 
            print ("BBarolo is running in automated mode. Check initial parameter estimate!")
            self._mod = libBB.Galfit_new(self.inp._cube)
        else:
            self._check_options()
            self._par.add_params(threads=threads,kwargs=self._opts)
            self._par.make_object()
            self._mod = libBB.Galfit_new_par(self.inp._cube,self._inri._rings,self._par._params)
            
        
        # Calculating the model
        self.modCalculated = libBB.Galfit_galfit(self._mod)
        if (self._opts['twostage'][0]): libBB.Galfit_secondStage(self._mod);
        
        # Write models
        libBB.Galfit_writeModel(self._mod,self._opts['norm'][0].encode('utf-8'),False)
        
        # Loading final rings
        try: self.bfit = np.genfromtxt(self._opts['outfolder'][0]+"/rings_final2.txt")
        except: self.bfit = np.genfromtxt(self._opts['outfolder'][0]+"/rings_final1.txt")
        
        # Loading final model
        for fname in os.listdir(self._opts['outfolder'][0]):
            if "mod_%s.fits"%(self._opts['norm'][0].lower()) in fname:
                filefits = self._opts['outfolder'][0]+"/"+fname
        self.outmodel = fits.open(filefits)[0]
        
        return (self.bfit, self.outmodel)
        
        
    def plot_model(self):
        """ Plot the model using BBarolo built-in output functions. """
        if not self.modCalculated:
            raise ValueError("3DFIT ERROR: you need to run fit() before plotting the model.")
        print (" Writing creative plots... ",end="")
        sys.stdout.flush()
        ret = libBB.Galfit_plotModel(self._mod) 
        if ret==0: print ("Done.")
        else: print(" Something went wrong! Check the plotscripts/ directory in your output directory.");
         
         
    def _check_options(self):
        """ Check if current options are ok for the C++ code """
        op = self._opts
        for key in op:
            op[key][0] = op[key][1](op[key][0])
            if key=='ltype' and op[key][0] not in range(1,6):
                raise ValueError("%s ERROR: %s can only be 1, 2, 3, 4 or 5."%(self.taskname,key))
            if key=='cdens' and op[key][0]<=0:
                raise ValueError("%s ERROR: %s can only be > 0."%(self.taskname,key))
            if key=='ftype' and op[key][0] not in range(1,4):
                raise ValueError("%s ERROR: %s can only be 1, 2 or 3."%(self.taskname,key))
            if key=='startrad' and op[key][0]<0:
                raise ValueError("%s ERROR: %s can only be positive."%(self.taskname,key))
            if key=='mask' and op[key][0] not in ['SMOOTH','SEARCH','THRESHOLD','NEGATIVE','SMOOTH&SEARCH','NONE']:
                if "FILE" not in op[key][0]:
                    raise ValueError("%s ERROR: %s can only 'SMOOTH','SEARCH','THRESHOLD','NEGATIVE', 'SMOOTH&SEARCH', 'NONE' or 'FILE(mask.fits)'."%(self.taskname,key))
            if key=='norm' and op[key][0].upper() not in ['LOCAL','AZIM','NONE']:
                raise ValueError("%s ERROR: %s can only 'LOCAL','AZIM' or 'NONE'."%(self.taskname,key))
            if key=='side' and op[key][0].upper() not in ['A','R','B']:
                raise ValueError("%s ERROR: %s can only 'A' (approaching),'R' (receding) or 'B' (both)."%(self.taskname,key))
            if key=='wfunc' and (op[key][0]<0 or op[key][0]>2):
                raise ValueError("%s ERROR: %s can only be 0, 1 or 2."%(self.taskname,key))

    def __str__ (self):
        self.show_options()
        return ""



class Search(Task):
    """3D Source finder 
    
    Args:
      fitsname (str): FITS file to search for sources
    
    """
    def __init__(self,fitsname):
        super(Search,self).__init__(fitsname=fitsname)
        self.taskname = 'SEARCH'
        self._opts = { "searchtype"  : ["spatial", str, "spectral or spatial search"],
                       "snrcut"      : [5, np.float32, "S/N cut for detection when sigma-clipping"],
                       "threshold"   : [0, np.float32, "Flux threshold for a detection"],
                       "adjacent"    : [True, bool, "Use the adjacent criterion for objects merger?" ],
                       "thrspatial"  : [-1, np.int32, "Maximum spatial separation between objects"],
                       "thrvelocity" : [2,  np.int32, "Maximum channels separation between objects"],
                       "minchannels" : [2,  np.int32, "Minimum channels to make an object"],
                       "minpixels"   : [-1, np.int32, "Minimum pixels required in an object"],
                       "minvoxels"   : [-1, np.int32, "Minimum voxels required in an object"],
                       "maxchannels" : [-1, np.int32, "Maximum channels to accept an object"],
                       "maxangsize"  : [-1, np.float32, "Maximum angular size in an object (arcmin)"],
                       "growth"      : [True, bool, "Growing objects once they are found?"],
                       "growthcut"   : [3, np.float32, "The SNR that we are growing objects down to"],
                       "growththresh": [0, np.float32, "The threshold for growing objects down to"],
                       "rejectbeforemerge" : [True, bool, "Whether to reject sources before merging"],
                       "twostagemerge" : [True, bool, "Whether to do a partial merge during search"]}


    def search(self,threads=1):
        return self._compute(threads)


    def _compute(self,threads=1):
        """ Perform the source finding """
        self._check_options()
        op = self._opts
        libBB.Search_search(self.inp._cube,op['searchtype'][0].encode('utf-8'),op['snrcut'][0],\
                            op['threshold'][0],op['adjacent'][0],op['thrspatial'][0],\
                            op['thrvelocity'][0],op['minpixels'][0],op['minchannels'][0],\
                            op['minvoxels'][0],op['maxchannels'][0],op['maxangsize'][0],\
                            op['growth'][0],op['growthcut'][0], op['growththresh'][0],\
                            op['rejectbeforemerge'][0],op['twostagemerge'][0],int(threads))
    
    
    def _check_options(self):
        """ Check if current options are ok for the C++ code """
        op = self._opts
        for key in op:
            op[key][0] = op[key][1](op[key][0])
            if key=='searchtype' and op[key][0] not in ['spatial','spectral']:
                raise ValueError("%s ERROR: %s can only be 'spatial' or 'spectral'"%(self.taskname,key))
            if key=='snrcut' and op[key][0]<=0:
                raise ValueError("%s ERROR: %s can only be > 0."%(self.taskname,key))
            if key=='threshold' and op[key][0]<0:
                raise ValueError("%s ERROR: %s can only be >= 0."%(self.taskname,key))
            if key=="growthcut" and op[key][0]>op['snrcut'][0]:
                raise ValueError("%s ERROR: %s can only be < snrcut"%(self.taskname,key))
            if key=="growththresh" and op[key][0]>op['threshold'][0]:
                raise ValueError("%s ERROR: %s can only be < threshold"%(self.taskname,key))



class FitMod2D(Task):
    """Classic 2D tilted-ring model
    
    Args:
      fitsname (str): FITS file to fit with a 2D tilted-ring model
    """
    def __init__(self,fitsname):
        super(FitMod2D,self).__init__(fitsname=fitsname)
        self.taskname = '2DFIT'
        # A pointer to the C++ model
        self._mod = None
        # Input rings
        self._inri = None
        # Output rings as a dictionary
        self.outrings = None
        # Options and arguments
        self._opts = { 'free'    : ['VROT', str, "Free parameters"],
                       'side'    : ['B', str, "Which side of the galaxy to fit"],
                       'wfunc'   : [2, np.int32, "Weighting function for major axis"],
                       'mask'    : ['SMOOTH', str, "Mask type"]}
        self._args = {'radii': [None, 'Radii of the model in arcsec (must be an array)'],
                      'xpos' : [None, 'X-center of the galaxy in pixels'],
                      'ypos' : [None, 'Y center of the galaxy in pixels'],
                      'vsys' : [None, 'Systemic velocity of the galaxy in km/s'],
                      'vrot' : [None, 'Rotation velocity in km/s'],
                      'vrad' : [None, 'Radial velocity in km/s'],
                      'inc'  : [None, 'Inclination angle in degrees'],
                      'phi'  : [None, 'Position angle of the receding part of the major axis (N->W)']}
        
    
    def __del__(self):
        if self._mod: libBB.Fit2D_delete(self._mod)
        
    
    def _input(self,radii,xpos,ypos,vsys,vrot,inc,phi,vrad=0):
        """ Initialize rings for the fit 
                
        Args:
          radii (list):        Radii of the model in arcsec
          xpos (float, list):  X center of the galaxy in pixels
          ypos (float, list):  Y center of the galaxy in pixels
          vsys (float, list):  Systemic velocity of the galaxy in km/s
          vrot (float, list):  Rotation velocity in km/s
          vrad (float, list):  Radial velocity in km/s
          inc (float, list):   Inclination angle in degrees
          phi (float, list):   Position angle of the receding part of the major axis (N->W)
    
        """
        if not isIterable(radii): raise ValueError("radii must be an array")
        self._inri = Rings(len(radii))
        self._inri.set_rings(radii,xpos,ypos,vsys,vrot,0.,vrad,0.,0.,0.,1.E20,0.,inc,phi)

    
    def _compute(self,threads=1):
        """ Fit the model.

        Run this function after having set initial parameters with :func:`init` and options with
        :func:`set_options`.

        Returns:
        a dictionary containing all the best fit rings
        """
        if self._inri is None: 
            raise ValueError("%s ERROR: you need to set the model with init(...) before calling compute()."%self.taskname)
        self.__del__()
        self._check_options()
        op = self._opts
        self._mod = libBB.Fit2D_new(self.inp._cube,self._inri._rings,op['mask'][0].encode('utf-8'),op['free'][0].encode('utf-8'),\
                                    op['side'][0].encode('utf-8'),op['wfunc'][0],int(threads))
        libBB.Fit2D_compute(self._mod)
        libBB.Fit2D_write(self._mod,self.inp._cube,'rings.txt'.encode('utf-8'))
        a = np.genfromtxt('rings.txt') 
        try: os.remove('rings.txt')
        except: pass
        self.outrings = {'rad':a[:,1],'vsys':a[:,2],'vrot':a[:,3],'vrad':a[:,4],\
                         'phi':a[:,5],'inc':a[:,6],'xpos':a[:,7],'ypos':a[:,8]}
        return self.outrings
        
    
    def _check_options(self):
        """ Check if current options are ok for the C++ code """
        op = self._opts
        for key in op:
            op[key][0] = op[key][1](op[key][0])
            if key=='mask' and op[key][0] not in ['SMOOTH','SEARCH','THRESHOLD','NEGATIVE','SMOOTH&SEARCH','NONE']:
                if "FILE" not in op[key][0]:
                    raise ValueError("%s ERROR: %s can only 'SMOOTH','SEARCH','THRESHOLD','NEGATIVE', 'SMOOTH&SEARCH', 'NONE' or 'FILE(mask.fits)'."%(self.taskname,key))
            if key=='side' and op[key][0].upper() not in ['A','R','B']:
                raise ValueError("%s ERROR: %s can only be 'A' (approaching),'R' (receding) or 'B' (both)."%(self.taskname,key))
            if key=='wfunc' and (op[key][0]<0 or op[key][0]>2):
                raise ValueError("%s ERROR: %s can only be 0, 1 or 2."%(self.taskname,key))
    
    
    def writeto(self,filename):
        """ Write the best-fit ring model to a text file """
        if self._inri is None or self._mod is None: 
            raise ValueError("%s ERROR: you need to compute the model with compute (...) before calling writeto()."%self.taskname)
        libBB.Fit2D_write(self._mod,self.inp._cube,filename.encode('utf-8'))



class Ellprof(Task):
    """Find statistics over rings
    
    Args:
      fitsname (str): FITS file to fit with a 2D tilted-ring model
    """
    def __init__(self,fitsname):
        super(Ellprof,self).__init__(fitsname=fitsname)
        self.taskname = 'ELLPROF'
        # Pointer to the C++ Ellprof class
        self._mod = None
        # Input rings
        self._inri = None
        # Output rings as a dictionary
        self.outrings = None
        # Options and arguments
        self._opts = {'side'    : ['B', str, "Which side of the galaxy to fit"],
                      'mask'    : ['SMOOTH', str, "Mask type"]}
        self._args = {'radii': [None, 'Radii of the model in arcsec (must be an array)'],
                      'xpos' : [None, 'X-center of the galaxy in pixels'],
                      'ypos' : [None, 'Y center of the galaxy in pixels'],
                      'inc'  : [None, 'Inclination angle in degrees'],
                      'phi'  : [None, 'Position angle of the receding part of the major axis (N->W)']}
        
    
    def __del__(self):
        if self._mod: libBB.Ellprof_delete(self._mod)
        
    
    def _input(self,radii,xpos,ypos,inc,phi):
        """ Initialize rings 
        
        Args:
          radii (list):        Radii of the model in arcsec
          xpos (float, list):  X center of the galaxy in pixels
          ypos (float, list):  Y center of the galaxy in pixels
          vsys (float, list):  Systemic velocity of the galaxy in km/s
          inc (float, list):   Inclination angle in degrees
          phi (float, list):   Position angle of the receding part of the major axis (N->W)
            
        """
        if not isIterable(radii): raise ValueError("radii must be an array")
        self._inri = Rings(len(radii))
        self._inri.set_rings(radii,xpos,ypos,0.,0.,0.,0.,0.,0.,0.,1.E20,0.,inc,phi)
        

    def _compute(self,threads=1):
        """ Extract profiles along the rings.

        Run this function after having set initial parameters with :func:`init` and options with
        :func:`set_options`.

        Returns:
        a dictionary containing all the statistics along the rings
        """
        if self._inri is None: 
            raise ValueError("%s ERROR: you need to set the model with init(...) before calling compute()."%self.taskname)
        self.__del__()
        self._check_options()
        op = self._opts
        self._mod = libBB.Ellprof_new(self.inp._cube,self._inri._rings,op['mask'][0].encode('utf-8'),op['side'][0].encode('utf-8'),int(threads))
        libBB.Ellprof_compute(self._mod)
        libBB.Ellprof_write(self._mod,'rings.txt'.encode('utf-8'))
        a = np.genfromtxt('rings.txt') 
        try: os.remove('rings.txt')
        except: pass
        self.outrings = {'rad':a[:,0],'sum':a[:,1],'mean':a[:,2],'median':a[:,3],'rms':a[:,4],'mad':a[:,5],'numpix':a[:,6],\
                         'surfdens':a[:,7],'surfdens_err':a[:,8],'surfdens_fo':a[:,9]}
        if (a.shape[1]>10): self.outrings.update({'msurfdens':a[:,10],'msurfdens':a[:,11]})
        
        return self.outrings
        
    
    def _check_options(self):
        """ Check if current options are ok for the C++ code """
        op = self._opts
        for key in op:
            op[key][0] = op[key][1](op[key][0])
            if key=='mask' and op[key][0] not in ['SMOOTH','SEARCH','THRESHOLD','NEGATIVE','SMOOTH&SEARCH','NONE']:
                if "FILE" not in op[key][0]:
                    raise ValueError("%s ERROR: %s can only 'SMOOTH','SEARCH','THRESHOLD','NEGATIVE', 'SMOOTH&SEARCH', 'NONE' or 'FILE(mask.fits)'."%(self.taskname,key))
            if key=='side' and op[key][0].upper() not in ['A','R','B']:
                raise ValueError("%s ERROR: %s can only be 'A' (approaching),'R' (receding) or 'B' (both)."%(self.taskname,key))

    
    def writeto(self,filename):
        """ Write the profiles to a text file """
        if self._inri is None or self._mod is None: 
            raise ValueError("%s ERROR: you need to compute the profiles with compute() before calling writeto()."%self.taskname)
        libBB.Ellprof_write(self._mod,filename.encode('utf-8'))



class SpectralSmoothing(Task):
    """ Spectral smoothing of a datacube with a filtering window.
    
    Args:
      fitsname (str): FITS file to fit with a 2D tilted-ring model
    """
    def __init__(self,fitsname):
        super(SpectralSmoothing,self).__init__(fitsname=fitsname)
        self.taskname = 'SPECTRALSMOOTHING'
        # Pointer to the C++ Hanning class
        self._ptr = None  
    
    
    def __del__(self):
        if self._ptr: libBB.SpectralSmooth3D_delete(self._ptr)
        

    def smooth(self,window_type='hanning',window_size=3,threads=1):
        """ Compute an smoothed 3D array 

        Args:
          window_type (str): type of filtering window (hanning,bartlett,welch,flattop,boxcar,blackman)
          window_size (int): width of the window (it has to be an odd number)
          threads (int):     number of cpus to use for computation.
        
        Returns:
        a 3D spectrally-smoothed array 
        """
        known_windows = ['HANNING','HANNING2','BOXCAR','BARTLETT','WELCH','BLACKMAN','FLATTOP']
        self.__del__()
        if window_size % 2 == 0:
            raise ValueError("%s ERROR: window_size must be an odd number."%self.taskname)
        if window_type.upper() not in known_windows:
            raise ValueError("%s ERROR: window_type could only be HANNING, HANNING2, \
                             BARTLETT, WELCH, BLACKMAN or FLATTOP"%(self.taskname))
            
        # Initializing SmoothSpectral smooth class
        self._ptr = libBB.SpectralSmooth3D_new(window_type.encode('utf-8'),int(window_size))
        # Performing Hanning smoothing class
        libBB.SpectralSmooth3D_compute(self._ptr,self.inp._cube,int(threads))
        return reshapePointer(libBB.SpectralSmooth3D_array(self._ptr),self.inp.dim[::-1])
        
        
    def writeto(self,filename="",average=False):
        """ Write the spectrally-smoothed 3D array to a FITS file 
        
        Args:
          filename (str): name of output FITS file
          average (bool): if true, average over window_size channels
        """
        if self._ptr is None: 
            raise ValueError("%s ERROR: you need to Hanning smooth the data with compute() before calling writeto()."%self.taskname)
        libBB.SpectralSmooth3D_write(self._ptr,self.inp._cube,filename.encode('utf-8'),average)