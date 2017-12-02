import os,sys
import numpy as np
import multiprocessing as mpr
from .BB_interface import libBB
from ctypes import *

def reshapePointer (p, shape):
    """ 
    This function take a POINTER to c_float and 
    reshape according to shape parameter
    \param p:       a POINTER(c_float) object
    \param shape:   tuple or list of axes sizes
    \return:        Reshaped Numpy array 
    """
    return np.ctypeslib.as_array(p, shape=tuple(shape))


class FitsCube(object):
    def __init__(self,fitsname):
        # File name of the FITS file to read
        self.fname = fitsname.encode('utf-8')
        # Pointer to the C++ Cube object
        self._cube = None
        if not os.path.isfile(fitsname):
            raise ValueError("FitsCube ERROR: file %s does not exist."%self.fname)
        self._cube  = libBB.Cube_new(self.fname)
        # Axis dimensions
        self.dim   = libBB.Cube_axisdim(self._cube)
        # Getting array data and reshaping in (z,y,x)
        self.data  = reshapePointer(libBB.Cube_array(self._cube),self.dim[::-1])
        
    def __del__(self):
        if self._cube: libBB.Cube_delete(self._cube)
        

class Rings(object):
    """ Wrapper class for C++ Rings struct (galmod.hh) """
    
    def __init__(self,nrings):
        # Number of rings
        self.nr     = nrings
        # Pointer to the C++ Rings object
        self._rings = libBB.Rings_new();
        self.rinDef = False
    
    def set_rings (self,radii,xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens,z0,inc,phi,nv):
        """ Define rings given the input parameters 
            
            \param radii: List or array
            \param other: float or array
        """
        
        typs, fl = (list,tuple,np.ndarray), np.float32
        
        if not isinstance(self.nr,(int,float)): raise ValueError("nrings must be a number")
            
        if isinstance(radii,typs): radii = np.array(radii,dtype=fl)
        else: raise ValueError("radii must be an array")
        
        xpos  = np.array(xpos,dtype=fl)  if isinstance(xpos,typs)  else np.full(self.nr,xpos,dtype=fl)
        ypos  = np.array(ypos,dtype=fl)  if isinstance(ypos,typs)  else np.full(self.nr,ypos,dtype=fl)
        vsys  = np.array(vsys,dtype=fl)  if isinstance(vsys,typs)  else np.full(self.nr,vsys,dtype=fl)
        vrot  = np.array(vrot,dtype=fl)  if isinstance(vrot,typs)  else np.full(self.nr,vrot,dtype=fl)
        vdisp = np.array(vdisp,dtype=fl) if isinstance(vdisp,typs) else np.full(self.nr,vdisp,dtype=fl)
        vrad  = np.array(vrad,dtype=fl)  if isinstance(vrad,typs)  else np.full(self.nr,vrad,dtype=fl)
        vvert = np.array(vvert,dtype=fl) if isinstance(vvert,typs) else np.full(self.nr,vvert,dtype=fl)
        dvdz  = np.array(dvdz,dtype=fl)  if isinstance(dvdz,typs)  else np.full(self.nr,dvdz,dtype=fl)
        zcyl  = np.array(zcyl,dtype=fl)  if isinstance(zcyl,typs)  else np.full(self.nr,zcyl,dtype=fl)
        dens  = np.array(dens,dtype=fl)  if isinstance(dens,typs)  else np.full(self.nr,dens,dtype=fl)
        z0    = np.array(z0,dtype=fl)    if isinstance(z0,typs)    else np.full(self.nr,z0,dtype=fl)
        inc   = np.array(inc,dtype=fl)   if isinstance(inc,typs)   else np.full(self.nr,inc,dtype=fl)
        phi   = np.array(phi,dtype=fl)   if isinstance(phi,typs)   else np.full(self.nr,phi,dtype=fl)
        nv    = np.array(nv,dtype=np.int) if isinstance(nv,typs)   else np.full(self.nr,nv,dtype=np.int)
        
        allr = (radii,xpos,ypos,vsys,vdisp,vrad,vvert,dvdz,zcyl,dens,z0,inc,phi,nv)
        
        for i in allr:
            if len(i)!=self.nr: raise ValueError("All quantities must have size = %i"%self.nr)
        
        libBB.Rings_set(self._rings,self.nr,radii,xpos,ypos,vsys,vrot,
                        vdisp,vrad,vvert, dvdz,zcyl,dens,z0,inc,phi,nv)
        
        self.rinDef = True
        
        
class GalWind(object):
    """ Wrapper class for C++ GalWind class (galwind.hh) """
    def __init__(self,fitsname):
        # Input datacube 
        self.inp = FitsCube(fitsname)
        # The C++ wind model
        self._mod = None
        # The output array model
        self.data = None
        
        
    def __del__(self):
        if self._mod: libBB.Galwind_delete(self._mod)   
            
            
    def input(self,x0,y0,pa,inc,disp,dens,vsys,vw,openang,htot,dtype=1,ntot=25,cdens=10,nv=10):
        """ Set the parameters of the model """
        if self._mod: self.__del__()
        self._mod = libBB.Galwind_new(self.inp._cube,x0,y0,pa,\
        inc,disp,dens,vsys,vw,openang,htot,dtype,ntot,cdens,nv)
    

    def compute(self):
        """ Compute the model """
        if self._mod:
            libBB.Galwind_compute(self._mod)
            self.data  = reshapePointer(libBB.Galwind_array(self._mod),self.inp.dim[::-1])
        else:
            raise ValueError("GALWIND ERROR: you need to set the model with input(...) before calling compute().")
        return self.data
        

    def smooth(self):
        """ Smooth the wind model and return the smoothed array """
        if self.data is not None:
            libBB.Galwind_smooth(self._mod)
            self.data  = reshapePointer(libBB.Galwind_array(self._mod),self.inp.dim[::-1])
        else:
            raise ValueError("GALWIND ERROR: you need to compute the model before calling smooth().")
        return self.data
    


class FitMod3D(object):
    
    def __init__(self,fitsname):
        # Input datacube 
        self.inp = FitsCube(fitsname)
        # A Rings object with initial conditions
        self._inri = None
        # The C++ Galfit object
        self._mod = None
        self._modCalculated = False
        # The output array model
        self.data = None
        # A dictionary with options and defaults values
        self._opts = {'CDENS'   : [10, np.int],
                      'NV'      : [-1, np.int],
                      'LTYPE'   : [2, np.int],
                      'SM'      : [True, np.bool],
                      'DELTAINC': [5., np.float32],
                      'DELTAPHI': [15., np.float32],
                      'FTYPE'   : [2, np.int],
                      'WFUNC'   : [2, np.int],
                      'BWEIGHT' : [1, np.int],
                      'TOL'     : [1E-03, np.float64],
                      'MASK'    : ['SMOOTH', str],
                      'NORM'    : ['LOCAL', str],
                      'FREE'    : ['VROT VDISP INC PA', str],
                      'SIDE'    : ['B', str],
                      'TWOSTAGE': [True, np.bool],
                      'POLYN'   : ['bezier', str],
                      'STARTRAD': [0, np.int],
                      'ERRORS'  : [False, np.bool],
                      'DISTANCE': [-1., np.float32],
                      'REDSHIFT': [-1., np.float64],
                      'RESTWAVE': [-1., np.float64],
                      'OUTFOLDER' : ['./', str]}
    
    
    def __del__(self):
        if self._mod: libBB.Galfit_delete(self._mod)
    
        
    def init(self,radii,xpos,ypos,vsys,vrot,vdisp,vrad,z0,inc,phi):
        """ Initialize rings for the fit (initial guesses) """
        if not isinstance(radii,(list,tuple,np.ndarray)): raise ValueError("radii must be an array")
        self._inri = Rings(len(radii))
        self._inri.set_rings(radii,xpos,ypos,vsys,vrot,vdisp,vrad,0,0,0,1.E20,z0,inc,phi,0)
      
      
    def show_options(self):
        print ("\nCurrent options for 3DFIT task:")
        for key, value in list(self._opts.items()):
            print (" %-8s = %s "%(key,value[0]))
    
    
    def set_options(self,**kwargs):
        """ Set options for the task. Keywords **kwargs are given in self.options """
        for key in kwargs:
            if key in self._opts:
                self._opts[key][0] = self._opts[key][1](kwargs[key]) 
            else:
                raise ValueError('Option %s unknown. Try show_options() for a list of keywords'%key)
        
        self._check_options()
        
        
    
    def fit(self):
        if self._inri is None: 
            print ("BBarolo is running in automated mode. Check initial parameter estimate!")
            self._mod = libBB.Galfit_new(self.inp._cube)
        else:
            self._check_options()
            op = self._opts
            self._mod = libBB.Galfit_new_all(self.inp._cube,self._inri._rings,op['DELTAINC'][0],op['DELTAPHI'][0],\
                                             op['LTYPE'][0],op['FTYPE'][0],op['WFUNC'][0],op['BWEIGHT'][0],op['NV'][0],\
                                             op['TOL'][0],op['CDENS'][0],op['STARTRAD'][0],op['MASK'][0].encode('utf-8'),\
                                             op['NORM'][0].encode('utf-8'),op['FREE'][0].encode('utf-8'),\
                                             op['SIDE'][0].encode('utf-8'),op['TWOSTAGE'][0],op['POLYN'][0].encode('utf-8'),\
                                             op['ERRORS'][0],op['SM'][0],op['DISTANCE'][0],op['REDSHIFT'][0],op['RESTWAVE'][0],\
                                             op['OUTFOLDER'][0].encode('utf-8'))
                                 
        self.modCalculated = libBB.Galfit_galfit(self._mod)        
        if (op['TWOSTAGE'][0]): libBB.Galfit_secondStage(self._mod);
        print (self.modCalculated)
        
    def write_model(self):
        if not self.modCalculated:
            raise ValueError("3DFIT ERROR: you need to run fit() before plotting the model.")
        libBB.Galfit_writeModel(self._mod,self._opts['NORM'][0].encode('utf-8'))
        
         
         
    def _check_options(self):
        """ Check if current options are ok for the C++ code """
        op = self._opts
        for key in op:
            op[key][0] = op[key][1](op[key][0])
            if key=='LTYPE' and op[key][0] not in range(1,6):
                raise ValueError("3DFIT ERROR: LTYPE can only be 1, 2, 3, 4 or 5.")
            if key=='FTYPE' and op[key][0] not in range(1,4):
                raise ValueError("3DFIT ERROR: FTYPE can only be 1, 2 or 3.")
            if key=='STARTRAD' and op[key][0]<0:
                raise ValueError("3DFIT ERROR: STARTRAD can only be positive.")
            if key=='MASK' and op[key][0] not in ['SMOOTH','SEARCH','THRESHOLD','NEGATIVE','NONE']:
                raise ValueError("3DFIT ERROR: MASK can only 'SMOOTH','SEARCH','THRESHOLD','NEGATIVE' or 'NONE'.")
            if key=='NORM' and op[key][0] not in ['LOCAL','AZIM','NONE']:
                raise ValueError("3DFIT ERROR: NORM can only 'LOCAL','AZIM' or 'NONE'.")

    def __str__ (self):
        self.show_options()
        return ""
            



