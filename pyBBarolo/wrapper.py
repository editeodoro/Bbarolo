"""
This module defines a few classes to call BBarolo C++ executable 
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

import os, io, stat, subprocess
import warnings as warn
from distutils.spawn import find_executable

# Global variable pointing to BBarolo executable
BBexe = None

# BBarolo should be in module directory. We make sure it is executable
moddir = os.path.dirname(os.path.realpath(__file__))
if os.path.isfile(moddir+'/BBarolo'): 
    BBexe = moddir+'/BBarolo'
    os.chmod(BBexe,  stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


class BBaroloWrapper(object):
    """ 
    A class to directly call BBarolo given some parameters
    
    Attributes
    ----------
    opts: (dict)
      A dictionary {param : value} with stored BBarolo's parameters.
    
    
    Methods
    -------
    add_options(**kwargs):
      Add new options to the opts attribute
    
    remove_option(toremove):
      Remove a previously added option from opts list
    
    run(exe=None,stdout=None):
      Run BBarolo after checking that the executable exists and stdout
    
    run_nochecks(exe='BBarolo',stdout=subprocess.DEVNULL):
      Run BBarolo without any checks
    
    write_parameterfile(fileout='param.par'):
      Write a parameter file with the parameters stored in the class
    
    reset():
      Reset BBarolo's options to an empty list
    
    """
    def __init__(self,params=None,**kwargs):
        """ Initialize a BBaroloWrapper instance.
            BBarolo's parameters can be given in three different ways:
            
            1) A list of strings (param=value) through the params variable
            2) A dictionary of (param,value) through the params variable
            3) Parameter=value set through the **kwargs
            
            Examples:
            1) BBaroloWrapper(["FITSFILE=mycube","FIT3D=true"])
            2) BBaroloWrapper(dict=(fitsfile='mycube',fit3d=True))
            3) BBaroloWrapper(fitsfile='mycube',fit3d=True)
            
        Args:
          params (list,dict):   A list of string or a dictionary with BB's parameters
          kwargs ():            Additional (or all) BB's parameters and values

        """
        
        if params is None and not kwargs:
            raise ValueError("BBaroloWrapper must be initialised with a list of BBarolo's parameters")
        
        self.opts = {}         # A dictionary with BB's parameters
        if params:
            # Parameters are given through params
            if isinstance(params,dict):
                self.opts = params
            elif isinstance(params,(list,tuple)):
                if not all(isinstance(k, str) for k in params):
                    raise TypeError("All elements of parameter list must be strings.")
                for s in params:
                    a = s.split("=")
                    self.opts.update({a[0].lower().strip() : a[1].strip()})
            else: 
                raise TypeError("Parameters must given as either a list of strings or a dictionary.")
        
        self.add_options(**kwargs)


    def add_options(self,**kwargs):
        """ Add new options to the parameter list """
        if kwargs: self.opts.update(kwargs)

    
    def remove_option(self,toremove):
        """ Remove an option from the parameter list """
        self.opts.pop(toremove, None)
    
    
    def run(self,exe=None,stdout=None):
        """ Run BBarolo with the given parameters after checking executable and stdout
    
        Args:
          exe (str):      Path to BBarolo executable. If None, search in $PATH
          stdout (str):   How to report BBarolo's messages. 
                          None=screen, str=file, 'null'=NULL
        """
        
        # Check if given BBarolo executable exists
        if exe is not None and not os.path.exists(exe):
            raise NameError(f"{exe} does not exist.")
        
        # If not given using global BBexe variable
        exe = BBexe if exe is None else exe
        if exe is None:
            # Try to find BBarolo executable exists in path
            exe = find_executable('BBarolo')
            if exe is None:
                raise NameError("Cannot find any BBarolo executable in $PATH. Please specify it with exe=yourBBarolo")
        
        # Deciding where to print BB's messages (screen,file or null)
        if isinstance(stdout,str): 
            if 'null' in stdout.lower(): stdout = subprocess.DEVNULL
            else: stdout = open(stdout,'w')

        # Running BBarolo
        if (self.run_nochecks(exe=exe,stdout=stdout)):
            warn.warn("BBarolo returned an unsuccessful exit code", RuntimeWarning)
        
        # Closing log file if open
        if isinstance(stdout,io.IOBase): stdout.close()


    def run_nochecks(self,exe=BBexe,stdout=subprocess.DEVNULL):
        """ Run BBarolo with no further checks """
        cmd = [f'{exe}', '-c']
        params = [f'{k}={self.opts[k]}' for k in self.opts]
        cmd.extend(params)
        return subprocess.call(cmd,stdout=stdout,stderr=stdout)

    
    def write_parameterfile(self,fileout='param.par'):
        """ Write a BBarolo's parameter file in fileout """
        with open(fileout,'w') as f:
            f.write(self.__str__())
    

    def reset(self):
        """ Reset options """
        self.opts.clear()


    def __getitem__(self,key):
       return self.opts[key]


    def __str__(self):
        s = "##### Input parameters for BBarolo #####\n"
        for k in self.opts:
            s += f'{k.upper():18s} {self.opts[k]} \n' 
        return s



class GenericTask(BBaroloWrapper):
    def __init__(self,task,fitsname,**kwargs):
        specopt = {'fitsfile' : fitsname, task : True }
        super(GenericTask,self).__init__(**specopt,**kwargs)

class GalMod(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(GalMod,self).__init__(task='galmod',fitsname=fitsname,**kwargs)

class GalWind(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(GalWind,self).__init__(task='galwind',fitsname=fitsname,**kwargs)

class FitMod3D(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(FitMod3D,self).__init__(task='fit3d',fitsname=fitsname,**kwargs)

class FitMod2D(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(FitMod2D,self).__init__(task='fit2d',fitsname=fitsname,**kwargs)
        
class Search(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(Search,self).__init__(task='search',fitsname=fitsname,**kwargs)
 
class EllProf(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(EllProf,self).__init__(task='ellprof',fitsname=fitsname,**kwargs)

class SpectralSmoothing(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(SpectralSmoothing,self).__init__(task='smoothspec',fitsname=fitsname,**kwargs)
     
class SpatialSmoothing(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(SpatialSmoothing,self).__init__(task='smooth',fitsname=fitsname,**kwargs)

class MakeMask(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(MakeMask,self).__init__(task='makemask',fitsname=fitsname,**kwargs)

class PVSlice(GenericTask):
    def __init__(self,fitsname,**kwargs):
        super(PVSlice,self).__init__(task='pvslice',fitsname=fitsname,**kwargs)
        

def read_parameter_file(filein):
    """ This function reads in a BBarolo's parameter file 
    
    Args:
      filein (str):   Parameter file name
    """
    if not os.path.exists(filein):
        raise FileNotFoundError(f"File {filein} does not exist")
        
    s = []
    with open(filein,'r') as fp:
       for line in fp:
           l = line.strip()
           if not l: continue
           if l[0:2] in ['#','//']: continue
           ssplit = l.split()
           s.append(f'{ssplit[0]}={" ".join(ssplit[1:])}')
           
    return s
