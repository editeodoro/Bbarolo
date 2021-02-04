import os, io, subprocess
import warnings as warn
from distutils.spawn import find_executable

class BBaroloWrapper(object):
    """ 
    A class to directly call BBarolo given some parameters
    
    Attributes
    ----------
    opts: (str)
      A list of strings (param=value) with stored BBarolo's parameters.
    
    
    Methods
    -------
    add_options(**kwargs):
      Add new options to the opts attribute
    
    run(exe=None,stdout=None):
      Run BBarolo after checking that the executable exists and stdout
    
    run_nochecks(exe='BBarolo',stdout=subprocess.DEVNULL):
      Run BBarolo without any checks
    
    write_parameterfile(fileout='param.par'):
      Write a parameter file with the parameters stored in the class
    
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
        
        # Creating a list of BBarolo's parameters
        self.opts = []
        if params:
            # Parameters are given through params
            if isinstance(params,dict):
                self.opts = [f'{k}={params[k]}' for k in params]
            elif isinstance(params,(list,tuple)):
                if not all(isinstance(k, str) for k in params):
                    raise TypeError("All elements of parameter list must be strings.")
                self.opts = [str(k.replace(" = ", "=")) for k in params]
            else: 
                raise TypeError("Parameters must given as either a list of strings or a dictionary.")
        
        self.add_options(**kwargs)


    def add_options(self,**kwargs):
        " Add new options to the parameter list"
        if kwargs:
            # Parameters are added through **kwargs 
            self.opts.extend([f'{k}={kwargs[k]}' for k in kwargs])


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
        
        # Try to find BBarolo executable exists somewhere
        exe = find_executable('BBarolo') if exe is None else exe
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


    def run_nochecks(self,exe='BBarolo',stdout=subprocess.DEVNULL):
        """ Run BBarolo with no further checks """
        cmd = [f'{exe}', '-c']
        cmd.extend(self.opts)
        return subprocess.call(cmd,stdout=stdout,stderr=stdout)
        

    def write_parameterfile(self,fileout='param.par'):
        """ Write a BBarolo's parameter file in fileout """
        with open(fileout,'w') as f:
            f.write(self.__str__())
    
    
    def __str__(self):
        s = "##### Input parameters for BBarolo #####\n"
        for ss in self.opts:
            a = ss.split("=")
            s += f'{a[0].upper():18s} {a[1]} \n' 
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
