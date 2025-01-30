"""
This module allows to use pyBBarolo in a Bayesian fashion.
It uses either dynesty or nautilus libraries for this.
"""

########################################################################
# Copyright (C) 2024 Enrico Di Teodoro
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

import ctypes, time, copy
import numpy as np
import scipy.stats
from .BB_interface import libBB
from .pyBBarolo import Param, Rings, FitMod3D, reshapePointer, vprint, isIterable

parallel = True
try: 
    from schwimmbad import MultiPool, MPIPool
    from mpi4py import MPI
except ImportError:
    print ("WARNING! Parallelization is disabled. Please install schwimmbad and mpi4py to enable it.")
    parallel = False

try: 
    import dynesty as dyn
    from dynesty.utils import resample_equal
except ImportError:
    raise ImportError("BayesianBB requires the package 'dynesty'. Please install it.")


class BayesianBBarolo(FitMod3D):
    
    """ A class to fit a galaxy model to a 3D datacube using a Bayesian framework
    
    Attributes
    ----------
    taskname : str
      A string with the name of the task (BAYESIAN3DFIT)
    
    ...
      
    Methods
    -------
    set_options(**kwargs):
        Set all options for BBarolo.

    TO BE COMPLETED
    """
    
    def __init__(self,fitsname,**kwargs):
        """ Initialize the BayesianBBarolo class.
    
        Parameters
        ----------
        fitsname : str
            The name of the fits file with the datacube to fit.
        **kwargs : dict
            Any other parameter to be passed to the BBarolo's library.
        """

        super(BayesianBBarolo,self).__init__(fitsname=fitsname)
        # Task name
        self.taskname = "BAYESIAN3DFIT"
        # Resetting FitMod3d._args. Not used here.
        self._args = {}
        # self._opts contains any other parameter given to the code
        defaults = dict(threads=1,outfolder='./output/',twostage=False)
        self._opts = Param(fitsfile=fitsname,**defaults)
        if kwargs:
            self._opts.add_params(kwargs)
        
        # A dictionary with the names of parameters to be fitted and their indexes
        self.freepar_idx = None
        self.ndim = 0
        # A list with the "names" of fitted parameters
        self.freepar_names = None
        # A dictionary with the boundaries for the priors of parameters
        self.bounds = None
        # A dictionary for prior probabilities
        self.priors = None
        # A pointer to the C++ Galfit object
        self._galfit = None
        self.modCalculated = False

    
     
    def set_options(self,**kwargs):
        """ Add options to be passed to BBarolo """
        self._opts.add_params(**kwargs)
    
    
    def show_options(self):
        """ Print the options to be passed to BBarolo """
        print(self._opts)
    
    
    def init(self,radii,xpos,ypos,vsys,vrot,vdisp,inc,phi,z0=0,vrad=0,dens=1,vvert=0,dvdz=0,zcyl=0):
        """ Initialize rings for the fit. Parameters that are not fitted will be fixed.
        
        Parameters
        ----------
        radii : list
            Radii of the model in arcsec
        xpos : float, list
            X center of the galaxy in pixels
        ypos : float, list
            Y center of the galaxy in pixels
        vsys : float, list
            Systemic velocity of the galaxy in km/s
        vrot : float, list
            Rotation velocity in km/s
        vdisp : float, list
            Velocity dispersion in km/s
        inc : float, list
            Inclination angle in degrees
        phi : float, list
            Position angle of the receding part of the major axis (N->W)
        z0 : float, list (optional)
            Disk scaleheight in arcsec. Default is 0.
        vrad : float, list (optional)
            Radial velocity in km/s. Default is 0.
        dens : float, list (optional)
            Face-on column density in units of 1E20 cm^-2. Default is 1.
        vvert : float, list (optional)
            Vertical velocity in km/s. Default is 0.
        dvdz : float, list (optional)
            Vertical velocity gradient in km/s/pc. Default is 0.
        zcyl : float, list (optional)
            Cylindrical height where the dvdz starts in arcsec. Default is 0.            
        """
        # Using the Fit3D initialization function (see pyBBarolo.py)
        # super(BayesianBBarolo,self)._input(**kwargs)

        # Initialize rings
        if not isIterable(radii): raise ValueError("radii must be an array")
        self._inri = Rings(len(radii))
        self._inri.set_rings(radii,xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens*1E20,z0,inc,phi)

        # Defining default boundaries
        rr = self._inri.r
        self.bounds = dict(xpos=[rr['xpos'][0]-5,rr['xpos'][0]+5],     #+- 5 pixels
                           ypos=[rr['ypos'][0]-5,rr['ypos'][0]+5],     #+- 5 pixels
                           vsys=[rr['vsys'][0]-50,rr['vsys'][0]+50],   #+- 50 km/s
                           vrot=[0,350],vdisp=[0,50],vrad=[0,50],
                           inc=[0,90],phi=[0,360],
                           vvert=[0,30],dvdz=[0,2],zcyl=[0,5],
                           dens=[0.1,200],z0=[0,10])
        self.priors = {key: None for key in self.bounds.keys()}
    

    
    def _setup(self,freepar,useBBres=True,**kwargs):
        """ Setup function for the fit 
            It sets the following class attributes:
            - self.ndim: number of parameters to fit
            - self.freepar_idx: indexes of the parameters to fit
            - self.freepar_names: names of the parameters to fit
            - self.priors: dictionary with prior probabilities
            - self._galfit: pointer to the C++ Galfit object
            - self.mask: mask from the input cube (3D array)
            - self._ellprof: a pointer to a C++ Ellprof object   
        """


        if self._inri is None: 
            print ("Error: you need to initialize rings first by calling the function init()")
            return
        
        self.useBBres = useBBres 

        # Determining the number of parameters to fit and indexes for theta
        self.freepar_idx = {}
        self.ndim = 0
        self.freepar_names = []

        for pname in freepar:
            s = pname.split('_')
            if s[0] not in self._inri.r:
                raise ValueError(f"ERROR! The requested free parameter is unknown: {s[0]}")
            if len(s)==2 and 'single' in pname:
                self.freepar_idx[s[0]] = np.array([self.ndim])
                self.ndim += 1
                self.freepar_names.append(s[0])
            elif len(s)==2 and 'func' in pname:
                # @TODO: add possibility to use functional forms, e.g. freepar=['dens_func']
                ...
            elif len(s)==1:
                self.freepar_idx[s[0]] = np.arange(self.ndim,self.ndim+self._inri.nr,dtype='int')
                self.ndim += self._inri.nr
                for i in range(self._inri.nr):
                    self.freepar_names.append(f'{s[0]}{i+1}')

        # Setting the priors
        for key in self.freepar_idx:
            if self.priors[key] is None:
                self.priors[key] = scipy.stats.uniform(loc=self.bounds[key][0],\
                                                       scale=self.bounds[key][1]-self.bounds[key][0])
            elif not isinstance(self.priors[key],scipy.stats._distn_infrastructure.rv_continuous_frozen):
                raise ValueError(f"ERROR! Prior for {key} should be None or a scipy.stats distribution.")


        # Making a Param C++ object and a Galfit object
        self._opts.add_params(verbose=False)
        self._opts.make_object()
        
        self._galfit = libBB.Galfit_new_par(self.inp._cube,self._inri._rings,self._opts._params)

        # Getting the mask from the input cube
        self.mask  = reshapePointer(libBB.Cube_getMask(self.inp._cube),self.inp.dim[::-1])
        self.data  = reshapePointer(libBB.Cube_array(self.inp._cube),self.inp.dim[::-1])

        # Checking whether the density is fitted or not. In case it is not, use a normalization
        self.useNorm = not any('dens' in sub for sub in self.freepar_names)

        # Setting up the Ellprof object if a normalization and not using BB residuals 
        self._ellprof = libBB.Ellprof_new_alt(self.inp._cube,self._inri._rings)
        if self.useNorm and not self.useBBres:
            #self._update_profile(self._inri)
            # Check if we need to update the profile in each fit iteration (= only if geometry is fitted)
            self.update_prof = any(sub in string for string in self.freepar_names for sub in ['inc','phi','xpos','ypos'])

        ################################################################################################################
        # TO DELETE! This is just for testing purposes
        self.update_prof = False
        ################################################################################################################


    def _calculate_model(self,rings):
        """ This function calculates a unnormalized model given a set of rings.
        
            Parameters
            ----------
            rings : Rings
                A Rings object with the parameters of the model.

            Returns
            -------
            mod : 3D array
                A 3D array with the model.
            bhi, blo : 1D arrays
                The boundaries of the model in the input cube.
            galmod : pointer 
                A pointer to the Galmod object. Need to be freed after use.
        """
        # Model is built in a smaller array (blo,bhi) that only contains the galaxy 
        # within the last ring. This is to calculate the residuals faster.
        bhi, blo = (ctypes.c_int * 2)(0), (ctypes.c_int * 2)(0)
        libBB.Galfit_getModelSize(self._galfit,rings._rings,bhi,blo)
        galmod = libBB.Galfit_getModel(self._galfit,rings._rings,bhi,blo)
        bhi, blo = np.array(bhi), np.array(blo)

        # Reshaping the model to the correct 3D shape
        mod_shape = (self.inp.dim[-1], bhi[1]-blo[1],bhi[0]-blo[0])
        mod = reshapePointer(libBB.Galmod_array(galmod),mod_shape)

        return mod, bhi, blo, galmod



    def _normalize_model(self,model,data):
        """ This is the default normalization function for the model. 
            It can be overwritten by the user if needed.
        """
        # Normalizing the model to the maximum value of the data
        nrm = np.nanmax(data)/np.nanmax(model)
        return nrm*model



    def _calculate_residuals(self,model,data,mask=None):
        """ This is the default residual calculation function. 
             It can be overwritten by the user if needed.

            @TODO: residuals could use a WFUNC and BWEIGHT, similar to BB

        """
        # Applying the mask to the data if requested
        data_4res  = data*mask if mask is not None else data
        model_4res = self._normalize_model(model,data_4res) if self.useNorm else model

        # Some weighting steps can be added here (see below for example)
        wi = 1.0

        # Calculating the residual
        res = np.nansum(wi*np.abs(data_4res-model_4res))

        #okpix = (self.mask==1) & (mod!=0)
            
        # Weighting the residual by the azimuthal angle
        #x0  = np.mean(rings.r['xpos'])
        #y0  = np.mean(rings.r['ypos'])
        #inc = np.radians(np.mean(rings.r['inc']))
        #phi = np.radians(np.mean(rings.r['phi']))

        #x, y = np.arange(self.mask.shape[2]), np.arange(self.mask.shape[1])
        #xr =  -(x[None,:]-x0)*np.sin(phi)+(y[:,None]-y0)*np.cos(phi);
        #yr = (-(x[None,:]-x0)*np.cos(phi)-(y[:,None]-y0)*np.sin(phi))/np.cos(inc);

        #wpow = 2.0
        #th = np.arctan2(yr, xr)
        #wi = np.abs(np.cos(th))**wpow[None,:,:]

        # Calculating number of blanks (= model present, but data not in the mask)
        #nblanks = np.sum((self.mask==0) & (mod!=0))
        #ntot    = np.sum(self.mask)
        #if nblanks>ntot: nblanks = ntot-1
        #print (ntot,nblanks)
        #bweight = 1.0
        #res *= (1+nblanks/ntot)**bweight/(ntot-nblanks)
        return res


    def _log_likelihood(self,theta):
        """ Likelihood function for the fit """
        
        rings = self._update_rings(self._inri,theta)

        if self.useBBres:
            # Calculating residuals through BBarolo directly
            res = libBB.Galfit_calcresiduals(self._galfit,rings._rings)
        else: 
            # Calculating residuals manually            

            # Recompute the density profile along the current rings and update the rings
            if self.useNorm and self.update_prof:
                self._update_profile(rings)

            # Calculate the model and the boundaries
            mod, bhi, blo, galmod = self._calculate_model(rings)
            
            # Calculate the residuals
            mask = self.mask[:,blo[1]:bhi[1],blo[0]:bhi[0]]
            data = self.data[:,blo[1]:bhi[1],blo[0]:bhi[0]]
            res  = self._calculate_residuals(mod,data,mask)

            libBB.Galmod_delete(galmod)
            
        return -1000*res
        # This 1000 factor arbitrary, but with small residual there is no convergence. 
        # We need to understand why...
        
    
    def _prior_transform(self,u):
        """ Prior default transform function for dynesty.
            It defines flat priors for all parameters with min/max values 
            given by self.bounds
        """
        p = np.zeros_like(u)
        for key in self.freepar_idx:
           #p_min,p_max = self.bounds[key]
           #p[self.freepar_idx[key]] = p_min + u[self.freepar_idx[key]]*(p_max-p_min)
           p[self.freepar_idx[key]] = self.priors[key].ppf(u[self.freepar_idx[key]])
        return p
    
    


    def _compute(self,threads=1,freepar=['vrot'],method='dynesty', useBBres = True, \
                 sampler_kwargs : dict = {}, run_kwargs : dict = {}, **kwargs):
        
        """ Front-end function to fit a model.

        Run this function after having set initial parameters with :func:`init` and options with
        :func:`set_options`.

        Args: TBD
        
        Returns: TBD

        Kwargs: verbose (bool), dynamic (bool), other BB options
          
        """

        verbose  = kwargs.get('verbose',True)

        # Setting up all the needed class attributes
        self._setup(freepar,useBBres,**kwargs)

        # These are needed for the parallelization
        global prior_transform
        global log_likelihood

        def prior_transform(u):
            return self._prior_transform(u)

        def log_likelihood(theta):
            return self._log_likelihood(theta)

        mpisize = MPI.COMM_WORLD.Get_size()
        threads = 1 if mpisize>1 else threads

        pool = None
        if parallel:
            if   mpisize>1: pool = MPIPool()
            elif threads>1: pool = MultiPool(processes=threads)
        else:
            print ("WARNING! Parallelization is disabled.")

        # Now running the sampling 
        toc = time.time()
        if method=='dynesty': 
            
            DynestySampler = dyn.DynamicNestedSampler if kwargs.get('dynamic',True) else dyn.NestedSampler
            self.sampler = DynestySampler(log_likelihood, prior_transform, ndim=self.ndim, \
                                                    bound='multi',pool=pool,**sampler_kwargs)
            
            self.sampler.run_nested(print_progress=verbose,**run_kwargs)


            # Saving the results into class attributes
            self.results = self.sampler.results
            self.samples = self.results.samples  # Posterior samples
            self.weights = np.exp(self.results.logwt - self.results.logz[-1])
                                
        elif method=='nautilus':
            
            try: 
                import nautilus
            except ImportError:
                raise ImportError("BayesianBB requires the package nautilus when method='nautilus'.")
            
            if 'dlogz' in run_kwargs and 'f_live' not in run_kwargs:
                run_kwargs['f_live'] = run_kwargs.pop('dlogz')
    
            self.sampler = nautilus.Sampler(prior_transform,log_likelihood,n_dim=self.ndim, \
                                            pool=pool,pass_dict=False,**sampler_kwargs)
            self.sampler.run(verbose=verbose,**run_kwargs)

            self.samples, weights, _ = self.sampler.posterior()
            self.weights = np.exp(weights-weights.max())

        else: 
            raise ValueError(f"ERROR! Unknown method '{method}'.")
        
        
        self.samples = resample_equal(self.samples,self.weights)
        self.params = np.median(self.samples,axis=0)

        dt = time.time()-toc
        dt = f'{dt:.2f} seconds' if dt<60.00 else f'{dt/60.00:.2f} minutes' 
        print(f'Sampling with {method} took {dt} to run.')
        
        if pool is not None:
            pool.close()
            pool.join()
        
        self.modCalculated = True
    
    
    def _update_rings(self,rings,theta):
        """ Update rings with the parameters theta """
        for key in self.freepar_idx:
            pvalue = theta[self.freepar_idx[key]]
            if key=='dens': pvalue *= 1E20
            if len(pvalue)==1: pvalue = pvalue[0]
            rings.modify_parameter(key,pvalue)
        rings.make_object()
        return rings


    def _update_profile(self,rings):
        """ Get the density profile from the rings using Ellprof"""
        
        libBB.Ellprof_update_rings(self._ellprof,rings._rings)
        libBB.Ellprof_compute(self._ellprof)
        dens = reshapePointer(libBB.Ellprof_dens_array(self._ellprof),(1,self._inri.nr))[0]
        # To avoid problems with Galmod, we normalize the profile such that the minimum value is 1
        dens *= 1./np.nanmin(dens[dens>0])
        rings.modify_parameter("dens",np.abs(dens)*1E20,makeobj=True)

        
    def write_bestmodel(self,plots=True,**kwargs):
        
        if not self.modCalculated:
            print ("Sampler has not been run yet. Please run compute() before running this function.")
            return
        
        verbose = kwargs.get('verbose',True)

        # Creating a new Rings object for the outputs
        self.outri = Rings(self._inri.nr)
        self.outri.set_rings_from_dict(self._inri.r)
        
        # Updating output rings with best parameters from the sampling
        self.outri = self._update_rings(self.outri,self.params)

        # Setting up the output rings in the Galfit object.
        # @TODO: Add support for errors in the rings (already in the C++ code)
        libBB.Galfit_setOutRings(self._galfit,self.outri._rings)

        vprint(verbose,"\nWriting the best model to output directory ...",end=' ',flush=True)
        
        if self.useBBres:
            # Writing the model and plots to the output directory
            libBB.Galfit_writeModel(self._galfit,"AZIM".encode('utf-8'),plots)
        
        else:
            
            # If we are not fitting the dens, we need to update the final profile
            if self.useNorm:
                self._update_profile(self.outri)

            # Deriving the last model
            _, ys, xs = self.data.shape
            bhi, blo = (ctypes.c_int * 2)(xs,ys), (ctypes.c_int * 2)(0)
            galmod = libBB.Galfit_getModel(self._galfit,self.outri._rings,bhi,blo)

            if self.useNorm:
                # Reshaping the model to the correct 3D shape
                mod = reshapePointer(libBB.Galmod_array(galmod),self.data.shape)
                # Normalizing and copying it back to the C++ Galmod object
                mod = self._normalize_model(mod,self.data)
                libBB.Galmod_set_array(galmod,np.ravel(mod).astype('float32'))
            
            # Writing all the outputs
            libBB.Galfit_writeOutputs(self._galfit,galmod,self._ellprof,plots)
        
        vprint(verbose,"Done!")


    def print_stats(self,burnin=None,percentiles=[15.865,50,84.135]):

        if not self.modCalculated:
            print ("Sampler has not been run yet. Please run compute() before running this function.")
            return
        
        if not isIterable(percentiles) and len(percentiles)!=3:
            raise ValueError("ERROR! Percentiles must be a list of three values.")

        results = copy.copy(self.sampler.results)
        if burnin is not None:
            results = results[burnin:]
        weights = np.exp(results.logwt - results.logz[-1])
        samples = resample_equal(results.samples,weights)

        # Printing central values and error of posterior distributions
        print("\nBest-fit parameters:")
        for i in range(len(self.freepar_names)):
            values = np.percentile(samples[:, i], percentiles)
            q = np.diff(values)
            txt = "  %10s = %10.3f %+10.3f %+10.3f"%(self.freepar_names[i], values[1], -q[0], q[1])
            print (txt)

        # Printing the maximum likelihood
        max_logl   = results.logl.max()
        max_index  = results.logl.argmax()
        max_params = samples[max_index]
        print(f"\nMaximum logl: {max_logl} at:")
        for i in range(len(self.freepar_names)):
            txt = "  %10s = %10.3f "%(self.freepar_names[i], max_params[i])
            print(txt)

        return results, samples







        