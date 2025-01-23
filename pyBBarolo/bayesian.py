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

import numpy as np
import scipy.stats
from .BB_interface import libBB
from .pyBBarolo import Param, Rings, FitMod3D, reshapePointer, vprint

import time

from schwimmbad import MultiPool, MPIPool
from mpi4py import MPI

try: 
    # We need global import for resample_equal. Is it really needed?
    import dynesty as dyn
    from dynesty.utils import resample_equal
except ImportError:
    raise ImportError("BayesianBB requires the package dynesty'")


class BayesianBBarolo(FitMod3D):
    
    """ A class to fit a galaxy model to a 3D datacube using a Bayesian framework
    
    Args:
      fitsname (str): FITS file of the galaxy to fit
    
    """
    
    def __init__(self,fitsname,**kwargs):
        super(BayesianBBarolo,self).__init__(fitsname=fitsname)
        # Task name
        self.taskname = "BAYESIAN3DFIT"
        # Resetting FitMod3d._args. Not used here.
        self._args = {}
        # self._opts contains any other parameter given to the code
        self._opts = Param(fitsfile=fitsname,outfolder='./output/')
        if kwargs:
            self._opts.add_params(kwargs)
        
        # A dictionary with the names of parameters to be fitted and their indexes
        self.freepar_idx = None
        self.ndim = 0
        # A list with the "names" of fitted parameters
        self.freepar_names = None
        # A dictionary with the boundaries for the priors of parameters
        self.bounds = None
        # A dicionary for prior probabilities
        self.priors = None
        # A pointer to the C++ Galfit object
        self._mod = None
        
        self.modCalculated = False
    
     
    def set_options(self,**kwargs):
        self._opts.add_params(**kwargs)
    
    
    def show_options(self):
        print(self._opts)
    
    
    def init(self,**kwargs):
        """ Initialize rings for the fit. Parameters are those of FitMod3D.init()
        
        Set initial guesses for the fit and fixed parameters. Parameters not set are estimated.
        @TODO: This could be easily modified to fit also dens, vvert, etc...
        
        """
        # Using the Fit3D initialization function (see pyBBarolo.py)
        super(BayesianBBarolo,self)._input(**kwargs)
        
        # Defining default boundaries
        rr = self._inri.r
        self.bounds = dict(xpos=[rr['xpos'][0]-5,rr['xpos'][0]+5],     #+- 5 pixels
                           ypos=[rr['ypos'][0]-5,rr['ypos'][0]+5],     #+- 5 pixels
                           vsys=[rr['vsys'][0]-50,rr['vsys'][0]+50],   #+- 50 km/s
                           vrot=[0,350],vdisp=[0,50],vrad=[0,50],
                           inc=[0,90],phi=[0,360],
                           vvert=[0,30],dvdz=[0,2],zcyl=[0,5],
                           dens=[0.01,200],z0=[0,10])
        self.priors = {key: None for key in self.bounds.keys()}
    
    
    def _log_likelihood(self,theta,useBBres=True):
        """ Likelihood function for the fit """
        
        # Interpreting theta based on free parameters fitted and update rings.
        rings = self._inri
        for key in self.freepar_idx:
            pvalue = theta[self.freepar_idx[key]]
            if len(pvalue)==1: pvalue = pvalue[0]
            rings.modify_parameter(key,pvalue)
        rings.make_object()
        
        if useBBres:
            # Calculating residuals through BBarolo directly
            res = libBB.Galfit_calcresiduals(self._mod,rings._rings)
        else: 
            # Calculating residuals manually
            # @TODO for having a residual calculation outside BB.
            # There are two possibilities:
            # 1) Fit the density profile (dens) as any other parameter. 
            #    - It could be functional density -> just build model with profile
            #    - Normalize to maximum or top 10% before residual calculation.
            # 2) Use azimuthal normalization. In this case:
            #   - Need to get a mask (BB or external function?)
            #   - Need a density profile given the ring geometry (BB or external function?) 
            #   - Build a model with derived density profile 
            #   - Normalize to maximum or top 10% before residual calculation.
            #
            #   In both cases, residual could use a WFUNC and BWEIGHT
            
            
            self._mod = libBB.Galmod_new_par(self.inp._cube,rings._rings,self._opts._params)
            libBB.Galmod_compute(self._mod)
            libBB.Galmod_smooth(self._mod)

            mod = reshapePointer(libBB.Galmod_array(self._mod),self.inp.dim[::-1])
            nrm = 1.00

            # @TODO: This is just a temporary solution for checking whether _log_likelihood is working.
            #        It should be integrated with a normalization step.
            res = np.sum(np.abs(self.inp._cube-nrm*mod))
        
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
    
    
    def _compute(self,threads=1,freepar=['vrot'],method='dynesty',\
                 sampler_kwargs : dict = {}, run_kwargs : dict = {}, **kwargs):
        
        """ Front-end function to fit a model.

        Run this function after having set initial parameters with :func:`init` and options with
        :func:`set_options`.

        Args: TBD
        
        Returns: TBD
          
        """

        if self._inri is None: 
            print ("Error: you need to initialize rings first by calling the function init()")
            return
        
        # Making a Param C++ object and a Galfit/Galmod object
        self._opts.add_params(verbose=False,twostage=False)
        self._opts.make_object()
        
        useBBres = kwargs.get('useBBres',True)
        libBBres = libBB.Galfit_new_par if useBBres else libBB.Galmod_new_par
        self._mod = libBBres(self.inp._cube,self._inri._rings,self._opts._params)

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
            
        # These are needed for the parallelization
        global prior_transform
        global log_likelihood

        def prior_transform(u):
            return self._prior_transform(u)

        def log_likelihood(theta):
            return self._log_likelihood(theta,useBBres)

        mpisize = MPI.COMM_WORLD.Get_size()
        threads = 1 if mpisize>1 else threads

        if   mpisize>1: pool = MPIPool()
        elif threads>1: pool = MultiPool(processes=threads)
        else: pool = None
    
        verbose = kwargs.get('verbose',True)
        
        ##############
        return 
        ##############
        
        # Now running the sampling 
        toc = time.time()
        if method=='dynesty': 
            
            DynestySampler = dyn.DynamicNestedSampler if kwargs.get('dynamic',True) else dyn.NestedSampler
            self.sampler = DynestySampler(log_likelihood, prior_transform, ndim=self.ndim, \
                                                    bound='multi',pool=pool,**sampler_kwargs)
            
            self.sampler.run_nested(print_progress=verbose,**run_kwargs)

            self.results = self.sampler.results

            # Extract the best-fit parameters
            self.samples = self.results.samples  # Posterior samples
            weights = np.exp(self.results.logwt - self.results.logz[-1])
                                
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
            weights = np.exp(weights-weights.max())
        else: 
            raise ValueError(f"ERROR! Unknown method {method}.")
        
        
        self.samples = resample_equal(self.samples,weights)
        self.params = np.median(self.samples,axis=0)

        dt = time.time()-toc
        dt = f'{dt:.2f} seconds' if dt<60.00 else f'{dt/60.00:.2f} minutes' 
        print(f'Sampling with {method} took {dt} to run.')
        
        if pool is not None:
            pool.close()
            pool.join()
        
        self.modCalculated = True
    
    
    def write_bestmodel(self,**kwargs):
        
        #if not self.modCalculated:
        #    print ("Best model not calculated yet. Please run compute() before running this function.")
        #    return
        
        verbose = kwargs.get('verbose',True)


        self.params = np.array([ 55.67056348, 108.24562833,  84.81053554,  95.2849607,    7.51556276,
                       10.1689289 ,  14.79084093, 12.08063944 , 65.48071344,  30.08482877])

        #self.outri = copy.deepcopy(self._inri)
        self.outri = Rings(self._inri.nr)
        self.outri.set_rings_from_dict(self._inri.r)
        
        # Updating output rings with best parameters from the sampling
        for key in self.freepar_idx:
            pvalue = self.params[self.freepar_idx[key]]
            if len(pvalue)==1: pvalue = pvalue[0]
            self.outri.modify_parameter(key,pvalue)
        self.outri.make_object()

        
        vprint(verbose,"Writing the best model to output directory ...",end=' ',flush=True)
        # Setting up the output rings in the Galfit object.
        #@TODO: Add support for errors in the rings (already in the C++ code)
        libBB.Galfit_setOutRings(self._mod,self.outri._rings)
        # Writing the model and plots to the output directory
        libBB.Galfit_writeModel(self._mod,"AZIM".encode('utf-8'),True)
        vprint(verbose,"Done!")




        