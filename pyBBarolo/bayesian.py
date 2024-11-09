"""
This module allows to use pyBBarolo in a Bayesian fashion.
It uses either dynesty or emcee libraries for this.
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

import os,sys
import numpy as np
from .BB_interface import libBB
from .pyBBarolo import Param, Rings, FitMod3D
from dynesty import DynamicNestedSampler

#import emcee
#from scipy.optimize import minimize

from schwimmbad import MultiPool, MPIPool
from mpi4py import MPI

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
        self._opts = Param(fitsfile=fitsname)
        if kwargs:
            self._opts.add_params(kwargs)
        
        # A dictionary with the names of parameters to be fitted and their indexes
        self.freepar_idx = None
        self.ndim = 0
        # A list with the "names" of fitted parameters
        self.freepar_names = None
        # A dictionary with the boundaries for the priors of parameters
        self.bounds = None
        # A pointer to the C++ Galfit object
        self._mod = None
    
     
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

    
    def _log_likelihood(self,theta):
        """ Likelihood function for the fit """
        
        # Interpreting theta based on free parameters fitted and update rings.
        rings = self._inri
        for key in self.freepar_idx:
            pvalue = theta[self.freepar_idx[key]]
            if len(pvalue)==1: pvalue = pvalue[0]
            rings.modify_parameter(key,pvalue)
        rings.make_object()
        
        # Calculating residuals through BBarolo directly
        # @TODO: This could be easily modified to just let BB making a model
        #        and then computing the residuals directly in python
        res = libBB.Galfit_calcresiduals(self._mod,rings._rings)
        
        ###### STUFF TO BE REMOVED #################
        self.count += 1
        if res<self.bestres:
            self.bestres = res
            self.bestfit = theta
        #print (theta, rings.r)
        #if self.count>10: exit()
        ############################################
        
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
            p_min,p_max = self.bounds[key]
            p[self.freepar_idx[key]] = p_min + u[self.freepar_idx[key]]*(p_max-p_min)
        
        return p
    
    
    def _compute(self,threads=1,freepar=['vrot'],method='dynesty',\
                 sampler_kwargs : dict = {}, run_kwargs : dict = {}):
        
        """ Front-end function to fit a model.

        Run this function after having set initial parameters with :func:`init` and options with
        :func:`set_options`.

        Args: TBD
        
        Returns: TBD
          
        """
        
        if self._inri is None: 
            print ("Error: you need to initialize rings first by calling the function init()")
            return
        
        # Making a Param C++ object and a Galfit object
        self._opts.add_params(verbose=False)
        self._opts.make_object()
        self._mod = libBB.Galfit_new_par(self.inp._cube,self._inri._rings,self._opts._params)
        
        
        ######### ALL STUFF THAT CAN BE REMOVED ############
        self.count = 0
        self.bestfit = 0
        self.bestres = 1E20
        ####################################################
        
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

        # These are needed for the parallelization
        global prior_transform
        global log_likelihood

        def prior_transform(u):
            return self._prior_transform(u)

        def log_likelihood(theta):
            return self._log_likelihood(theta)


        mpisize = MPI.COMM_WORLD.Get_size()
        threads = 1 if mpisize>1 else threads

        if   mpisize>1: pool = MPIPool()
        elif threads>1: pool = MultiPool()
        else: pool = None
        
        # Now fitting
        if method=='dynesty':
            self.sampler = DynamicNestedSampler(log_likelihood, prior_transform, ndim=self.ndim, \
                                                bound='multi',pool=pool,**sampler_kwargs)
            self.sampler.run_nested(**run_kwargs)
            self.results = self.sampler.results

            # Extract the best-fit parameters
            samples = self.results.samples  # Posterior samples
            weights = np.exp(self.results.logwt - self.results.logz[-1])
            params = np.average(samples, axis=0, weights=weights)
            
            print (params)

        ''' We could support emcee as well...
        elif method=='emcee':
        
            # Log-prior function
            def log_prior(theta):
                m, b = theta
                if 0 < m < 200 and 0 < b < 15:
                    return 0.0
                return -np.inf

            # Log-probability function
            def log_probability(theta):
                lp = log_prior(theta)
                if not np.isfinite(lp):
                    return -np.inf
                return lp + self.log_likelihood(theta)
            
            n_walkers=50 
            n_steps=3000
            burn_in=1000
            # Initial guess and setting up the sampler
            initial = np.array([120, 10])  # Initial guess for slope and intercept
            pos = initial + 10 * np.random.randn(n_walkers, 2)
            sampler = emcee.EnsembleSampler(n_walkers, 2, log_probability)

            # Run the MCMC chain
            sampler.run_mcmc(pos, n_steps, progress=True)
            samples = sampler.get_chain(discard=burn_in, flat=True)  # Discard burn-in samples

            # Extract best-fit values
            m_median, b_median = np.median(samples, axis=0)
            print(f"Best-fit slope (m): {m_median}")
            print(f"Best-fit intercept (b): {b_median}")
            
        
        elif method=='simplex':
            
            def funcmin(theta):
                if np.any(theta<0):
                    return 1E10
                else:
                    return -self.log_likelihood(theta)
            
            # Initial guess for the parameters
            initial_guess = [100, 100, 100, 100, 10,10,10,10]

            # Minimize the chi-squared function using the Nelder-Mead (downhill simplex) method
            result = minimize(funcmin, initial_guess, method='Nelder-Mead',tol=1E-10)
            bf = result.x  # Extract best-fit parameters

            print (bf)
        '''
        else: 
            raise ValueError(f"ERROR! Unknown method {method}.")

        
        if pool is not None:
            pool.close()
            pool.join()
            