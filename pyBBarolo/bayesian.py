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

import ctypes, time, copy, os
from matplotlib import axes
import numpy as np
from .BB_interface import libBB
from .pyBBarolo import Param, Rings, FitMod3D, reshapePointer, vprint, isIterable

try: 
    import dynesty as dyn
    from dynesty import utils as dyfunc
    from scipy import stats
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("BayesianBB requires the packages 'scipy', 'matplotlib' and 'dynesty'. Please install them.")

parallel = True
try: 
    from schwimmbad import MultiPool, MPIPool
except ImportError:
    print ("WARNING! Parallelization is disabled. Please install schwimmbad and mpi4py to enable it.")
    parallel = False

# If we are running in parallel with MPI, we need to import mpi4py
isMPI = "OMPI_COMM_WORLD_SIZE" in os.environ or "PMI_RANK" in os.environ
if isMPI:
    try:
        from mpi4py import MPI
    except ImportError:
        raise ValueError ("ERROR! To run in parallel with MPI, you need to install mpi4py.")


def get_distribution(distr_name, *params, **kwargs):
    """Create a frozen scipy.stats distribution given its name and parameters."""
    try:
        # Get the distribution class from scipy.stats
        distr = getattr(stats, distr_name)
        # Return the frozen distribution with given parameters
        return distr(*params, **kwargs)
    except AttributeError:
        raise ValueError(f"ERROR: Distribution '{distr_name}' not found in scipy.stats.")


class BayesianBBarolo(FitMod3D):
    
    """ A class to fit a galaxy model to a 3D datacube using a Bayesian framework
    
    Attributes
    ----------
    fitsname : str
      A string with the name of the fits file to be used.

    ...
      
    Methods
    -------
    set_options(**kwargs):
        Set all options for BBarolo.

    show_options():
        Print the options set for BBarolo.
    
    init(radii,xpos,ypos,vsys,vrot,vdisp,inc,phi,z0,vrad,dens,vvert,dvdz,zcyl):
        Initialize the rings of the model.
    
    compute(threads,freepar,method,useBBres,sampler_kwargs,run_kwargs,like_kwargs, **kwargs):
        Fit the model to the data using a Bayesian sampler.

    write_bestmodel(plots,**kwargs):
        Write the best model and plots using the BBarolo library.

    corner_plot(samples,labels,truths,quantiles,plot_dynesty,**kwargs):
        Create a corner plot of the posterior distributions.

    plot_profiles(**kwargs):
        Plot the best-fit profiles of the fitted parameters.

    print_priors():
        Print the priors set for the fit.
    
    print_stats():
        Compute and print statistics of the fit.

    load_results(inputfile):
        Load previously saved results from a file.
    
    save_results(outfile):
        Save results to a file.
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
        # Calling the parent class to initialize the input cube
        super(BayesianBBarolo,self).__init__(fitsname=fitsname, astropy_pointer=False)
        # Task name
        self.taskname = "BAYESIAN3DFIT"
        # Resetting FitMod3d._args. Not used here.
        self._args = {}
        # self._opts contains any other parameter given to the code
        defaults = dict(threads=1,outfolder='./output/',twostage=False)
        self._opts = Param(fitsfile=fitsname,**defaults)
        if kwargs:
            self._opts.add_params(**kwargs)
        
        # A dictionary with the names of parameters to be fitted and their indexes
        self.freepar_idx = None
        self.ndim = 0
        # A list with the "names" of fitted parameters
        self.freepar_names = None
        # A dictionary with the priors of parameters
        self.priors = {}
        # A dictionary for prior probability distributions
        self.prior_distr = None
        # A dictionary with functions for parametric fits
        self.funcs = {}
        # A pointer to the C++ Galfit and Ellprof objects
        self._galfit = self._ellprof = None
        self.modCalculated = False
        # Input and output rings
        self._inri = self.outri = None

    
    def set_options(self,**kwargs):
        """ Add options to be passed to BBarolo """
        self._opts.add_params(**kwargs)
    
    
    def show_options(self):
        """ Print the options to be passed to BBarolo """
        print(self._opts)
    
    
    def init(self,radii,xpos,ypos,vsys=0,vrot=100,vdisp=10,inc=0,phi=0,z0=0,vrad=0,dens=1,vvert=0,dvdz=0,zcyl=0):
        """ Initialize rings for the fit. Parameters that are not fitted will be fixed.
        
        Parameters
        ----------
        radii : list
            Radii of the model in arcsec
        xpos : float, list (optional)
            X center of the galaxy in pixels
        ypos : float, list (optional)
            Y center of the galaxy in pixels
        vsys : float, list (optional)
            Systemic velocity of the galaxy in km/s
        vrot : float, list (optional)
            Rotation velocity in km/s
        vdisp : float, list (optional)
            Velocity dispersion in km/s
        inc : float, list (optional)
            Inclination angle in degrees
        phi : float, list (optional)
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
        if self._inri is not None:
            print ("Error: rings can not be re-initialized. Please create a new BayesianBBarolo object.")
            return

        # Initialize rings
        if not isIterable(radii): raise ValueError("radii must be an array")
        self._inri = Rings(len(radii))
        self._inri.set_rings(radii,xpos,ypos,vsys,vrot,vdisp,vrad,vvert,dvdz,zcyl,dens,z0,inc,phi)
        
        # Defining default uniform priors 
        rr = self._inri.r
        self.priors['xpos']  = dict(name='uniform',loc=rr['xpos'][0]-5,scale=10)    # +- 5 pixels
        self.priors['ypos']  = dict(name='uniform',loc=rr['ypos'][0]-5,scale=10)    # +- 5 pixels
        self.priors['vsys']  = dict(name='uniform',loc=rr['vsys'][0]-20,scale=40)   # +- 20 km/s
        self.priors['vrot']  = dict(name='uniform',loc=0,scale=350)                 # 0-350 km/s
        self.priors['vdisp'] = dict(name='uniform',loc=0,scale=50)                  # 0-50 km/s    
        self.priors['vrad']  = dict(name='uniform',loc=0,scale=50)                  # 0-50 km/s
        self.priors['inc']   = dict(name='uniform',loc=0,scale=90)                  # 0-90 degrees
        self.priors['phi']   = dict(name='uniform',loc=0,scale=360)                 # 0-360 degrees
        self.priors['vvert'] = dict(name='uniform',loc=0,scale=30)                  # 0-30 km/s
        self.priors['dvdz']  = dict(name='uniform',loc=0,scale=2)                   # 0-2 km/s/pc
        self.priors['zcyl']  = dict(name='uniform',loc=0,scale=5)                   # 0-5 arcsec
        self.priors['dens']  = dict(name='uniform',loc=0.1,scale=200)               # 0.1-200 1E20 cm^-2
        self.priors['z0']    = dict(name='uniform',loc=0,scale=10)                  # 0-10 arcsec
        self.priors['sigma'] = dict(name='halfcauchy',scale=0.5)
        
        # Priors for maximum radius
        radsep = np.abs(rr['radii'][-1]-rr['radii'][-2])
        self.priors['radmax']= dict(name='uniform',loc=rr['radii'][-1]-radsep/2.,scale=radsep) # +- radsep

        # Defining the parametric functions. Default only for dens, vrot and vdisp
        self.funcs['vrot']   = lambda R,p1,p2: 2./np.pi*p1*np.arctan(R/p2)          # Arctan function for vrot
        self.funcs['vdisp']  = lambda R,p1,p2: p1 + p2*R                            # Line for vdisp
        self.funcs['dens']   = lambda R,p1,p2: p1*np.exp(-R/p2)                     # Exponential for dens

        # Defining default priors for the functional forms
        self.priors['vrot_p1']  = dict(name='uniform',loc=1,scale=350)  
        self.priors['vrot_p2']  = dict(name='uniform',loc=1,scale=1000) 
        self.priors['vdisp_p1'] = dict(name='uniform',loc=0,scale=10)  
        self.priors['vdisp_p2'] = dict(name='uniform',loc=-10,scale=20) 
        self.priors['dens_p1']  = dict(name='uniform',loc=1,scale=20)  
        self.priors['dens_p2']  = dict(name='uniform',loc=1,scale=200) 

        # Initializing the prior distributions to None
        self.prior_distr = {key: None for key in self.priors.keys()}

        # Lambda parameters for smoothness penalty in the likelihood
        self.lambdas = {key: 0.0 for key in ['vrot','vdisp','dens','vrad','phi','inc','xpos','ypos','z0']}   
    
    
    def _autoinit():
        """ This function initializes the rings in case the user does not call init() before compute(). 
            The rings are estimated with BBarolo's ParamGuess class.
        """
        raise NotImplementedError("Automatic initialization of rings is not implemented yet. Please call the init() function before compute().")
        

    def _setup(self,freepar,useBBres=False,**kwargs):
        """ Setup function for the fit 
            It sets the following class attributes:
            - self.ndim: number of parameters to fit
            - self.freepar_idx: indexes of the parameters to fit
            - self.freepar_names: names of the parameters to fit
            - self._freepar_inp: input freepar list
            - self.prior_distr: dictionary with prior probabilities
            - self._galfit: pointer to the C++ Galfit object
            - self.mask: mask from the input cube (3D array)
            - self._ellprof: a pointer to a C++ Ellprof object
        """
        

        if self._inri is None:
            self._autoinit()

        self.useBBres = useBBres 

        # Determining the number of parameters to fit and indexes for theta
        self.freepar_idx = {}
        self.ndim = 0
        self.freepar_names = []
        self._freepar_inp = freepar

        dup = {key : False for key in self.priors.keys()}

        okpars = ['sigma','radmax']

        for pname in freepar:
            s = pname.split('_')
            if s[0] not in self._inri.r and s[0] not in okpars:
                raise ValueError(f"ERROR! The requested free parameter is unknown: {s[0]}")
            if dup[s[0]]:
                raise ValueError(f"ERROR! Multiple entries for {s[0]} parameter. Please choose one!")
            dup[s[0]] = True

            if len(s)==2 and 'single' in pname:
                self.freepar_idx[s[0]] = np.array([self.ndim])
                self.ndim += 1
                self.freepar_names.append(s[0])
            elif len(s)==2 and 'func' in pname:
                if s[0] not in self.funcs:
                    raise ValueError(f"ERROR! The functional form for {s[0]} is unknown. Please define it first.")
                nparams = self.funcs[s[0]].__code__.co_argcount-1
                for i in range(nparams):
                    self.freepar_names.append(f'{s[0]}_p{i+1}')
                    self.freepar_idx[f'{s[0]}_p{i+1}'] = np.arange(self.ndim,self.ndim+1,dtype='int')
                    self.ndim += 1
            elif len(s)==1:
                if s[0] in okpars:
                    self.freepar_idx[s[0]] = np.array([self.ndim])
                    self.ndim += 1
                    self.freepar_names.append(s[0])
                else:
                    self.freepar_idx[s[0]] = np.arange(self.ndim,self.ndim+self._inri.nr,dtype='int')
                    self.ndim += self._inri.nr
                    for i in range(self._inri.nr):
                        self.freepar_names.append(f'{s[0]}{i+1}')
        
        # Setting the priors
        for key in self.freepar_idx:
            if self.prior_distr[key] is None:
                # Prior distribution is not set externally, so we use the defined self.priors
                try:
                    distr_kw = {key: value for key, value in self.priors[key].items() if key != "name"}
                    self.prior_distr[key] = get_distribution(self.priors[key]['name'],**distr_kw)
                except:
                    raise ValueError(f"Something is wrong in the prior for parameter '{key}'. Please double check!")
                    
            elif not isinstance(self.prior_distr[key], stats._distn_infrastructure.rv_continuous_frozen):
                raise ValueError(f"ERROR! Prior for {key} should be None or a scipy.stats distribution.")
        
        # Making a Param C++ object and a Galfit object
        self._opts.add_params(verbose=False,**kwargs)
        self._opts.make_object()

        if self._galfit is None: 
            self._galfit = libBB.Galfit_new_par(self.inp._cube,self._inri._rings,self._opts._params)
            # Getting the mask from the input cube
            self.mask  = reshapePointer(libBB.Cube_getMask(self.inp._cube),self.inp.dim[::-1])
            self.data  = reshapePointer(libBB.Cube_array(self.inp._cube),self.inp.dim[::-1])
            self.data_mom0 = np.nansum(self.data*self.mask,axis=0)


        ##########################################################################################
        # Normalizing data to a maximum of 10. This is to avoid numerical issues in the fit.
        # This is just for testing purposes
        
        # Calculate the noise from the median absolute deviation of the data in the masked region
        self.noiserms = 1.4826*np.median(np.abs(self.data[self.mask==0]-np.median(self.data[self.mask==0])))
        # The below is a pointer, so it modifies also the data array in the Galfit object!
        #self.scaling_factor = 1./np.nanmax(self.data)
        #self.data *= self.scaling_factor
        #self.data /= self.noiserms
        self.map_rms   = self.noiserms*np.sqrt(np.nansum(self.mask, axis=0))
        ##########################################################################################

        # Checking whether the density is fitted or not. In case it is not, use a normalization
        self.useNorm = not any('dens' in sub for sub in self.freepar_names)
   
        if self._ellprof is None:
            self._ellprof = libBB.Ellprof_new_alt(self.inp._cube,self._inri._rings)
        
        # Calculating the initial profile if we are using normalization. This is updated in each fit iteration if geometry is fitted.
        if self.useNorm and not self.useBBres:
            self._update_profile(self._inri)
            # Check if we need to update the profile in each fit iteration (= only if geometry is fitted)
            self.update_prof = any(sub in string for string in self.freepar_names for sub in ['inc','phi','xpos','ypos'])


    def _calculate_model(self,rings,fullcube=False):
        """ This function calculates a unnormalized model given a set of rings.
        
            Parameters
            ----------
            rings : Rings
                A Rings object with the parameters of the model.
            fullcube : bool (False by default)
                Whether to calculate the model in the full input cube (True) or only in a smaller cube containing the galaxy (False). 

            Returns
            -------
            mod : 3D array
                A 3D array with the model.
            bhi, blo : 1D arrays
                The boundaries of the model in the input cube.
            galmod : pointer 
                A pointer to the Galmod object. Need to be freed after use.
        """
        
        _, ys, xs = self.data.shape
        bhi, blo = (ctypes.c_int * 2)(xs,ys), (ctypes.c_int * 2)(0)

        if not fullcube:
            # Model is built in a smaller array (blo,bhi) that only contains the galaxy 
            # within the last ring. This is to calculate the residuals faster.
            bhi, blo = (ctypes.c_int * 2)(0), (ctypes.c_int * 2)(0)
            libBB.Galfit_getModelSize(self._galfit,rings._rings,bhi,blo)

        galmod = libBB.Galfit_getModel(self._galfit,rings._rings,bhi,blo,True)
        bhi, blo = np.array(bhi), np.array(blo)

        # Reshaping the model to the correct 3D shape
        mod_shape = (self.inp.dim[2], bhi[1]-blo[1],bhi[0]-blo[0])
        mod = reshapePointer(libBB.Galmod_array(galmod),mod_shape)

        return mod, bhi, blo, galmod


    def _calculate_residuals(self,model,data,mask=None,**kwargs):
        """ This is the default residual calculation function. 
             It can be overwritten by the user if needed.
        """
        # Applying the mask to the data if requested
        data_4res  = data*mask if mask is not None else data
        sigma = kwargs.get('sigma',self.noiserms)
        if sigma<0:
            return np.inf

        # Calculating the residual
        res  = 0.5*np.nansum( (data_4res-model)**2 / sigma**2)
        res += 0.5*np.log(2.*np.pi*sigma**2)

        return res


    def _log_likelihood(self,theta,**kwargs):
        """ Likelihood function for the fit 
        
            Accepted kwargs:
            - factor: a multiplicative factor for the log-likelihood. Default is 1.
            - lambdas: a dictionary with the lambda parameters for the smoothness penalty. Default is self.lambdas.
        """
        
        # Checking for non-acceptable negative values
        a = ('vrot','vdisp','inc','xpos','ypos','z0','dens','radmax')
        for k in self.freepar_idx:
            if k.startswith(a) and np.any(theta[self.freepar_idx[k]]<0): 
                return -np.inf
    
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
            model, bhi, blo, galmod = self._calculate_model(rings)
            
            # Calculate the residuals
            mask = self.mask[:,blo[1]:bhi[1],blo[0]:bhi[0]]
            data = self.data[:,blo[1]:bhi[1],blo[0]:bhi[0]]

            kwargs['sigma'] = theta[self.freepar_idx['sigma']][0] if 'sigma' in self.freepar_names else 1

            logL  = self._calculate_residuals(model,data,mask,**kwargs)

            libBB.Galmod_delete(galmod)

            # Adding a smoothness penalty if requested in the ring-by-ring fit
            lambdas = kwargs.get('lambdas', self.lambdas)
            smoothness = any(lambdas[key]>0 for key in lambdas.keys())
            if smoothness:
                for key in self.freepar_idx:
                    if key in lambdas and lambdas[key]>0 and len(self.freepar_idx[key])>1:
                        v = rings.r[key]
                        
                        if key=='vrot':
                            # Adding a zero at the beginning to force the rotation curve to start from zero
                            v = np.concatenate([[0.0], v])

                        dr = rings.r['radii'][1] - rings.r['radii'][0]
                        d2v = (v[2:] - 2*v[1:-1] + v[:-2]) / dr**2

                        penalty = np.sum(d2v**2)
                        #penalty = np.sum((d2v / (1 + r[1:-1]))**2)
                        logL += 0.5 * lambdas[key] * penalty
            
            factor = kwargs.get('factor',1.0)

        return -factor*logL
    

    def _log_likelihood_geometry(self,theta,**kwargs):
        """ Likelihood function for the fit """
        
        # Checking for non-acceptable negative values
        a = ('inc','xpos','ypos','dens','z0','radmax')
        for k in self.freepar_idx:
            if k.startswith(a) and np.any(theta[self.freepar_idx[k]]<0): 
                return -np.inf
    
        rings = self._update_rings(self._inri,theta)

        # Recompute the density profile along the current rings and update the rings
        if self.useNorm and self.update_prof:
            self._update_profile(rings)

        # Calculate the model and the boundaries
        model, bhi, blo, galmod = self._calculate_model(rings,fullcube=False)
        
        # Calculate the residuals
        map_mod  = np.nansum(model,axis=0)
        map_data = self.data_mom0[blo[1]:bhi[1],blo[0]:bhi[0]]
        map_rms  = self.map_rms[blo[1]:bhi[1],blo[0]:bhi[0]]

        map_data = map_data[map_rms!=0]
        map_mod  = map_mod[map_rms!=0]
        map_rms  = 1#np.nanmedian(map_rms[map_rms!=0])

        # Some normalization (arbitrary) to help convergence
        maxdata = 10.#*map_rms #np.nanmax(map_data)
        map_data *= 1/maxdata
        map_mod  *= 1/maxdata
        
        res  = 0.5*np.nansum((map_data-map_mod)**2 / map_rms**2)
        res += 0.5*np.log(2.*np.pi*map_rms**2)

        libBB.Galmod_delete(galmod)

        return -res
        
    
    def _prior_transform(self,u):
        """ Prior default transform function for sampler.
            It defines priors based on the user-defined distributions.
        """
        p = np.zeros_like(u)
        for key in self.freepar_idx:
           p[self.freepar_idx[key]] = self.prior_distr[key].ppf(u[self.freepar_idx[key]])
        return p
    
    
    def _compute(self,threads=1,freepar=['vrot'],method='dynesty', useBBres=False, \
                 log_like = None, like_kwargs : dict = {},
                 sampler_kwargs : dict = {}, run_kwargs : dict = {}, **kwargs):
        
        """ Front-end function to fit a model.

        Run this function after having set initial parameters with :func:`init` and options with
        :func:`set_options`.

        Args: TBD
        
        Returns: TBD

        Kwargs: verbose (bool), dynamic (bool), other BB options
          
        """
        # Setting up all the needed class attributes
        self._setup(freepar,useBBres,**kwargs)
        
        log_like = log_like if log_like is not None else self._log_likelihood
        verbose  = kwargs.get('verbose',True)

        if verbose:
            print (f"\nFitting {self.ndim} parameters: {self.freepar_names}")
            self.print_priors()


        # These are needed for the parallelization
        global prior_transform
        global log_likelihood

        def prior_transform(u):
            return self._prior_transform(u)

        def log_likelihood(theta):
            return log_like(theta,**like_kwargs)


        # Setting up the parallelization
        pool = None
        if isMPI:
            mpisize = MPI.COMM_WORLD.Get_size()
            pool = MPIPool() if mpisize else None
        elif threads>1:
            if parallel:
                pool = MultiPool(processes=threads)
            else:
                print ("WARNING: Parallelization is disabled!")

        # Now running the sampling 
        toc = time.time()
        if method=='dynesty': 
            
            DynestySampler = dyn.DynamicNestedSampler if kwargs.get('dynamic',True) else dyn.NestedSampler
            self.sampler = DynestySampler(log_likelihood, prior_transform, ndim=self.ndim, \
                                                    bound='multi',pool=pool,**sampler_kwargs)
            
            self.sampler.run_nested(print_progress=verbose,**run_kwargs)

            self.results = self.sampler.results
            samples = self.results.samples
            weights = np.exp(self.results.logwt - self.results.logz[-1])

            #self.quantiles = np.array([dyfunc.quantile(s, [0.16, 0.5, 0.84], weights=self.weights) for s in self.samples.T]).T
                        
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
            
            samples, weights, _ = self.sampler.posterior()
            weights = np.exp(weights-weights.max())
            self.results = None
            
        else: 
            raise ValueError(f"ERROR! Unknown method '{method}'.")
        
        dt = time.time()-toc
        dt = f'{dt:.2f} seconds' if dt<60.00 else f'{dt/60.00:.2f} minutes' 

        vprint(verbose,f'\nSampling with {method} took {dt} to run.')
        
        if pool is not None:
            pool.close()
            if not isMPI:
                pool.join()
        
        self.modCalculated = True

        # Calculating sampling statistics and setting class attributes
        self._set_sampling_stats(samples,weights)

        # Creating a new Rings object for the best-fit rings.
        self.outri = Rings(self._inri.nr)
        self.outri.set_rings_from_dict(self._inri.r)
        
        # Updating output rings with best parameters from the sampling
        self.outri = self._update_rings(self.outri,self.params)

        # Setting up the output rings in the Galfit object.
        # @TODO: Add support for errors in the rings (already in the C++ code)
        libBB.Galfit_setOutRings(self._galfit,self.outri._rings)
 
    
    def compute_geometry(self,freepar=['xpos','ypos','inc','phi'],**kwargs):
        """ A wrapper to compute only the geometry of the model using the 0th moment map (Cannubi-style) """

        allowed = {'xpos_single','ypos_single','inc_single','phi_single','z0_single','dens','dens_single','dens_func','radmax'}
        if not (freepar and set(freepar).issubset(allowed)): 
            raise ValueError("ERROR! The function compute_geometry only accepts the following free parameters: ", allowed)

        self._compute(freepar=freepar,log_like=self._log_likelihood_geometry,**kwargs)


    def _set_sampling_stats(self,samples,weights):
        """ Calculate statistics of the sampling results and set class attributes. """
        self.samples = samples
        self.weights = weights / np.nansum(weights)
        self.samples_equal = dyfunc.resample_equal(self.samples, self.weights)
        self.modCalculated = True
        
        # Below only ok with Numpy > 2.0, but nautilus currently does not support it. So using equal weighting samples.
        #self.quantiles = np.quantile(self.samples, [0.16,0.50,0.84], weights=self.weights, axis=0, method='inverted_cdf')        
        self.quantiles = np.quantile(self.samples_equal, [0.16,0.50,0.84], axis=0)
        self.params = self.quantiles[1]


    def _update_rings(self,rings,theta):
        """ Update rings with the parameters theta """

        freepar_idx = copy.copy(self.freepar_idx)

        # Treating the radmax parameter
        if 'radmax' in freepar_idx:
            radmax = theta[freepar_idx['radmax']][0]
            #I only adjust the latest ring
            rings.r['radii'][-1] = radmax
            #rings.modify_parameter('radii',radii)

        # Treating cases with function forms
        whatFunctional = np.unique([key.split('_')[0] for key in freepar_idx if '_p' in key])
        if len(whatFunctional)>0:
            for p in whatFunctional:
                keys = [key for key in freepar_idx if p in key]
                fvalue = [theta[freepar_idx[key]] for key in keys]
                pvalue = self.funcs[p](rings.r['radii'],*fvalue)
                rings.modify_parameter(p,pvalue)
                for key in keys: freepar_idx.pop(key)

        # Now treating the rest of the parameters with regular rings
        for key in freepar_idx:
            pvalue = theta[freepar_idx[key]]
            if key=='sigma' or key=='radmax': continue
            if len(pvalue)==1: pvalue = pvalue[0]
            rings.modify_parameter(key,pvalue)
        
        rings.make_object()
        return rings


    def _update_profile(self,rings):
        """ Get the density profile from the rings using Ellprof"""
        libBB.Ellprof_update_rings(self._ellprof,rings._rings)
        libBB.Ellprof_compute(self._ellprof)
        dens = reshapePointer(libBB.Ellprof_dens_array(self._ellprof),(1,self._inri.nr))[0]

        #mindens = np.nanmin(dens[dens>1E-10])
        #dens *= 1./mindens
        rings.modify_parameter("dens",np.abs(dens),makeobj=True)

        
    def write_bestmodel(self,plots=True,**kwargs):
        
        if not self.modCalculated:
            print ("Sampler has not been run yet. Please run compute() before running this function.")
            return
        
        verbose = kwargs.get('verbose',True)
        vprint(verbose,"\nWriting the best model to output directory ...",end=' ',flush=True)
        
        if self.useBBres:
            # Writing the model and plots to the output directory
            libBB.Galfit_writeModel(self._galfit,"AZIM".encode('utf-8'),plots)
        
        else:
            
            # If we are not fitting the dens, we need to update the final profile
            if self.useNorm:
                self._update_profile(self.outri)

            # Deriving the last model
            model, _, _, galmod = self._calculate_model(self.outri,fullcube=True)

            if self.useNorm:
                # Normalizing and copying it back to the C++ Galmod object
                model = self._normalize_model(model,self.data,self.mask)
                libBB.Galmod_set_array(galmod,np.ravel(model).astype('float32'))
                model = reshapePointer(libBB.Galmod_array(galmod),self.data.shape)
            
            # Writing all the outputs
            libBB.Galfit_writeOutputs(self._galfit,galmod,self._ellprof,plots)
            #libBB.Galmod_delete(galmod)

        vprint(verbose,"Done!")


    def corner_plot(self,samples=None,labels=None,truths=None,quantiles=[0.15865,0.50,0.84135],plot_dynesty=False,**kwargs):
        """ Create a corner plot of the posterior distributions """
        
        if samples is None and not self.modCalculated:
            print ("Sampler has not been run yet. Please run compute() before running this function or load samples with load_samples().")
            return

        defaults  = dict(max_n_ticks=5, color='darkblue', show_titles=True)

        if self.results is None:
            plot_dynesty = False

        if plot_dynesty:
            try: 
                from dynesty import plotting as dyplot
            except ImportError:
                raise ImportError("BayesianBB requires the package 'dynesty' to produce dynesty plots.")

            # Dynesty corner plot needs the results object
            samples = self.results
            cornerplot = dyplot.cornerplot

        else:
            try: 
                import corner
            except ImportError:
                raise ImportError("BayesianBB requires the package 'corner' to produce corner plots.")

            cornerplot = corner.corner
            defaults['plot_datapoints'] = False

        if samples is None:
            samples = self.samples_equal
        if labels is None:
            labels = self.freepar_names
        
        mergedkey = {**defaults, **kwargs}

        c = cornerplot(samples, title_quantiles=quantiles, quantiles=quantiles, \
                       labels=labels, label_kwargs=dict(fontsize=20), **mergedkey)
        
        return c[0] if plot_dynesty else c


    def plot_profiles(self,**kwargs):
        
        if not self.modCalculated:
            print ("Sampler has not been run yet. Please run compute() before running this function or load samples with load_samples().")
            return

        rings = self.outri.r
        toplot = ['vrot','vdisp','dens','inc','phi','z0','xpos','ypos','vsys']
        labels = [r'$V_\mathrm{rot}$ (km/s)',r'$\sigma$ (km/s)','Density','Inclination (deg)','Position angle (deg)',\
                  r'$z_0$ (arcsec)','Xpos (pix)','Ypos (pix)',r'$V_\mathrm{sys}$ (km/s)']

        # Creating the figure and axes
        fig = plt.figure(figsize=(12,12))
        nrows, ncols = 3, 3
        xlen, ylen, x_sep, y_sep = 0.27, 0.13, 0.07, 0.015
        for i in range(nrows):
            yl = ylen * (1.6 if i == 0 else 1.0)
            y = 0.9 - i * (ylen + y_sep)
            for j in range(ncols):
                x = 0.1 + j * (xlen + x_sep)
                fig.add_axes([x, y, xlen, yl])
                # Some formatting of axes
                isbottom = i == nrows - 1
                fig.axes[-1].tick_params(direction='in',top=True,right=True,labelbottom=isbottom)
                if isbottom:
                    fig.axes[-1].set_xlabel('Radius (arcsec)', fontsize=14)
                fig.axes[-1].set_xlim(0, rings['radii'][-1]*1.05)

        ax = np.ravel(fig.axes) 
        for i in range(len(toplot)):
            a = toplot[i]
            q = self.quantiles[:,self.freepar_idx[a]]
            err_low, err_up = q[1]-q[0], q[2]-q[1]
            ax[i].errorbar(rings['radii'],rings[a],yerr=[err_low, err_up],fmt='o',c='#B22222')    
            ax[i].set_ylabel(labels[i])
            ax[i].set_ylim(0.9*np.nanmin(rings[a]-err_low), 1.1*np.nanmax(rings[a]+err_up))

            if a=='vrot':
                ax[i].set_ylim(0, None)
        
        return fig


    def print_stats(self,quantiles=[0.15865,0.50,0.84135]):

        if not self.modCalculated:
            print ("Sampler has not been run yet. Please run compute() before running this function.")
            return
        
        if not isIterable(quantiles) and len(quantiles)!=3:
            raise ValueError("ERROR! Percentiles must be a list of three values.")

        results = self.results
        samples = self.samples_equal

        # Printing central values and error of posterior distributions
        logl = self._log_likelihood(self.params)
        print(f"\nBest-fit parameters (log_likelihood = {logl:.2f}):") 
        for i in range(len(self.freepar_names)):
            values = self.quantiles.T[i]
            q = np.diff(values)
            txt = "  %10s = %10.3f %+10.3f %+10.3f"%(self.freepar_names[i], values[1], -q[0], q[1])
            print (txt)

        # Printing the maximum likelihood
        if results is not None:
            max_logl   = results.logl.max()
            max_index  = results.logl.argmax()
            max_params = samples[max_index]
            print(f"\nMaximum logl: {max_logl} at:")
            for i in range(len(self.freepar_names)):
                txt = "  %10s = %10.3f "%(self.freepar_names[i], max_params[i])
                print(txt)


    def print_priors(self):
        """ Print the priors for each parameter """
        print ("\nPrior distributions for free parameters:")
        for key in self.freepar_idx:
            print(f"{key:>10s} = {self.prior_distr[key].dist.name:10s}", self.prior_distr[key].kwds)
        print ("\n")
    

    def save_results(self,outfile):
        """ Save the results as a dictionary in a pickle file."""

        tosave = dict()
        tosave['results'] = self.results if self.results is not None else None
        tosave['samples'] = self.samples
        tosave['freepar_inp'] = self._freepar_inp
        tosave['freepar_names'] = self.freepar_names
        tosave['freepar_idx']   = self.freepar_idx
        tosave['inrings'] = self._inri.r
        tosave['options'] = self._opts.opts
        np.save(outfile,tosave)

 
    def load_results(self,inputfile,sampler='dynesty'):
        """ Load previously saved results from a file. Input file must be a numpy .npy file."""
        
        r = np.load(inputfile,allow_pickle=True).item()
        
        # Initializing rings
        rr = r['inrings']
        self.init(radii=rr['radii'],xpos=rr['xpos'],ypos=rr['ypos'],vsys=rr['vsys'],vrot=rr['vrot'],\
                  vdisp=rr['vdisp'],inc=rr['inc'],phi=rr['phi'],z0=rr['z0'],vrad=rr['vrad'],\
                  dens=rr['dens'],vvert=rr['vvert'],dvdz=rr['dvdz'],zcyl=rr['zcyl'])
        
        # Setting options and all necessary attributes
        self.set_options(**r['options'])
        self._setup(r['freepar_inp'])

        # Loading sampling results
        if sampler=='dynesty':
            self.results = r['results']
            weights = np.exp(self.results.logwt - self.results.logz[-1])
            self._set_sampling_stats(self.results.samples, weights)
            self.outri = Rings(self._inri.nr)
            self.outri.set_rings_from_dict(self._inri.r)
            self.outri = self._update_rings(self.outri,self.params)

        elif sampler=='nautilus':
            raise NotImplementedError()
        else:
            raise ValueError(f"Error: sampler {sampler} not known")
