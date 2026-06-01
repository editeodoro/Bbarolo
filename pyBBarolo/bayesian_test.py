from .bayesian import *

class BayesianBBarolo_test(BayesianBBarolo):
    
    def __init__(self, fitsname):
        super(BayesianBBarolo_test,self).__init__(fitsname=fitsname)


    def _calculate_model(self,rings):
        # Model is built in a smaller array (blo,bhi) that only contains the galaxy 
        # within the last ring. This is to calculate the residuals faster.
        bhi, blo = (ctypes.c_int * 2)(0), (ctypes.c_int * 2)(0)
        libBB.Galfit_getModelSize(self._galfit,rings._rings,bhi,blo)

        ############################################################
        # Below is the call to the new function that returns unconvolved and convolved model
        galmod = libBB.Galfit_getModel_BBB(self._galfit,rings._rings,bhi,blo,-1)
        mod_shape = (2, self.inp.dim[2], bhi[1]-blo[1],bhi[0]-blo[0])
        mod = reshapePointer(galmod,mod_shape)
        ############################################################
        
        return mod, bhi, blo
    

    def _log_likelihood(self,theta,**kwargs):
        """ Likelihood function for the fit """
        
        rings = self._update_rings(self._inri,theta)

        # Recompute the density profile along the current rings and update the rings
        if self.useNorm and self.update_prof:
            self._update_profile(rings)

        # Calculate the model and the boundaries
        model, bhi, blo = self._calculate_model(rings)
            
        # Calculate the residuals
        mask = self.mask[:,blo[1]:bhi[1],blo[0]:bhi[0]]
        data = self.data[:,blo[1]:bhi[1],blo[0]:bhi[0]]

        kwargs['sigma'] = theta[self.freepar_idx['sigma']][0] if 'sigma' in self.freepar_names else 1
            
        return -self._calculate_residuals(model,data,mask,**kwargs)