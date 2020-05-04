import os
import numpy as np
from pyBBarolo import *

if __name__ == '__main__':

    fitsname = "./examples/ngc2403.fits"

    #'''
    ### 3D FIT TUTORIAL ############################################################
    f3d = FitMod3D(fitsname)
    f3d.init(radii=np.arange(15,450,30),xpos=77,ypos=77,vsys=132.8,\
             vrot=120,vdisp=8,vrad=0,z0=10,inc=60,phi=123.7)
    f3d.show_options()
    f3d.set_options(mask="SEARCH",free="VROT VDISP",wfunc=2,distance=3.2,ltype=2)
    f3d.set_options(outfolder="output/ngc2403")
    f3d.set_beam(bmaj=60,bmin=60,bpa=-80)
    bfrings, bestmod = f3d.compute(threads=4) 
    f3d.plot_model()
    ################################################################################
    #'''

    #'''
    ### GALMOD TUTORIAL ############################################################
    gm = GalMod(fitsname)
    gm.show_arguments()
    gm.init(radii=np.arange(15,1200,30),xpos=74,ypos=74,vsys=132.8,\
             vrot=120,vdisp=8,z0=10,inc=60,phi=123.7)
    gm.show_options()
    gm.set_options(ltype=1)
    mymodel = gm.compute()
    mymodel = gm.smooth()
    mymodel.writeto("awesome_model.fits",overwrite=True)
    ################################################################################
    #'''

    #'''
    ### ELLPROF TUTORIAL #############################################################
    el = Ellprof(fitsname)
    el.show_arguments()
    el.init(radii=np.arange(15,1200,30),xpos=74,ypos=74,inc=60,phi=123.7)
    el.show_options()
    el.set_options(mask="SMOOTH")
    rings = el.compute()   
    print (rings['rad'],rings['msurfdens']) 
    el.writeto("rings_ellprof.txt")
    ################################################################################
    #'''

    #'''
    ### 2DFIT TUTORIAL #############################################################
    rm = FitMod2D(fitsname)
    rm.show_arguments()
    rm.init(radii=np.arange(15,1200,30),xpos=74,ypos=74,vsys=132.8,vrot=120,inc=60,phi=123.7)
    rm.show_options()
    rm.set_options(wfunc=1, free="VROT INC")
    rings = rm.compute(threads=1)   
    print (rings['rad'],rings['vrot']) 
    rm.writeto("rings_2dfit.txt")
    ################################################################################
    #'''
    
    #'''
    ### SEARCH TUTORIAL ############################################################
    s = Search(fitsname)
    s.set_options(growth=True)
    s.show_options()
    s.search(threads=4)    
    ################################################################################
    #'''
    
    #'''
    ### SPECTRAL SMOOTHING TUTORIAL ################################################
    s = SpectralSmoothing(fitsname)
    a = s.smooth(window_type='hanning',window_size=11,threads=8)
    s.writeto("smoothedcube.fits",average=True)
    ################################################################################
    #'''
    