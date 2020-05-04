
Quickstart
#################

Running a task
===============

In pyBBarolo, BBarolo's task are wrapped as python classes. The generic procedure to run a task is as follow:

1. Import the task from pyBBarolo module::

    from pyBBarolo import Task
    
  where Task is one of pyBBarolo :ref:`tasks <pybbtasks>`.


2. Create an object of the task class::
    
    bb = Task(fitsname)
    
  where *fitsname* is the input FITS file. All tasks need an initial FITS file.
  
3. Initialize the task::

    bb.init(args)
    
  where *args* are required arguments that depend on the various tasks. Required arguments can be printed with::
  
    bb.show_arguments()

4. Set options::

    bb.set_options(opts)
    
  where *opts* are the task-dependent available options. All options have default values, so this step is not mandatory. To see a list of available options::
  
    bb.show_options()
    
5. Run the task::
    
    bb.compute()


.. _pybbtasks:

Available tasks
================

* **FitMod3D()**: wrapped class for BBarolo's :ref:`3DFIT <3dfit>` task.

* **GalMod()**: wrapped class for BBarolo's :ref:`GALMOD <galmodtask>` task.

* **Search()**: wrapped class for BBarolo's :ref:`SEARCH <searchtask>` task.

* **FitMod2D()**: wrapped class for BBarolo's :ref:`2DFIT <2dfit>` task.

* **Ellprof()**: wrapped class for BBarolo's :ref:`ELLPROF <ellproftask>` task.

* **SpectralSmooth()**: wrapped class for BBarolo's :ref:`SMOOTHSPEC <spectralsmoothtask>` task.

 
Example 1: 3D fit of a galaxy
=============================

Suppose you have an astonishing observation of your favorite galaxy and you want to fit a 3D kinematic model to your emission line datacube. 

Let's fit the HI datacube of NGC 2403, which is available as a part of BBarolo's working `examples <http://editeodoro.github.io/Bbarolo/downloads/examples/>`_. 
We just have to set initial conditions, options and then run the task. 

First of all, we import and start FitMod3D::

    from pyBBarolo import FitMod3D

    # FITS file of the galaxy to model
    filen = "./examples/ngc2403.fits"
    # Initializing a 3DFIT object
    f3d = FitMod3D(filen)
    
Secondly, we initialize rings with initial guesses for the fit::

    # Initializing rings. Parameters can be values or arrays
    f3d.init(radii=np.arange(15,1200,30),xpos=77,ypos=77,vsys=132.8,vrot=120,vdisp=8,vrad=0,z0=10,inc=60,phi=123.7)

A list of needed arguments can be printed with ``f3d.show_arguments()``.

Thirdly, we can change some default options for the fit. For a list of available options: ``f3d.show_options()``.

For instance, we can set a mask made through the source-finding algorithm (``mask="SEARCH"``), parameters to fit (``free="VROT VDISP"``), the distance of the galaxy in Mpc (``distance=3.2``) and the directory for outputs (``outfolder='output/ngc2403``)::

    f3d.set_options(mask="SEARCH",free="VROT VDISP",wfunc=2,distance=3.2,outfolder='output/ngc2403')

If the beam information is not available in the FITS header, it is fundamental to set the size of the beam::

    f3d.set_beam(bmaj=60,bmin=60,bpa=-80)
    
It is now time to run the fit::
    
    bfrings, bestmod = f3d.compute(threads=4) 

This function performs the fit and writes relevant FITS files in the output directory. The function returns a *n* x *m* matrix containing the best-fit rings (``bfrings``), where *n* = number of rings and *m* = number of parameters, and a FITS astropy object (``bestmod``) containing the best-fit model. These are also written in ``bfit`` and ``outmodel`` methods of FitMod3D class. 

Finally, we can use BBarolo built-in routines to write plots of data and model, like channel maps, moment maps, position-velocity diagrams and best-fit parameters::

    f3d.plot_model()




Example 2: 3D model of a galaxy
===============================

It is also possible to simply build a 3D model datacube from given parameters. This is accomplished with the GalMod task. The procedure is similar to the one above::

    from pyBBarolo import GalMod

    # FITS file of the galaxy to model
    filen = "./examples/ngc2403.fits"
    # Initializing a GalMod object
    gm = GalMod(filen)
    # Initializing rings. Parameters can be values or arrays
    gm.init(radii=np.arange(15,1200,30),xpos=74,ypos=74,vsys=132.8,vrot=120,vrad=10,vvert=5,vdisp=8,z0=10,inc=60,phi=123.7)
    # Now, let's take a look to the default options (see BB documentation)
    gm.show_options()
    # Changing some options
    gm.set_options(ltype=1)
    # Compute the model
    mymodel = gm.compute()
    # Smooth to the same resolution of data
    mymodel = gm.smooth()
    # mymodel is an astropy cube and we can do whatever we like with it.
    mymodel.writeto("awesome_model.fits",overwrite=True)
        


Example 3: All tasks 
================================
This :download:`script <pyBB_test.py>` shows how to use all wrapped classes:

.. literalinclude:: pyBB_test.py
   :language: python


.. Example 3: 3D model of an outflow
   =================================





