.. _galmodtask:

GALMOD task
###########

GALMOD is the routine underlying the 3DFIT task. It builds a 3D simulated datacube of a disk galaxy starting from a set of concentric rings with given column density and kinematics. The routine is an updated version of the namesake routine in GIPSY (see also GIPSY `GALMOD <https://www.astro.rug.nl/~gipsy/tsk/galmod.dc1>`_). 

Parameters
==========

Parameters for rings are the same of the :ref:`3DFIT task <ringio>`. Options are LTYPE, CDENS, NV and VELDEF (see :ref:`3DFIT options <3dfitopt_add>`).

Additional GALMOD-specific parameters are:

* **GALMOD** [false]. This flag enables the 3D disk modelling. Can be *true* or *false*.

* **VVERT** [0]. Vertical velocity in *km/s*. 

* **DVDZ** [0]. Gradient of rotation velocity as we move away from the disk plane. This is in *km/s/arcs*.

* **ZCYL** [0]. Height in *arcsec* from the disk plane where the gradient DVDZ begins.

* **SM** [true]. Whether to smooth the model to the same spatial resolution of data.

Outputs
========

The task produces the following outputs. Here *NAME* is the name of the galaxy and *NORM* is the kind normalization used.

* A FITS file *NAMEmod_NORM.fits*, containing the model datacube.

* A FITS file *mask.fits*, containing the mask used. 

* FITS files of position-velocity cuts taken along the average major and minor axes for the input datacube and the model. In particular:

  * *NAME_pv_a.fits*: P-V of the data along the major axis.
  * *NAME_pv_b.fits*: P-V of the data along the minor axis.
  * *NAMEmod_pv_a_NORM.fits*: P-V of the model along the major axis.
  * *NAMEmod_pv_b_NORM.fits*: P-V of the model along the minor axis.

* FITS files of the moment maps for the input data and the model. These can be found in the *maps* subdirectory:
  
  * *NAME_0mom.fits*, *NAME_1mom.fits*, *NAME_2mom.fits:* 0th, 1st and 2nd moment maps of the data.
  * *NAME_NORM_0mom.fits*, *NAME_NORM_1mom.fits*, *NAME_NORM_2mom.fits:* 0th, 1st and 2nd moment maps of the model.

* A text file *densprof.txt*, with the radial intensity profiles along the chosen rings.


Example
========
Above outputs can be obtained with the following :download:`parameter <examples/n2403_galmod.par>` file and the usual example :download:`datacube <examples/ngc2403.fits>`.

.. literalinclude:: examples/n2403_galmod.par
   :language: c