.. _ellproftask:

ELLPROF task
#############
This task can be used to calculate the radial density profile of a galaxy.  

Parameters
==========

* **ELLPROF** [false]. This flag enables the radial profile task.

Parameters for the task are: **RADII**, **NRADII**, **RADSEP**, **XPOS**, **YPOS**, **PA**, **INC**, **SIDE** (see :ref:`3DFIT <3dfit>`). If **FITSFILE** is a datacube, the profile is calculated from the column density map calculated after masking the cube accordingly to the **MASK** parameter. If **FITSFILE** is a 2D intensity map, this is used to extract the profile.

Outputs
========

The task produces the following outputs. *NAME* is the name of the galaxy.

* A FITS file *NAME_densmap.fits* with the map used to extract the profile.

* A text file *NAME_densprof.txt*, containing with the radial profiles.


Example
========
An example :download:`parameter <examples/n2403_ellprof.par>` file to extract the radial profile from the usual NGC2403 :download:`datacube <examples/ngc2403.fits>`, given a set of rings. 

.. literalinclude:: examples/n2403_ellprof.par
   :language: c