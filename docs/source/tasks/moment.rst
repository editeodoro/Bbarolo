
.. _momentstask:

Moment maps and position-velocity cuts 
######################################

BBarolo can be used to extract global profiles, moment maps and position velocity diagrams. For moment maps and profile, the input datacube can be masked using the MASK parameter (see :ref:`3DFIT <3dfitopt>`).

Parameters for maps
===================

* **GLOBALPROFILE** [false]. It *true*, calculate the total line profile from a datacube and write it to a text file.

* **TOTALMAP** [false]. It *true*, calculate the total intensity map from a datacube and write it to a FITS file. 

* **VELOCITYMAP** [false]. It *true*, calculate the velocity field from a datacube and write it to a FITS file. 

* **DISPERSIONMAP** [false]. It *true*, calculate the velocity dispersion field from a datacube and write it to a FITS file. 

* **MAPTYPE** [MOMENT]. It specifies the way the kinematic maps are derived. Can be either *MOMENT* (classical moments) or *GAUSSIAN* (gaussian fit).

* **RMSMAP** [false]. It *true*, calculate the RMS map, i.e. the RMS in each spectrum, from a datacube and write it to a FITS file. The RMS is calculated in an iterative way. RMS is the standard deviation for normal statistics and MADFM/0.6745 for robust statistics (**FLAGROBUSTSTATS** parameter).

* **MASSDENSMAP** [false]. It *true*, calculate a mass surface-density map in units of Msun/pc^2 from a datacube and write it to a FITS file. This is just for HI data and the input datacube is required to have JY/BEAM flux density units.

Parameters for PV slices
========================

* **PVSLICE** [false]. If *true*, extract a position-velocity image from a datacube and write it to a FITS file. The slice can be defined by either a point and an angle (see **XPOS_PV**, **YPOS_PV**, **PA_PV** below), or by two points (see **P1_PV**, **P2_PV** below). The former has priority over the latter.

* **XPOS_PV** [none]. Reference X of the slice. Can be a pixel or a coordinate (see also **XPOS** for :ref:`3DFIT <ringio>`).

* **YPOS_PV** [none]. Reference X of the slice. Can be a pixel or a coordinate (see also **YPOS** for :ref:`3DFIT <ringio>`).

* **P1_PV** [none]. X and Y coordinates of the first point defining a slice. For example: '10 30' for a pixel at X=10 and Y=30. 

* **P2_PV** [none]. X and Y coordinates of the second point defining a slice. 

* **WIDTH_PV** [0]. Width of the slice in arcsec. The PV cut will be averaged over a rectangular window extending from -width/2 to +width/2 from the slice defined above.

* **ANTIALIAS** [0.5]. It defines how many pixels to use for antialiasing algorithm. It can be an integer or 0.5. The number of pixels used for antialiasing will be (1+2*ANTIALIAS)^2. Set it to 0 for no antialiasing.  


Outputs
========

The required map/P-V is written in a FITS file.

Example
=======
The following :download:`parameter <examples/n2403_moment.par>` file extract maps and PV diagrams from the usual :download:`datacube <examples/ngc2403.fits>`.

.. literalinclude:: examples/n2403_moment.par
   :language: c


.. figure:: examples/moment.png
   :alt: Maps 
   :scale: 80%
   :align: center
   
   Moment maps, RMS map and P-Vs along the major and minor axis


