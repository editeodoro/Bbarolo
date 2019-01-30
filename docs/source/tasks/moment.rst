
.. _momentstask:

Moment maps and position-velocity cuts 
######################################

BBarolo can be used to extract global profiles, moment maps and position velocity diagrams. For moment maps and profile, the input datacube can be masked using the MASK parameter (see :ref:`3DFIT <3dfitopt>`).

Parameters
===========

* **GLOBALPROFILE** [false]. It *true*, calculate the total line profile from a datacube and write it to a text file.

* **TOTALMAP** [false]. It *true*, calculate the total intensity map from a datacube and write it to a FITS file. 

* **VELOCITYMAP** [false]. It *true*, calculate the velocity field from a datacube and write it to a FITS file. 

* **DISPERSIONMAP** [false]. It *true*, calculate the velocity dispersion field from a datacube and write it to a FITS file. 

* **MAPTYPE** [MOMENT]. Available only for v1.4.2+. It specifies the way the kinematic maps are derived. Can be either *MOMENT* (classical moments) or *GAUSSIAN* (gaussian fit).

* **RMSMAP** [false]. It *true*, calculate the RMS map, i.e. the RMS in each spectrum, from a datacube and write it to a FITS file. The RMS is calculated in an iterative way. RMS is the standard deviation for normal statistics and MADFM/0.6745 for robust statistics (**FLAGROBUSTSTATS** parameter).

* **MASSDENSMAP** [false]. It *true*, calculate a mass surface-density map in units of Msun/pc^2 from a datacube and write it to a FITS file. This is just for HI data and the input datacube is required to have JY/BEAM flux density units.

* **FLAGPV** [false]. If *true*, extract position-velocity image from a datacube and write it to a FITS file. The cut is defined by a point and an angle, as set with the following parameters.

* **XPOS_PV** [none]. Reference X pixel of the cut.

* **YPOS_PV** [none]. Reference Y pixel of the cut.

* **PA_PV** [none]. Position angle of the cut, defined anti-clockwise from the X direction. 


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


