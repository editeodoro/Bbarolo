.. _general:

General parameters
=======================

In the following, a list of general parameters (e.g.\, not task-specific).

Input/output parameters
^^^^^^^^^^^^^^^^^^^^^^^

* **FITSFILE** [none]. The name of the input FITS file. This is a mandatory parameter for **all tasks**.

* **OUTFOLDER** [./output]. The directory where the output files will be written. 

* **VERBOSE** [true]. Enable all the output messages.

* **THREADS** [1]. Number of CPUs to use for task execution. All BBarolo's tasks have shared-memory parallelization. The code needs to be compiled with OpenMP support. If you encounter any problem with multi-thread execution, switch back to single-thread mode and signal the problem.

* **PLOTS** [true]. If true, output plots will be produced (Python/Gnuplot needed).

* **SHOWBAR** [true]. Whether to show progress bars.

.. _beam:

Beam parameters
^^^^^^^^^^^^^^^^^^^^^^^

Following parameters can be used to specify the size and shape of the Point Spread Function (PSF or beam). These parameters are ignored if beam information is written in the header of the input FITS, either through BMAJ, BMIN and BPA keywords or in the HISTORY. The code defines the beam following the priority order: header -> bmaj,bmin,bpa params -> beamfwhm param -> default to 30 arcsec.

* **BMAJ** [none]. The FWHM of the major axis of the elliptical Gaussian beam in *arcsec*. 

* **BMIN** [none]. The FWHM of the minor minor axis of the elliptical Gaussian beam in *arcsec*.

* **BPA** [none]. The position angle of the major axis of the elliptical Gaussian beam in *degrees*, counter-clock from the North direction.

* **BEAMFWHM** [none]. The FWHM of a circular Gaussian beam in *arcsec* . 
