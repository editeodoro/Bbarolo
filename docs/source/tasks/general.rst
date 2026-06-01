.. _general:

General parameters
=======================

In the following, a list of general parameters (e.g.\, not task-specific).

Input/output parameters
^^^^^^^^^^^^^^^^^^^^^^^

* **FITSFILE** [none]. The name of the input FITS file. This is a mandatory parameter for **all tasks**.

* **OUTFOLDER** [none]. The directory where the output files will be written. If not set, the output directory will default to "./output/OBJECT", where OBJECT is the corresponding keyword in the FITS header of the **FITSFILE**. 

* **OUTPREFIX** [none]. Prefix for output files written in the **OUTFOLDER**. If not set, the output prefix will default to "OBJECT", where OBJECT is the corresponding keyword in the FITS header of the **FITSFILE**. 

* **VERBOSE** [true]. Enable all the output messages.

* **THREADS** [max CPUs]. Number of CPUs to use for task execution. All BBarolo's tasks have shared-memory parallelization. The code needs to be compiled with OpenMP support. If you encounter any problem with multi-thread execution, switch back to single-thread mode and signal the problem.

* **PLOTS** [true]. If true, output plots will be produced (Python/Gnuplot needed). 

* **SHOWBAR** [true]. Whether to show progress bars.

* **STATS** [false]. If true, calculate and print statistics of input FITS file.

* **CHECKCUBE** [0]. Check the datacubes for bad data in either channels, columns or rows. Replace the CHECKCHANNELS parameter in previous versions. Allowed values are 0 = no check performed; 1 = channels and spatial rows and columns; 2 = only channels; 3 = rows and columns; 4 = only columns; 5 = only rows.


.. _beam:

Beam parameters
^^^^^^^^^^^^^^^^^^^^^^^

Following parameters can be used to specify the size and shape of the Point Spread Function (PSF or beam). These parameters are ignored if beam information is written in the header of the input FITS, either through BMAJ, BMIN and BPA keywords or in the HISTORY. The code defines the beam following the priority order: header -> bmaj,bmin,bpa params -> beamfwhm param -> default to 30 arcsec.

* **BMAJ** [none]. The FWHM of the major axis of the elliptical Gaussian beam in *arcsec*. 

* **BMIN** [none]. The FWHM of the minor minor axis of the elliptical Gaussian beam in *arcsec*.

* **BPA** [none]. The position angle of the major axis of the elliptical Gaussian beam in *degrees*, counter-clock from the North direction.

* **BEAMFWHM** [none]. The FWHM of a circular Gaussian beam in *arcsec* . 

* **FLUXCONVERT** [true]. Whether to try to convert all flux units to either Jy or Jy/beam in moment maps and in the source finder. 

