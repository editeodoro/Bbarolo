.. _smoothtask:

SMOOTH task
############

This task convolves each channel map in a datacube with a given elliptical Gaussian. 

Parameters
==========

* **SMOOTH** [false]. This flag enables the smooth algorithm. Can be *true* or *false*.

* **OBMAJ** [none]. Major axis of the initial beam in *arcsec*. Do not set if you want to use the beam information in the input FITS file (the parameter overrides it).

* **OBMIN** [none]. Minor axis of the initial beam in *arcsec*. Do not set if you want to use the beam information in the input FITS file (the parameter overrides it).

* **OBPA** [none]. Position angle of the major axis of the initial beam in *degrees*. Do not set if you want to use the beam information in the input FITS file (the parameter overrides it).

* **BMAJ** [none]. Major axis of the final beam in *arcsec*.

* **BMIN** [none]. Minor axis of the final beam in *arcsec*.

* **BPA** [none]. Position angle of the major axis of the final beam in *degrees*.

* **FACTOR** [2]. If set, the beam of the output cube is [FACTOR\*OBMAJ,FACTOR\*OBMIN,OBPA]. Ignored if BMAJ, BMIN, BPA are specified.

* **SCALEFACTOR** [none]. Scaling factor for output datacube. BBarolo will calculate an appropriate one if left unset.

* **FFT** [true]. Whether to convolve by using Fast Fourier Transform or not.

* **REDUCE** [false]. If *true*, BBarolo repixels the output datacube to preserve the number of pixels in a beam.

* **SMOOTHOUTPUT** [none]. Output smoothed FITS file. Default is input file name with a suffix indicating the new beam size.

Outputs
========

The task writes the smoothed datacubes in the FITS file *NAME_sN.fits*, where *NAME* is the name of the galaxy and *N* is the new beam size.


Example
========
Below, an example :download:`parameter <examples/n2403_smooth.par>` file to smooth the usual :download:`datacube <examples/ngc2403.fits>` to a coarser spatial resolution.

.. literalinclude:: examples/n2403_smooth.par
   :language: c

|

|

SMOOTHSPEC task (v1.5.2+)
#########################
This task convolves each spectrum in a datacube with a given window, i.e. it performs spectral smoothing. 

Parameters
==========
* **SMOOTHSPEC** [false]. This flag enables the spectral smoothing algorithm. Can be *true* or *false*.

* **WINDOW_TYPE** [HANNING]. Type of the smoothing window. Implemented windows include *HANNING*, *BOXCAR*, *BARTLETT*, *WELCH*, *BLACKMAN*, *FLATTOP*. 

* **WINDOW_SIZE** [3]. Size of the smoothing window in channels. 

* **REDUCE** [false]. If *true*, output data will be averaged over WINDOW_SIZE channels.

Outputs
========

The task writes the smoothed datacube in the FITS file *NAME_hN.fits*, where *NAME* is the name of the galaxy and *N* is the size of the Hanning window.

Example
========
Below, a :download:`parameter <examples/n2403_hanning.par>` file to Hanning smooth the usual :download:`datacube <examples/ngc2403.fits>`.

.. literalinclude:: examples/n2403_hanning.par
   :language: c