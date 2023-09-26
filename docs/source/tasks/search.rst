.. _searchtask:

SEARCH task
############

BBarolo's search algorithm is derived from `Duchamp <https://www.atnf.csiro.au/people/Matthew.Whiting/Duchamp>`_, a 3D source finder for spectral-line data developed by `Matthew Whiting <https://www.atnf.csiro.au/people/Matthew.Whiting/>`_. BBarolo adds a few new functionalities and parallelization. For a comprehensive description of the algorithm and the input parameters, see Duchamp's `main paper <http://adsabs.harvard.edu/abs/2012MNRAS.421.3242W>`_ and `user guide <http://www.atnf.csiro.au/people/Matthew.Whiting/Duchamp/downloads/UserGuide-1.6.1.pdf>`_.

Parameters
==========

* **SEARCH** [false]. This flag enables the source finding algorithm. Can be *true* or *false*.

* **FLAGROBUSTSTATS** [true]. Whether to use to robust estimators (median and MADFM) instead of normal estimators (mean and standard deviation) when calculating cube statistics.

* **ITERNOISE** [false]. Whether to use an iterative algorithm to estimate the noise level. If true, it will reiterate over the array, masking pixels above 3sigma and re-calculating noise statistics until convergence. 

* **SORTSOURCES** [NVOX]. This parameter specifies how to sort detections. Accepted values are *XVALUE*, *YVALUE*, *ZVALUE*, *RA*, *DEC*, *VEL*, *IFLUX* (integrated flux), *PFLUX* (peak flux), *NPIX* (number of pixels), *NVOX* (number of voxels), *W20* (width at 20% peak flux), *W50* (width at 50% peak flux) and *SNR* (average signal-to-noise). By default, sorting is descending, but it can be done ascending by adding a *-* (e.g. *-VEL*).

* **CUBELETS** [false]. If true, it produces individual cubelets and sub-images for each detected source.

* **SEARCHTYPE** [spatial]. How the search is performed. Accepted values are *spatial* and *spectral*. Spatial search is done in 2D channel maps, spectral search along 1D spectra.

* **SNRCUT** [5]. The primary S/N cut (number of σ above the mean/median). 

* **THRESHOLD** [none]. Alternatively to SNRCUT, the primary threshold can be given in the same flux units of the input datacube. This overrides SNRCUT.

* **FLAGGROWTH** [true]. Whether to grow detected sources to a secondary threshold.

* **GROWTHCUT** [3]. Secondary S/N cut used when growing objects (number of σ above the mean/median).

* **GROWTHTHRESHOLD** [none]. Alternatively to GROWTHCUT, the secondary threshold can be given in the same flux units of the input datacube. This overrides GROWTHCUT.

* **MINPIX** [beam area]. The minimum number of spatial pixels for a detection to be accepted. Default is the area covered by the observational beam.

* **MINCHANNELS** [2]. The minimum number of channels for a detection to be accepted.

* **MINVOXELS** [none].  The minimum number of voxels for a detection to be accepted. If not set, MINVOXELS = MINPIX*MINCHANNELS.

* **MAXCHANNELS** [none]. The maximum number of channels for a detection to be accepted. Default is no limits.

* **MAXANGSIZE** [none]. The maximum angular size of a detection to be accepted in *arcmin*. Default is no limits.

* **FLAGADJACENT** [true]. Whether to use the adjacent criterion to merge objects. If *false*, the next two parameters are used to determine whether objects are to be merged.

* **THRESHSPATIAL** [2]. The maximum minimum spatial separation in *pixels* for two objects to be merged into a single one. Ignored if FLAGADJACENT is *true*.

* **THRESHVELOCITY** [3]. The maximum minimum channel separation in *channels* for two objects to be merged into a single one. Ignored if FLAGADJACENT is *true*.

* **REJECTBEFOREMERGE** [true]. Whether to reject sources before merging them.

* **TWOSTAGEMERGING** [true]. Whether to do a partial merge during search.

Outputs
========

The task produces the following outputs. Here *NAME* is the name of the galaxy.

* A FITS file *detections.fits*, a datacube with just the detections.

* A FITS file *DetectMap.fits*, containing a 2D map telling how many channels are detected in each spaxel. 

* FITS files of the moment maps for the detections (*NAME_mom0th.fits*, *NAME_mom1st.fits*, *NAME_mom2nd.fits*).
  
* A text file *detections.txt*, containinig a list of detections and their properties.


Example
========
With the following :download:`parameter <examples/n2403_search.par>` file, the :download:`datacube <examples/ngc2403.fits>` of NGC2403 is searched. Unsurprisingly, the galaxy is detected...

.. literalinclude:: examples/n2403_search.par
   :language: c