.. _alltasks:

List of tasks and parameters
############################

BBarolo's main algorithm for fitting 3D kinematic models to emission line data (:ref:`3DFIT <3dfit>`) makes use of a number of utilities. These tasks include, for example, the disk modeling (:ref:`GALMOD <galmodtask>`), the source finder (:ref:`SEARCH <searchtask>`) and the smoothing utility (:ref:`SMOOTH <smoothtask>`), and can be conveniently used outside the main algorithm as well.

In this page, I list the main tasks and related input parameters available in BBarolo. Parameter names are in **boldface**, default values are in brackets. The names of parameters are not case-sensitive.

General parameters
=======================

In the following, a list of general parameters (e.g.\, not task-specific).

Input/output parameters
^^^^^^^^^^^^^^^^^^^^^^^

* **FITSFILE** [none]. The name of the input FITS file. This is a mandatory parameter for **all tasks**.

* **OUTFOLDER** [./output]. The directory where the output files will be written. 

* **VERBOSE** [true]. Enable all the output messages.

* **THREADS** [1]. Number of CPUs to use for parallelized tasks.

* **SHOWBAR** [true]. Whether to show progress bars.

Beam parameters
^^^^^^^^^^^^^^^^^^^^^^^

Following parameters can be used to specify the size and shape of the Point Spread Function (PSF or beam). These parameters are ignored if beam information is written in the header of the input FITS, either through BMAJ, BMIN and BPA keywords or in the HISTORY. The code defines the beam following the priority order: header -> bmaj,bmin,bpa params -> beamfwhm param -> default to 30 arcsec.

* **BMAJ** [none]. The FWHM of the major axis of the elliptical Gaussian beam in *arcsec*. 

* **BMIN** [none]. The FWHM of the minor minor axis of the elliptical Gaussian beam in *arcsec*.

* **BPA** [none]. The position angle of the major axis of the elliptical Gaussian beam in *degrees*, counter-clock from the North direction.

* **BEAMFWHM** [none]. The FWHM of a circular Gaussian beam in *arcsec* . 


.. _3dfit:

3DFIT task
==========

3DFIT is the main BBarolo's routine: it fits a 3D tilted-ring model to an emission-line data-cube. Algorithms used are described in `this paper <http://adsabs.harvard.edu/abs/2015MNRAS.451.3021D>`_.


* **3DFIT** [false]. This flag enables the 3D fitting algorithm. Can be *true* or *false*. The old flag GALFIT is now deprecated and will be no more supported in future BBarolo's releases.

.. _ringio:

Rings IO
^^^^^^^^^^^^^^^^

Following parameters are used to define the set of rings used for the fit. All parameters are allowed to vary ring-by-ring or to be fixed. In the first case, given values represent initial guesses for the fit.

All parameters listed below (except NRADII and RADSEP) can be given in the form of a single value valid for all rings or through a text file containing values at different radii. In this second case, the sintax to be used is *file(filename,N,M)*, where *filename* is the name of the file with values, *N* is the column number (counting from 1) and *M* is the starting row (all rows if omitted).


* **NRADII** [none]. The number of rings to be used and fitted. If not given, BBarolo tries to guess it from the size of the galaxy.

* **RADSEP** [none]. The separation between rings in *arcsec*. If N radii have been requested, the rings will be placed at N*RADSEP + RADSEP/2. If not given, BBarolo assumes the FWHM of the beam major axis as radius separation. 

* **RADII** [none]. This parameter can be used as an alternative to NRADII and RADSEP. This needs to be a text file (see above).

* **XPOS** [none]. X-center of rings. Accepted format are in *pixels* (starting from 0, unlike GIPSY) or in WCS coordinates in the format +000.0000d (*degrees*) or +00:00:00.00 (*sexagesimal*). If not specified, BBarolo tries to guess it using the 3D source finder.

* **YPOS** [none]. Like XPOS, but for the y-axis.

* **VSYS** [none]. Systemic velocity in *km/s*. If not given, BBarolo tries to guess it from the global line profile.

* **VROT** [none]. Rotation velocity in *km/s*. If not given, BBarolo tries to guess it.

* **VDISP** [8]. Velocity dispersion in *km/s*. 

* **VRAD** [0]. Radial velocity in *km/s*. 

* **INC** [none]. Inclination in *degrees*. If not given, BBarolo tries to guess it from the column density map.

* **PA** [none]. Position angle in *degrees*, measured anti-clockwise form the North direction. If not given, BBarolo tries to guess it from the velocity field.

* **Z0** [0]. Scale-height of the disc in *arcsec*. 

* **DENS** [1]. Gas surface density in units of *1E20 atoms/cm2*. Fit of this parameter is not currently implemented and its value is not relevant if a normalization is used. 

.. _3dfitopt:

Additional options
^^^^^^^^^^^^^^^^^^

Additional parameters to control and refine the fit. All following parameters have default values and are therefore optional.

* **DELTAINC** [5]. This parameter fixes the boundaries of parameter space at [INC-DELTAINC, INC+DELTAINC]. It is not advisable to let the inclination varying over the whole range [0,90].

* **DELTAPA** [15]. This parameter fixes the boundaries of parameter space at [PA-DELTAINC, PA+DELTAPA]. It is not advisable to let the position angle varying over the whole range [0,360].

* **FREE** [VROT VDISP INC PA]. The list of parameters to fit.

* **FTYPE** [2]. Function to be minimized. Accepted values are: 1 = chi-squared, 2 = \|mod-obs\|, (default) and 3 = \|mod-obs\|/(mod+obs)).

* **WFUNC** [2]. Weighting function to be used in the fit. Accepted values are: 0 = uniform weight, 1 = \|cos(θ)\| and 2 = cos(θ)^2, default), where θ is the azimuthal angle (= 0 for galaxy major axis).

* **LTYPE** [1]. Layer type along z. Accepted values are: 1 = Gaussian (default), 2  = sech^2, 3 = exponential, 4 = Lorentzian and 5 = box.

* **CDENS** [10]. Surface density of clouds in the plane of the rings per area of a pixel in units of *1E20 atoms/cm^2* (see also GIPSY `GALMOD <https://www.astro.rug.nl/~gipsy/tsk/galmod.dc1>`_).

* **NV** [nchan]. Number of subclouds in the velocity profile of a single cloud (see also GIPSY `GALMOD <https://www.astro.rug.nl/~gipsy/tsk/galmod.dc1>`_). Default is the number of channels in the datacube.

* **SIDE** [B]: Side of the galaxy to be fitted. Accepted values are: A = approaching, R = receding and B = both (default)

* **MASK** [SMOOTH]. This parameter tells the code how to build a mask to identify the regions of genuine galaxy emission. Accepted values are *SMOOTH*, *SEARCH*, *THRESHOLD*, *NONE* or a FITS mask file:

  * *SMOOTH*: the input cube is smoothed according to the :ref:`smooth parameters <smoothtask>` and the mask built from the region at S/N>BLANKCUT, where **BLANKCUT** is a parameter representing the S/N cut to apply in the smoothed datacube. Defaults are to smooth by a FACTOR = 2 and cut at BLANKCUT = 3.
  
  * *SEARCH*: the source finding is run and the largest detection used to determine the mask. The :ref:`source finding parameters <searchtask>` can be set to change the default values. 
  
  * *THRESHOLD*: blank all pixels with flux < THRESHOLD. A **THRESHOLD** parameter must be specified in the same flux units of the input datacube. 
  
  *  *NONE*: all regions with flux > 0 are used. 
  
  * *file(fitsname.fits)*: A mask FITS file (i.e. filled with 0,1).

* **NORM** [LOCAL]. Type of normalization of the model. Accepted values are: *LOCAL* (pixel by pixel), *AZIM* (azimuthal) or *NONE*.

* **TWOSTAGE** [true]. This flag enables the second fitting stage after parameter regularization. This is relevant just if the user wishes to fit parameters other than VROT, VDISP, VRAD and VVERT. The inclination and the position angle are regularized by polynomials of degree POLYN or a Bezier function (default), while the other parameters by constant functions.

* **POLYN** [-1]. Degree of polynomials for the regularization of inclination and position angles. -1 enables the Bezier function.

* **BWEIGHT** [2]. Exponent of weight for blank pixels. See Section 2.4 of reference paper for details.

* **FLAGERRORS** [false]. Whether the code has to estimate the errors. This usually heavily slows down the run.

* **STARTRAD** [0]. This parameter allows the user to start the fit from the given ring.

* **LINEAR** [0.85]. This parameter controls the spectral broadening of the instrument. It is in units of channel and it represents the standard deviation, not the FWHM. The default is for data that has been Hanning smoothed, so that FWHM = 2 channels and σ = FWHM/2.355.


Additional parameters for high-z galaxies (BBarolo > 1.2.1)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For high-z galaxies you need to set two additional parameters.

* **RESTWAVE** [none]. The rest wavelength of the line you want to fit. Units must be the same of the spectral axis of the cube. For example, if we want fit H-alpha and CUNIT3 = "angstrom", set 6563.

* **REDSHIFT** [none]. The redshift of the galaxy.

These two parameters are used to calculate the conversion from wavelengths to velocities. The velocity reference is set to 0 at RESTWAVE*(REDSHIFT+1). VSYS has to be set to 0, but can be also used to fine-tune the redshift. Finally, if these two parameters are not set, BBarolo will use the CRPIX3 as velocity reference and the proper VSYS has to be set based on that.


.. _galmodtask:

GALMOD task
============

GALMOD is the routine underlying the 3DFIT task. It builds a 3D simulated datacube of a disk galaxy starting from the a set of concentric rings with given column density and kinematics. The routine is an updated version of the namesake routine in GIPSY (see also GIPSY `GALMOD <https://www.astro.rug.nl/~gipsy/tsk/galmod.dc1>`_). 

Parameters for rings are the same of the :ref:`3DFIT task <ringio>`. Options are LTYPE, CDENS and NV (see :ref:`3DFIT options <3dfitopt>`).

Additional GALMOD-specific parameters are:

* **GALMOD** [false]. This flag enables the 3D disk modelling. Can be *true* or *false*.

* **VVERT** [0]. Vertical velocity in *km/s*. 

* **DVDZ** [0]. Gradient of rotation velocity as we move away from the disk plane. This is in *km/s/arcs*.

* **ZCYL** [0]. Height in *arcsec* from the disk plane where the gradient DVDZ begins.

* **SM** [true]. Whether to smooth the model to the same spatial resolution of data.


.. _searchtask:

SEARCH task
============

BBarolo's search algorithm is derived from `Duchamp <https://www.atnf.csiro.au/people/Matthew.Whiting/Duchamp>`_, a 3D source finder for spectral-line data developed by `Matthew Whiting <https://www.atnf.csiro.au/people/Matthew.Whiting/>`_. BBarolo adds a few functionalities and a (mild) parallelization. For a comprehensive description of the algorithm and the input parameters, see Duchamp's `main paper <http://adsabs.harvard.edu/abs/2012MNRAS.421.3242W>`_ and `user guide <http://www.atnf.csiro.au/people/Matthew.Whiting/Duchamp/downloads/UserGuide-1.6.1.pdf>`_.

Main parameters to control the source finder are as follows.

* **SEARCH** [false]. This flag enables the source finding algorithm. Can be *true* or *false*.

* **FLAGROBUSTSTATS** [true]. Whether to use to robust estimators (median and MADFM) instead of normal estimators (mean and standard deviation) when calculating cube statistics.

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


.. _smoothtask:

SMOOTH task
============

This task convolves each channel map in a datacube with a given elliptical Gaussian. 

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


2DFIT task
==========
The classical 2D tilted-ring modelling of a galaxy: a model velocity field is fitted to the observed velocity field (see, e.g., `Begeman 1987 <http://adsabs.harvard.edu/abs/1987PhDT.......199B>`_). This technique is fast and good for high spatial resolution data, but completely unreliable for low resolution data (no beam smearing correction).

* **2DFIT** [false]. This flag enables the 2D fitting of the velocity field.

Parameters and options that control the task are in common with :ref:`3DFIT <3dfitopt>`. In particular, 2DFIT supports the following parameters: **NRADII**, **RADSEP**, **XPOS**, **YPOS**, **VSYS**, **VROT**, **VRAD**, **PA**, **INC**, **FREE**, **SIDE**, **WFUNC**. If **FITSFILE** is a datacube, the velocity field to fit is extracted as 1st moment using a mask for the input datacube defined by the **MASK** parameter (written in the output directory). If **FITSFILE** is a 2D velocity map, this is used to fit the tilted-ring model.


.. _ellproftask:

ELLPROF task
==========
This task can be used to calculate the radial density profile of a galaxy.  

* **ELLPROF** [false]. This flag enables the radial profile task.

Parameters for the task are: **RADII**, **NRADII**, **RADSEP**, **XPOS**, **YPOS**, **PA**, **INC**, **SIDE** (see :ref:`3DFIT <3dfit>`). If **FITSFILE** is a datacube, the profile is calculated from the column density map calculated after masking the cube accordingly to the **MASK** parameter. If **FITSFILE** is a 2D intensity map, this is used to extract the profile.


.. _momentstask:

Moment maps and position-velocity cuts 
======================================

BBarolo can be used to extract global profiles, moment maps and position velocity diagrams. For moment maps and profile, the input datacube can be masked using the MASK parameter (see :ref:`3DFIT <3dfitopt>`).

* **GLOBALPROFILE** [false]. It *true*, calculate the total line profile from a datacube and write it to a text file.

* **TOTALMAP** [false]. It *true*, calculate the total intensity map from a datacube and write it to a FITS file. 

* **VELOCITYMAP** [false]. It *true*, calculate the velocity field from a datacube and write it to a FITS file. 

* **DISPERSIONMAP** [false]. It *true*, calculate the velocity dispersion field from a datacube and write it to a FITS file. 

* **MASSDENSMAP** [false]. It *true*, calculate a mass surface-density map in units of Msun/pc^2 from a datacube and write it to a FITS file. This is just for HI data and the input datacube is required to have JY/BEAM flux density units.

* **FLAGPV** [false]. If *true*, extract position-velocity image from a datacube and write it to a FITS file. The cut is defined by a point and an angle, as set with the following parameters.

* **XPOS_PV** [none]. Reference X pixel of the cut.

* **YPOS_PV** [none]. Reference Y pixel of the cut.

* **PA_PV** [none]. Position angle of the cut, defined anti-clockwise from the X direction. 
