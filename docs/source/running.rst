 

Running BBarolo
################



Graphical User Interface
^^^^^^^^^^^^^^^^^^^^^^^^^

BBarolo comes with a Graphical User Interface (GUI) that can help the user in setting up the input parameters. 
Running BBarolo through the GUI should be quite straightforward: you do not have to learn the annoying list of :ref:`available parameters <alltasks>`: the user needs just to fill the required fields in the GUI and this will create a text file with the correct parameters and run BBarolo.

**N.B.:** Although the GUI allows the user to set the main parameters, many options can only be enabled through the command line tool. Moreover, I stress that the GUI has not been updated together with the last releases of BBarolo. If you experience any issues with the GUI or want to have full control of the code, I recommend you to use the command line.

Command line
^^^^^^^^^^^^
BBarolo is mainly meant to be run from the command line. For a very quick guide and to appreciate the biggest achievement of my PhD, just type ``BBarolo`` on your keyboard in a terminal window. 


**Execution with a parameter file**: BBarolo takes input parameters specified through a parameter file, provided at the runtime. This is a text file containing a list of parameter names and values::

     PARAM1     VALUE1
     PARAM2     VALUE2
     PARAM3     VALUE3
     ...        ...
     
All available parameters are described in the :ref:`task documentation <alltasks>`. In the input file, parameter names are not case-sensitive and lines starting with \# or \/\/ are not read in. The order in which parameters are listed is unimportant, but, if a parameter is listed more than once, only the last value is considered.

Some parameters are mandatory, some others are optional and have default values which are assumed when not explicitly set. 
A template parameter file for the 3DFIT task can be obtained with the command ``> BBarolo -t``. A list of all parameters with their default values can be printed with ``> BBarolo -d``. An example of parameter file can be found `here <http://editeodoro.github.io/Bbarolo/resources/param.par>`_, full runnable instances can be downloaded from `this page <http://editeodoro.github.io/Bbarolo/downloads/examples>`_. The command ``> BBarolo -v`` will return information about the code version and compiler flags.

After your parameter file is ready, BBarolo can be run with the following::

    > BBarolo -p paramfile
    
where *paramfile* is the name of the user-defined input file. Since version v1.6, individual parameter values can also be overridden from the command line directly, without changing the parameter file. For example:: 

    > BBarolo -p paramfile INC=60 PA=120
    
will set a INC of 60 degrees and a PA of 120 degrees and ignore the values listed in the *paramfile*.

|

**Execution with command-line arguments**: BBarolo can alternatively read all parameters directly from the command line using the ``-c`` option. Parameters are given in the form ``PARAM=VALUE`` with no blank spaces::

    > BBarolo -c PARAM1=VALUE1 PARAM2=VALUE2 PARAM3=VALUE3

If ``VALUE`` must contain white spaces, just include it in between quotation marks (e.g. ``FREE="VROT VDISP"``). Names of parameters are not case-sensitive. This way is very convenient when the user only needs to run a task with few parameters. For example, to run the source finder with default parameters::

    > BBarolo -c fitsfile=yourfits.fits search=true
    
|

**Automated execution**: BBarolo can otherwise be run in a completely automated way, i.e. providing no parameters but the input FITS datacube. In this case, the code uses the source finder to identify the galaxy in the datacube, tries to guess initial values for the rings and fits a 3D model to the data. Although this procedure might work and return nice best-fit models if used with high resolution and high S/N data, it should be used carefully. In particular, the algorithm for guessing initial values for the fit is still quite coarse. Wrong initial guesses may lead to completely inappropriate models.

If you still want to try the automated execution::

    > BBarolo -f fitsfile
    
where *fitsfile* is the name of the FITS file of the galaxy to analyse.



