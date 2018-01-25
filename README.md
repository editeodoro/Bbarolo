----------------- BBAROLO SHORT GUIDE -----------------------

BBarolo is a 3D fitting tool to derive the kinematics of galaxies from emission-line observations.

Needed dependencies: FFTW3, WCSLIB, CFITSIO. 
Optional: GNUPLOT, PYTHON with ASTROPY (for output plots), QT Kit (for the GUI) 


1) Compiling the code.

 > ./configure
 > make
 > make install (optional)
 
For the GUI (optional): 
 
 > make gui



2) Running examples (to be separetely downloaded): 

 > ./BBarolo -p examples/n2403.par


3) A few commands:

 - Template input file: 
	> ./BBarolo -t
	
 - List of available parameters:
	> ./BBarolo -d


A deeper description of parameters can be found at BBarolo's website (http://editeodoro.github.io/Bbarolo/). 

