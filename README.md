<p align="center">
  <img src="http://editeodoro.github.io/Bbarolo/files/bbarolo.jpg" alt="BBlogo"/>
</p>


BBarolo is a 3D fitting tool to derive the kinematics of galaxies from emission-line observations.

Full documentation is hosted at http://bbarolo.readthedocs.io/en/latest. 


## Dependencies

Needed: FFTW3, WCSLIB, CFITSIO. 

Optional: GNUPLOT, PYTHON with ASTROPY (for output plots), QT Kit (for the GUI) 


## Installing BBarolo

Clone the repository:

````
git clone -b master --single-branch https://github.com/editeodoro/Bbarolo 
cd Bbarolo
````

Compile:
````
./configure
make
make install
 ````
Last command is optional to install BBarolo executable in a given path (default /usr/local/bin).


To compile the GUI (optional): 
 ````
 make gui
````

To install Python wrapper pyBBarolo:
 ````
 make pybbinst
````

