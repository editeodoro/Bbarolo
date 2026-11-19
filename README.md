<p align="center">
  <img src="http://editeodoro.github.io/Bbarolo/files/bbarolo.jpg" alt="BBlogo"/>
</p>

> [!WARNING]
> ## BBarolo has moved to version 2.0: [github.com/editeodoro/bbarolo2](https://github.com/editeodoro/bbarolo2)
> ### This repository (BBarolo 1.x) is kept for legacy purposes only.

---


BBarolo is a 3D fitting tool to derive the kinematics of galaxies from emission-line observations.


## Dependencies (BBarolo 1.8)

Needed: FFTW3, WCSLIB, CFITSIO.

Optional: GNUPLOT, PYTHON with ASTROPY (for output plots), QT Kit (for the GUI)


## Installing BBarolo 1.8

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
