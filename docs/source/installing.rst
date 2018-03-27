

Installing BBarolo
##################


Pre-compiled binaries
=====================

Pre-compiled executable files for the GUI and the command line utility are available `at this page <http://editeodoro.github.io/Bbarolo/downloads/binaries/>`_. Binaries are for Linux x86_64 and Mac OS X (> 10.6). 

The command line executable (``BBarolo``) is also included in the GUI binary packages and can be found at:

* **Linux**: same directory of BBaroloGUI
* **MacOS**: BBaroloGUI.app/Contents/MacOS/BBarolo

The pre-compiled executables should work on most systems. If they don't, please compile BBarolo from source following the instructions below.


Compiling from source
=====================

The best way to exploit BBarolo functionalities is to compile the code directly on your computer.
BBarolo compiles and runs on Unix machines only.


.. _requirements:

Requirements
^^^^^^^^^^^^

To compile the code, all you need is:

- a C++ compiler supporting C++11 standard, like the `GNU <https://gcc.gnu.org/>`_ compiler. OpenMP support is required for multi-threading.
- `CFITSIO <http://heasarc.gsfc.nasa.gov/fitsio/>`_ library.
- `FFTW <http://www.fftw.org/>`_ library.
- `WCS <http://www.atnf.csiro.au/people/mcalabre/WCS/>`_ library.
- `QT toolkit <http://www.qt.io/developers/>`_ (> 4.0), only if you would like to use the GUI.
- `Gnuplot <http://www.gnuplot.info/>`_ and `Python <https://www.python.org/>`_ (> 2.6) with the `Astropy <http://www.astropy.org/>`_ package. 

Most of these libraries and packages should already be installed on scientific machines. Otherwise, you can easily install them through the terminal commands of the various package managers, i.e. *apt-get* on Ubuntu-based, *pacman* on Arch-based, *yum* on RPM-based distros, *brew* or *port* on Mac OS X. Note that the QT toolkit is only needed to compile the GUI, which is not mandatory. Gnuplot and Python are not needed to successfully compile the code, but without them BBarolo will not produce any outputs. 

.. _compiling:

Compiling
^^^^^^^^^^^^

If your machine satisfies the above requirements, compiling BBarolo will hopefully be a piece of cake. 

1. **Download** the `latest <https://github.com/editeodoro/Bbarolo/archive/1.4.tar.gz>`_ stable release. From a terminal::

    > wget https://github.com/editeodoro/Bbarolo/archive/X.Y.tar.gz .


  where *X.Y* is the release version. If you are brave, you can also try the latest (non stable) source code from `here <https://github.com/editeodoro/Bbarolo>`_.

2. **Uncompress** it and enter the BBarolo directory::

    > tar -xvf Bbarolo-X.Y.tar.gz
    > cd Bbarolo-X.Y


3. **Configure** running autoconfigure script::

    > ./configure 
     
  If the script is not executable: ``> chmod +x configure``. 
  The configure script takes a number of optional arguments. For instance, it is possible to specify installation and library directories, or the compiler to use::

    > ./configure CXX=icpc --prefix=/dir/to/install --with-cfitsio=/dir/of/cfitsio --with-wcslib=/dir/to/wcslib

  If the configuration fails, follow the suggestions given by the script to manually set the path of libraries.

4. **Compile** the source::

    > make
    
  To compile in parallel: ``> make -j N``, where N is the number of processors. If the compilation   succeeds, the executable **BBarolo** will appear in the current directory. 



Optional steps
""""""""""""""

5. Install, e.g. copy the executable in the installation path (default is /usr/local/bin)::
     
    > make install
     
6. Compile the GUI (QT > 4 needed)::
    
    > make gui 

  This can fail for a number of reasons. 

7. Compile BBarolo as a library::

    > make lib

8. Clean up unnecessary files::

    > make clean


