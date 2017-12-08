
################
Installing
################

=====================
Pre-compiled binaries
=====================




=====================
Compiling from source
=====================

The best way to exploit BBarolo functionalities is to compile the code directly on your computer.
BBarolo compiles and runs on Unix machines.

^^^^^^^^^^^^
Requirements
^^^^^^^^^^^^

To compile the code, all you need is:

- a C++ compiler supporting C++11 standard, like the `GNU <https://gcc.gnu.org/>`_ compiler.
- `CFITSIO <http://heasarc.gsfc.nasa.gov/fitsio/>`_ library.
- `FFTW <http://www.fftw.org/>`_ library.
- `WCS <http://www.atnf.csiro.au/people/mcalabre/WCS/>`_ library.
- `QT toolkit <http://www.qt.io/developers/>`_ (> 4.0), only if you would like to use the GUI.
- `Gnuplot <http://www.gnuplot.info/>`_ and `Python <https://www.python.org/>`_ (> 2.6) with the `Astropy <http://www.astropy.org/>`_ package. 

Most of these libraries and packages should already be installed on scientific machines. Otherwise, you can easily install them through the terminal commands of the various package managers, i.e. *apt-get* on Ubuntu-based, *pacman* on Arch-based, *yum* on RPM-based distros, *brew* or *port* on Mac OS X. Note that the QT toolkit is only needed to compile the GUI, which is not mandatory. Gnuplot and Python are not needed to successfully compile the code, but without them BBarolo will not produce any outputs. 

^^^^^^^^^^^^
Compiling
^^^^^^^^^^^^

If your machine satisfies the above requirements, compiling BBarolo will be a piece of cake. 

1. **Download** the `latest <https://github.com/editeodoro/Bbarolo/archive/1.3.tar.gz>`_ stable release. From a terminal::

    > wget https://github.com/editeodoro/Bbarolo/archive/1.3.tar.gz .


  If you are brave, you can also try the last (non stable) source code from `here <https://github.com/editeodoro/Bbarolo>`_.

2. **Uncompress** it and enter the BBarolo directory::

    > tar -xvzf Bbarolo-X.Y.tar
    > cd Bbarolo-X.Y
    
  where *X.Y* is the release version.


3. **Configure** running autoconfigure script::

    > ./configure 
     
  If the script is not executable: ``> chmod +x configure``. 
  The configure script takes a number of optional arguments. For instance, it is possible to specify installation and library directories, or the compiler to use::

    > ./configure CXX=icpc --prefix=/dir/to/install --with-cfitsio=/dir/of/cfitsio --with-wcslib=/dir/to/wcslib

  If the configuration fails, follow the suggestions given by the script for manually set the path of libraries.

4. **Compile** the source::

    > make
    
  To compile in parallel: ``> make -j N``, where N is the number of processors. If the compilation   succeeds, the executable **BBarolo** will appear in the current directory. 


""""""""""""""
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


