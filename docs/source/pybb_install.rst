
Installing pyBBarolo
#####################


From pip
=====================
pyBBarolo is available as a package in the Python Package Index (`PyPI <https://pypi.python.org/pypi>`_). The easiest way of installing it is through ``pip``::

    > pip install pyBBarolo
    
This will download the package, compile BBarolo source code and install pyBBarolo in the python library path. Make sure you have permissions to write in the installing directory. Depending on your computer setup, you may need to run pip with superuser privileges (e.g.: ``> sudo pip install pyBBarolo``).

**N.B.:** The above command will compile BBarolo, which means that your machine needs to have pre-installed all the libraries which the C++ code depends from (see :ref:`requirements <requirements>`). If pip fails during compilation, please follow the procedure below.


Build and install
=====================
The python package can be alternatively installed from the main repository. You'll need to compile BBarolo as a library and then install pyBBarolo:

1. Follow steps 1-4 of procedure to :ref:`compile BBarolo from source <compiling>`::

    > wget https://github.com/editeodoro/Bbarolo/archive/X.Y.tar.gz .
    > tar -xvzf Bbarolo-X.Y.tar && cd Bbarolo-X.Y
    > ./configure
    > make

  where *X.Y* is the software release. PyBBarolo is only available for BBarolo's release 1.4 and over.


2. Install the python package::

    > python setup.py install
    
    
If either compilation or installation fail, refer to BBarolo :ref:`compiling <compiling>` and :ref:`troubleshooting <troubleshooting>` pages. 