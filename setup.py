#!/usr/bin/env python
from __future__ import print_function
import os, sys, pip, subprocess
from distutils.core import setup
from distutils.dir_util import remove_tree, mkpath
import multiprocessing as mpr
#from pyBBarolo import __version__ as version
version = "1.1.0"
logfile = "setup.log"
try: os.remove(logfile)
except: pass
finally: f = open("setup.log", "a")

def checkModule(module):
    print('Checking %s... '%module,end='')
    try:
        __import__(module)
        print ("OK.")
    except ImportError:
        print("Module '%s' is not present, I will try to install it."%module)
        pip.main(['install',module])
  
def compileBB():
    if not os.path.isfile("Makefile"):
        print ("Running BBarolo configure... ",end="")
        sys.stdout.flush()
        ret = subprocess.call(["./configure"], shell=True, stdout=f)
        if ret!=0: sys.exit("\nConfiguration script failed. Check %s for errors.\n"%logfile)
        print ("OK.")    

    print ("Compiling BBarolo... ",end="")
    sys.stdout.flush()
    ret = subprocess.call("make -j%i lib"%mpr.cpu_count(), shell=True, stdout=f)
    if ret!=0: sys.exit("\nCompilation failed. Check %s for errors.\n"%logfile)
    print ("OK.")
    
      
# First: check if dependencies are available
#modules = ['numpy','astropy']
#for m in modules: checkModule(m)
    

# If installing, the additional data are the compiled libraries
package_data = {'pyBBarolo': ['libBB*']}


# Installing pyBBarolo package
setup(name='pyBBarolo',
      version=version,
      description='a Python wrapper to BBarolo code',
      author=['Enrico Di Teodoro'],
      author_email=['enrico.diteodoro@gmail.com'], 
      url='https://github.com/editeodoro/Bbarolo',
      download_url="https://github.com/editeodoro/Bbarolo",
      packages=['pyBBarolo'],
      package_dir={'pyBBarolo':'pyBBarolo'}, 
      package_data=package_data,
          classifiers=[
                   "Development Status :: 3 - Alpha",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
          ],
    )
    


