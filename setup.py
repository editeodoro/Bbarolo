#!/usr/bin/env python
from __future__ import print_function
import os, sys, pip, subprocess
from distutils.core import setup
from distutils.dir_util import remove_tree, mkpath
import multiprocessing as mpr
#from pyBBarolo import __version__ as version
version = "1.0.4"
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
    
      
if sys.argv[1]=='sdist':  
    # If we are creating the sdist package, make a tar with BB source
    try: remove_tree("pyBBarolo/BBarolo")
    except: pass
    try: remove_tree("pyBBarolo/BBarolo.tar.gz")
    except: pass
    mkpath("pyBBarolo/BBarolo")    
    subprocess.call("cp -r src/ pyBBarolo/BBarolo/src", shell=True,stdout=f)
    subprocess.call("cp -r config/ pyBBarolo/BBarolo/config", shell=True,stdout=f)
    subprocess.call("cp -r configure pyBBarolo/BBarolo/", shell=True,stdout=f)
    subprocess.call("cp -r Makefile.in pyBBarolo/BBarolo/", shell=True,stdout=f)
    subprocess.call("cp -r Makefile.in pyBBarolo/BBarolo/", shell=True,stdout=f)
    subprocess.call("rm -rf pyBBarolo/BBarolo/src/Build", shell=True,stdout=f)
    subprocess.call("cd pyBBarolo && tar -czvf BBarolo.tar.gz BBarolo", shell=True,stdout=f)
    remove_tree("pyBBarolo/BBarolo")
    
    # If we creating the dist, additional file is the tar just created
    package_data = {'pyBBarolo': ['*.tar.gz']}
    
else:
    # First: check if dependencies are available
    modules = ['numpy','astropy']
    for m in modules: checkModule(m)
    
    # Second: compile BB library
    if os.path.isdir("./src"):
        # If we are in BB root, compile from here
        compileBB()
        subprocess.call("mv libBB* pyBBarolo/", shell=True, stdout=f)
    else:
        # Enter pyBB dir and look for tar file
        os.chdir("pyBBarolo")
        if not os.path.isfile("BBarolo.tar.gz"): raise ValueError("BBarolo.tar.gz not found")
        subprocess.call(["tar xvzf BBarolo.tar.gz"], shell=True, stdout=f)
        os.chdir("BBarolo")
        compileBB()
        os.chdir("../")
        subprocess.call("mv BBarolo/libBB* .", shell=True, stdout=f)
        #remove_tree("BBarolo")
        os.chdir("../")
    
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
    


