#!/usr/bin/env python
from __future__ import print_function
import os, sys, pip, subprocess
from distutils.core import setup
from distutils.dir_util import remove_tree, mkpath
import multiprocessing as mpr
version=open("pyBBarolo/_version.py").readlines()[-1].split()[-1].strip("\"'")
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
    ret = subprocess.call("make -j%i"%mpr.cpu_count(), shell=True, stdout=f)
    if ret!=0: sys.exit("\nCompilation failed. Check %s for errors.\n"%logfile)
    ret = subprocess.call("make lib", shell=True, stdout=f)
    if ret!=0: sys.exit("\nLibrary linking failed. Check %s for errors.\n"%logfile)
    print ("OK.")
    

if sys.argv[1]=='sdist':  
    # If we are creating the sdist package, make a tar with BB source
    try: remove_tree("pyBBarolo/BB")
    except: pass
    try: os.remove("pyBBarolo/BBarolo.tar.gz")
    except: pass
    mkpath("pyBBarolo/BB")
    subprocess.call("cp -r src/ pyBBarolo/BB/src", shell=True,stdout=f)
    subprocess.call("cp -r config/ pyBBarolo/BB/config", shell=True,stdout=f)
    subprocess.call("cp -r configure pyBBarolo/BB/", shell=True,stdout=f)
    subprocess.call("cp -r Makefile.in pyBBarolo/BB/", shell=True,stdout=f)
    subprocess.call("cp -r Makefile.in pyBBarolo/BB/", shell=True,stdout=f)
    subprocess.call("rm -rf pyBBarolo/BB/src/Build", shell=True,stdout=f)
    subprocess.call("tar -czf ./pyBBarolo/BBarolo.tar.gz -C ./pyBBarolo BB", shell=True,stdout=f)
    remove_tree("pyBBarolo/BB/")

    # If we creating the dist, additional file is the tar just created
    package_data = {'pyBBarolo': ['*.tar.gz']}

else:
    print ("------ Installing pyBBarolo v%s ------"%version)
    
    # First: check if dependencies are available
    modules = ['numpy','astropy']
    for m in modules: checkModule(m)
    
    # Second: compile BB library
    if os.path.isdir("./src"):
        # If we are in BB root, compile from here
        compileBB()
        subprocess.call("cp BBarolo pyBBarolo/", shell=True, stdout=f)
        subprocess.call("mv libBB* pyBBarolo/", shell=True, stdout=f)
    else:
        # Enter pyBB dir and look for tar file
        os.chdir("pyBBarolo")
        if not os.path.isfile("BBarolo.tar.gz"): raise ValueError("BBarolo.tar.gz not found")
        subprocess.call(["tar xvzf BBarolo.tar.gz"], shell=True, stdout=f)
        os.chdir("BB")
        compileBB()
        os.chdir("../")
        subprocess.call("mv BB/libBB* .", shell=True, stdout=f)
        subprocess.call("mv BB/BBarolo .", shell=True, stdout=f)
        #remove_tree("BB")
        os.chdir("../")
    
    # If installing, the additional data are the compiled libraries
    package_data = {'pyBBarolo': ['libBB*','BBarolo']}


# Installing pyBBarolo package
setup(name='pyBBarolo',
      version=version,
      description='a Python wrapper to BBarolo code',
      author='Enrico Di Teodoro',
      author_email='enrico.diteodoro@gmail.com', 
      url='https://github.com/editeodoro/Bbarolo',
      download_url="https://github.com/editeodoro/Bbarolo",
      packages=['pyBBarolo'],
      package_dir={'pyBBarolo':'pyBBarolo'}, 
      package_data=package_data,
      classifiers=["Development Status :: 3 - Alpha",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
          ],
    )


