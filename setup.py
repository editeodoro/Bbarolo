import os, sys, re
import subprocess
import multiprocessing as mpr
from setuptools import setup

def read_version():
    """Extracts the version number from the first line using Regex"""
    version_file = os.path.join(os.path.dirname(__file__), "pyBBarolo", "_version.py")
    try:
        with open(version_file, "r") as f:
            first_line = f.readline().strip()
        match = re.search(r"(\d+\.\d+\.\d+)", first_line)
        if match:
            return match.group(1)
        else:
            sys.exit(f"\nError: Could not find a valid version number (X.Y.Z) in the first line: '{first_line}'\n")
    except IOError:
        sys.exit(f"\nError: Could not read version from {version_file}\n")

def compileBB():
    logfile = "setup.log"
    if os.path.exists(logfile):
        os.remove(logfile)

    with open(logfile, "a") as f:
        if not os.path.isfile("Makefile"):
            print("Running BBarolo configure... ")
            ret = subprocess.call("./configure", shell=True, stdout=f)
            if ret != 0:
                sys.exit(f"\nConfiguration script failed. Check {logfile} for errors.\n")

        print("Compiling BBarolo... ")
        ret = subprocess.call(f"make -j{mpr.cpu_count()}", shell=True, stdout=f)
        if ret != 0:
            sys.exit(f"\nCompilation failed. Check {logfile} for errors.\n")
        
        ret = subprocess.call("make lib", shell=True, stdout=f)
        if ret != 0:
            sys.exit(f"\nLibrary linking failed. Check {logfile} for errors.\n")
        
        subprocess.call("mv BBarolo pyBBarolo/", shell=True, stdout=f)
        subprocess.call("mv libBB* pyBBarolo/", shell=True, stdout=f)

# 1. Run compilation linearly before setuptools executes
current_version = read_version()
print(f"Building pyBBarolo version: {current_version}")

compileBB()

# 2. Run the actual setup
setup(
    name="pyBBarolo",
    version=current_version,
    packages=['pyBBarolo'],
    package_dir={'pyBBarolo': 'pyBBarolo'},
    package_data={'pyBBarolo': ['libBB*', 'BBarolo']},
    include_package_data=False, 
)