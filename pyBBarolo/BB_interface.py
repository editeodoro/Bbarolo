import os
from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer

# Load the library, compiling it if required
libdir = os.path.dirname(os.path.realpath(__file__))
try:
    libBB = np.ctypeslib.load_library("libBBarolo",libdir)
except OSError:
    print ("Failed to load libBBarolo... trying to compile it")
    os.system("cd "+libdir+"; make lib")
    libBB = np.ctypeslib.load_library("libBBarolo", libdir)


# Type definitions
array_1d_int    = ndpointer(dtype=np.int, ndim=1,flags="CONTIGUOUS")
array_1d_float  = ndpointer(dtype=np.float32, ndim=1,flags="CONTIGUOUS")
array_1d_double = ndpointer(dtype=np.double, ndim=1,flags="CONTIGUOUS")


# Class Cube interface #################################################################
libBB.Cube_new.restype = c_void_p
libBB.Cube_new.argtypes = [c_char_p]
libBB.Cube_delete.restype = None
libBB.Cube_delete.argtypes = [c_void_p]
libBB.Cube_axisdim.restype = ndpointer(dtype=c_int, shape=(3))
libBB.Cube_axisdim.argtypes = [c_void_p]
libBB.Cube_array.restype = POINTER(c_float)
libBB.Cube_array.argtypes = [c_void_p]
########################################################################################

# Struct Rings interface################################################################
libBB.Rings_new.restype = c_void_p
libBB.Rings_new.argtypes = [ ]
libBB.Rings_set.restype = c_void_p
libBB.Rings_set.argtypes = [c_void_p,c_int,array_1d_float,array_1d_float,array_1d_float,\
                            array_1d_float,array_1d_float,array_1d_float,array_1d_float,\
                            array_1d_float,array_1d_float,array_1d_float,array_1d_float,\
                            array_1d_float,array_1d_float,array_1d_float,array_1d_int]
########################################################################################

# Class Galfit  interface ##############################################################
libBB.Galfit_new.restype = c_void_p
libBB.Galfit_new.argtypes = [c_void_p]
libBB.Galfit_new_all.restype = c_void_p
libBB.Galfit_new_all.argtypes = [c_void_p,c_void_p,c_float,c_float,c_int,c_int,c_int,c_int,c_int,\
                                 c_double,c_int,c_int,c_char_p,c_char_p,c_char_p,c_char_p,\
                                 c_bool,c_char_p,c_bool,c_bool,c_float,c_double,c_double,c_char_p]
libBB.Galfit_delete.restype = None
libBB.Galfit_delete.argtypes = [c_void_p]
libBB.Galfit_galfit.restype = c_bool
libBB.Galfit_galfit.argtypes = [c_void_p]
libBB.Galfit_secondStage.restype = c_bool
libBB.Galfit_secondStage.argtypes = [c_void_p]
libBB.Galfit_writeModel.restype = None
libBB.Galfit_writeModel.argtypes = [c_void_p,c_char_p]
########################################################################################

# Class Galwind interface ##############################################################
libBB.Galwind_new.restype = c_void_p
libBB.Galwind_new.argtypes = [c_void_p,c_double,c_double,c_double,c_double,c_double,c_double,\
                              c_float,c_double,c_double,c_double,c_int,c_int,c_int,c_int]
libBB.Galwind_delete.restype = None
libBB.Galwind_delete.argtypes = [c_void_p]
libBB.Galwind_array.restype = POINTER(c_float)
libBB.Galwind_array.argtypes = [c_void_p]
libBB.Galwind_compute.restype = c_bool
libBB.Galwind_compute.argtypes = [c_void_p]
libBB.Galwind_smooth.restype = c_bool
libBB.Galwind_smooth.argtypes = [c_void_p]
libBB.Galwind_writeFITS.restype = c_bool
libBB.Galwind_writeFITS.argtypes = [c_void_p]
libBB.Galwind_writeMomentMaps.restype = c_bool
libBB.Galwind_writeMomentMaps.argtypes = [c_void_p]

########################################################################################

