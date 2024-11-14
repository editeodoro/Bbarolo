"""
Interface between C++ and Python functions.
"""

########################################################################
# Copyright (C) 2017 Enrico Di Teodoro
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################

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
array_1d_int    = ndpointer(dtype=np.int64, ndim=1,flags="CONTIGUOUS")
array_1d_float  = ndpointer(dtype=np.float32, ndim=1,flags="CONTIGUOUS")
array_1d_double = ndpointer(dtype=np.double, ndim=1,flags="CONTIGUOUS")


# Class Param interface ###############################################################
libBB.Param_new.restype = c_void_p
libBB.Param_new.argtypes = [ ]
libBB.Param_setfromfile.restype = None
libBB.Param_setfromfile.argtypes = [c_void_p,c_char_p]
libBB.Param_setfromstr.restype = None
libBB.Param_setfromstr.argtypes = [c_void_p,c_char_p]
libBB.Param_delete.restype = None
libBB.Param_delete.argtypes = [c_void_p]
########################################################################################


# Class Cube interface #################################################################
libBB.Cube_new.restype = c_void_p
libBB.Cube_new.argtypes = [c_char_p]
libBB.Cube_delete.restype = None
libBB.Cube_delete.argtypes = [c_void_p]
libBB.Cube_axisdim.restype = ndpointer(dtype=c_int, shape=(3))
libBB.Cube_axisdim.argtypes = [c_void_p]
libBB.Cube_array.restype = POINTER(c_float)
libBB.Cube_array.argtypes = [c_void_p]
libBB.Cube_setBeam.restype = None
libBB.Cube_setBeam.argtypes = [c_void_p, c_float, c_float, c_float]
libBB.Cube_getBeam.restype = POINTER(c_float)
libBB.Cube_getBeam.argtypes = [c_void_p]
########################################################################################


# Struct Rings interface ###############################################################
libBB.Rings_new.restype = c_void_p
libBB.Rings_new.argtypes = [ ]
libBB.Rings_delete.restype = None
libBB.Rings_delete.argtypes = [c_void_p]
libBB.Rings_set.restype = None
libBB.Rings_set.argtypes = [c_void_p,c_int,array_1d_float,array_1d_float,array_1d_float,\
                            array_1d_float,array_1d_float,array_1d_float,array_1d_float,\
                            array_1d_float,array_1d_float,array_1d_float,array_1d_float,\
                            array_1d_float,array_1d_float,array_1d_float]
########################################################################################


# Class Galmod interface ##############################################################
libBB.Galmod_new.restype = c_void_p
libBB.Galmod_new.argtypes = [c_void_p,c_void_p,c_int,c_int,c_int,c_float,c_int]
libBB.Galmod_new_par.restype = c_void_p
libBB.Galmod_new_par.argtypes = [c_void_p,c_void_p,c_void_p]
libBB.Galmod_delete.restype = None
libBB.Galmod_delete.argtypes = [c_void_p]
libBB.Galmod_array.restype = POINTER(c_float)
libBB.Galmod_array.argtypes = [c_void_p]
libBB.Galmod_compute.restype = c_bool
libBB.Galmod_compute.argtypes = [c_void_p]
libBB.Galmod_smooth.restype = c_bool
libBB.Galmod_smooth.argtypes = [c_void_p]
########################################################################################


# Class Galfit interface ##############################################################
libBB.Galfit_new.restype = c_void_p
libBB.Galfit_new.argtypes = [c_void_p]
libBB.Galfit_new_par.restype = c_void_p
libBB.Galfit_new_par.argtypes = [c_void_p,c_void_p,c_void_p]
libBB.Galfit_delete.restype = None
libBB.Galfit_delete.argtypes = [c_void_p]
libBB.Galfit_initialGuesses.restype = POINTER(c_float)
libBB.Galfit_initialGuesses.argtypes = [c_void_p,c_char_p,c_char_p,c_char_p,c_char_p]
libBB.Galfit_galfit.restype = c_bool
libBB.Galfit_galfit.argtypes = [c_void_p]
libBB.Galfit_secondStage.restype = c_bool
libBB.Galfit_secondStage.argtypes = [c_void_p]
libBB.Galfit_writeModel.restype = None
libBB.Galfit_writeModel.argtypes = [c_void_p,c_char_p,c_bool]
libBB.Galfit_plotModel.restype = c_int
libBB.Galfit_plotModel.argtypes = [c_void_p]
libBB.Galfit_calcresiduals.restype = c_float
libBB.Galfit_calcresiduals.argtypes = [c_void_p,c_void_p]
########################################################################################


# Class Galwind interface ##############################################################
libBB.Galwind_new.restype = c_void_p
libBB.Galwind_new.argtypes = [c_void_p,c_float,c_float,c_float,c_float,c_float,c_float,\
                              c_float,c_float,c_float,c_float,c_int,c_int,c_int,c_int,c_int]
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


# Class Search interface ##############################################################
libBB.Search_search.restype = None
libBB.Search_search.argtypes = [c_void_p,c_char_p,c_float,c_float,c_bool,c_int,c_int,\
                                c_int,c_int,c_int,c_int,c_float,c_bool,c_float,c_float,\
                                c_bool,c_bool,c_int]
########################################################################################


# Class Ringmodel interface ############################################################
libBB.Fit2D_new.restype = c_void_p
libBB.Fit2D_new.argtypes = [c_void_p,c_void_p,c_char_p,c_char_p,c_char_p,c_int,c_int]
libBB.Fit2D_delete.restype = None
libBB.Fit2D_delete.argtypes = [c_void_p]
libBB.Fit2D_compute.restype = None
libBB.Fit2D_compute.argtypes = [c_void_p]
libBB.Fit2D_write.restype = None
libBB.Fit2D_write.argtypes = [c_void_p,c_void_p,c_char_p]
########################################################################################


# Class Ellprof interface ############################################################
libBB.Ellprof_new.restype = c_void_p
libBB.Ellprof_new.argtypes = [c_void_p,c_void_p,c_char_p,c_char_p,c_int]
libBB.Ellprof_delete.restype = None
libBB.Ellprof_delete.argtypes = [c_void_p]
libBB.Ellprof_compute.restype = None
libBB.Ellprof_compute.argtypes = [c_void_p]
libBB.Ellprof_write.restype = None
libBB.Ellprof_write.argtypes = [c_void_p,c_char_p]
########################################################################################


# Class SpectralSmooth3D interface ############################################################
libBB.SpectralSmooth3D_new.restype = c_void_p
libBB.SpectralSmooth3D_new.argtypes = [c_char_p,c_size_t]
libBB.SpectralSmooth3D_delete.restype = None
libBB.SpectralSmooth3D_delete.argtypes = [c_void_p]
libBB.SpectralSmooth3D_compute.restype = None
libBB.SpectralSmooth3D_compute.argtypes = [c_void_p,c_void_p,c_int]
libBB.SpectralSmooth3D_write.restype = None
libBB.SpectralSmooth3D_write.argtypes = [c_void_p,c_void_p,c_char_p,c_bool]
libBB.SpectralSmooth3D_array.restype = POINTER(c_float)
libBB.SpectralSmooth3D_array.argtypes = [c_void_p]
########################################################################################