#TEMPLATE = app
#CONFIG += console c++11
#CONFIG -= app_bundle
#CONFIG -= qt

#QT += core gui
#greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

SOURCES += \
    Map/detection.cpp \
    Map/object2D.cpp \
    Map/object3D.cpp \
    Map/objectgrower.cpp \
    Map/scan.cpp \
    Map/voxel.cpp \
    Arrays/cube.cpp \
    Arrays/header.cpp \
    Arrays/image.cpp \
    Arrays/param.cpp \
    Arrays/stats.cpp \
    Tasks/ellprof.cpp \
    Tasks/galfit_errors.cpp \
    Tasks/galfit_min.cpp \
    Tasks/galfit_out.cpp \
    Tasks/galfit.cpp \
    Tasks/galmod.cpp \
    Tasks/galwind.cpp \
    Tasks/moment.cpp \
    Tasks/ringmodel.cpp \
    Tasks/search.cpp \
    Tasks/slitfit.cpp \
    Tasks/smooth3D.cpp \
    Tasks/spacepar.cpp \
    Utilities/conv2D.cpp \
    Utilities/converter.cpp \
    Utilities/fitsUtils.cpp \
    Utilities/interpolation.cpp \
    Utilities/lsqfit.cpp \
    Utilities/paramguess.cpp \
    Utilities/progressbar.cpp \
    Utilities/statistics.cpp \
    Utilities/string.cpp \
    Utilities/utils.cpp \
    Utilities/wcsUtils.cpp \
    BB_interface.cpp \
    bbarolo_mpi.cpp \
    bbarolo.cpp

HEADERS += \
    Map/detection.hh \
    Map/object2D.hh \
    Map/object3D.hh \
    Map/objectgrower.hh \
    Map/scan.hh \
    Map/voxel.hh \
    Arrays/cube.hh \
    Arrays/header.hh \
    Arrays/image.hh \
    Arrays/param.hh \
    Arrays/rings.hh \
    Arrays/stats.hh \
    Tasks/ellprof.hh \
    Tasks/galfit.hh \
    Tasks/galmod.hh \
    Tasks/galwind.hh \
    Tasks/moment.hh \
    Tasks/ringmodel.hh \
    Tasks/search.hh \
    Tasks/smooth3D.hh \
    Tasks/spacepar.hh \
    Utilities/allocator.hpp \
    Utilities/conv2D.hh \
    Utilities/converter.hh \
    Utilities/gnuplot.hh \
    Utilities/lsqfit.hh \
    Utilities/optimization.hh \
    Utilities/paramguess.hh \
    Utilities/progressbar.hh \
    Utilities/utils.hh \
    Tasks/rendering3D.hh

FORMS += \
    GUI/bbarolowindow.ui
