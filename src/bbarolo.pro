TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QT += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

SOURCES +=     bbarolo.cpp \
    Arrays/cube.cpp \
    Arrays/header.cpp \
    Arrays/image.cpp \
    Arrays/param.cpp \
    Arrays/search.cpp \
    Arrays/stats.cpp \
    GUI/bbarolowindow.cpp \
    GUI/bbarolowindow_run.cpp \
    GUI/main.cpp \
    GUI/qcustomplot.cpp \
    Map/detection.cpp \
    Map/object2D.cpp \
    Map/object3D.cpp \
    Map/objectgrower.cpp \
    Map/scan.cpp \
    Map/voxel.cpp \
    Utilities/allocator.cpp \
    Utilities/conv2D.cpp \
    Utilities/converter.cpp \
    Tasks/galfit.cpp \
    Tasks/galfit_errors.cpp \
    Tasks/galfit_min.cpp \
    Tasks/galmod.cpp \
    Tasks/galwind.cpp \
    Utilities/interpolation.cpp \
    Utilities/lsqfit.cpp \
    Utilities/progressbar.cpp \
    Tasks/ringmodel.cpp \
    Tasks/smooth3D.cpp \
    Utilities/statistics.cpp \
    Utilities/string.cpp \
    Utilities/utils.cpp \
    Utilities/paramguess.cpp \
    Tasks/galfit_out.cpp \
    Tasks/slitfit.cpp \
    Utilities/wcsUtils.cpp \
    Utilities/fitsUtils.cpp \
    Tasks/ellprof.cpp \
    Tasks/moment.cpp

HEADERS += Arrays/cube.hh \
    Arrays/header.hh \
    Arrays/image.hh \
    Arrays/param.hh \
    Arrays/stats.hh \
    GUI/bbarolowindow.h \
    GUI/qcustomplot.h \
    GUI/ui_bbarolowindow.h \
    Map/detection.hh \
    Map/object2D.hh \
    Map/object3D.hh \
    Map/objectgrower.hh \
    Map/scan.hh \
    Map/voxel.hh \
    Utilities/conv2D.hh \
    Utilities/converter.hh \
    Tasks/galfit.hh \
    Tasks/galmod.hh \
    Utilities/gnuplot.hh \
    Utilities/lsqfit.hh \
    Tasks/moment.hh \
    Utilities/paramguess.hh \
    Utilities/progressbar.hh \
    Tasks/ringmodel.hh \
    Tasks/smooth3D.hh \
    Utilities/spacepar.hh \
    Utilities/utils.hh \
    Tasks/ellprof.hh \
    Tasks/galwind.hh \
    GUI/consolestream.h \
    GUI/q_streamdebug.h

FORMS += \
    GUI/bbarolowindow.ui
