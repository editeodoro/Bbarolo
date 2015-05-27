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
    Utilities/galfit.cpp \
    Utilities/galfit_errors.cpp \
    Utilities/galfit_min.cpp \
    Utilities/galmod.cpp \
    Utilities/interpolation.cpp \
    Utilities/lsqfit.cpp \
    Utilities/progressbar.cpp \
    Utilities/ringmodel.cpp \
    Utilities/smooth3D.cpp \
    Utilities/statistics.cpp \
    Utilities/string.cpp \
    Utilities/utils.cpp \
    Utilities/galfit_out.cpp \
    Utilities/slitfit.cpp \
    Utilities/wcsUtils.cpp

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
    Utilities/galfit.hh \
    Utilities/galmod.hh \
    Utilities/gnuplot.hh \
    Utilities/lsqfit.hh \
    Utilities/moment.hh \
    Utilities/paramguess.hh \
    Utilities/progressbar.hh \
    Utilities/ringmodel.hh \
    Utilities/smooth3D.hh \
    Utilities/spacepar.hh \
    Utilities/utils.hh \
    GUI/consolestream.h \
    GUI/q_streamdebug.h

FORMS += \
    GUI/bbarolowindow.ui
