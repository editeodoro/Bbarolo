QT += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = BBaroloGUI
TEMPLATE = app

CONFIG += qt

DESTDIR= ../../
OBJECTS_DIR= ../Build
MOC_DIR= ../Build

SOURCES += main.cpp\
        bbarolowindow.cpp \
    bbarolowindow_run.cpp \
    ../Arrays/param.cpp \
    ../Utilities/utils.cpp \
    ../Utilities/wcsUtils.cpp \
    ../Utilities/string.cpp \
    ../Arrays/header.cpp \
#    ../Arrays/stats.cpp \
#    ../Arrays/search.cpp \
#    ../Arrays/image.cpp \
#    ../Arrays/cube.cpp \
#    ../Map/voxel.cpp \
#    ../Map/scan.cpp \
#    ../Map/objectgrower.cpp \
#    ../Map/object3D.cpp \
#    ../Map/object2D.cpp \
#    ../Map/detection.cpp \
#    ../Utilities/statistics.cpp \
#    ../Utilities/smooth3D.cpp \
#    ../Utilities/progressbar.cpp \
#    ../Utilities/lsqfit.cpp \
    ../Utilities/interpolation.cpp \
#    ../Utilities/galmod.cpp \
#    ../Utilities/galfit_min.cpp \
#    ../Utilities/galfit_errors.cpp \
#    ../Utilities/galfit_out.cpp \
#    ../Utilities/galfit.cpp \
#    ../Utilities/conv2D.cpp \
    ../Utilities/fitsUtils.cpp \
    qcustomplot.cpp

HEADERS  += bbarolowindow.h \
#    ../Arrays/cube.hh \
    ../Arrays/header.hh \
#    ../Arrays/image.hh \
    ../Arrays/param.hh \
#    ../Arrays/stats.hh \
#    ../Map/voxel.hh \
#    ../Map/scan.hh \
#    ../Map/objectgrower.hh \
#    ../Map/object3D.hh \
#    ../Map/object2D.hh \
#    ../Map/detection.hh \
    ../Utilities/utils.hh \
#    ../Utilities/spacepar.hh \
#    ../Utilities/ringmodel.hh \
#    ../Utilities/progressbar.hh \
#    ../Utilities/paramguess.hh \
#    ../Utilities/moment.hh \
#    ../Utilities/lsqfit.hh \
#    ../Utilities/gnuplot.hh \
#    ../Utilities/galmod.hh \
#    ../Utilities/galfit.hh \
#    ../Utilities/conv2D.hh \
	qcustomplot.h

FORMS    += bbarolowindow.ui

ICON = resources/Bbarolo.icns

RESOURCES += \
    resources.qrc

INCLUDEPATH += /usr/local/include ../
