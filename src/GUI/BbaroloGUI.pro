QT += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = BBaroloGUI
TEMPLATE = app

CONFIG += qt

BBBASE= $$PWD/../../
DESTDIR= ../../
OBJECTS_DIR= ../Build
MOC_DIR= ../Build

SOURCES += main.cpp\
        bbarolowindow.cpp \
        bbarolowindow_run.cpp \
        qcustomplot.cpp 
        
HEADERS  += bbarolowindow.h \
            qcustomplot.h

FORMS += bbarolowindow.ui

ICON = resources/Bbarolo.icns

RESOURCES += resources.qrc

INCLUDEPATH += /usr/local/include /opt/homebrew/include ../ $$PWD/. ../Arrays/ ../Utilities

#QMAKE_CC = gcc
#QMAKE_CXX = g++

unix: LIBS += -lwcs -lcfitsio -lfftw3
unix: LIBS += -L$$PWD/../.. -lBBarolo-1.7
unix: PRE_TARGETDEPS += $$PWD/../../libBBarolo-1.7.a
