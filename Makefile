CFLAGS = -O2 -ftree-vectorize -fPIC -Wuninitialized 
#CFLAGS = -O0 -g -fbounds-check -Wall -Wextra 
#CFLAGS = -fast -openmp -Wall

OS = LINUX
CC =    gcc $(CFLAGS)
CXX =   g++ $(CFLAGS)
LINK=   g++ $(CFLAGS) -flto

BASE = ./src

INSTALLDIR = /usr/local/bin
LIBDIR = /usr/local/lib
INCDIR = /usr/local/include/barolo

V = 1.0
EXEC = Bbarolo-$(V)
EXEC-STUB = Bbarolo

EXEC-STUB = Bbarolo


AR = ar cq
LIB = libbarolo-$(V).a
LIB_LN = libbarolo.a
RANLIB = ranlib
SHAREDLIB = $(LINK) -fPIC -shared -Wl,-soname,libbarolo.so.1.0
#SHAREDFLAGS = -shared -Wl,-soname,libbarolo.so.1.0
LIBSO = libbarolo.so.1.0
LIBSO_LN = libbarolo.so

INSTALL = install
LN_S = ln -s

CFITSIOINC = -I/usr/local/include
CFITSIOLIB = -L/usr/local/lib  -lcfitsio

FFTW3INC = -I/usr/include -I/usr/include/fftw3lib 
FFTW3LIB = -L/usr/local/lib  -lfftw3 

WCSINC = -I/usr/include -I/usr/include/wcslib 
WCSLIB = -L/usr/local/lib  -lwcs

OPT = -DHAVE_GNUPLOT -DHAVE_FFTW3 -DHAVE_PYTHON

ifeq ($(OS),MACOSX)
	OPT += -DMACOSX
endif

OTHERLIB = -lm 

CINC = -I$(BASE) $(FFTW3INC) $(CFITSIOINC) $(WCSINC)

LIBS = $(FFTW3LIB) $(CFITSIOLIB) $(WCSLIB) $(OTHERLIB)
STATICLIBS = $(LIBDIR)/libcfitsio.a $(LIBDIR)/libfftw3.a $(LIBDIR)/libwcs.a 

ARRAYDIR = $(BASE)/Arrays
MAPDIR = $(BASE)/Map
UTILDIR = $(BASE)/Utilities
OBJDIR = $(BASE)/Build

HEADS = $(BASE)/bbarolo.hh \
  	$(ARRAYDIR)/param.hh\
	$(ARRAYDIR)/header.hh\
	$(ARRAYDIR)/cube.hh\
	$(ARRAYDIR)/image.hh\
	$(ARRAYDIR)/stats.hh\
	$(MAPDIR)/voxel.hh\
	$(MAPDIR)/scan.hh\
	$(MAPDIR)/object2D.hh\
	$(MAPDIR)/object3D.hh\
	$(MAPDIR)/detection.hh\
	$(MAPDIR)/objectgrower.hh\
	$(UTILDIR)/gnuplot.hh\
	$(UTILDIR)/progressbar.hh\
	$(UTILDIR)/lsqfit.hh\
	$(UTILDIR)/moment.hh\
	$(UTILDIR)/paramguess.hh\
	$(UTILDIR)/conv2D.hh\
	$(UTILDIR)/smooth3D.hh\
	$(UTILDIR)/galmod.hh\
	$(UTILDIR)/galfit.hh\
	$(UTILDIR)/utils.hh\
	$(UTILDIR)/converter.hh

OBJECTS = $(OBJDIR)/bbarolo.o\
	 $(OBJDIR)/header.o\
	 $(OBJDIR)/param.o\
	 $(OBJDIR)/stats.o\
	 $(OBJDIR)/cube.o\
	 $(OBJDIR)/search.o\
	 $(OBJDIR)/image.o\
	 $(OBJDIR)/object2D.o\
	 $(OBJDIR)/object3D.o\
	 $(OBJDIR)/objectgrower.o\
	 $(OBJDIR)/utils.o\
	 $(OBJDIR)/conv2D.o\
	 $(OBJDIR)/smooth3D.o\
	 $(OBJDIR)/galmod.o\
	 $(OBJDIR)/galfit.o\
	 $(OBJDIR)/galfit_min.o\
	 $(OBJDIR)/galfit_errors.o\
	 $(OBJDIR)/galfit_out.o\
	 $(OBJDIR)/slitfit.o\
	 $(OBJDIR)/progressbar.o\
	 $(OBJDIR)/allocator.o\
	 $(OBJDIR)/string.o\
	 $(OBJDIR)/statistics.o\
	 $(OBJDIR)/interpolation.o\
	 $(OBJDIR)/converter.o\
	 $(OBJDIR)/wcsUtils.o\

bbarolo : build $(OBJECTS)
	$(LINK) -o $(EXEC-STUB) $(OBJECTS) $(LIBS) $(OPT)

lib     : $(OBJECTS)
	$(AR) $(LIB) $(OBJECTS)
	$(RANLIB) $(LIB)
	$(SHAREDLIB) $(CINC) $(LIBS) -o $(LIBSO) $(OBJECTS)
#	$(CXX) $(SHAREDFLAGS) -o $(LIBSO) $(OBJECTS)

linux : $(OBJECTS)
	$(LINK) -o Bbarololinux $(OBJECTS) $(LIBS)

install : 
	$(INSTALL) -d -m 2755 $(INSTALLDIR)
	$(INSTALL) -m 755 $(EXEC) $(INSTALLDIR)
	$(RM) $(INSTALLDIR)/$(EXEcC-STUB)
	cd $(INSTALLDIR) && $(LN_S) $(EXEC) $(EXEC-STUB)
	-test ! -f $(LIB) || $(INSTALL) -d -m 2755 $(LIBDIR)
	-test ! -f $(LIB) || $(INSTALL) -m 644 $(LIB) $(LIBDIR)
	-test ! -f $(LIB) || cd $(LIBDIR) && $(RM) $(LIB_LN) && $(LN_S) $(LIB) $(LIB_LN)
	-test ! -f $(LIBSO) || $(INSTALL) -m 755 $(LIBSO) $(LIBDIR) && \
		if [ "libbarolo.so" != "" ] ; then \
			cd $(LIBDIR) && $(RM) $(LIBSO_LN) && $(LN_S) $(LIBSO) $(LIBSO_LN); \
		fi
	$(INSTALL) -d -m 2755 $(INCDIR)
	$(INSTALL) -m 644 $(BASE)/*.hh $(INCDIR)
	$(INSTALL) -m 644 $(BASE)/*.h $(INCDIR)
	$(INSTALL) -d -m 2755 $(INCDIR)/Arrays
	$(INSTALL) -m 644 $(ARRAYDIR)/*.hh $(INCDIR)/Arrays
	$(INSTALL) -d -m 2755 $(INCDIR)/Map
	$(INSTALL) -m 644 $(MAPDIR)/*.hh $(INCDIR)/Map
	$(INSTALL) -d -m 2755 $(INCDIR)/Utilities
	$(INSTALL) -m 644 $(UTILDIR)/*.hh $(INCDIR)/Utilities


$(OBJECTS) : $(HEADS)

$(OBJDIR)/%.o: $(BASE)/%.cpp
	$(CXX) $(OPT) -c $< -o $@ $(CINC) 

$(OBJDIR)/%.o: $(ARRAYDIR)/%.cpp
	$(CXX) $(OPT) -c $< -o $@ $(CINC) 

$(OBJDIR)/%.o: $(MAPDIR)/%.cpp
	$(CXX) $(OPT) -c $< -o $@ $(CINC) 

$(OBJDIR)/%.o: $(UTILDIR)/%.cpp
	$(CXX) $(OPT) -c $< -o $@ $(CINC) 

.cc.o:
	$(CXX) -c $< $(CINC) -o $@ $(OPT)

.c.o:
	$(CC) -c $< $(CINC) -o $@ $(OPT)
	

build: 
	mkdir -p $(OBJDIR)
	
clean : 
	rm -rf $(OBJDIR) $(LIB) $(LIBSO) && \
	cd ./src/GUI && make clean

cleanest: clean
	rm -rf Makefile autom4te.cache config.log config.status src/config.h $(EXEC)

cleangui:
	cd ./src/GUI && make clean

gui :
	cd ./src/GUI  && mv BbaroloGUI.pro BbaroloGUI.or && \
	cp BbaroloGUI.or BbaroloGUI.pro && \
	echo "LIBS += $(LIBS) -static-libgcc -static-libstdc++" >> BbaroloGUI.pro && \
	echo "QMAKE_CXXFLAGS += $(OPT) $(CINC)" >> BbaroloGUI.pro && \
	echo "QMAKE_CXXFLAGS -= -stdlib=libc++" >> BbaroloGUI.pro && \
	echo "QMAKE_CC = $(CC)" >> BbaroloGUI.pro &&\
	echo "QMAKE_CXX = $(CXX)" >> BbaroloGUI.pro &&\
	echo "QMAKE_LINK = $(CXX)" >> BbaroloGUI.pro &&\
	qmake && \
	rm BbaroloGUI.pro && mv BbaroloGUI.or BbaroloGUI.pro && \
	make
ifeq ($(OS),MACOSX)
	cp ../../Bbarolo ../../BbaroloGUI.app/Contents/MacOS/
	cp $(BASE)/GUI/Bbarolo.tiff ./BbaroloGUI.app/Contents/Resources/
endif

all: 
	make && make gui

static: build $(OBJECTS)
		$(LINK) -o $(EXEC-STUB) $(OBJECTS) $(STATICLIBS) $(OPT)
		make guistatic

guistatic :
		rm -rf ./BbaroloGUI.app/Contents/Resources/qt.conf
		macdeployqt ./BbaroloGUI.app/
		cp $(BASE)/GUI/Bbarolo.tiff ./BbaroloGUI.app/Contents/Resources/
		$(LINK) -o ./BbaroloGUI.app/Contents/MacOS/Bbarolo $(OBJECTS) $(STATICLIBS) -static-libgcc -static-libstdc++

guistaticlinux :
	cd ./src/GUI  && mv BbaroloGUI.pro BbaroloGUI.or && \
	cp BbaroloGUI.or BbaroloGUI.pro && \
	echo "LIBS += $(STATICLIBS) -static-libgcc -static-libstdc++ " >> BbaroloGUI.pro && \
	echo "QMAKE_CXXFLAGS += $(OPT) $(CINC)" >> BbaroloGUI.pro && \
	echo "QMAKE_CXXFLAGS -= -stdlib=libc++" >> BbaroloGUI.pro && \
	echo "QMAKE_CC = $(CC)" >> BbaroloGUI.pro &&\
	echo "QMAKE_CXX = $(CXX)" >> BbaroloGUI.pro &&\
	echo "QMAKE_LINK = $(CXX)" >> BbaroloGUI.pro &&\
	export PATH=/home/edt/Downloads/qt/bin/:$(PATH) &&\
	qmake -config release &&\
	rm BbaroloGUI.pro && mv BbaroloGUI.or BbaroloGUI.pro && \
	make && \
	cd ../.. &&\
	$(LINK) -o Bbarolo $(OBJECTS) $(STATICLIBS) -static-libgcc -static-libstdc++
