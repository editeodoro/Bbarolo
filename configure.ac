# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([BAROLO], [1.0], [enrico.diteodoro@unibo.it])
#AC_CONFIG_HEADER([src/config.h])
#AC_PREFIX_DEFAULT(/usr/local)

##### ------ Macro that checks for specific modules in python ----- ######
AC_DEFUN([AC_PYTHON_MODULE],
[AC_MSG_CHECKING(for module $1 in python)
echo "import $1" | python -
if test $? -ne 0 ; then
	AC_MSG_RESULT(no)
	AC_MSG_WARN([
    -------------------------------------------------------------------
    WARNING: Cannot find python module '$1' in ordinary PYTHONPATH. If 
             you have '$1' please permanently add its directory to the 
             PYTHONPATH variable in your .profile or .bash_profile with:

              export PYTHONPATH=\$PYTHONPATH:/moduledir     (BASH shell)
              setenv PYTHONPATH \$PYTHONPATH:/moduledir        (C-shell) 

             Then run again ./configure. 
             Otherwise, BBarolo won't produce all the output files.
    -------------------------------------------------------------------
])
	AC_SUBST([PYTHON], "")
else
	AC_MSG_RESULT(yes)
fi
])
##########################################################################

AC_CONFIG_AUX_DIR([config])

# Get the system type.
AC_CANONICAL_BUILD
ARCH="${build_cpu}-$build_os"
AC_SUBST([ARCH])
# Check for which host we are on and setup a few things
# specifically based on the host
case $build_os in
  darwin* )
        # Do something specific for mac
        OP="MACOSX"
				;;
  linux*)
        # Do something specific for linux
        OP="LINUX"
				;;
   *BSD*)
        # Do something specific for BSD
        OP="BSD"
				;;
    *)
        #Default Case
        AC_MSG_ERROR([Your platform is not currently supported])
        ;;
esac
AC_SUBST([OP])


# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_LN_S
AC_CHECK_PROGS(INSTALL, install)

# Check for maths library -- will add -lm to $LIBS
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([unistd.h ctime cmath iostream fstream])
if test "x$ac_cv_header_stdc" = xno; then
  AC_MSG_ERROR([
    -------------------------------------------------------------------
    An ANSI standard C library is required to build BBarolo.  
    One of the ANSI C header files it requires is missing or unusable.

    ERROR: BBarolo configuration failure.
    -------------------------------------------------------------------], [1])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

#Check for openMP and set the flag
AC_OPENMP
CFLAGS=$OPENMP_CXXFLAGS

BASEDIR=./src
AC_SUBST([BASEDIR])

# Library and installation utilities.
#------------------------------------
# Static library generation.
AC_PROG_RANLIB
# Installation utilities.
#AC_PROG_LN_S
#AC_PROG_INSTALL


# Shared library generation.
if test "x$ac_cv_c_compiler_gnu" = xyes ; then
  # Note that -fPIC is on by default for Macs, this just makes it obvious.
  SHRFLAGS="-fPIC"
  SHRLD="\$(LINK) $SHRFLAGS"

  case "$build_os" in
  darwin*)
    AC_MSG_NOTICE([ Using darwin build option])
    SHRLIB="libbarolo.$PACKAGE_VERSION.dylib"
    SHRLD="$SHRLD -dynamiclib -single_module"
    SHRLD="$SHRLD -compatibility_version $PACKAGE_VERSION -current_version $PACKAGE_VERSION"
    SHRLN=

    case "$build_cpu" in
    powerpc*)
      AC_MSG_NOTICE([ Using darwin build option with powerpc subcase])		
      # Switch off -fPIC (not applicable for Intel Macs).
      CFLAGS="$CFLAGS -mdynamic-no-pic"
      ;;
    esac
    ;;
  *)
    # Covers Linux and Solaris at least.
    AC_MSG_NOTICE([ Using alternative (linux) build option])
    SHRLIB="libbarolo.so.$PACKAGE_VERSION"
    SHRLD="$SHRLD -shared -Wl,-soname,$SHRLIB"
    SHRLN="libbarolo.so"
    ;;
  esac

else
    AC_MSG_NOTICE([ Not setting shared library options])
  SHRLIB=
  SHRFLAGS=
  SHRLD=
  SHRSFX=
  SHRLN=
fi

AC_SUBST([SHRLIB])
AC_SUBST([SHRFLAGS])
AC_SUBST([SHRLD])
AC_SUBST([SHRSFX])
AC_SUBST([SHRLN])


# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor cos lround])

# Extra places to look for third-party include files and libraries.
INCDIRS="$INCDIRS           \
         /usr/include       \
         /usr/local/include \
         /usr/local/fftw3   \
         /usr/local/cfitsio \
         /local/cfitsio     \
				 /usr/local/wcslib	\
         /opt/local/include"

LIBDIRS="$LIBDIRS           \
         /usr/lib           \
         /usr/local/lib     \
         /usr/local/cfitsio \
         /local/lib         \
         /local/cfitsio     \
				 /usr/local/wcslib	\
	     	/opt/local/lib"

AC_MSG_NOTICE(LIBDIRS)

for LIBDIR in $LIBDIRS ; do
  AC_CHECK_FILE([$LIBDIR], [LDFLAGS="$LDFLAGS -L$LIBDIR"], [continue])
done

##########################################################################
# Search for CFITSIO.
# use --with-cfitsio=directory to specify a particular CFITSIO directory.
# else do the normal searching for libraries.

AC_MSG_NOTICE([  -- CFITSIO -- ])
AC_ARG_WITH(cfitsio,
	AC_HELP_STRING([--with-cfitsio=<dir>],
		       [The CFITSIO library and include files are in directory <dir>.]),
[ #if we specify a directory...
if test "x$withval" = xno; then
  AC_MSG_NOTICE([
    The CFITSIO library is required.
    Ignoring the "--with-cfitsio=no" and searching for the library.], [1])
  CFITSIOINCDIRS=$INCDIRS
  CFITSIOLIBDIRS=$LIBDIRS
else
  AC_MSG_NOTICE([Using directory "$withval"])
  CFITSIOINCDIRS="$withval \
              $withval/include"
  CFITSIOLIBDIRS="$withval \
              $withval/lib"
  for LIBDIR in $CFITSIOLIBDIRS ; do
    AC_CHECK_FILE([$LIBDIR], [LDFLAGS="$LDFLAGS -L$LIBDIR"], [continue])
  done
  for INCDIR in $CFITSIOINCDIRS ; do
    AC_CHECK_FILE([$INCDIR], [CPPFLAGS="$CPPFLAGS -I$INCDIR"], [continue])
  done
fi
],
[ # if directory not specified, go searching.
CFITSIOINCDIRS=$INCDIRS
CFITSIOLIBDIRS=$LIBDIRS
])

for INCDIR in $CFITSIOINCDIRS ; do
  AC_CHECK_FILE([$INCDIR/fitsio.h], [CFITSIOINC=-I$INCDIR; break])
  INCDIR=$INCDIR/cfitsio
  AC_CHECK_FILE([$INCDIR/fitsio.h], [CFITSIOINC=-I$INCDIR; break])
done

AC_CHECK_LIB([socket],  [recv],   [CFITSIOLIB="$CFITSIOLIB -lsocket"], [], [$LIBS])
AC_CHECK_LIB([cfitsio], [ffopen], [CFITSIOLIB="$CFITSIOLIB -lcfitsio"], [],
             [$CFITSIOLIB $LIBS])

# CFITSIO is essential for BBarolo, so exit with error message if it
#   is not present.
if test "x$CFITSIOINC" = x -o "x$CFITSIOLIB" = x; then
  AC_MSG_ERROR([
    -------------------------------------------------------------
    Could not find the CFITSIO library.

    ERROR: BBarolo configuration failure. If you have the libra-
           ry in non ordinary path, try to run ./configure with 
           --with-cfitsio=<dir> option
    -------------------------------------------------------------], [1])
else
  AC_MSG_NOTICE([CFITSIO appears to be available.])
  AC_DEFINE([HAVE_CFITSIO], [1], [Define to 1 if CFITSIO is available.])
for LIBDIR in $CFITSIOLIBDIRS ; do
  AC_CHECK_FILE([$LIBDIR/libcfitsio.a], [CFITSIOLIB="-L$LIBDIR $CFITSIOLIB"; break])
  LIBDIR=$LIBDIR/cfitsio
  AC_CHECK_FILE([$LIBDIR/libcfitsio.a], [CFITSIOLIB="-L$LIBDIR $CFITSIOLIB"; break])
done

fi

AC_SUBST([CFITSIOINC])
AC_SUBST([CFITSIOLIB])
##########################################################################

##########################################################################
# Search for FFTW3#.
# use --with-fftw3=directory to specify a particular FFTW3 directory.
# else do the normal searching for libraries.

AC_MSG_NOTICE([  -- FFTW3 -- ])
FFTW3INCDIR=""
AC_ARG_WITH(fftw3,
	AC_HELP_STRING([--with-fftw3=<dir>],
		       [The FFTW3 library and include files are in directory <dir>.]),
[ #if we specify a directory...
if test "x$withval" = xno; then
  AC_MSG_NOTICE([
    The FFTW3 library is required.
    Ignoring the "--with-fftw3lib=no" and searching for the library.], [1])
  FFTW3INCDIRS=$INCDIRS
  FFTW3LIBDIRS=$LIBDIRS
else
  AC_MSG_NOTICE([Using directory "$withval"])
  FFTW3INCDIRS="$withval \
              $withval/include"
  FFTW3LIBDIRS="$withval \
              $withval/lib"
  for LIBDIR in $FFTW3LIBDIRS ; do
    AC_CHECK_FILE([$LIBDIR], [LDFLAGS="$LDFLAGS -L$LIBDIR"], [continue])
  done
  for INCDIR in $FFTW3INCDIRS ; do
    AC_CHECK_FILE([$INCDIR/fftw3lib], 
                  [CPPFLAGS="$CPPFLAGS -I$INCDIR"; FFTW3INCDIR=$INCDIR; break], 
		  [continue])
  done
fi
],
[ # if directory not specified, go searching.
FFTW3INCDIRS=$INCDIRS
FFTW3LIBDIRS=$LIBDIRS
])


if test "x$FFTW3INCDIR" = x; then
  for INCBASE in $FFTW3INCDIRS ; do
    AC_CHECK_FILE([$INCBASE/fftw3.h], [FFTW3INC="-I$INCBASE -I$INCBASE/fftw3lib $FFTW3INC"; break])
  done
else
  AC_CHECK_FILE([$FFTW3INCDIR/fftw3lib/fftw3.h], [FFTW3INC="-I$FFTW3INCDIR -I$FFTW3INCDIR/fftw3lib $FFTW3INC"])
fi


# Look for libfftw3.a
AC_CHECK_LIB([fftw3], [fftw_execute], [FFTW3LIB="$FFTW3LIB -lfftw3"], [],
             [$FFTW3LIB $LIBS $CFITSIOLIB])
             
#AC_CHECK_LIB([fftw3_threads], [fftw_init_threads], [FFTW3LIB="$FFTW3LIB -lfftw3_omp"], [],
#             [$FFTW3LIB $LIBS $CFITSIOLIB])


# FFTW3 is essential for BBarolo, so exit with error message if it
#   is not present. Need at least libfftw3 to be available.
if test "x$FFTW3INC" = x -o "x$FFTW3LIB" = x; then
  AC_MSG_ERROR([
    ------------------------------------------------------------
    Could not find the FFTW3 library.

    ERROR: BBarolo configuration failure. If you have the libra-
           ry in non ordinary path, try to run ./configure with 
           --with-fftw3=<dir> option.
    ------------------------------------------------------------], [1])
else
  AC_MSG_NOTICE([FFTW3LIB appears to be available.])
  AC_DEFINE([HAVE_FFTW3LIB], [1], [Define to 1 if FFTW3LIB is available.])
# Now add the -L statements to start of FFTW3LIB.
# Search for libfftw3.a.

  DIRLIST=". \
           fftw3 \
           fftw3lib"

  for LIBBASE in $FFTW3LIBDIRS ; do
    for DIR in $DIRLIST; do
      if test "x$DIR" = x.; then
       LIBDIR=$LIBBASE
      else
       LIBDIR="$LIBBASE/$DIR"
      fi
      AC_CHECK_FILE([$LIBDIR/libfftw3.a], 
             [FFTW3LIBFRONT="-L$LIBDIR"; break])
    done
  done
  FFTW3LIB="$FFTW3LIBFRONT $FFTW3LIB"
fi


AC_SUBST([FFTW3INC])
AC_SUBST([FFTW3LIB])
AC_SUBST([FFTW3], "-DHAVE_FFTW3")
##########################################################################


##########################################################################
# Search for WCSLIB.
# use --with-wcslib=directory to specify a particular WCSLIB directory.
# else do the normal searching for libraries.

AC_MSG_NOTICE([  -- WCSLIB -- ])
WCSINCDIR=""
AC_ARG_WITH(wcslib,
	AC_HELP_STRING([--with-wcslib=<dir>],
		       [The WCSLIB library and include files are in directory <dir>.]),
[ #if we specify a directory...
if test "x$withval" = xno; then
  AC_MSG_NOTICE([
    The WCSLIB library is required.
    Ignoring the "--with-wcslib=no" and searching for the library.], [1])
  WCSINCDIRS=$INCDIRS
  WCSLIBDIRS=$LIBDIRS
else
  AC_MSG_NOTICE([Using directory "$withval"])
  WCSINCDIRS="$withval \
              $withval/include"
  WCSLIBDIRS="$withval \
              $withval/lib"
  for LIBDIR in $WCSLIBDIRS ; do
    AC_CHECK_FILE([$LIBDIR], [LDFLAGS="$LDFLAGS -L$LIBDIR"], [continue])
  done
  for INCDIR in $WCSINCDIRS ; do
    AC_CHECK_FILE([$INCDIR/wcslib], 
                  [CPPFLAGS="$CPPFLAGS -I$INCDIR"; WCSINCDIR=$INCDIR; break], 
		  [continue])
  done
fi
],
[ # if directory not specified, go searching.
WCSINCDIRS=$INCDIRS
WCSLIBDIRS=$LIBDIRS
])


if test "x$WCSINCDIR" = x; then
  for INCBASE in $WCSINCDIRS ; do
    AC_CHECK_FILE([$INCBASE/wcslib/wcs.h], [WCSINC="-I$INCBASE -I$INCBASE/wcslib $WCSINC"; break])
  done
else
  AC_CHECK_FILE([$WCSINCDIR/wcslib/wcs.h], [WCSINC="-I$WCSINCDIR -I$WCSINCDIR/wcslib $WCSINC"])
fi


# Look for libwcs.a, and libpgsbox.a
AC_CHECK_LIB([wcs], [wcss2p], [WCSLIB="$WCSLIB -lwcs"], [],
             [$WCSLIB $LIBS $CFITSIOLIB $PGPLOTLIB])

# Only use libpgsbox if we have pgplot enabled - determine this by testing the PGPLOTINCDIRS variable
if test "x$PGPLOTINCDIRS" != x; then 
     AC_CHECK_LIB([pgsbox], [cpgsbox], [WCSLIB="$WCSLIB -lpgsbox"], [],
                 [$WCSLIB $LIBS $CFITSIOLIB $PGPLOTLIB])
fi

if test "x$WCSINC" = x -o "x$WCSLIB" = x; then
  AC_MSG_ERROR([
    ------------------------------------------------------------
    Could not find the WCSLIB library.

    ERROR: BBarolo configuration failure. If you have the libra-
           ry in non ordinary path, try to run ./configure with 
           --with-wcslib=<dir> option.
    ------------------------------------------------------------], [1])
else
  AC_MSG_NOTICE([WCSLIB appears to be available.])
  AC_DEFINE([HAVE_WCSLIB], [1], [Define to 1 if WCSLIB is available.])
# Now add the -L statements to start of WCSLIB.
# Search for libwcs.a, and, if we find it, search for libpgsbox in same directory and
#   then in other directories, keeping the appropriate order of -L calls.

  DIRLIST=". \
           wcs \
           wcslib"

  for LIBBASE in $WCSLIBDIRS ; do
    for DIR in $DIRLIST; do
      if test "x$DIR" = x.; then
       LIBDIR=$LIBBASE
      else
       LIBDIR="$LIBBASE/$DIR"
      fi
      AC_CHECK_FILE([$LIBDIR/libwcs.a], 
             [WCSLIBFRONT="-L$LIBDIR"; 
             AC_CHECK_FILE([$LIBDIR/libpgsbox.a], [], [
               for LIBBASE2 in $WCSLIBDIRS ; do
                for DIR2 in $DIRLIST; do
                 if test "x$DIR" = x.; then
                  LIBDIR=$LIBBASE
                 else
                  LIBDIR="$LIBBASE/$DIR"
                 fi
                 AC_CHECK_FILE([$LIBDIR2/libpgsbox.a],
                               [WCSLIBFRONT="$WCSLIB -L$LIBDIR2"; break])
                done
               done])
              break])
    done
  done
  WCSLIB="$WCSLIBFRONT $WCSLIB"
fi

AC_SUBST([WCSINC])
AC_SUBST([WCSLIB])
##########################################################################

################### CHECKING FOR PYTHON ##################################
AC_MSG_NOTICE([  -- PYTHON -- ])
AC_CHECK_PROG(PYTHON, python, yes,no)
if test "x$PYTHON" = "xyes"; then
  AC_SUBST([PYTHON], "-DHAVE_PYTHON")
else 
  AC_SUBST([PYTHON], "")
	AC_MSG_WARN([
    -------------------------------------------------------------------
    WARNING: Cannot find python in ordinary PATH. If you have python,
             please permanently add its directory to the PATH variable 
             in your .profile or .bash_profile with:

              export PATH=\$PATH:/yourpythondir            (BASH shell)
              setenv PATH \$PATH:/yourpythondir               (C-shell) 

             Then run again ./configure. 
             Otherwise, BBarolo won't produce all the output files.
    -------------------------------------------------------------------
])
fi
AC_PYTHON_MODULE(matplotlib)
AC_PYTHON_MODULE(astropy)
##########################################################################


################### CHECKING FOR GNUPLOT ################################
AC_MSG_NOTICE([  -- GNUPLOT -- ])
AC_CHECK_PROG(GNUPLOT, gnuplot, yes,no)
if test "x$GNUPLOT" = "xyes"; then
  AC_SUBST([GNUPLOT], "-DHAVE_GNUPLOT")
else 
  AC_SUBST([GNUPLOT], "")
	AC_MSG_WARN([

    -------------------------------------------------------------------
    WARNING: Cannot find gnuplot in ordinary PATH. If you have gnuplot,
             please permanently add its directory to the PATH variable 
             in your .profile or .bash_profile with:

              export PATH=\$PATH:/yourgnuplotdir           (BASH shell)
              setenv PATH \$PATH:/yourgnuplotdir              (C-shell) 
             
             or create a link to ordinary path with:

              sudo ln -s /yourgnuplotdir/gnuplot /usr/local/bin/gnuplot

             Then run again ./configure. 
             Gnuplot is however unnecessary if python is installed.
    -------------------------------------------------------------------
])
fi
##########################################################################





AC_SUBST([LDFLAGS])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

AC_MSG_NOTICE([

   Configuration Summary:
      Operating System:  ${build_os}
      Compilator:        ${CXX}
      Source directory:  ${BASEDIR}
      Target directory:  ${prefix}
      CFITSIO library:   ${CFITSIOLIB}
      FFTW3 library:     ${FFTW3LIB}
      WCS library:       ${WCSLIB}
      Have Gnuplot:      ${GNUPLOT}
      Have Python:       ${PYTHON}
])

AC_MSG_NOTICE([
   --------------------------------------------------------
   CONGRATULATIONS!

   It seems your computer has the needed to compile BBarolo! 
   Now type "make" and pray your favoured God.  
   --------------------------------------------------------])

