# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# If autoconf fails with 'possibly undefined macro', run autoreconf -i --force
# Macro can be installed with the autoconf-archive. After that, run aclocal before autoconf
# to pull the required macros into aclocal.m4

AC_PREREQ(2.59)
AC_INIT([BBAROLO], [1.5], [enrico.diteodoro@gmail.com])
#AC_CONFIG_HEADER([src/config.h])
#AC_PREFIX_DEFAULT(/usr/local)

##### ------ Macro that checks for specific modules in python ----- ######
AC_DEFUN([AC_PYTHON_MODULE],
[AC_MSG_CHECKING(for module $1 in python)
echo "import $1" | ${PYTHON} -
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
###########################################################################

######### ------ Macro that checks for MPI implementation ----- ###########
AC_DEFUN([AX_PROG_CXX_MPI], [
AC_PREREQ(2.50)

# Check for compiler
# Needs to be split off into an extra macro to ensure right expansion
# order.
AC_REQUIRE([_AX_PROG_CXX_MPI],[_AX_PROG_CXX_MPI([$1])])

AS_IF([test x"$_ax_prog_cxx_mpi_mpi_wanted" = xno],
  [ _ax_prog_cxx_mpi_mpi_found=no ],
  [
    AC_LANG_PUSH([C++])

    # test whether MPI_Init() is available
    # We do not use AC_SEARCH_LIBS here, as it caches its outcome and
    # thus disallows corresponding calls in the other AX_PROG_*_MPI
    # macros.
    for lib in NONE mpi mpich; do
      save_LIBS=$LIBS
      if test x"$lib" = xNONE; then
        AC_MSG_CHECKING([for function MPI_Init])
      else
        AC_MSG_CHECKING([for function MPI_Init in -l$lib])
        LIBS="-l$lib $LIBS"
      fi
      AC_LINK_IFELSE([
        AC_LANG_PROGRAM([
extern "C" { void MPI_Init(); }
],[MPI_Init();])],
        [ _ax_prog_cxx_mpi_mpi_found=yes ],
        [ _ax_prog_cxx_mpi_mpi_found=no ])
      AC_MSG_RESULT($_ax_prog_cxx_mpi_mpi_found)
      if test "x$_ax_prog_cxx_mpi_mpi_found" = "xyes"; then
        break;
      fi
      LIBS=$save_LIBS
    done

    # Check for header
    AS_IF([test x"$_ax_prog_cxx_mpi_mpi_found" = xyes], [
      AC_MSG_CHECKING([for mpi.h])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <mpi.h>])],
        [ AC_MSG_RESULT(yes)],
        [ AC_MSG_RESULT(no)
         _ax_prog_cxx_mpi_mpi_found=no
      ])
    ])
    AC_LANG_POP([C++])
])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test x"$_ax_prog_cxx_mpi_mpi_found" = xyes], [
        ifelse([$2],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$2])
        :
],[
        $3
        :
])

])dnl AX_PROG_CXX_MPI

dnl _AX_PROG_CXX_MPI is an internal macro required by AX_PROG_CXX_MPI.
dnl To ensure the right expansion order, the main function AX_PROG_CXX_MPI
dnl has to be split into two parts.
dnl
dnl Known MPI C++ compilers:
dnl  mpic++
dnl  mpicxx
dnl  mpiCC
dnl  sxmpic++     NEC SX
dnl  hcp
dnl  mpxlC_r
dnl  mpxlC
dnl  mpixlcxx_r
dnl  mpixlcxx
dnl  mpg++
dnl  mpc++
dnl  mpCC
dnl  cmpic++
dnl  mpiFCC       Fujitsu
dnl  CC
dnl
AC_DEFUN([_AX_PROG_CXX_MPI], [
  AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
  ifelse([$1],,[_ax_prog_cxx_mpi_mpi_wanted=yes],[
    AC_MSG_CHECKING([whether to compile using MPI])
    if $1; then
      _ax_prog_cxx_mpi_mpi_wanted=yes
    else
      _ax_prog_cxx_mpi_mpi_wanted=no
    fi
    AC_MSG_RESULT($_ax_prog_cxx_mpi_mpi_wanted)
  ])
  if test x"$_ax_prog_cxx_mpi_mpi_wanted" = xyes; then
    if test -z "$CXX" && test -n "$MPICXX"; then
      CXX="$MPICXX"
    else
      AC_CHECK_TOOLS([CXX], [mpic++ mpicxx mpiCC sxmpic++ hcp mpxlC_r mpxlC mpixlcxx_r mpixlcxx mpg++ mpc++ mpCC cmpic++ mpiFCC CCicpc pgCC pathCC sxc++ xlC_r xlC bgxlC_r bgxlC openCC sunCC crayCC g++ c++ gpp aCC CC cxx cc++ cl.exe FCC KCC RCC])
    fi
  fi
  AC_PROG_CXX
])dnl _AX_PROG_CXX_MPI
###########################################################################



# Here starts configure

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

# If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard C compiler instead.
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
      [compile with MPI (parallelization) support. If none is found,
      MPI is not used. Default: no])
  ],,[with_mpi=no])
   AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[AC_SUBST([HAVEMPI], "-DHAVE_MPI")],[
    if test x"$with_mpi" = xyes; then
      AC_MSG_WARN([MPI compiler requested, but couldn't use MPI.])
    else
      AC_MSG_WARN([No MPI compiler found, won't use MPI.])
    fi
])

# Setting CFLAGS and CXXFLAGS to user-defined variable for C++11
CFLAGS="$CFLAGS"

AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_LN_S
AC_CHECK_PROGS(INSTALL, install)

# Check for maths library -- will add -lm to $LIBS
#AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_CHECK_HEADERS([iostream fstream])
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

#Check for openMP and set the flag in $OPENMP_CXXFLAGS variable
AC_OPENMP

BASEDIR=`pwd`/src
AC_SUBST([BASEDIR])

# Library and installation utilities.
#------------------------------------
# Static library generation.
AC_PROG_RANLIB

# Shared library generation.
if test "x$ac_cv_c_compiler_gnu" = xyes ; then
  # Note that -fPIC is on by default for Macs, this just makes it obvious.
  SHRFLAGS="-fPIC"
  SHRLD="\$(LINK) $SHRFLAGS"

  case "$build_os" in darwin*)
    AC_MSG_NOTICE([ Using darwin build option])
    SHRLIB="libBBarolo.dylib.$PACKAGE_VERSION"
    SHRLD="$SHRLD -dynamiclib -single_module"
    SHRLD="$SHRLD -compatibility_version $PACKAGE_VERSION -current_version $PACKAGE_VERSION"
    SHRLN="libBBarolo.dylib"

    case "$build_cpu" in powerpc*)
      AC_MSG_NOTICE([ Using darwin build option with powerpc subcase])		
      # Switch off -fPIC (not applicable for Intel Macs).
      CFLAGS="$CFLAGS -mdynamic-no-pic"
      ;;
    esac
    ;;
  *)
    # Covers Linux and Solaris at least.
    AC_MSG_NOTICE([ Using alternative (linux) build option])
    SHRLIB="libBBarolo.so.$PACKAGE_VERSION"
    SHRLD="$SHRLD -shared -Wl,-soname,$SHRLIB"
    SHRLN="libBBarolo.so"
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


# Extra places to look for third-party include files and libraries.
INCDIRS="$INCDIRS           \
         /usr/include       \
         /usr/local/include \
         /opt/local/include"

LIBDIRS="$LIBDIRS           \
         /usr/lib           \
         /usr/local/lib     \
         /local/lib         \
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
AC_CHECK_LIB([cfitsio], [ffopen], [CFITSIOLIB="$CFITSIOLIB -lcfitsio"], [], [$CFITSIOLIB $LIBS])

# CFITSIO is essential for BBarolo, so exit with error message if not present.
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
  # Now add the -L statements to start of CFITSIOLIB.
  # Search for libcfitsio.a.

  DIRLIST=". \
           cfitsio"

  for LIBBASE in $CFITSIOLIBDIRS ; do
    for DIR in $DIRLIST; do
      if test "x$DIR" = x.; then
        LIBDIR=$LIBBASE
      else
        LIBDIR="$LIBBASE/$DIR"
      fi
      AC_CHECK_FILE([$LIBDIR/libcfitsio.a], [CFITSIOLIBFRONT="-L$LIBDIR"; \
                     STATICLIBS+="$LIBDIR/libcfitsio.a "; break 2])
     done
  done
  CFITSIOLIB="$CFITSIOLIBFRONT $CFITSIOLIB"
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
    Ignoring the "--with-fftw3=no" and searching for the library.], [1])
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
                  [CPPFLAGS="$CPPFLAGS -I$INCDIR"; FFTW3INCDIR=$INCDIR; break;],[continue])
  done
fi
],
[ # if directory not specified, go searching.
FFTW3INCDIRS=$INCDIRS
FFTW3LIBDIRS=$LIBDIRS
])


if test "x$FFTW3INCDIR" = x; then
  for INCBASE in $FFTW3INCDIRS ; do
    AC_CHECK_FILE([$INCBASE/fftw3.h], [FFTW3INC="-I$INCBASE $FFTW3INC"; break])
  done
else
  AC_CHECK_FILE([$FFTW3INCDIR/fftw3lib/fftw3.h], [FFTW3INC="-I$FFTW3INCDIR/fftw3lib $FFTW3INC"])
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
             [FFTW3LIBFRONT="-L$LIBDIR"; STATICLIBS+="$LIBDIR/libfftw3.a "; break 2])
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
    AC_CHECK_FILE([$INCBASE/wcslib/wcs.h], [WCSINC="-I$INCBASE $WCSINC"; break])
  done
else
  AC_CHECK_FILE([$WCSINCDIR/wcslib/wcs.h], [WCSINC="-I$WCSINCDIR/wcslib $WCSINC"])
fi


# Look for libwcs.a
AC_CHECK_LIB([wcs], [wcss2p], [WCSLIB="$WCSLIB -lwcs"], [], [$WCSLIB $LIBS $CFITSIOLIB])

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
# Search for libwcs.a

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
             [WCSLIBFRONT="-L$LIBDIR"; STATICLIBS+="$LIBDIR/libwcs.a "; break 2])
    done
  done
  WCSLIB="$WCSLIBFRONT $WCSLIB"
fi

AC_SUBST([WCSINC])
AC_SUBST([WCSLIB])
##########################################################################

AC_SUBST([STATICLIBS])

################### CHECKING FOR PYTHON ##################################
AC_MSG_NOTICE([  -- PYTHON -- ])
AM_PATH_PYTHON
AC_ARG_VAR([PYTHON_INCLUDE], [Include flags for python, bypassing python-config])
AC_ARG_VAR([PYTHON_LIBS], [Library flags for python, bypassing python-config])
AC_ARG_VAR([PYTHON_CONFIG], [Path to python-config])
AS_IF([test -z "$PYTHON_INCLUDE"], [
  AS_IF([test -z "$PYTHON_CONFIG"], [
    AC_PATH_PROGS([PYTHON_CONFIG],
                  [python$PYTHON_VERSION-config python-config],
                  [no],
                  [`dirname $PYTHON`])
    AS_IF([test "$PYTHON_CONFIG" = no], [AC_MSG_WARN([cannot find python-config for $PYTHON.])])
  ])
  AC_MSG_CHECKING([python include flags])
  PYTHON_INCLUDE=`$PYTHON_CONFIG --includes`
  AC_MSG_RESULT([$PYTHON_INCLUDE])
  AC_MSG_CHECKING([python lib flags])
  PYTHON_LIBS=`$PYTHON_CONFIG --ldflags`
  AC_MSG_RESULT([$PYTHON_LIBS])
  AC_SUBST([HAVEPYTHON], "-DHAVE_PYTHON")], [
  AC_SUBST([HAVEPYTHON], "")
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
])
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
      Compiler:          ${CXX}
      Source directory:  ${BASEDIR}
      Target directory:  ${prefix}
      CFITSIO library:   ${CFITSIOLIB}
      CFITSIO headers:   ${CFITSIOINC}
      FFTW3 library:     ${FFTW3LIB}
      FFTW3 headers:     ${FFTW3INC}
      WCS library:       ${WCSLIB}
      WCS headers:       ${WCSINC}
      Have Gnuplot:      ${GNUPLOT}
      Have Python:       ${HAVEPYTHON}
])

AC_MSG_NOTICE([
   --------------------------------------------------------
   CONGRATULATIONS!

   It seems your computer has the needed to compile BBarolo! 
   Now type "make" and pray to your favoured God.  
   --------------------------------------------------------])

