dnl Process this file with autoconf to produce a configure script.
dnl $Id$

dnl Initialize.
AC_PREREQ(2.57)
AC_INIT([UNURAN], [0.4.5], [unuran@statistik.wu-wien.ac.at])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/unuran_config.h])
AC_CONFIG_HEADER([config.h])

dnl Library versioning (current:revision:age)
dnl See the libtool manual for an explanation of the numbers
dnl 
dnl unuran-0.4.0   0:0:0
dnl unuran-0.4.2   1:0:1
dnl unuran-0.4.4   1:0:2

UNURAN_LT_VERSION="1:0:2"
AC_SUBST([UNURAN_LT_VERSION])

dnl Add dmalloc support.
AM_WITH_DMALLOC

dnl Check for which system.
AC_CANONICAL_HOST

dnl Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_CHECK_TOOL(AR, ar, :)
# AC_CHECK_TOOL(RANLIB, ranlib, :)
AC_PROG_RANLIB
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

dnl Checks for libraries.
dnl libm: the standard math library
AC_CHECK_LIB([m], [sqrt])
dnl prng: Otmar Lendl's library for pseudo random number generators
AC_CHECK_LIB([prng], [prng_new])
AC_CHECK_HEADERS([prng.h])
dnl RngStream: Pierre L'Ecuyers multiple random stream generator
AC_CHECK_LIB([rngstream], [RngStream_CreateStream])
AC_CHECK_HEADERS([RngStream.h])
dnl GSL: GNU Scientific Library
dnl AC_CHECK_LIB([gslcblas], [gsl_blas_ddot])
AC_CHECK_LIB([gslcblas], [main])
AC_CHECK_LIB([gsl], [gsl_poly_eval])
AC_CHECK_HEADERS([gsl/gsl_rng.h])

dnl cephes: Stephen L. Moshier's library of mathematical functions 
dnl However, at the moment we do not use it any more.
dnl AC_CHECK_LIB([md], [lgam])

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h sys/time.h unistd.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME

dnl Checks for library functions.
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strcasecmp strchr strtol strtoul])

dnl Checks for GNU extensions
AC_CHECK_DECLS([getopt, vsnprintf])

dnl We need Mathematica(R) for some tests
AC_CHECK_PROG(have_mathematica,math,true,false) 
AM_CONDITIONAL(HAVE_MATHEMATICA, test x$have_mathematica = xtrue)
dnl AM_CONDITIONAL(HAVE_MATHEMATICA,true)
dnl AM_CONDITIONAL(HAVE_MATHEMATICA,false)

dnl Set flags for compiler
if test X"$GCC" = Xyes ; then
	AC_SUBST(AM_CFLAGS,-Wall)
fi


AC_CONFIG_FILES([\
	Makefile \
        scripts/Makefile \
	src/Makefile \
	src/distr/Makefile \
	src/distributions/Makefile \
	src/methods/Makefile \
	src/parser/Makefile \
	src/specfunct/Makefile \
	src/tests/Makefile \
	src/uniform/Makefile \
	src/utils/Makefile \
	doc/Makefile \
        tests/Makefile \
	examples/Makefile \
	experiments/Makefile ])
AC_OUTPUT

dnl -----------------------------------------




