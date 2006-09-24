dnl Process this file with autoconf to produce a configure script.
dnl $Id$

dnl Initialize.
AC_PREREQ(2.59)
AC_INIT([UNURAN],[0.7.3],[unuran@statistik.wu-wien.ac.at])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/unuran_config.h])
AC_CONFIG_HEADER([config.h])

dnl Library versioning (current:revision:age)
dnl See the libtool manual for an explanation of the numbers
dnl 
dnl unuran-0.4.0   0:0:0
dnl unuran-0.4.2   1:0:1
dnl unuran-0.4.4   1:0:2
dnl unuran-0.5.0   2:0:0
dnl unuran-0.6.0   3:0:0
dnl unuran-0.7.0   4:0:0
dnl unuran-0.7.1   4:0:1
dnl unuran-0.7.2   4:0:2

UNURAN_LT_VERSION="4:0:3"
AC_SUBST([UNURAN_LT_VERSION])

dnl Whether experimental code should be included
dnl (set to false before releasing a new version)
AM_CONDITIONAL(INCLUDE_EXPERIMENTAL,true)
dnl AM_CONDITIONAL(INCLUDE_EXPERIMENTAL,false)

dnl Add dmalloc support.
AM_WITH_DMALLOC

dnl Check for which system.
AC_CANONICAL_HOST

dnl Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_CC
AC_PROG_CXX
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
AC_CHECK_LIB([rngstreams], [RngStream_CreateStream])
AC_CHECK_HEADERS([RngStream.h])
dnl GSL: GNU Scientific Library
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

dnl For the following function replacements exist
AH_BOTTOM([/* Substitute UNU.RAN functions for missing system functions */])

AC_CHECK_DECLS(log1p,,,[#include <math.h>])
AH_BOTTOM([#if !HAVE_DECL_LOG1P
#define log1p _unur_log1p
#endif])

dnl Checks for GNU extensions
AC_CHECK_DECLS([getopt, vsnprintf])

dnl Check for C99 functions */
AC_CHECK_DECLS(isfinite,,,[#include <math.h>])
AC_CHECK_DECLS(isinf,,,[#include <math.h>])
AC_CHECK_DECLS(isnan,,,[#include <math.h>])
AC_CHECK_DECLS(INFINITY,,,[#include <math.h>])
AC_CHECK_DECLS(HUGE_VAL,,,[#include <math.h>])

dnl Check IEEE comparisons, whether "x != x" is true for NaNs
AC_CACHE_CHECK([for IEEE comparisons], ac_cv_c_ieee_comparisons,
[AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <math.h>
int main (void) 
{ 
   int status; double inf, nan;
   /* inf = exp(1.0e10); */
   /* nan = inf / inf ; */
   nan = 0.0 / 0.0;
   status = (nan == nan);
   exit (status);
}]])],
[ac_cv_c_ieee_comparisons="yes"],
[ac_cv_c_ieee_comparisons="no"],
[ac_cv_c_ieee_comparisons="yes"])
])

if test "$ac_cv_c_ieee_comparisons" != no ; then
  AC_DEFINE([HAVE_IEEE_COMPARISONS], [1], 
      [Define to 1 if "x != x" is true for NaNs])
  AC_SUBST(HAVE_IEEE_COMPARISONS)
fi


dnl We need Mathematica(R) for some tests
AC_CHECK_PROG(have_mathematica,math,true,false) 
AM_CONDITIONAL(HAVE_MATHEMATICA, test x$have_mathematica = xtrue)
dnl AM_CONDITIONAL(HAVE_MATHEMATICA,true)
dnl AM_CONDITIONAL(HAVE_MATHEMATICA,false)

dnl check for Valgrind
dnl (Valgrind is a GPL'd system for debugging and profiling x86-Linux programs)
AC_CHECK_PROG(have_valgrind,valgrind,true,false) 
AM_CONDITIONAL(HAVE_VALGRIND, test x$have_valgrind = xtrue)

dnl Set flags for compiler
if test X"$GCC" = Xyes ; then
	AC_SUBST(AM_CFLAGS,-Wall)
	AC_SUBST(AM_CPPFLAGS,-Wall)
fi

dnl create Makefiles
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
	src/urng/Makefile \
	src/utils/Makefile \
	doc/Makefile \
	doc/src/Makefile \
	doc/figures/Makefile \
        tests/Makefile \
        tests/testdistributions/Makefile \
	examples/Makefile \
	experiments/Makefile ])

dnl Makefiles for experimental code
AC_CHECK_FILE(${srcdir}/src/experimental/Makefile.am, 
	AC_CONFIG_FILES(src/experimental/Makefile))
AC_CHECK_FILE(${srcdir}/run/Makefile.am, 
	AC_CONFIG_FILES(run/Makefile))

AC_OUTPUT

dnl -----------------------------------------




