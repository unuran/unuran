dnl Process this file with autoconf to produce a configure script.
dnl $Id$

dnl Initialize.
AC_PREREQ(2.60)
dnl AC_INIT([UNURAN],[0.9.0],[unuran@statistik.wu-wien.ac.at])
AC_INIT([UNURAN],[devel],[unuran@statistik.wu-wien.ac.at])
AC_CONFIG_AUX_DIR(autoconf)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/unuran_config.h])
AC_CONFIG_HEADER([config.h])

dnl Library versioning (current:revision:age)
dnl See the libtool manual  for an explanation of the numbers
dnl 
dnl unuran-0.4.0   0:0:0
dnl unuran-0.4.2   1:1:0
dnl unuran-0.4.4   1:2:0
dnl unuran-0.5.0   2:0:0
dnl unuran-0.6.0   3:0:0
dnl unuran-0.7.0   4:0:0
dnl unuran-0.7.1   4:1:0
dnl unuran-0.7.2   4:2:0
dnl unuran-0.8.1   5:0:0
dnl unuran-0.9.x   6:0:0

UNURAN_LT_VERSION="6:0:0"
AC_SUBST([UNURAN_LT_VERSION])

dnl Whether experimental code should be included
dnl (set to false before releasing a new version)
AM_CONDITIONAL(INCLUDE_EXPERIMENTAL,true)
dnl AM_CONDITIONAL(INCLUDE_EXPERIMENTAL,false)

dnl Whether we distribute source code only
dnl (without documentation, tests, and examples)
dnl DO NOT CHANGE THIS HERE!
AM_CONDITIONAL(INCLUDE_SOURCE_ONLY,false)

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

dnl Optional support for external uniform random number generators

dnl RngStream: Pierre L'Ecuyer's multiple random stream generator
AC_ARG_WITH(urng_rngstream,
	[AS_HELP_STRING([--with-urng-rngstream],
	   [URNG: use Pierre L'Ecuyer's RNGSTREAM library @<:@default=no@:>@])],
	[case "${withval}" in
	   yes) ;; no) ;; 
	   *) AC_MSG_ERROR(bad value '${withval}' for --with-urng-rngstream) ;;
	esac],
	[with_urng_rngstream=no])
AS_IF([test "x$with_urng_rngstream" != xno],
	[AC_CHECK_LIB([rngstreams], [RngStream_CreateStream], [], [])
	 AC_CHECK_HEADERS([RngStream.h],
	   [AC_SUBST(UNURAN_SUPPORTS_RNGSTREAM,-DUNURAN_SUPPORTS_RNGSTREAM)
	    AC_DEFINE([UNURAN_HAS_RNGSTREAM], [1], 
             [Define to 1 if you use Pierre L'Ecuyer's RNGSTREAM library.])],
	   [AC_MSG_FAILURE([
****************************************************************
*                                                              *
*   --with-urng-rngstream given, but 'RngStream.h' not found   *
*                                                              *
****************************************************************
	   ]) ]) ], [])


dnl prng: Otmar Lendl's library for pseudo random number generators
AC_ARG_WITH(urng_prng,
	[AS_HELP_STRING([--with-urng-prng],
	   [URNG: use Otmar Lendl's PRNG library @<:@default=no@:>@])],
	[case "${withval}" in
	   yes) ;; no) ;; 
	   *) AC_MSG_ERROR(bad value '${withval}' for --with-urng-prng) ;;
	esac],
	[with_urng_prng=no])
AS_IF([test "x$with_urng_prng" != xno],
	[AC_CHECK_LIB([prng], [prng_new], [], [])
	 AC_CHECK_HEADERS([prng.h],
	   [AC_SUBST(UNURAN_SUPPORTS_PRNG,-DUNURAN_SUPPORTS_PRNG)
	    AC_DEFINE([UNURAN_HAS_PRNG], [1], 
             [Define to 1 if you use Otmar Lendl's PRNG library.])],
	   [AC_MSG_FAILURE([
******************************************************
*                                                    *
*   --with-urng-prng given, but 'prng.h' not found   *
*                                                    *
******************************************************
	   ]) ]) ], [])


dnl GSL: GNU Scientific Library
AC_ARG_WITH(urng_gsl,
	[AS_HELP_STRING([--with-urng-gsl],
	   [URNG: use random number generators from GNU Scientific Library @<:@default=no@:>@])],
	[case "${withval}" in
	   yes) ;; no) ;; 
	   *) AC_MSG_ERROR(bad value '${withval}' for --with-urng-gsl) ;;
	esac],
	[with_urng_gsl=no])
AS_IF([test "x$with_urng_gsl" != xno],
	[AC_CHECK_LIB([gslcblas], [main])
	 AC_CHECK_LIB([gsl], [gsl_rng_alloc], [], [])
	 AC_CHECK_HEADERS([gsl/gsl_rng.h gsl/gsl_qrng.h],
	   [AC_SUBST(UNURAN_SUPPORTS_GSL,-DUNURAN_SUPPORTS_GSL)
	    AC_DEFINE([UNURAN_HAS_GSL], [1], 
             [Define to 1 if you use uniform random number generators from GSL.])],
	   [AC_MSG_FAILURE([
********************************************************
*                                                      *
*   --with-urng-gsl given, but 'gsl_rng.h' not found   *
*                                                      *
********************************************************
	   ]) ]) ], [])


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
   inf = exp(1.0e10);
   nan = inf / inf;
   /* nan = 0.0 / 0.0; */
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

dnl Some routines are deprecated.
dnl Compiling and linking these into the library can be controlled
dnl by a configure flag. 
dnl default is disabled.

AC_ARG_ENABLE(deprecated,
	[AS_HELP_STRING([--enable-deprecated],
		[enable support for deprecated routines @<:@default=use both@:>@])],
	[],
	[enable_deprecated=no])
AM_CONDITIONAL(ENABLE_DEPRECATED, test x$enable_deprecated = xyes)

dnl Some tests compare the return values of PDFs, CDFs and of
dnl the function parser with the output of corresponding
dnl Mathematica(R) scripts.
dnl Thus for these tests we need Mathematica(R).
dnl However, simply checking for the existance of programm 'math'
dnl breaks the configure script when 'math' is available but
dnl ask for a password. 
dnl Hence these tests are enabled by means of a configure option.
AC_ARG_WITH(mathematica,
	[AS_HELP_STRING([--with-mathematica],
		[Use Mathematica(R) for some checks @<:@default=no@:>@])],
	[case "${withval}" in
	  yes) mathematica=true ;;
	  no)  mathematica=false ;;
	  *) AC_MSG_ERROR(bad value '${withval}' for --with-mathematica) ;;
	esac],
	[mathematica=false])
AM_CONDITIONAL(HAVE_MATHEMATICA, test x$mathematica = xtrue)

dnl check for Valgrind
dnl (Valgrind is a GPL'd system for debugging and profiling x86-Linux programs)
AC_CHECK_PROG(have_valgrind,valgrind,true,false) 
AM_CONDITIONAL(HAVE_VALGRIND, test x$have_valgrind = xtrue)

dnl Set flags for compiler
if test X"$GCC" = Xyes ; then
	AC_SUBST(AM_CFLAGS,-Wall)
fi

dnl create Makefiles for source
AC_CONFIG_FILES([\
	Makefile \
	src/Makefile \
	src/distr/Makefile \
	src/distributions/Makefile \
	src/methods/Makefile \
	src/parser/Makefile \
	src/specfunct/Makefile \
	src/tests/Makefile \
	src/uniform/Makefile \
	src/urng/Makefile \
	src/utils/Makefile ] )

dnl create Makefiles for auxiliary scripts (optional)
AC_CHECK_FILE(${srcdir}/scripts/Makefile.am,
	AC_CONFIG_FILES( [ \
	   scripts/Makefile ] ))

dnl create Makefiles for documentation (optional)
AC_CHECK_FILE(${srcdir}/doc/Makefile.am,
	AC_CONFIG_FILES( [ \
	   doc/Makefile \
	   doc/src/Makefile \
	   doc/figures/Makefile ] ))

dnl create Makefiles for tests (optional)
AC_CHECK_FILE(${srcdir}/tests/Makefile.am,
	AC_CONFIG_FILES( [ \
	   tests/Makefile \
	   tests/testdistributions/Makefile ] ))

dnl create Makefiles for examples (optional)
AC_CHECK_FILE(${srcdir}/examples/Makefile.am,
	AC_CONFIG_FILES( [ \
	   examples/Makefile ] ))
AC_CHECK_FILE(${srcdir}/experiments/Makefile.am,
	AC_CONFIG_FILES( [ \
	   experiments/Makefile ] ))

dnl Makefiles for experimental code (optional)
AC_CHECK_FILE(${srcdir}/src/experimental/Makefile.am,
	AC_CONFIG_FILES( [ \
	   src/experimental/Makefile ] ))
AC_CHECK_FILE(${srcdir}/run/Makefile.am,
	AC_CONFIG_FILES( [ \
	   run/Makefile ] ))

AC_OUTPUT

dnl -----------------------------------------




