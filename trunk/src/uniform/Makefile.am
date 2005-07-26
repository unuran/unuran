## Process this file with automake to produce Makefile.in
## $Id$

INCLUDES = -I$(top_srcdir)/src

noinst_LTLIBRARIES = libuniform.la

libuniform_la_SOURCES = \
 urng_struct.h urng_source.h urng_typedefs.h \
 urng_default.c urng_set.c urng_API.c \
 mrg31k3p.c fish.c mstd.c \
 urng_fvoid.c urng_fvoid.h \
 urng_gsl.c urng_gsl.h \
 urng_gslqrng.c urng_gslqrng.h \
 urng_prng.c urng_prng.h \
 urng_randomshift.c urng_randomshift.h \
 urng_rngstreams.c urng_rngstreams.h

noinst_HEADERS = \
 urng.h \
 unur_uniform.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
