## Process this file with automake to produce Makefile.in
## $Id$

INCLUDES = -I$(top_srcdir)/src

noinst_LTLIBRARIES = libuniform.la

libuniform_la_SOURCES = \
 urng_struct.h urng_source.h urng_typedefs.h \
 urng_default.c urng_set.c urng_API.c \
 mrg31k3p.c fish.c mstd.c \
 urng_fvoid.c \
 urng_gsl.c urng_gslqmc.c \
 urng_prng.c \
 urng_randomshift.c \
 urng_rngstreams.c

noinst_HEADERS = \
 urng.h \
 unur_uniform.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
