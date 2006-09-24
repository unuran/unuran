## Process this file with automake to produce Makefile.in
## $Id$

INCLUDES = -I$(top_srcdir)/src

noinst_LTLIBRARIES = liburng.la

liburng_la_SOURCES = \
 urng_source.h urng_struct.h urng_typedefs.h \
 urng.h urng.c urng_set.c

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
