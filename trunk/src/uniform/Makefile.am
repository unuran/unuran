## Process this file with automake to produce Makefile.in
## $Id$

INCLUDES = -I$(top_srcdir)/src

noinst_LTLIBRARIES = libuniform.la

libuniform_la_SOURCES = \
  urng_default.c \
  mrg31k3p.c fish.c mstd.c

noinst_HEADERS = \
 urng.h \
 unur_uniform.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
