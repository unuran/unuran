## Process this file with automake to produce Makefile.in
## $Id$

INCLUDES = -I$(top_srcdir)/src

noinst_LTLIBRARIES = libdistr.la

libdistr_la_SOURCES = \
 distr.c  distr.h  distr_struct.h distr_source.h \
 cont.c   cont.h \
 corder.c corder.h \
 cemp.c   cemp.h \
 cvec.c   cvec.h \
 cvemp.c  cvemp.h\
 discr.c  discr.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

