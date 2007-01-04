## Process this file with automake to produce Makefile.in
## $Id$

AM_CPPFLAGS = \
 -I$(top_srcdir)/src

noinst_LTLIBRARIES = libspecfunct.la

# routines taken and slightly modified from 
# Cephes Math Library Release 2.3:  June, 1995
# by Stephen L. Moshier
CEPHES = \
 cephes_source.h \
 mtherr.c \
 gamma.c \
 igam.c \
 incbet.c \
 ndtr.c \
 ndtri.c \
 polevl.c

OTHER = \
 log1p.c

libspecfunct_la_SOURCES = $(CEPHES) $(OTHER)

noinst_HEADERS = \
 unur_specfunct_source.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
