## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/src

noinst_LIBRARIES = libspecfunct.a

libspecfunct_a_SOURCES = \
 mconf_source.h \
 mtherr.c \
 gamma.c \
 igam.c \
 incbet.c \
 ndtr.c \
 ndtri.c \
 polevl.c

noinst_HEADERS = \
 unur_specfunct_source.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
