## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/distributions

noinst_LIBRARIES = libspecfunct.a

libspecfunct_a_SOURCES = \
 mtherr.c \
 gamma.c \
 igam.c \
 incbet.c \
 ndtr.c \
 ndtri.c \
 polevl.c

noinst_HEADERS = 

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
