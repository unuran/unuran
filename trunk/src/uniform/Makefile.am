## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/src

noinst_LIBRARIES = libuniform.a

libuniform_a_SOURCES = mrg31k3p.c fish.c mstd.c

noinst_HEADERS = unur_uniform.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
