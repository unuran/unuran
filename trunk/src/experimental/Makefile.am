## Process this file with automake to produce Makefile.in
## $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/methods

lib_LTLIBRARIES = libunurexperimental.la

libunurexperimental_la_SOURCES = \
 expect.c          expect_source.h \
 hitrou.c hitrou.h hitrou_struct.h \
 varou.c  varou.h  varou_struct.h 

noinst_HEADERS = 

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

