## Process this file with automake to produce Makefile.in
## $Id$

AM_CPPFLAGS = \
 -I$(top_srcdir)/src \
 -I$(unuran_srcdir)/src \
 -I$(unuran_srcdir)/src/methods

lib_LTLIBRARIES = libunurexperimental.la

libunurexperimental_la_SOURCES = \
 cext.c    cext.h    cext_struct.h \
 mvtdr.c   mvtdr.h   mvtdr_struct.h \
   mvtdr_newset.ch mvtdr_init.ch mvtdr_sample.ch mvtdr_debug.ch

#libunurexperimental_la_SOURCES = \
# ball.c    ball.h    ball_struct.h \
# varou.c   varou.h   varou_struct.h \
# walk.c    walk.h    walk_struct.h \
# mvtdr.c   mvtdr.h   mvtdr_struct.h \
#   mvtdr_newset.ch mvtdr_init.ch mvtdr_sample.ch mvtdr_debug.ch

noinst_HEADERS = 

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

