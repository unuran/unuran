## Process this file with automake to produce Makefile.in
## $Id$

AM_CPPFLAGS = \
 -I$(top_srcdir)/src \
 -I$(unuran_srcdir)/src \
 -I$(unuran_srcdir)/src/methods

lib_LTLIBRARIES = libunurexperimental.la

libunurexperimental_la_SOURCES = \
 pinv.c  pinv.h  pinv_struct.h

#libunurexperimental_la_SOURCES = \
# ball.c    ball.h    ball_struct.h \
# varou.c   varou.h   varou_struct.h \
# walk.c    walk.h    walk_struct.h

noinst_HEADERS = 

noinst_PROGRAMS = pinvwh
pinvwh_SOURCES = pinvwh_main.c

pinvwh_LDADD = \
 $(top_builddir)/src/libunurexperimental.la \
 $(unuran_builddir)/src/libunuran.la



# clean backup files
CLEANFILES = \
 unuran.log \
 valgrind-* \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

