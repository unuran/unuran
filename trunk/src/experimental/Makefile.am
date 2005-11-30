## Process this file with automake to produce Makefile.in
## $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/methods

lib_LTLIBRARIES = libunurexperimental.la

libunurexperimental_la_SOURCES = \
 ball.c    ball.h    ball_struct.h \
 hitro.c  hitro.h  hitro_struct.h \
 hitrou.c  hitrou.h  hitrou_struct.h \
 varou.c   varou.h   varou_struct.h \
 walk.c    walk.h    walk_struct.h 

noinst_HEADERS = 

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

