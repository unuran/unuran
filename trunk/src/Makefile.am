## Process this file with automake to produce Makefile.in
## $Id$

SUBDIRS = distributions methods specfunct tests uniform

AM_CFLAGS = -Wall

lib_LIBRARIES = libunuran.a

libunuran_a_SOURCES = 

libunuran_a_LIBADD = \
 methods/*.o \
 distributions/*.o \
 specfunct/*.o \
 tests/*.o \
 uniform/*.o


include_HEADERS = \
 unuran.h \
 unuran_config.h

unuran.h: $(top_srcdir)/src/methods/*.h $(top_srcdir)/src/distributions/*.h 
	$(top_srcdir)/scripts/merge_h.pl $(top_srcdir)/src/methods/in_unuran.h > unuran.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 unuran.h
