## Process this file with automake to produce Makefile.in
## $Id$

SUBDIRS = distributions methods tests uniform

AM_CFLAGS = -Wall

pkglib_LIBRARIES = libunuran.a

libunuran_a_SOURCES = 

libunuran_a_LIBADD = \
 methods/*.o \
 distributions/*.o   \
 tests/*.o   \
 uniform/*.o

pkginclude_HEADERS = \
 unur_defs.h \
 unuran_config.h

# clean backup files too
clean-local:
	rm -f *~
