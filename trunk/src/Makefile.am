## Process this file with automake to produce Makefile.in
## $Id$

SUBDIRS = distributions methods tests uniform utils

AM_CFLAGS = -Wall

pkglib_LIBRARIES = libunuran.a

libunuran_a_SOURCES = 

libunuran_a_LIBADD = \
 methods/*.o \
 distributions/*.o   \
 utils/*.o   \
 tests/*.o   \
 uniform/*.o

pkginclude_HEADERS = \
 unuran_config.h \
 unur_defs.h


# clean backup files too
clean-local:
	rm -f *~
