## Process this file with automake to produce Makefile.in
## $Id$

SUBDIRS = distributions methods tests uniform

AM_CFLAGS = -Wall

pkglib_LIBRARIES = libunuran.a

libunuran_a_SOURCES = 

libunuran_a_LIBADD = \
 methods/*.o \
 distributions/*.o \
 tests/*.o \
 uniform/*.o

pkginclude_HEADERS = \
 unuran.h \
 unuran_config.h

unuran.h:	$(top_srcdir)/src/distributions/libdistributions.a \
		$(top_srcdir)/src/methods/libmethods.a \
		$(top_srcdir)/src/tests/libtests.a \
		$(top_srcdir)/src/uniform/libuniform.a
	$(top_srcdir)/scripts/merge_h.pl $(top_srcdir)/src/methods/in_unuran.h > unuran.h

# clean backup files too
clean-local:
	rm -f *~

maintainer-clean-local:
	rm -f unuran.h