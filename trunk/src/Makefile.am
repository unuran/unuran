## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

SUBDIRS = utils methods specfunct distr distributions parser tests uniform .

SUBLIBS = \
 methods/libmethods.la \
 distr/libdistr.la \
 distributions/libdistributions.la \
 specfunct/libspecfunct.la \
 parser/libparser.la \
 tests/libtests.la \
 uniform/libuniform.la \
 utils/libutils.la

lib_LTLIBRARIES = libunuran.la
libunuran_la_SOURCES = 
libunuran_la_LIBADD = $(SUBLIBS)
libunuran_la_LDFLAGS = --version-info $(UNURAN_LT_VERSION)

include_HEADERS = \
 unuran.h \
 unuran_config.h

noinst_HEADERS = \
 unuran.h.in \
 unur_source.h \
 unur_cookies.h \
 unur_struct.h \
 unur_typedefs.h

##unuran.h: $(top_srcdir)/*.h $(top_srcdir)/src/*/*.h 
unuran.h:
	$(top_srcdir)/scripts/merge_h.pl $(top_srcdir)/src/unuran.h.in > unuran.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 unuran.h
