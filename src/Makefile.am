## Process this file with automake to produce Makefile.in
## $Id$

# Experimental code
if INCLUDE_EXPERIMENTAL
EXPERIMENTAL_DIRS = experimental
else
EXPERIMENTAL_DIRS = 
endif


SUBDIRS = \
 utils \
 methods \
 specfunct \
 distr \
 distributions \
 parser \
 tests \
 uniform \
 urng \
 $(EXPERIMENTAL_DIRS) \
 .

DIST_SUBDIRS = $(SUBDIRS)

SUBLIBS = \
 methods/libmethods.la \
 distr/libdistr.la \
 distributions/libdistributions.la \
 specfunct/libspecfunct.la \
 parser/libparser.la \
 tests/libtests.la \
 uniform/libuniform.la \
 urng/liburng.la \
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


# make unuran.h (process file unuran.h.in)
sinclude .dep-unuran_h
unuran.h: unuran.h.in 
	$(top_srcdir)/scripts/merge_h.pl $(top_srcdir)/src/unuran.h.in > unuran.h
	sed -e "s#^unuran\.h\:#stringparser\.c\:#" .dep-unuran_h | sed -e "s#\./#\.\./#g" > ./parser/.dep-stringparser_c

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 .dep-unuran_h \
 unuran.h
