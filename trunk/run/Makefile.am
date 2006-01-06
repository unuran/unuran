## Process this file with automake to produce Makefile.in

AM_LDFLAGS = -L$(top_builddir)/src

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/experimental \
 -I$(top_srcdir)/src/tests

LDADD = \
 $(top_builddir)/src/libunuran.la \
 $(top_builddir)/tests/testdistributions/libtestdistributions.la \
 $(top_builddir)/src/experimental/libunurexperimental.la

noinst_PROGRAMS = \
 run \
 something

## make_correlation_table

# clean backup files
CLEANFILES = \
 unuran.log \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

