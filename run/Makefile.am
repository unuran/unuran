## Process this file with automake to produce Makefile.in

AM_LDFLAGS = -L$(top_builddir)/src

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests \
 -I$(top_srcdir)/tests/testdistributions

LDADD = \
 $(top_builddir)/tests/testdistributions/libtestdistributions.la \
 $(top_builddir)/src/experimental/libunurexperimental.la \
 $(top_builddir)/src/libunuran.la

noinst_PROGRAMS = \
 run \
 something

noinst_HEADERS = meanvarcor.c

## make_correlation_table

# clean log files and backup files
CLEANFILES = *~ *.log valgrind-*

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

