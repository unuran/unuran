## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

AM_LDFLAGS = \
 -L$(top_builddir)/src

LDADD = \
 $(top_builddir)/src/libunuran.la

noinst_PROGRAMS = \
 make_timing_table

## make_correlation_table

# clean backup files
CLEANFILES = \
 unuran.log \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
