## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall\
 -L$(top_builddir)/src

INCLUDES = -I$(top_srcdir)/src

LDADD = -lunuran

noinst_PROGRAMS = make_timing_table

make_timing_table_SOURCE = make_timing_table

# clean backup files
CLEANFILES = \
 unuran.log \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
