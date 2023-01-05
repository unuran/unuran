## Process this file with automake to produce Makefile.in

AM_LDFLAGS = -L$(unuran_builddir)/src

AM_CPPFLAGS = \
 -I$(top_srcdir)/src \
 -I$(unuran_srcdir)/src \
 -I$(unuran_srcdir)/src/tests

LDADD = \
 $(top_builddir)/src/libunuran.la
 
noinst_PROGRAMS = run 

# clean log files and backup files
CLEANFILES = *~ *.log valgrind-*

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 \
 $(top_srcdir)/aclocal.m4 \
 $(top_srcdir)/config.h.in \
 $(top_srcdir)/configure

# Remark: the last three lines are an ugly hack
# to remove these files in the top source directory.
# It cannot be done in $(top_srcdir)/Makefile.am as
# this breaks the cleaning process.
