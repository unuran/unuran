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
 Makefile.in

