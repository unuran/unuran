## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/stdgen \
 -I$(top_srcdir)/src/uniform \
 -I$(top_srcdir)/src/utils \
 -I$(top_srcdir)/src/tests

noinst_PROGRAMS = run
run_SOURCES = run.c
run_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml 

# clean backup files too
clean-local:
	rm -f *~
	rm -f *.scatterplot
