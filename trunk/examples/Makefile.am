## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall\
 -L$(top_builddir)/src/methods\
 -L$(top_builddir)/src/distributions\
 -L$(top_builddir)/src/tests\
 -L$(top_builddir)/src/uniform

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/uniform \
 -I$(top_srcdir)/src/tests

noinst_PROGRAMS = run
run_SOURCES = run.c

run_LDADD = -lunuran_tests -lunuran_distributions -lunuran -luniform -lprng -lml 

# clean backup files too
clean-local:
	rm -f *~
	rm -f *.scatterplot




