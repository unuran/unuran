## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall\
 -L$(top_builddir)/src

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/uniform \
 -I$(top_srcdir)/src/tests

noinst_PROGRAMS = run test

run_SOURCES = run.c
example_SOURCE = test.c

run_LDADD = -lunuran -lprng -lml -lm
test_LDADD = -lunuran -lprng -lml -lm

# clean backup files too
clean-local:
	rm -f *~
	rm -f *.scatterplot




