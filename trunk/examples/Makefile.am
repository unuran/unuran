## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall\
 -L$(top_builddir)/src

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

noinst_PROGRAMS = run example1 example2 example3 example4 test

run_SOURCES = run.c
example1_SOURCE = example1.c
example2_SOURCE = example2.c
example3_SOURCE = example3.c
example4_SOURCE = example4.c
test_SOURCE = test.c

run_LDADD = -lunuran
example1_LDADD = -lunuran
example2_LDADD = -lunuran
example3_LDADD = -lunuran
example4_LDADD = -lunuran
test_LDADD = -lunuran

# clean backup files
CLEANFILES = \
 *~ \
 *.scatterplot

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in



