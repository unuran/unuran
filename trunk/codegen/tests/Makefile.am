## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/distributions


EXTRA_DIST = \
 test_acg.pl \
 test_acg_conf.pl \
 test_urand.pl


TESTS = \
 test_urand.pl \
 test_acg.pl


# clean backup files
CLEANFILES = \
 run_test_* \
 *.java *.class \
 unuran.log \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

