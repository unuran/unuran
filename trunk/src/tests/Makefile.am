## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/methods

noinst_LIBRARIES = libtests.a

libtests_a_SOURCES = \
 chi2test.c \
 counturn.c \
 printsample.c \
 scatter.c \
 tests.c \
 timing.c

pkginclude_HEADERS = unur_tests.h

# clean backup files too
clean-local:
	rm -f *~
