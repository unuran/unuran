## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distr \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/stdgen \
 -I$(top_srcdir)/src/utils

noinst_LIBRARIES = libtests.a

libtests_a_SOURCES = \
 chi2test.c \
 counturn.c \
 printsample.c \
 scatter.c \
 tests.c \
 timing.c \
 transformed.c 

pkginclude_HEADERS = unur_tests.h

# clean backup files too
clean-local:
	rm -f *~
