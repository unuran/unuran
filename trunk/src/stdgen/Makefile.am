## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/utils

noinst_LIBRARIES = libstdgen.a

libstdgen_a_SOURCES = gamma.c normal.c

pkginclude_HEADERS = unur_stdgen.h

# clean backup files too
clean-local:
	rm -f *~
