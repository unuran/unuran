## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/parser

noinst_LIBRARIES = libtests.a

libtests_a_SOURCES = \
 chi2test.c \
 correlation.c\
 counturn.c \
 moments.c \
 printsample.c \
 quantiles.c \
 scatter.c \
 tests.c \
 timing.c

# scatter.c \ removed from source

include_HEADERS = \
 unuran_tests.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
