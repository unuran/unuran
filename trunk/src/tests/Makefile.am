## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/src

noinst_LIBRARIES = libtests.a

libtests_a_SOURCES = \
 chi2test.c \
 correlation.c\
 counturn.c \
 moments.c \
 printsample.c \
 quantiles.c \
 tests.c \
 timing.c

include_HEADERS = \
 unuran_tests.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
