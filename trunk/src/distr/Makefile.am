## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/utils

noinst_LIBRARIES = libdistr.a

libdistr_a_SOURCES = \
 beta.c \
 cauchy.c \
 chi2.c \
 exponential.c \
 gamma.c \
 laplace.c \
 lognormal.c \
 lomax.c \
 normal.c \
 pareto.c \
 powerexpon.c \
 rayleigh.c \
 student.c \
 uniform.c

pkginclude_HEADERS = unur_distr.h

# clean backup files too
clean-local:
	rm -f *~
