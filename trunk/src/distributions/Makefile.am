## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/methods

noinst_LIBRARIES = libdistr.a

libdistr_a_SOURCES = \
 beta.c        beta_gen.c \
 cauchy.c \
 chisquare.c \
 exponential.c exponential_gen.c \
 gamma.c       gamma_gen.c \
 laplace.c \
 lognormal.c \
 lomax.c \
 normal.c      normal_gen.c \
 pareto.c \
 powerexponential.c \
 rayleigh.c \
 student.c \
 uniform.c

pkginclude_HEADERS = unur_distribution.h

noinst_HEADERS = unur_distribution_lib.h

# clean backup files too
clean-local:
	rm -f *~
