## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/methods

lib_LIBRARIES = libunuran_distributions.a

libunuran_distributions_a_SOURCES = \
 beta.c        beta_gen.c \
 cauchy.c      cauchy_gen.c \
 chi.c         chi_gen.c \
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
 student.c     student_gen.c \
 uniform.c

pkginclude_HEADERS = \
 unuran_distributions.h

noinst_HEADERS = \
 source_distributions.h \
 source_specfunct.h \
 source_stddistr.h

# clean backup files too
clean-local:
	rm -f *~
