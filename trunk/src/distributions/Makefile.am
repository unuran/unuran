## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/methods

lib_LIBRARIES = libunuran_distributions.a

libunuran_distributions_a_SOURCES = \
 c_beta.c        c_beta_gen.c \
 c_burr.c        c_burr_gen.c \
 c_cauchy.c      c_cauchy_gen.c \
 c_chi.c         c_chi_gen.c \
 c_chisquare.c \
 c_exponential.c c_exponential_gen.c \
 c_extremeI.c    c_extremeI_gen.c \
 c_extremeII.c   c_extremeII_gen.c \
 c_gamma.c       c_gamma_gen.c \
 c_gig.c         c_gig_gen.c \
 c_laplace.c     c_laplace_gen.c \
 c_logistic.c    c_logistic_gen.c \
 c_lognormal.c \
 c_lomax.c \
 c_normal.c      c_normal_gen.c \
 c_pareto.c \
 c_powerexponential.c c_powerexponential_gen.c \
 c_rayleigh.c \
 c_slash.c       c_slash_gen.c \
 c_student.c     c_student_gen.c \
 c_triangular.c  c_triangular_gen.c \
 c_uniform.c     c_uniform_gen.c \
 c_weibull.c     c_weibull_gen.c \
 d_zipf.c        d_zipf_gen.c

pkginclude_HEADERS = \
 unuran_distributions.h

noinst_HEADERS = \
 source_distributions.h \
 source_specfunct.h \
 source_stddistr.h

# clean backup files too
clean-local:
	rm -f *~
