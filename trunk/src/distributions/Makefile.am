## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/parser \
 -I$(top_srcdir)/src/uniform

noinst_LIBRARIES = libdistributions.a

libdistributions_a_SOURCES = \
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
 c_lomax.c       c_lomax_gen.c \
 c_normal.c      c_normal_gen.c \
 c_pareto.c      c_pareto_gen.c \
 c_powerexponential.c c_powerexponential_gen.c \
 c_rayleigh.c \
 c_slash.c       c_slash_gen.c \
 c_student.c     c_student_gen.c \
 c_triangular.c  c_triangular_gen.c \
 c_uniform.c     c_uniform_gen.c \
 c_weibull.c     c_weibull_gen.c \
 d_binomial.c    d_binomial_gen.c \
 d_geometric.c   d_geometric_gen.c \
 d_hypergeometric.c d_hypergeometric_gen.c \
 d_logarithmic.c d_logarithmic_gen.c \
 d_negativebinomial.c \
 d_poisson.c     d_poisson_gen.c \
 d_zipf.c        d_zipf_gen.c \
 vc_multinormal.c 

noinst_HEADERS = \
 unuran_distributions.h \
 unur_stddistr.h \
 source_distributions.h \
 source_specfunct.h

EXTRA_DIST = \
 README

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
