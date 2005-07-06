## Process this file with automake to produce Makefile.in
## $Id$

INCLUDES = -I$(top_srcdir)/src

noinst_LTLIBRARIES = libtestdistributions.la

libtestdistributions_la_SOURCES = \
 c_wo_logpdf.c \
 c_sawtooth_discpdf.c \
 c_sawtooth_contpdf.c \
 vc_special_correlations_ar1.c \
 vc_special_correlations_constantrho.c

noinst_HEADERS = testdistributions.h 

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
