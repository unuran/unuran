## Process this file with automake to produce Makefile.in
## $Id$

# AUTOMAKE_OPTIONS = 

AM_CFLAGS = -Wall

noinst_PROGRAMS = compute_machine_constants

compute_machine_constants_SOURCE = compute_machine_constants.c

EXTRA_DIST = \
 read_PDF.pl \
 merge_h.pl

noinst_HEADERS = \
 unur_fp_const_source.h

# create files with floating point constants
unur_fp_const_source.h: compute_machine_constants
	./compute_machine_constants > unur_fp_const_source.h	
	cp unur_fp_const_source.h $(top_srcdir)/src/utils/

# we do not put unur_fp_const_source.h into the distribution
dist-hook:
	rm -f $(distdir)/unur_fp_const_source.h

# clean backup files
CLEANFILES = \
 *~ \
 source_fp_const.h


# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
