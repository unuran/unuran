## Process this file with automake to produce Makefile.in
## $Id$

# AUTOMAKE_OPTIONS = 

AM_CFLAGS = -Wall

noinst_PROGRAMS = compute_machine_constants

compute_machine_constants_SOURCE = compute_machine_constants.c

EXTRA_DIST = \
 merge_h.pl

noinst_HEADERS = \
 source_fp_const.h

# create files with floating point constants
source_fp_const.h: compute_machine_constants
	./compute_machine_constants > source_fp_const.h	
	cp source_fp_const.h $(top_srcdir)/src/methods/

# we do not put source_fp_const.h into the distribution
dist-hook:
	rm -f $(distdir)/source_fp_const.h

# clean backup files
CLEANFILES = \
 *~ \
 source_fp_const.h


# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
