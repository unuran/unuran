## Process this file with automake to produce Makefile.in
## $Id$

noinst_PROGRAMS = compute_machine_constants

EXTRA_DIST = \
 read_PDF.pl \
 merge_h.pl

noinst_HEADERS = \
  $(top_srcdir)/src/utils/unur_fp_const_source.h

# create files with floating point constants
$(top_srcdir)/src/utils/unur_fp_const_source.h: compute_machine_constants
	test -f $(top_srcdir)/src/utils/unur_fp_const_source.h \
	|| ./compute_machine_constants > $(top_srcdir)/src/utils/unur_fp_const_source.h

# clean backup files
CLEANFILES = \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
