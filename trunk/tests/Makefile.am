## Process this file with automake to produce Makefile.in

SUBDIRS = .

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

noinst_LTLIBRARIES = libtestroutines.la
libtestroutines_la_SOURCES = testroutines.c

LDADD = libtestroutines.la $(top_builddir)/src/libunuran.la

# Special files for tests with Mathematica(TM)
if HAVE_MATHEMATICA
MATHEMATICA_PROGS = test_StdDistr  test_functionparser
MATHEMATICA_TESTS = test_StdDistr  test_functionparser
else
MATHEMATICA_PROGS = 
MATHEMATICA_TESTS = 
endif


check_PROGRAMS = \
 t_distr \
 t_distr_cemp \
 t_distr_cont \
 t_distr_corder \
 t_distr_cvec \
 t_distr_cvemp \
 t_distr_discr \
 t_distr_matr \
 t_auto \
 t_arou \
 t_corder \
 t_cstd \
 t_dari \
 t_dau \
 t_dgt \
 t_dsrou \
 t_dss \
 t_dstd \
 t_empk \
 t_empl \
 t_hinv \
 t_hrb \
 t_hrd \
 t_hri \
 t_mcorr \
 t_ninv \
 t_nrou \
 t_srou \
 t_ssr \
 t_tabl \
 t_tdr \
 t_unif \
 t_utdr \
 t_vempk \
 t_vmt \
 t_vnrou \
 t_stringparser \
 test_special_matrix \
 $(MATHEMATICA_PROGS)

## Troubles with dmalloc:
##   t_dstd
##   t_hinv

## Disabled:
##   t_timing_dsrou
##   t_timing_ninv

TESTS = $(check_PROGRAMS)

# This is necessary until we fix the bug in method AROU
XFAIL_TESTS = t_arou 

noinst_HEADERS = \
 testunuran.h

EXTRA_DIST = \
 test_StdDistr.m \
 test_functionparser.m \
 make_test_files.pl \
 README.conf \
 t_distr.conf \
 t_distr_cemp.conf \
 t_distr_cont.conf \
 t_distr_corder.conf \
 t_distr_cvec.conf \
 t_distr_cvemp.conf \
 t_distr_discr.conf \
 t_distr_matr.conf \
 t_stringparser.conf \
 t_corder.conf \
 t_auto.conf \
 t_arou.conf \
 t_cstd.conf \
 t_dari.conf \
 t_dau.conf \
 t_dgt.conf \
 t_dsrou.conf \
 t_dss.conf \
 t_dstd.conf \
 t_empk.conf \
 t_empl.conf \
 t_hinv.conf \
 t_hrb.conf \
 t_hrd.conf \
 t_hri.conf \
 t_mcorr.conf \
 t_ninv.conf \
 t_srou.conf \
 t_ssr.conf \
 t_tabl.conf \
 t_tdr.conf \
 t_unif.conf \
 t_utdr.conf \
 t_vempk.conf \
 t_vmt.conf \
 t_vnrou.conf \
 t_timing_dsrou.conf \
 t_timing_ninv.conf


# Compare results of CDF, PDF, and dPDF of standard distributions
# and results of function parser to those produced with 
# Mathematica(TM) ...

if HAVE_MATHEMATICA

test_StdDistr_SOURCES = test_StdDistr.c
test_StdDistr_DEPENDENCIES = t_StdDistr.data

t_StdDistr.data: test_StdDistr.m
	math -initfile $(top_srcdir)/tests/test_StdDistr.m


test_functionparser_SOURCES = test_functionparser.c
test_functionparser_DEPENDENCIES = t_functionparser.data

t_functionparser.data: test_functionparser.m
	math -initfile $(top_srcdir)/tests/test_functionparser.m

endif


# Generate C sources for tests
SUFFIXES = .c
.conf.c: make_test_files.pl
	./make_test_files.pl $< > $@

t_*.c: make_test_files.pl

# run tests with valgrind
if HAVE_VALGRIND

VALGRINDOPTIONS = -v \
 --leak-check=yes \
 --leak-resolution=high \
 --num-callers=40 \
 --show-reachable=yes 
valgrind: $(check_PROGRAMS)
	@list='$(check_PROGRAMS)'; for p in $$list; do \
		valgrind --logfile=valgrind-$$p $(VALGRINDOPTIONS) ./$$p; \
	done

endif

# clean log files and backup files
CLEANFILES = *~ *.log valgrind-*pid*

# clean generated files
DISTCLEANFILES = \
 rm -f t_*.data

# clean generated files
MAINTAINERCLEANFILES = \
 rm -f t_*.c \
 Makefile.in
