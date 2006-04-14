## Process this file with automake to produce Makefile.in

SUBDIRS = testdistributions .

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

# Experimental code
if INCLUDE_EXPERIMENTAL
EXPERIMENTAL_LDADD = $(top_builddir)/src/experimental/libunurexperimental.la
# EXPERIMENTAL_PROGS = t_varou
EXPERIMENTAL_PROGS = t_mvtdr 
EXPERIMENTAL_XFAIL = 
else
EXPERIMENTAL_LDADD = 
EXPERIMENTAL_PROGS =
EXPERIMENTAL_XFAIL = 
endif


noinst_LTLIBRARIES = libtestroutines.la
libtestroutines_la_SOURCES = testroutines.c

LDADD = \
 libtestroutines.la \
 $(EXPERIMENTAL_LDADD) \
 $(top_builddir)/tests/testdistributions/libtestdistributions.la \
 $(top_builddir)/src/libunuran.la

# Special files for tests with Mathematica(TM)
if HAVE_MATHEMATICA
MATHEMATICA_PROGS = test_StdDistr  test_functionparser
else
MATHEMATICA_PROGS = 
endif


check_PROGRAMS = \
 t_unur_tests \
 t_urng \
 t_distr \
 t_distr_cemp \
 t_distr_cont \
 t_distr_corder \
 t_distr_cvec \
 t_distr_condi \
 t_distr_cvemp \
 t_distr_cxtrans \
 t_distr_discr \
 t_distr_matr \
 t_auto \
 t_arou \
 t_cstd \
 t_dari \
 t_dau \
 t_dgt \
 t_dsrou \
 t_dss \
 t_dstd \
 t_empk \
 t_empl \
 t_gibbs \
 t_hinv \
 t_hitro \
 t_hrb \
 t_hrd \
 t_hri \
 t_mcorr \
 t_ninv \
 t_norta \
 t_nrou \
 t_srou \
 t_ssr \
 t_tabl \
 t_tdr t_tdr_gw t_tdr_ps t_tdr_ia \
 t_tdrgw \
 t_unif \
 t_utdr \
 t_vempk \
 t_vmt \
 t_vnrou \
 t_stringparser \
 test_special_matrix \
 $(EXPERIMENTAL_PROGS) \
 $(MATHEMATICA_PROGS)


## Troubles with dmalloc:
##   t_dstd
##   t_hinv

TESTS = $(check_PROGRAMS)

# This is necessary until we fix the bug in method AROU
##XFAIL_TESTS = t_arou 
##XFAIL_TESTS = $(EXPERIMENTAL_XFAIL)

noinst_HEADERS = \
 testunuran.h

EXTRA_DIST = \
 test_StdDistr.m \
 test_functionparser.m \
 make_test_files.pl \
 README.conf \
 t_stringparser.conf \
 t_unur_tests.conf \
 t_urng.conf \
 t_distr.conf \
 t_distr_cemp.conf \
 t_distr_cont.conf \
 t_distr_corder.conf \
 t_distr_cvec.conf \
 t_distr_condi.conf \
 t_distr_cvemp.conf \
 t_distr_cxtrans.conf \
 t_distr_discr.conf \
 t_distr_matr.conf \
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
 t_gibbs \
 t_hinv.conf \
 t_hitro.conf \
 t_hrb.conf \
 t_hrd.conf \
 t_hri.conf \
 t_mcorr.conf \
 t_ninv.conf \
 t_norta.conf \
 t_nrou.conf \
 t_srou.conf \
 t_ssr.conf \
 t_tabl.conf \
 t_tdr.conf t_tdr_gw.conf t_tdr_ps.conf t_tdr_ia.conf \
 t_tdrgw.conf \
 t_unif.conf \
 t_utdr.conf \
 t_vempk.conf \
 t_vmt.conf \
 t_vnrou.conf


# Compare results of CDF, PDF, and dPDF of standard distributions
# and results of function parser to those produced with 
# Mathematica(TM) ...

if HAVE_MATHEMATICA

test_StdDistr_SOURCES = test_StdDistr.c
test_StdDistr_DEPENDENCIES = t_StdDistr.data

test_StdDistr.o: t_StdDistr.data

t_StdDistr.data: test_StdDistr.m
	math -initfile $(top_srcdir)/tests/test_StdDistr.m


test_functionparser_SOURCES = test_functionparser.c
test_functionparser_DEPENDENCIES = t_functionparser.data

test_functionparser.o: t_functionparser.data

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
valgrind: $(check_PROGRAMS)
	@list='$(check_PROGRAMS)'; for p in $$list; do \
		$(top_srcdir)/scripts/run_valgrind.sh ./$$p; \
	done

endif

# clean log files and backup files
CLEANFILES = *~ *.log valgrind-*

# clean generated files
DISTCLEANFILES = \
 rm -f t_*.data

# clean generated files
MAINTAINERCLEANFILES = \
 rm -f t_*.c \
 Makefile.in
