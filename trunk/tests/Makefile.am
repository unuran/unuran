## Process this file with automake to produce Makefile.in

SUBDIRS = .

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

# Special files for tests with Mathematica
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
 t_stringparser \
 t_corder \
 t_auto \
 t_arou \
 t_cstd \
 t_dari \
 t_dau \
 t_dgt \
 t_dsrou \
 t_dstd \
 t_empk \
 t_hinv \
 t_ninv \
 t_srou \
 t_ssr \
 t_tabl \
 t_tdr \
 t_unif \
 t_utdr \
 t_vmt \
 t_timing_dsrou \
 t_timing_ninv \
 $(MATHEMATICA_PROGS)

TESTS = \
 t_distr \
 t_distr_cemp \
 t_distr_cont \
 t_distr_corder \
 t_distr_cvec \
 t_distr_cvemp \
 t_distr_discr \
 t_stringparser \
 t_corder \
 t_auto \
 t_arou \
 t_cstd \
 t_dari \
 t_dau \
 t_dgt \
 t_dsrou \
 t_dstd \
 t_empk \
 t_hinv \
 t_ninv \
 t_srou \
 t_ssr \
 t_tabl \
 t_tdr \
 t_unif \
 t_utdr \
 t_vmt \
 $(MATHEMATICA_TESTS)


noinst_HEADERS = \
 testunuran.h

EXTRA_DIST = \
 test_StdDistr.m \
 test_functionparser.m \
 make_test_files.pl \
 README.conf \
 distr.conf \
 distr_cemp.conf \
 distr_cont.conf \
 distr_corder.conf \
 distr_cvec.conf \
 distr_cvemp.conf \
 distr_discr.conf \
 stringparser.conf \
 corder.conf \
 auto.conf \
 arou.conf \
 cstd.conf \
 dari.conf \
 dau.conf \
 dgt.conf \
 dsrou.conf \
 dstd.conf \
 empk.conf \
 hinv.conf \
 ninv.conf \
 srou.conf \
 ssr.conf \
 tabl.conf \
 tdr.conf \
 unif.conf \
 utdr.conf \
 vmt.conf \
 timing_dsrou.conf \
 timing_ninv.conf


# compare results of CDF, PDF, and dPDF of standard distributions
# and results of function parser
# to those produced with Mathematica ...

if HAVE_MATHEMATICA

test_StdDistr_SOURCES = test_StdDistr.c
test_StdDistr_LDADD = $(top_builddir)/src/libunuran.a
test_StdDistr_DEPENDENCIES = $(top_builddir)/src/libunuran.a t_StdDistr.data

t_StdDistr.data: test_StdDistr.m
	math -initfile $(top_srcdir)/tests/test_StdDistr.m


test_functionparser_SOURCES = test_functionparser.c
test_functionparser_LDADD = $(top_builddir)/src/libunuran.a
test_functionparser_DEPENDENCIES = $(top_builddir)/src/libunuran.a t_functionparser.data

t_functionparser.data: test_functionparser.m
	math -initfile $(top_srcdir)/tests/test_functionparser.m

endif


# run tests for methods ...

t_distr_SOURCES = t_distr.c testroutines.c
t_distr_LDADD = $(top_builddir)/src/libunuran.a

t_distr_cemp_SOURCES = t_distr_cemp.c testroutines.c
t_distr_cemp_LDADD = $(top_builddir)/src/libunuran.a

t_distr_cont_SOURCES = t_distr_cont.c testroutines.c
t_distr_cont_LDADD = $(top_builddir)/src/libunuran.a

t_distr_corder_SOURCES = t_distr_corder.c testroutines.c
t_distr_corder_LDADD = $(top_builddir)/src/libunuran.a

t_distr_cvec_SOURCES = t_distr_cvec.c testroutines.c
t_distr_cvec_LDADD = $(top_builddir)/src/libunuran.a

t_distr_cvemp_SOURCES = t_distr_cvemp.c testroutines.c
t_distr_cvemp_LDADD = $(top_builddir)/src/libunuran.a

t_distr_discr_SOURCES = t_distr_discr.c testroutines.c
t_distr_discr_LDADD = $(top_builddir)/src/libunuran.a

t_stringparser_SOURCES = t_stringparser.c testroutines.c
t_stringparser_LDADD = $(top_builddir)/src/libunuran.a

t_corder_SOURCES = t_corder.c testroutines.c
t_corder_LDADD = $(top_builddir)/src/libunuran.a

t_auto_SOURCES = t_auto.c testroutines.c
t_auto_LDADD = $(top_builddir)/src/libunuran.a

t_arou_SOURCES = t_arou.c testroutines.c
t_arou_LDADD = $(top_builddir)/src/libunuran.a

t_cstd_SOURCES = t_cstd.c testroutines.c
t_cstd_LDADD = $(top_builddir)/src/libunuran.a

t_dari_SOURCES = t_dari.c testroutines.c
t_dari_LDADD = $(top_builddir)/src/libunuran.a

t_dau_SOURCES = t_dau.c testroutines.c
t_dau_LDADD = $(top_builddir)/src/libunuran.a

t_dgt_SOURCES = t_dgt.c testroutines.c
t_dgt_LDADD = $(top_builddir)/src/libunuran.a

t_dsrou_SOURCES = t_dsrou.c testroutines.c
t_dsrou_LDADD = $(top_builddir)/src/libunuran.a

t_dstd_SOURCES = t_dstd.c testroutines.c
t_dstd_LDADD = $(top_builddir)/src/libunuran.a

t_empk_SOURCES = t_empk.c testroutines.c
t_empk_LDADD = $(top_builddir)/src/libunuran.a

t_hinv_SOURCES = t_hinv.c testroutines.c
t_hinv_LDADD = $(top_builddir)/src/libunuran.a

t_ninv_SOURCES = t_ninv.c testroutines.c
t_ninv_LDADD = $(top_builddir)/src/libunuran.a

t_srou_SOURCES = t_srou.c testroutines.c
t_srou_LDADD = $(top_builddir)/src/libunuran.a

t_ssr_SOURCES = t_ssr.c testroutines.c
t_ssr_LDADD = $(top_builddir)/src/libunuran.a

t_tabl_SOURCES = t_tabl.c testroutines.c
t_tabl_LDADD = $(top_builddir)/src/libunuran.a

t_tdr_SOURCES = t_tdr.c testroutines.c
t_tdr_LDADD = $(top_builddir)/src/libunuran.a

t_unif_SOURCES = t_unif.c  testroutines.c
t_unif_LDADD = $(top_builddir)/src/libunuran.a

t_utdr_SOURCES = t_utdr.c testroutines.c
t_utdr_LDADD = $(top_builddir)/src/libunuran.a

t_vmt_SOURCES = t_vmt.c testroutines.c
t_vmt_LDADD = $(top_builddir)/src/libunuran.a

t_timing_dsrou_SOURCES = t_timing_dsrou.c testroutines.c
t_timing_dsrou_LDADD = $(top_builddir)/src/libunuran.a

t_timing_ninv_SOURCES = t_timing_ninv.c testroutines.c
t_timing_ninv_LDADD = $(top_builddir)/src/libunuran.a


t_distr.c: distr.conf make_test_files.pl
	./make_test_files.pl distr.conf          > t_distr.c

t_distr_cemp.c: distr_cemp.conf make_test_files.pl
	./make_test_files.pl distr_cemp.conf     > t_distr_cemp.c

t_distr_cont.c: distr_cont.conf make_test_files.pl
	./make_test_files.pl distr_cont.conf     > t_distr_cont.c

t_distr_corder.c: distr_corder.conf make_test_files.pl
	./make_test_files.pl distr_corder.conf   > t_distr_corder.c

t_distr_cvec.c: distr_cvec.conf make_test_files.pl
	./make_test_files.pl distr_cvec.conf     > t_distr_cvec.c

t_distr_cvemp.c: distr_cvemp.conf make_test_files.pl
	./make_test_files.pl distr_cvemp.conf    > t_distr_cvemp.c

t_distr_discr.c: distr_discr.conf make_test_files.pl
	./make_test_files.pl distr_discr.conf    > t_distr_discr.c

t_stringparser.c: stringparser.conf make_test_files.pl
	./make_test_files.pl stringparser.conf   > t_stringparser.c

t_corder.c: corder.conf make_test_files.pl
	./make_test_files.pl corder.conf         > t_corder.c

t_auto.c: auto.conf make_test_files.pl
	./make_test_files.pl auto.conf           > t_auto.c

t_arou.c: arou.conf make_test_files.pl
	./make_test_files.pl arou.conf           > t_arou.c

t_cstd.c: cstd.conf make_test_files.pl
	./make_test_files.pl cstd.conf           > t_cstd.c

t_dari.c: dari.conf make_test_files.pl
	./make_test_files.pl dari.conf           > t_dari.c

t_dau.c: dau.conf make_test_files.pl
	./make_test_files.pl dau.conf            > t_dau.c

t_dgt.c: dgt.conf make_test_files.pl
	./make_test_files.pl dgt.conf            > t_dgt.c

t_dsrou.c: dsrou.conf make_test_files.pl
	./make_test_files.pl dsrou.conf          > t_dsrou.c

t_dstd.c: dstd.conf make_test_files.pl
	./make_test_files.pl dstd.conf           > t_dstd.c

t_empk.c: empk.conf make_test_files.pl
	./make_test_files.pl empk.conf           > t_empk.c

t_hinv.c: hinv.conf make_test_files.pl
	./make_test_files.pl hinv.conf           > t_hinv.c

t_ninv.c: ninv.conf make_test_files.pl
	./make_test_files.pl ninv.conf           > t_ninv.c

t_srou.c: srou.conf make_test_files.pl
	./make_test_files.pl srou.conf           > t_srou.c

t_ssr.c: ssr.conf make_test_files.pl
	./make_test_files.pl ssr.conf            > t_ssr.c

t_tabl.c: tabl.conf make_test_files.pl
	./make_test_files.pl tabl.conf           > t_tabl.c

t_tdr.c: tdr.conf make_test_files.pl
	./make_test_files.pl tdr.conf            > t_tdr.c

t_unif.c: unif.conf make_test_files.pl
	./make_test_files.pl unif.conf           > t_unif.c

t_utdr.c: utdr.conf make_test_files.pl
	./make_test_files.pl utdr.conf           > t_utdr.c

t_vmt.c: vmt.conf make_test_files.pl
	./make_test_files.pl vmt.conf            > t_vmt.c

t_timing_dsrou.c: timing_dsrou.conf make_test_files.pl
	./make_test_files.pl timing_dsrou.conf   > t_timing_dsrou.c

t_timing_ninv.c: timing_ninv.conf make_test_files.pl
	./make_test_files.pl timing_ninv.conf    > t_timing_ninv.c


# clean backup files
CLEANFILES = *~ *.log

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 t_*

