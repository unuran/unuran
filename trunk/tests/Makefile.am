## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

check_PROGRAMS = \
 t_distr \
 t_stringparser \
 t_corder \
 t_timing \
 t_arou \
 t_cstd \
 t_dari \
 t_dau \
 t_dgt \
 t_dstd \
 t_empk \
 t_ninv \
 t_srou \
 t_ssr \
 t_tabl \
 t_tdr \
 t_unif \
 t_utdr

## commented out
## run_test_with_Mathematica
 
TESTS = \
 t_distr \
 t_stringparser \
 t_corder \
 t_arou \
 t_cstd \
 t_dari \
 t_dau \
 t_dgt \
 t_dstd \
 t_empk \
 t_ninv \
 t_srou \
 t_ssr \
 t_tabl \
 t_tdr \
 t_unif \
 t_utdr

## commented out
## run_test_with_Mathematica
## t_timing

noinst_HEADERS = \
 testunuran.h \
 test_with_Mathematica.h

EXTRA_DIST = \
 test_with_Mathematica.m \
 make_test_files.pl \
 README.conf \
 distr.conf \
 stringparser.conf \
 corder.conf \
 timing.conf \
 arou.conf \
 cstd.conf \
 dari.conf \
 dau.conf \
 dgt.conf \
 dstd.conf \
 empk.conf \
 ninv.conf \
 srou.conf \
 ssr.conf \
 tabl.conf \
 tdr.conf \
 unif.conf \
 utdr.conf

# compare results of CDF, PDF, and dPDF to those produced with Mathematica ...

data:
	math -initfile test_with_Mathematica.m

run_test_with_Mathematica.c: test_with_Mathematica.m
	math -initfile test_with_Mathematica.m

## run_test_with_Mathematica_SOURCES = run_test_with_Mathematica.c test_with_Mathematica.c 
## run_test_with_Mathematica_LDADD = $(top_builddir)/src/libunuran.a

# run tests for methods ...

t_distr_SOURCES = t_distr.c testroutines.c
t_distr_LDADD = $(top_builddir)/src/libunuran.a

t_stringparser_SOURCES = t_stringparser.c testroutines.c
t_stringparser_LDADD = $(top_builddir)/src/libunuran.a

t_corder_SOURCES = t_corder.c testroutines.c
t_corder_LDADD = $(top_builddir)/src/libunuran.a

t_timing_SOURCES = t_timing.c testroutines.c
t_timing_LDADD = $(top_builddir)/src/libunuran.a

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

t_dstd_SOURCES = t_dstd.c testroutines.c
t_dstd_LDADD = $(top_builddir)/src/libunuran.a

t_empk_SOURCES = t_empk.c testroutines.c
t_empk_LDADD = $(top_builddir)/src/libunuran.a

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


t_distr.c: distr.conf make_test_files.pl
	./make_test_files.pl distr.conf

t_stringparser.c: stringparser.conf make_test_files.pl
	./make_test_files.pl stringparser.conf

t_corder.c: corder.conf make_test_files.pl
	./make_test_files.pl corder.conf

t_timing.c: timing.conf make_test_files.pl
	./make_test_files.pl timing.conf

t_arou.c: arou.conf make_test_files.pl
	./make_test_files.pl arou.conf

t_cstd.c: cstd.conf make_test_files.pl
	./make_test_files.pl cstd.conf

t_dari.c: dari.conf make_test_files.pl
	./make_test_files.pl dari.conf

t_dau.c: dau.conf make_test_files.pl
	./make_test_files.pl dau.conf

t_dgt.c: dgt.conf make_test_files.pl
	./make_test_files.pl dgt.conf

t_dstd.c: dstd.conf make_test_files.pl
	./make_test_files.pl dstd.conf

t_empk.c: empk.conf make_test_files.pl
	./make_test_files.pl empk.conf

t_ninv.c: ninv.conf make_test_files.pl
	./make_test_files.pl ninv.conf

t_srou.c: srou.conf make_test_files.pl
	./make_test_files.pl srou.conf

t_ssr.c: ssr.conf make_test_files.pl
	./make_test_files.pl ssr.conf

t_tabl.c: tabl.conf make_test_files.pl
	./make_test_files.pl tabl.conf

t_tdr.c: tdr.conf make_test_files.pl
	./make_test_files.pl tdr.conf

t_utdr.c: utdr.conf make_test_files.pl
	./make_test_files.pl utdr.conf

t_unif.c: unif.conf make_test_files.pl
	./make_test_files.pl unif.conf

# clean backup files
CLEANFILES = *~ *.log

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 run_test_with_Mathematica.c \
 t_*

