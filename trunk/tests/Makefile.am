## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

check_PROGRAMS = \
 t_distr \
 t_corder \
 t_arou \
 t_cstd \
 t_dau \
 t_dgt \
 t_srou \
 t_ssr \
 t_tabl \
 t_tdr \
 t_unif \
 t_utdr

TESTS = \
 t_distr \
 t_corder \
 t_arou \
 t_cstd \
 t_dau \
 t_dgt \
 t_srou \
 t_ssr \
 t_tabl \
 t_tdr \
 t_unif \
 t_utdr

noinst_HEADERS = testunuran.h

EXTRA_DIST = \
 make_test_files.pl \
 README.conf \
 distr.conf \
 corder.conf \
 arou.conf \
 cstd.conf \
 dau.conf \
 dgt.conf \
 srou.conf \
 ssr.conf \
 tabl.conf \
 tdr.conf \
 unif.conf \
 utdr.conf


t_distr_SOURCES = t_distr.c testroutines.c
t_distr_LDADD = $(top_builddir)/src/libunuran.a

t_corder_SOURCES = t_corder.c testroutines.c
t_corder_LDADD = $(top_builddir)/src/libunuran.a

t_arou_SOURCES = t_arou.c testroutines.c
t_arou_LDADD = $(top_builddir)/src/libunuran.a

t_cstd_SOURCES = t_cstd.c testroutines.c
t_cstd_LDADD = $(top_builddir)/src/libunuran.a

t_dau_SOURCES = t_dau.c testroutines.c
t_dau_LDADD = $(top_builddir)/src/libunuran.a

t_dgt_SOURCES = t_dgt.c testroutines.c
t_dgt_LDADD = $(top_builddir)/src/libunuran.a

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

t_corder.c: corder.conf make_test_files.pl
	./make_test_files.pl corder.conf

t_arou.c: arou.conf make_test_files.pl
	./make_test_files.pl arou.conf

t_cstd.c: cstd.conf make_test_files.pl
	./make_test_files.pl cstd.conf

t_dau.c: dau.conf make_test_files.pl
	./make_test_files.pl dau.conf

t_dgt.c: dgt.conf make_test_files.pl
	./make_test_files.pl dgt.conf

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
 t_*

