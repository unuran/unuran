## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

check_PROGRAMS = \
 t_dau \
 t_dgt \
 t_ssr \
 t_srou \
 t_unif 

TESTS = \
 t_dau \
 t_dgt \
 t_ssr \
 t_srou \
 t_unif 

EXTRA_DIST = \
 make_test_files.pl \
 README.conf \
 dau.conf \
 dgt.conf \
 ssr.conf \
 srou.conf \
 unif.conf

t_dau_SOURCES = t_dau.c testroutines.c list_of_distributions.c
t_dau_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

t_dgt_SOURCES = t_dgt.c testroutines.c list_of_distributions.c
t_dgt_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

t_srou_SOURCES = t_srou.c testroutines.c list_of_distributions.c
t_srou_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

t_ssr_SOURCES = t_ssr.c testroutines.c list_of_distributions.c
t_ssr_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

t_unif_SOURCES = t_unif.c  testroutines.c list_of_distributions.c
t_unif_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm


t_dau.c: dau.conf make_test_files.pl
	./make_test_files.pl dau.conf

t_dgt.c: dgt.conf make_test_files.pl
	./make_test_files.pl dgt.conf

t_srou.c: srou.conf make_test_files.pl
	./make_test_files.pl srou.conf

t_ssr.c: ssr.conf make_test_files.pl
	./make_test_files.pl ssr.conf

t_unif.c: unif.conf make_test_files.pl
	./make_test_files.pl unif.conf


noinst_HEADERS = testunuran.h

# clean backup files too
clean-local:
	rm -f *~
	rm -f *.log

maintainer-clean-local:
	rm -f t_*
