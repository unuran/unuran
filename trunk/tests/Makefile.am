## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

check_PROGRAMS = \
 t_tdr

TESTS = \
 t_tdr

EXTRA_DIST = \
 make_test_files.pl \
 README.conf \
 distr.conf \
 dau.conf \
 dgt.conf \
 srou.conf \
 ssr.conf \
 tdr.conf \
 unif.conf


t_tdr_SOURCES = t_tdr.c testroutines.c list_of_distributions.c
t_tdr_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

t_tdr.c: tdr.conf make_test_files.pl
	./make_test_files.pl tdr.conf

noinst_HEADERS = testunuran.h

# clean backup files too
clean-local:
	rm -f *~
	rm -f *.log

maintainer-clean-local:
	rm -f t_*
