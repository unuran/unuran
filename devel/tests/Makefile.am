## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = \
 $(UNURAN_SUPPORTS_GSL) \
 $(UNURAN_SUPPORTS_PRNG) \
 $(UNURAN_SUPPORTS_RNGSTREAM) \
 -I$(unuran_srcdir) \
 -I$(unuran_srcdir)/src \
 -I$(unuran_srcdir)/src/tests \
 -I$(unuran_srcdir)/tests

LDADD = \
 $(unuran_builddir)/tests/libtestroutines.la \
 $(unuran_builddir)/tests/testdistributions/libtestdistributions.la \
 $(top_builddir)/src/libunurexperimental.la \
 $(unuran_builddir)/src/libunuran.la

check_PROGRAMS = 

TESTS = $(check_PROGRAMS)

# These tests fail until we fix some bugs
#XFAIL_TESTS = 

# Generate C sources for tests
SUFFIXES = .c
.conf.c: $(unuran_srcdir)/tests/make_test_files.pl
	$(unuran_srcdir)/tests/make_test_files.pl $(unuran_srcdir) $< > $@

t_*.c: $(unuran_srcdir)/tests/make_test_files.pl

# run tests with valgrind
if HAVE_VALGRIND
valgrind: $(check_PROGRAMS)
	@list='$(check_PROGRAMS)'; for p in $$list; do \
		$(unuran_srcdir)/scripts/run_valgrind.sh ./$$p; \
	done
endif

# clean log files and backup files
CLEANFILES = \
 *.log \
 valgrind-* \
 *~

# clean generated files
DISTCLEANFILES = \
 rm -f t_*.data

# clean generated files
MAINTAINERCLEANFILES = \
 rm -f t_*.c \
 Makefile.in
