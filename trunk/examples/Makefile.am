## Process this file with automake to produce Makefile.in

AM_LDFLAGS = -L$(top_builddir)/src

INCLUDES = -I$(top_srcdir)/src

LDADD = $(top_builddir)/src/libunuran.la

EXAMPLES = \
 example0 example0_str \
 example1 example1_str \
 example2 example2_str \
 example3 example3_str \
 example4 example4_str \
 example_anti  example_anti_str \
 example_discr example_discr_str \
 example_cont  example_cont_str \
 example_emp   example_emp_str \
 example_vemp  example_FuncStr

noinst_PROGRAMS = $(EXAMPLES)

TESTS = $(EXAMPLES)

# clean backup files
CLEANFILES = \
 unuran.log \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
