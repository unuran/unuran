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
 example_vemp  \
 example_FuncStr \
 example_mcorr \
 example_vmt

check_PROGRAMS = $(EXAMPLES)

TESTS = $(EXAMPLES)

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
CLEANFILES = \
 unuran.log \
 valgrind-*pid* \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
