## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall\
 -L$(top_builddir)/src

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

noinst_PROGRAMS = example1 example2 example3 example4 \
 example_anti example_discr example_cont example_emp  example_vemp

example1_SOURCE      = example1.c
example2_SOURCE      = example2.c
example3_SOURCE      = example3.c
example4_SOURCE      = example4.c
example_anti_SOURCE  = example_anti.c
example_discr_SOURCE = example_discr.c
example_discr_SOURCE = example_cont.c
example_emp_SOURCE   = example_emp.c
example_vemp_SOURCE  = example_vemp.c

example1_LDADD      = -lunuran
example2_LDADD      = -lunuran
example3_LDADD      = -lunuran
example4_LDADD      = -lunuran
example_anti_LDADD  = -lunuran
example_discr_LDADD = -lunuran
example_cont_LDADD  = -lunuran
example_emp_LDADD   = -lunuran
example_vemp_LDADD  = -lunuran

# clean backup files
CLEANFILES = \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

