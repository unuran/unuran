## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

check_PROGRAMS = \
 t_dau \
 t_dgt \
 t_srou \
 t_stdr

TESTS = \
 t_dau \
 t_dgt \
 t_srou \
 t_stdr

t_dau_SOURCES = t_dau.c support.c list_of_distributions.c
t_dau_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

t_dgt_SOURCES = t_dgt.c support.c list_of_distributions.c
t_dgt_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

t_srou_SOURCES = t_srou.c support.c list_of_distributions.c
t_srou_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

t_stdr_SOURCES = t_stdr.c support.c list_of_distributions.c
t_stdr_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

noinst_HEADERS = t_unuran.h

# clean backup files too
clean-local:
	rm -f *~
	rm -f *.log




