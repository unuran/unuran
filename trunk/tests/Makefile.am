## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/tests

check_PROGRAMS = \
 t_srou

TESTS = \
 t_srou

t_srou_SOURCES = t_srou.c support.c
t_srou_LDADD = $(top_builddir)/src/libunuran.a -lprng -lml -lm

noinst_HEADERS = t_unuran.h

# clean backup files too
clean-local:
	rm -f *~
	rm -f *.log




