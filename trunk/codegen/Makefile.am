## Process this file with automake to produce Makefile.in
## $Id$

SUBDIRS = . tests

AM_CFLAGS = -Wall

# acg does not build with -ansi -pedantic
# (there should be a better fix)
CFLAGS = -O2

INCLUDES = -I$(top_srcdir)/src

bin_PROGRAMS = acg
acg_SOURCES = acg.c
acg_LDADD = libunuracg.la $(top_builddir)/src/libunuran.la


lib_LTLIBRARIES = libunuracg.la

include_HEADERS = \
 unuran_acg.h

libunuracg_la_SOURCES = \
 PDFgen.c \
 headergen.c \
 codegen_C.c \
 codegen_FORTRAN.c \
 codegen_JAVA.c \
 codegen_UNURAN.c \
 codegen_source.h \
 tdr_codegen_C.c \
 tdr_codegen_FORTRAN.c \
 tdr_codegen_JAVA.c \
 tdr_codegen_UNURAN.c
 

noinst_HEADERS = 

EXTRA_DIST = \
 make_PDFgen.pl    PDFgen.c.in \
 make_headergen.pl headergen.c.in \
 make_anuran.pl    anuran.pl.in 


PDFgen.c: make_PDFgen.pl PDFgen.c.in $(top_srcdir)/scripts/read_PDF.pl
	./make_PDFgen.pl < PDFgen.c.in > PDFgen.c

headergen.c: make_headergen.pl headergen.c.in $(top_srcdir)/scripts/read_PDF.pl
	./make_headergen.pl < headergen.c.in > headergen.c

anuran.pl: make_anuran.pl anuran.pl.in $(top_srcdir)/scripts/read_PDF.pl
	./make_anuran.pl < anuran.pl.in > anuran.pl


# clean backup files
CLEANFILES = \
 *.log \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 PDFgen.c \
 headergen.c \
 anuran.pl
