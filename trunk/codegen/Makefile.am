## Process this file with automake to produce Makefile.in
## $Id$

SUBDIRS = . tests

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/distributions


bin_PROGRAMS = acg
acg_SOURCES = acg.c
acg_LDADD = libunuracg.a $(top_builddir)/src/libunuran.a



lib_LIBRARIES = libunuracg.a

include_HEADERS = \
 unuran_acg.h

libunuracg_a_SOURCES = \
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
 read_PDF.pl \
 make_PDFgen.pl \
 make_headergen.pl \
 make_acg.pl \
 make_anuran.pl


PDFgen.c: make_PDFgen.pl ../scripts/read_PDF.pl
	./make_PDFgen.pl

headergen.c: make_headergen.pl ../scripts/read_PDF.pl
	./make_headergen.pl

acg.c: make_acg.pl  ../scripts/read_PDF.pl
	./make_acg.pl > acg.c

anuran.pl: make_anuran.pl ../scripts/read_PDF.pl
	./make_anuran.pl > anuran.pl


# clean backup files
CLEANFILES = \
 *.log \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 PDFgen.c \
 headergen.c \
 acg.c \
 anuran.pl

