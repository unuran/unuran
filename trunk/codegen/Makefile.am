## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/distributions

lib_LIBRARIES = libunuracg.a

include_HEADERS = \
 unuran_acg.h

libunuracg_a_SOURCES = \
 PDFgen_source.h \
 PDFgen.c \
 codegen.c \
 codegen_source.h \
 tdr_codegen.c
 
noinst_PROGRAMS = \
 make_test_PDFgen \
 make_test_codegen


check_PROGRAMS = \
 test_PDFgen \
 test_codegen

TESTS = \
 test_PDFgen \
 test_codegen


noinst_HEADERS = 

EXTRA_DIST = \
 read_PDF.pl \
 make_PDFgen.pl \
 make_test_PDFgen.pl \
 make_test_PDFgen.c \
 make_test_codegen.pl \
 make_test_codegen.c


PDFgen.c: make_PDFgen.pl read_PDF.pl
	./make_PDFgen.pl

PDFgen_source.h: make_PDFgen.pl read_PDF.pl
	./make_PDFgen.pl


make_test_PDFgen.c: $(top_srcdir)/codegen/test.conf make_test_PDFgen.pl read_PDF.pl
	./make_test_PDFgen.pl $(top_srcdir)/codegen/test.conf

make_test_PDFgen_SOURCES = make_test_PDFgen.c
make_test_PDFgen_LDADD = libunuracg.a $(top_builddir)/src/libunuran.a

test_PDFgen.c: make_test_PDFgen
	./make_test_PDFgen

test_PDFgen_SOURCES = test_PDFgen.c
test_PDFgen_LDADD = $(top_builddir)/src/libunuran.a


make_test_codegen.c: $(top_srcdir)/codegen/test.conf make_test_codegen.pl read_PDF.pl
	./make_test_codegen.pl $(top_srcdir)/codegen/test.conf

make_test_codegen_SOURCES = make_test_codegen.c
make_test_codegen_LDADD = libunuracg.a  $(top_builddir)/src/libunuran.a

test_codegen.c: make_test_codegen
	./make_test_codegen

test_codegen_SOURCES = test_codegen.c
test_codegen_LDADD = $(top_builddir)/src/libunuran.a


# clean backup files
CLEANFILES = \
 PDFgen.c \
 PDFgen_source.h \
 make_test_PDFgen.c \
 test_PDFgen.c \
 make_test_codegen.c \
 test_codegen.c \
 *.log \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 PDFgen.c \
 PDFgen_source.h \
 make_test_PDFgen.c \
 test_PDFgen.c \
 make_test_codegen.c \
 test_codegen.c \
 Makefile.in

