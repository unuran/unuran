## Process this file with automake to produce Makefile.in
## $Id$

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
 make_headergen.pl \
 make_acg.pl \
 make_test_PDFgen.pl \
 make_test_PDFgen.c \
 make_test_codegen.pl \
 make_test_codegen.c


PDFgen.c: make_PDFgen.pl read_PDF.pl
	./make_PDFgen.pl

headergen.c: make_headergen.pl read_PDF.pl
	./make_headergen.pl

acg.c: make_acg.pl read_PDF.pl
	./make_acg.pl > acg.c

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
 headergen.c \
 acg.c \
 make_test_PDFgen.c \
 test_PDFgen.c \
 make_test_codegen.c \
 test_codegen.c \
 *.log \
 *~

# clean generated files
MAINTAINERCLEANFILES = \
 PDFgen.c \
 headergen.c \
 make_test_PDFgen.c \
 test_PDFgen.c \
 make_test_codegen.c \
 test_codegen.c \
 Makefile.in

