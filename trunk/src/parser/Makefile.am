## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/methods \
 -I$(top_srcdir)/src/uniform

noinst_LIBRARIES = libparser.a

libparser_a_SOURCES = \
 parser.c  parser.h  parser_source.h \
 stringparser.c \
 functparser.c \
 functparser_source.h \
 functparser_struct.h \
 functparser_symbols.h \
 functparser_stringgen.ch \
 functparser_codegen_C.ch \
 functparser_codegen_FORTRAN.ch \
 functparser_codegen_JAVA.ch \
 functparser_debug.ch \
 functparser_deriv.ch \
 functparser_eval.ch \
 functparser_init.ch \
 functparser_parser.ch \
 functparser_scanner.ch

EXTRA_DIST = \
 make_stringparser.pl \
 functparser_doc.h \
 stringparser_doc.h \
 stringparser.c.in

parser_srcdir = $(top_srcdir)/src/parser

stringparser.c: make_stringparser.pl stringparser.c.in
	$(parser_srcdir)/make_stringparser.pl < $(parser_srcdir)/stringparser.c.in > $(parser_srcdir)/stringparser.c


# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 stringparser.c \
 stringparser_doc.h \
 Makefile.in
