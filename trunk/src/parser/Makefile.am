## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/methods

noinst_LIBRARIES = libparser.a

libparser_a_SOURCES = \
 parser.c \
 functparser.c \
 stringparser.c

include_HEADERS =

noinst_HEADERS = \
 parser.h \
 parser_source.h \
 functparser_source.h \
 functparser_struct.h \
 functparser_symbols.h \
 functparser_stringgen.c \
 functparser_codegen_C.c \
 functparser_codegen_FORTRAN.c \
 functparser_debug.c \
 functparser_deriv.c \
 functparser_eval.c \
 functparser_init.c \
 functparser_parser.c \
 functparser_scanner.c


EXTRA_DIST = \
 make_stringparser.pl \
 parser_doc.h \
 stringparser.c.in

parser_srcdir = $(top_srcdir)/src/parser

stringparser.c: make_stringparser.pl stringparser.c.in
	$(parser_srcdir)/make_stringparser.pl < $(parser_srcdir)/stringparser.c.in > $(parser_srcdir)/stringparser.c


# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 stringparser.c \
 parser_doc.h \
 Makefile.in
