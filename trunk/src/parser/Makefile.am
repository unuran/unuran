## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/methods

noinst_LIBRARIES = libparser.a

libparser_a_SOURCES = stringparser.c

include_HEADERS =

EXTRA_DIST = \
 make_stringparser.pl \
 stringparser.c.in

parser_srcdir = $(top_srcdir)/src/parser

stringparser.c: make_stringparser.pl stringparser.c.in
	$(parser_srcdir)/make_stringparser.pl < $(parser_srcdir)/stringparser.c.in > $(parser_srcdir)/stringparser.c


# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 stringparser.c \
 Makefile.in
