## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/stdgen \
 -I$(top_srcdir)/src/utils

noinst_LIBRARIES = libmethods.a

libmethods_a_SOURCES = \
 distribution.c \
 misc.c \
 arou.c \
 cstd.c \
 dau.c  \
 dis.c  \
 rect.c \
 tabl.c \
 tdr.c  \
 unif.c \
 utdr.c

pkginclude_HEADERS = \
 unur_methods.h \
 unur_distribution.h \
 unur_misc.h  \
 unur_arou.h \
 unur_cstd.h \
 unur_dau.h  \
 unur_dis.h  \
 unur_rect.h \
 unur_tabl.h \
 unur_tdr.h  \
 unur_unif.h \
 unur_utdr.h

# clean backup files too
clean-local:
	rm -f *~
