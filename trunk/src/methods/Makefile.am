## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/utils

noinst_LIBRARIES = libmethods.a

libmethods_a_SOURCES = \
 distr.c \
 misc.c \
 arou.c \
 srou.c \
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
 unur_distr.h \
 unur_misc.h  \
 unur_arou.h \
 unur_srou.h \
 unur_cstd.h \
 unur_dau.h  \
 unur_dis.h  \
 unur_rect.h \
 unur_tabl.h \
 unur_tdr.h  \
 unur_unif.h \
 unur_utdr.h

noinst_HEADERS =  unur_methods_lib.h 

# clean backup files too
clean-local:
	rm -f *~
