## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions

noinst_LIBRARIES = libmethods.a

libmethods_a_SOURCES = \
 distr.c \
 x_error.c \
 x_misc.c \
 x_stream.c \
 x_umalloc.c \
 x_urng.c \
 arou.c \
 cstd.c \
 dau.c  \
 dis.c  \
 rect.c \
 srou.c \
 stdr.c \
 tabl.c \
 tdr.c  \
 unif.c \
 utdr.c

pkginclude_HEADERS = \
 unuran.h \
 unur_umalloc.h \
 unur_urng.h \
 unur_distr.h \
 unur_misc.h  \
 unur_arou.h \
 unur_cstd.h \
 unur_dau.h  \
 unur_dis.h  \
 unur_rect.h \
 unur_srou.h \
 unur_stdr.h \
 unur_tabl.h \
 unur_tdr.h  \
 unur_unif.h \
 unur_utdr.h

noinst_HEADERS = \
 source_unuran.h \
 source_cookies.h \
 source_errno.h \
 source_math.h

# clean backup files too
clean-local:
	rm -f *~



