## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions

lib_LIBRARIES = libunuran.a

libunuran_a_SOURCES = \
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
 dstd.c \
 ninv.c \
 rect.c \
 srou.c \
 stdr.c \
 tabl.c \
 tdr.c  \
 unif.c \
 utdr.c

pkginclude_HEADERS = \
 unuran.h \
 unuran_errno.h \
 unur_urng.h \
 unur_misc.h  \
 distr.h \
 arou.h \
 cstd.h \
 dau.h  \
 dis.h  \
 dstd.h \
 ninv.h  \
 rect.h \
 srou.h \
 stdr.h \
 tabl.h \
 tdr.h  \
 unif.h \
 utdr.h

noinst_HEADERS = \
 source_unuran.h \
 source_cookies.h \
 source_distr.h \
 source_math.h \
 source_methods.h \
 source_struct.h \
 distr_struct.h \
 x_umalloc_struct.h \
 arou_struct.h \
 cstd_struct.h \
 dau_struct.h  \
 dis_struct.h  \
 dstd_struct.h \
 ninv_struct.h  \
 rect_struct.h \
 srou_struct.h \
 stdr_struct.h \
 tabl_struct.h \
 tdr_struct.h  \
 unif_struct.h \
 utdr_struct.h

# clean backup files too
clean-local:
	rm -f *~
