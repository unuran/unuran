## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions

noinst_LIBRARIES = libmethods.a

libmethods_a_SOURCES = \
 distr.c \
 distr_cont.c \
 distr_corder.c \
 distr_cemp.c \
 distr_cvec.c \
 distr_demp.c \
 distr_discr.c \
 x_debug.c \
 x_error.c \
 x_gen.c \
 x_math.c \
 x_stream.c \
 x_umalloc.c \
 x_urng.c \
 arou.c \
 empk.c \
 cstd.c \
 dau.c  \
 dgt.c  \
 dstd.c \
 ninv.c \
 rect.c \
 srou.c \
 ssr.c \
 tabl.c \
 tdr.c  \
 unif.c \
 utdr.c

noinst_HEADERS = \
 in_unuran.h \
 source_unuran.h \
 source_cookies.h \
 source_debug.h \
 source_distr.h \
 source_fp.h \
 source_gen.h \
 source_math.h \
 source_misc.h \
 source_methods.h \
 source_struct.h \
 x_debug.h \
 x_errno.h \
 x_gen.h \
 x_math.h \
 x_urng.h \
 distr.h distr_struct.h \
 arou.h  arou_struct.h \
 empk.h  empk_struct.h \
 cstd.h  cstd_struct.h \
 dau.h   dau_struct.h \
 dgt.h   dgt_struct.h \
 dstd.h  dstd_struct.h \
 ninv.h  ninv_struct.h \
 rect.h  rect_struct.h \
 srou.h  srou_struct.h \
 ssr.h   ssr_struct.h \
 tabl.h  tabl_struct.h \
 tdr.h   tdr_struct.h tdr_newset.c tdr_init.c tdr_sample.c tdr_debug.c \
 unif.h  unif_struct.h \
 utdr.h  utdr_struct.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
