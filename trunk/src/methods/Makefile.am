## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/distributions \
 -I$(top_srcdir)/src/parser \
 -I$(top_srcdir)/src/uniform

noinst_LIBRARIES = libmethods.a

libmethods_a_SOURCES = \
 in_unuran.h \
 source_unuran.h source_cookies.h source_fp.h source_methods.h source_struct.h \
\
 distr.c        distr.h  distr_struct.h distr_source.h \
 distr_cont.c   distr_cont.h \
 distr_corder.c distr_corder.h \
 distr_cemp.c   distr_cemp.h \
 distr_cvec.c   distr_cvec.h \
 distr_cvemp.c  distr_cvemp.h\
 distr_discr.c  distr_discr.h\
\
 x_debug.c  x_debug.h x_debug_source.h \
 x_error.c  x_errno.h \
 x_gen.c    x_gen.h   x_gen_source.h \
 x_math.c   x_math.h  x_math_source.h \
 x_misc.c             x_misc_source.h \
 x_slist.c  x_slist.h x_slist_struct.h \
 x_stream.c \
 x_umalloc.c          x_umalloc_source.h \
 x_urng.c   x_urng.h \
\
 auto.c  auto.h  auto_struct.h \
 arou.c  arou.h  arou_struct.h \
 empk.c  empk.h  empk_struct.h \
 cstd.c  cstd.h  cstd_struct.h \
 dari.c  dari.h  dari_struct.h \
 dau.c   dau.h   dau_struct.h \
 dgt.c   dgt.h   dgt_struct.h \
 dsrou.c dsrou.h dsrou_struct.h \
 dstd.c  dstd.h  dstd_struct.h \
 hinv.c  hinv.h  hinv_struct.h \
 ninv.c  ninv.h  ninv_struct.h \
 srou.c  srou.h  srou_struct.h \
 ssr.c   ssr.h   ssr_struct.h \
 tabl.c  tabl.h  tabl_struct.h \
 tdr.c   tdr.h   tdr_struct.h tdr_source.h tdr_newset.ch tdr_init.ch tdr_sample.ch tdr_debug.ch \
 unif.c  unif.h  unif_struct.h \
 utdr.c  utdr.h  utdr_struct.h \
 vempk.c vempk.h vempk_struct.h \
 vmt.c   vmt.h   vmt_struct.h

# clean backup files
CLEANFILES = \
 *~ \
 source_fp_const.h

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

