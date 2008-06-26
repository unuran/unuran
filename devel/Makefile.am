## Process this file with automake to produce Makefile.in
## $Id$

AUTOMAKE_OPTIONS = foreign

SUBDIRS = . src run tests

# compile with almost all GCC warnings enabled
setgccWflags:
	export CFLAGS="`cat $(unuran_srcdir)/scripts/gccWflags`"; ./configure

setgccWflags-ansi:
	export CFLAGS="-ansi `cat $(unuran_srcdir)/scripts/gccWflags`"; ./configure

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 aclocal.m4 \
 configure \
 config.h.in \
 COPYING \
 INSTALL \
 autoconf/compile \
 autoconf/config.guess \
 autoconf/config.sub \
 autoconf/depcomp \
 autoconf/install-sh \
 autoconf/ltmain.sh \
 autoconf/mdate-sh \
 autoconf/missing \
 autoconf/mkinstalldirs

