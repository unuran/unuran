## Process this file with automake to produce Makefile.in
## $Id$

AUTOMAKE_OPTIONS = foreign

SUBDIRS = . run

## SUBDIRS = . src run tests

# compile with almost all GCC warnings enabled
setgccWflags:
	export CFLAGS="`cat $(unuran_srcdir)/scripts/gccWflags`"; ./configure

setgccWflags-ansi:
	export CFLAGS="-ansi `cat $(unuran_srcdir)/scripts/gccWflags`"; ./configure

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 $(srcdir)/Makefile.in \
 $(srcdir)/autoconf/compile \
 $(srcdir)/autoconf/config.guess \
 $(srcdir)/autoconf/config.sub \
 $(srcdir)/autoconf/depcomp \
 $(srcdir)/autoconf/install-sh \
 $(srcdir)/autoconf/ltmain.sh \
 $(srcdir)/autoconf/mdate-sh \
 $(srcdir)/autoconf/missing \
 $(srcdir)/autoconf/test-driver

# Remark: the three files
#   $(srcdir)/aclocal.m4
#   $(srcdir)/config.h.in
#   $(srcdir)/configure
# are removed in run/Makefile.am.
# It cannot be done here as this breaks the cleaning process.



#MAINTAINERCLEANFILES = \
# $(srcdir)/aclocal.m4 \
# $(srcdir)/configure \
# $(srcdir)/config.h.in

