## Process this file with automake to produce Makefile.in
## $Id$

SUBDIRS = utils methods specfunct distr distributions parser tests uniform .

AM_CFLAGS = -Wall

lib_LIBRARIES = libunuran.a

libunuran_a_SOURCES = 

libunuran_a_LIBADD = \
 methods/*.o \
 distr/*.o \
 distributions/*.o \
 specfunct/*.o \
 parser/*.o \
 tests/*.o \
 utils/*.o \
 uniform/*.o


include_HEADERS = \
 unuran.h \
 unuran_config.h

noinst_HEADERS = \
 unuran.h.in \
 unur_source.h \
 unur_cookies.h \
 unur_struct.h \
 unur_typedefs.h

##unuran.h: $(top_srcdir)/src/methods/*.h $(top_srcdir)/src/distributions/*.h 
unuran.h:
	$(top_srcdir)/scripts/merge_h.pl $(top_srcdir)/src/unuran.h.in > unuran.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 unuran.h
