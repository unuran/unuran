## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src

noinst_LIBRARIES = libutils.a

libutils_a_SOURCES = \
 error.c \
 stream.c \
 umalloc.c \
 urng.c \
 utils.c

pkginclude_HEADERS = \
 unur_cookies.h \
 unur_errno.h \
 unur_math.h \
 unur_umalloc.h \
 unur_urng.h \
 unur_utils.h

# clean backup files too
clean-local:
	rm -f *~
