## Process this file with automake to produce Makefile.in
## $Id$

# AUTOMAKE_OPTIONS = 

# Experimental code
if INCLUDE_EXPERIMENTAL
EXPERIMENTAL_DIRS = run
else
EXPERIMENTAL_DIRS = 
endif

SUBDIRS = . scripts src examples experiments tests doc $(EXPERIMENTAL_DIRS)

EXTRA_DIST = \
 Announcement-0.1.0 \
 KNOWN-PROBLEMS \
 autogen.sh

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 configure \
 aclocal.m4 \
 config.guess \
 config.h.in \
 config.sub \
 depcomp \
 install-sh \
 ltmain.sh \
 missing \
 mkinstalldirs

# remove windows directory
maintainer-clean-local:
	rm -rf win
