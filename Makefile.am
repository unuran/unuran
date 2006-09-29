## Process this file with automake to produce Makefile.in
## $Id$

# AUTOMAKE_OPTIONS = 

# Experimental code
if INCLUDE_EXPERIMENTAL
EXPERIMENTAL_DIRS = run
else
EXPERIMENTAL_DIRS = 
endif

if INCLUDE_SOURCE_ONLY
SUBDIRS = . src
else
SUBDIRS = . scripts src examples experiments tests $(EXPERIMENTAL_DIRS) doc
endif

DIST_SUBDIRS = $(SUBDIRS)

EXTRA_DIST = \
 Announcement-0.1.0 \
 KNOWN-PROBLEMS \
 autogen.sh


# make special distribution for ROOT
rootdir = tmp-unuran-root
distrootfile = $(PACKAGE)-$(VERSION)-root.tar.gz

remove_rootdir = \
  { test ! -d $(rootdir) \
    || { find $(rootdir) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(rootdir); }; }

dist-root: $(PACKAGE)-$(VERSION).tar.gz
	@echo "Create tar ball for ROOT" | \
          sed 'h;s/./=/g;p;x;p;x'
	@$(remove_rootdir)
	@make dist
	mkdir $(rootdir)
	tar zxf $(PACKAGE)-$(VERSION).tar.gz -C $(rootdir) 
	(cd $(rootdir)/$(PACKAGE)-$(VERSION); \
	   sed -r -e 's/(AM_CONDITIONAL\(\s*INCLUDE_SOURCE_ONLY)\s*,.*/\1,true)/' \
	   < configure.ac > configure.ac.tmp.1; \
	   sed -r -e 's/(AM_CONDITIONAL\(\s*INCLUDE_EXPERIMENTAL)\s*,.*/\1,false)/' \
	   < configure.ac.tmp.1 > configure.ac.tmp.2; \
	   sed -r -e 's/(AC_INIT\(\s*\[UNURAN\]\s*,\s*\[)(.*)\]\s*,\s*\[/\1\2-root\],\[/' \
	   < configure.ac.tmp.2 > configure.ac; \
	   rm -v configure.ac.tmp*; \
	   autoreconf -v;  \
	   ./configure; \
	   make dist; \
	)
	cp $(rootdir)/$(PACKAGE)-$(VERSION)/$(PACKAGE)-*.tar.gz .
	$(remove_rootdir)
	@echo "$(distrootfile) is ready for distribution" | \
          sed 'h;s/./=/g;p;x;p;x'


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
	@$(remove_rootdir)

