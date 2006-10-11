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

dist-root:
	@echo "Create tar ball for ROOT" | \
          sed 'h;s/./=/g;p;x;p;x'
	$(remove_rootdir)
	make
	make dist
	mkdir $(rootdir)
	tar zxf $(PACKAGE)-$(VERSION).tar.gz -C $(rootdir) 
	(cd $(rootdir)/$(PACKAGE)-$(VERSION); \
	   sed -e 's/\(AM_CONDITIONAL(\s*INCLUDE_SOURCE_ONLY\)\s*,.*/\1,true)/' \
		-e 's/\(AM_CONDITIONAL(\s*INCLUDE_EXPERIMENTAL\)\s*,.*/\1,false)/' \
		-e 's/\(AC_INIT(\s*\[UNURAN\]\s*,\s*\[\)\(.*\)\]\s*,\s*\[/\1\2-root\],\[/' \
	   < configure.ac > configure.ac.tmp; \
	   mv -v configure.ac.tmp configure.ac; \
	   for f in `find . -type f -name '*.[ch]'`; do \
		./scripts/remove_comments.pl $$f; \
	   done; \
	   for f in `find . -type f -name '*.[cd]h'`; do \
		./scripts/remove_comments.pl $$f; \
	   done; \
	   echo "AC_MSG_NOTICE([ " >> configure.ac; \
	   echo "*===================================================" >> configure.ac; \
	   echo "*" >> configure.ac; \
	   echo "*   This is a pure source distribution!" >> configure.ac; \
	   echo "*" >> configure.ac; \
	   echo "*   For a manual, documented sources, tests," >> configure.ac; \
	   echo "*   and examples please download" >> configure.ac; \
	   echo "*" >> configure.ac; \
	   echo "*       $(PACKAGE)-$(VERSION).tar.gz" >> configure.ac; \
	   echo "*" >> configure.ac; \
	   echo "*===================================================" >> configure.ac; \
	   echo "])" >> configure.ac; \
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

