## Process this file with automake to produce Makefile.in
## $Id$


info_TEXINFOS = unuran.texi

# these example files are generated from the corresponding
# files $(top_srcdir)/examples/example?.c
unuran_examples = \
 ref_example1.texi \
 ref_example2.texi \
 ref_example3.texi \
 ref_example4.texi \
 ref_example_cont.texi \
 ref_example_discr.texi \
 ref_example_emp.texi \
 ref_example_vemp.texi

unuran_TEXINFOS = \
 unuran.texi \
 unuran_src.texi \
 $(unuran_examples) \
 version.texi

noinst_DATA = 

scripts = \
 make_texi.pl

EXTRA_DIST = \
 examples.h \
 glossary.h \
 index.h \
 intro.h \
 methods.h \
 misc.h \
 references.h \
 stddist.h \
 top.h \
 unuran.html \
 unuran.txt \
 unuran.dvi \
 unuran.ps \
 unuran.pdf \
 texinfo.tex \
 $(scripts)

# generated texi files ...
unuran_src.texi: $(sripts) unuran.texi
	./make_texi.pl $(top_srcdir) > unuran_src.texi

example_dir  = $(top_srcdir)/examples
$(unuran_examples): $(example_dir)/example*.c
	sed 's/{/@{/g;s/}/@}/g' < $(example_dir)/`echo $* | sed 's/ref_//'`.c > $@


SUFFIXES = .html .txt .pdf

# make html file ...
.texi.html:
	$(MAKEINFO) --html $< -o $@

# make plain text file ...
.texi.txt:
	$(MAKEINFO) --fill-column=80 --no-headers --number-sections $< -o $@

# make pdf file ...
.texi.pdf:
	texi2pdf $<

# abbreviations
.PHONY: html txt pdf ps
html: version.texi unuran.html
txt:  version.texi unuran.txt
pdf:  version.texi unuran.pdf
ps:   version.texi unuran.ps

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 unuran.html \
 unuran.txt \
 unuran.pdf \
 mdate-sh \
 $(unuran_examples) \
 unuran_src.texi

