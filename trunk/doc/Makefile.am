## Process this file with automake to produce Makefile.in
## $Id$

SUBDIRS = src . 

info_TEXINFOS = unuran.texi

# these example files are generated from the corresponding
# files $(top_srcdir)/examples/example?.c
unuran_examples = \
 ref_example0.texi ref_example0_str.texi \
 ref_example1.texi ref_example1_str.texi \
 ref_example2.texi ref_example2_str.texi \
 ref_example3.texi ref_example3_str.texi \
 ref_example4.texi ref_example4_str.texi \
 ref_example_anti.texi  ref_example_anti_str.texi \
 ref_example_cont.texi  ref_example_cont_str.texi \
 ref_example_discr.texi ref_example_discr_str.texi \
 ref_example_emp.texi   ref_example_emp_str.texi \
 ref_example_vemp.texi  ref_example_FuncStr.texi


unuran_TEXINFOS = \
 unuran.texi \
 unuran_src.texi \
 unuran_method_requirements.texi \
 $(unuran_examples) \
 version.texi

noinst_DATA = 

scripts = \
 make_texi.pl

EXTRA_DIST = \
 unuran.html \
 unuran.txt \
 unuran.dvi \
 unuran.ps \
 unuran.pdf \
 $(scripts)

# location of same source files 
TEXINFO_TEX = ./src/texinfo.tex

# generated texi files ...
unuran_src.texi: $(scripts) unuran.texi $(top_srcdir)/src/unuran.h  
	./make_texi.pl $(top_srcdir) > unuran_src.texi

example_dir  = $(top_srcdir)/examples
$(unuran_examples): $(example_dir)/example*.c
	sed 's/{/@{/g;s/}/@}/g' < $(example_dir)/`echo $* | sed 's/ref_//'`.c > $@


# rules for creating info, dvi, ps, and pdf files are provided by automake.
# add new rules for html and txt.
SUFFIXES = .html .txt

.texi.html:	unuran.info
	$(MAKEINFO) --html --no-split $< -o $@

.texi.txt:
	$(MAKEINFO) --fill-column=80 --no-headers --number-sections $< -o $@

# dependencies
unuran.html: unuran.texi version.texi $(unuran_TEXINFOS)
unuran.txt: unuran.texi version.texi $(unuran_TEXINFOS)

# abbreviations
.PHONY: html txt
html: unuran.html
txt: unuran.txt

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 unuran.html \
 unuran.txt \
 unuran.pdf \
 $(unuran_examples) \
 unuran_src.texi

