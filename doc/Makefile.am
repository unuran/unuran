## Process this file with automake to produce Makefile.in
## $Id$

SUBDIRS = src figures . 

info_TEXINFOS = unuran.texi

# Location of same source files 
TEXINFO_TEX = ./src/texinfo.tex
AM_MAKEINFOFLAGS = -I $(top_srcdir)/doc -I $(top_srcdir)/doc/src

# Automatically created part of manual
unuran_src = ./src/unuran_src.texi

# Files with examples.
# These example files are generated from the corresponding
# files $(top_srcdir)/examples/example?.c
unuran_examples = \
 ./src/ref_example0.texi	./src/ref_example0_str.texi \
 ./src/ref_example1.texi	./src/ref_example1_str.texi \
 ./src/ref_example2.texi	./src/ref_example2_str.texi \
 ./src/ref_example3.texi	./src/ref_example3_str.texi \
 ./src/ref_example_gsl.texi	\
 ./src/ref_example_rngstreams.texi	\
 ./src/ref_example_anti.texi	./src/ref_example_anti_str.texi \
 ./src/ref_example_cont.texi	./src/ref_example_cont_str.texi \
 ./src/ref_example_discr.texi	./src/ref_example_discr_str.texi \
 ./src/ref_example_emp.texi	./src/ref_example_emp_str.texi \
 ./src/ref_example_vemp.texi	./src/ref_example_FuncStr.texi

# Sources of texinfo file
unuran_TEXINFOS = \
 unuran.texi \
 ./src/unuran_method_requirements.texi \
 $(unuran_src) \
 $(unuran_examples) \
 $(top_srcdir)/doc/figures/*.eps \
 $(top_srcdir)/doc/figures/*.pdf \
 $(top_srcdir)/doc/figures/*.png \
 $(top_srcdir)/doc/figures/*.txt

noinst_DATA = 

EXTRA_DIST = \
 $(unuran_src) \
 $(unuran_examples) \
 unuran.html \
 unuran.txt \
 unuran.dvi \
 unuran.ps \
 unuran.pdf


# Generated texi files ...
sinclude ./src/.dep-unuran_src_texi
$(unuran_src): $(top_srcdir)/doc/src/make_texi.pl
	$(top_srcdir)/doc/src/make_texi.pl $(top_srcdir) > $(unuran_src)

# Generated texi files for examples ...
$(unuran_examples): $(top_srcdir)/examples/example*.c
	sed 's/{/@{/g;s/}/@}/g' < $(top_srcdir)/examples/`echo $* | sed 's/src\/ref_//'`.c > $@

# Rules for creating info, dvi, ps, and pdf files are provided by automake.
# Add new rules for html and txt.
SUFFIXES = .html .txt

.texi.html:
	$(MAKEINFO) $(AM_MAKEINFOFLAGS) --html --no-split $< -o tmp-$@
	sed -e 's/\&OUML\;/\&Ouml\;/g' tmp-$@ > $@
	rm -f tmp-$@

.texi.txt:
	$(MAKEINFO) $(AM_MAKEINFOFLAGS) --fill-column=80 --no-headers --number-sections $< -o $@

# Dependencies
unuran.html: unuran.texi version.texi $(unuran_TEXINFOS)
unuran.txt: unuran.texi version.texi $(unuran_TEXINFOS)

# Abbreviations
.PHONY: html txt
html: unuran.html
txt: unuran.txt

# Clean backup and temp files
CLEANFILES = \
 *~ \
 tmp-unuran.html

# Clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 mdate-sh \
 $(unuran_examples) \
 $(unuran_src) \
 unuran.html \
 unuran.txt

