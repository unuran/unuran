## Process this file with automake to produce Makefile.in
## $Id$


info_TEXINFOS = unuran.texi

# these example files are generated from the corresponding
# files $(top_srcdir)/examples/example?.c
unuran_examples = \
 ref_example1.texi \
 ref_example2.texi \
 ref_example3.texi \
 ref_example4.texi

unuran_TEXINFOS = \
 unuran.texi \
 $(unuran_examples) \
 unuran_src.texi

noinst_DATA = 

scripts = \
 make_texi.pl \
 DISTRIBUTION_to_texi.pl \
 ERROR_to_texi.pl \
 METHOD_to_texi.pl \
 URNG_to_texi.pl

EXTRA_DIST = \
 unuran.html \
 unuran.txt \
 unuran.dvi \
 unuran.ps \
 unuran.pdf \
 texinfo.tex \
 $(scripts)

# generated texi files ...
unuran_src.texi: $(sripts)
	./make_texi.pl $(top_srcdir)/src/ > unuran_src.texi

example_dir  = $(top_srcdir)/examples
$(unuran_examples): $(example_dir)/example*.c
	sed 's/{/@{/g;s/}/@}/g' < $(example_dir)/`echo $* | sed 's/ref_//'`.c > $@


SUFFIXES = .html .txt .pdf

# make html file ...
.texi.html: unuran.info
	$(MAKEINFO) --html $< -o $@

# make plain text file ...
.texi.txt: unuran.info
	$(MAKEINFO) --fill-column=80 --no-headers --number-sections $< -o $@

# make pdf file ...
.texi.pdf: unuran.info
	texi2pdf $<

# abbreviations
.PHONY: html txt pdf ps
html: version.texi unuran.html
txt: version.texi unuran.txt
pdf: version.texi unuran.pdf
ps: version.texi unuran.ps

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




# to DELETE

# paths ...
#method_dir   = $(top_srcdir)/src/methods
#distr_dir    = $(top_srcdir)/src/distributions
#example_dir  = $(top_srcdir)/examples

#ref_geninfo.texi: $(method_dir)/x_gen.h
#	./geninfo.perl $(method_dir)/x_gen.h      > ref_geninfo.texi

#ref_distr.texi: $(method_dir)/distr.h
#	./distr.perl $(method_dir)/distr.h        > ref_distr.texi

#ref_distr_discr.texi: $(distr_dir)/d_*.c
#	./distr_discr.perl $(distr_dir)/d_*.c     > ref_distr_discr.texi

#ref_distr_cont.texi: $(distr_dir)/c_*.c
#	./distr_cont.perl $(distr_dir)/c_*.c      > ref_distr_cont.texi

#ref_distr_lib.texi: $(distr_dir)/unuran_distributions.h
#	./distr_lib.perl $(distr_dir)/unuran_distributions.h > ref_distr_lib.texi

#ref_methods.texi: $(method_dir)/*.h
#	./methods.perl $(method_dir)/*.h          > ref_methods.texi

#ref_method_functions.texi: $(method_dir)/*.h
#	./method_functions.perl $(method_dir)/*.h > ref_method_functions.texi

#ref_internals.texi: $(method_dir)/*.h
#	./internals.perl $(method_dir)/*.h        > ref_internals.texi


#texi_generated = \
# ref_geninfo.texi \
# ref_distr.texi \
# ref_distr_discr.texi \
# ref_distr_cont.texi \
# ref_distr_lib.texi \
# ref_methods.texi \
# ref_method_functions.texi \
# ref_internals.texi \
# $(unuran_examples)

