## Process this file with automake to produce Makefile.in


info_TEXINFOS = unuran-ref.texi

unuran_ref_TEXINFOS = \
 unuran-ref.texi \
 distr.texi \
 intro.texi \
 routines.texi \
 title.texi

noinst_DATA = 

EXTRA_DIST = \
 tex_macros.tex \
 fig-pdf-beta.png \
 fig-pdf-beta.tex \
 fig-pdf-cauchy.png \
 fig-pdf-cauchy.tex \
 fig-pdf-chisquare.png \
 fig-pdf-chisquare.tex \
 fig-pdf-exponential.png \
 fig-pdf-exponential.tex \
 fig-pdf-gamma.png \
 fig-pdf-gamma.tex \
 fig-pdf-laplace.png \
 fig-pdf-laplace.tex \
 fig-pdf-lognormal.png \
 fig-pdf-lognormal.tex \
 fig-pdf-normal.png \
 fig-pdf-normal.tex \
 fig-pdf-pareto.png \
 fig-pdf-pareto.tex \
 fig-pdf-powerexponential.png \
 fig-pdf-powerexponential.tex \
 fig-pdf-rayleigh.png \
 fig-pdf-rayleigh.tex \
 fig-pdf-student.png \
 fig-pdf-student.tex \
 fig-pdf-uniform.tex

method_PATH  = ../src/methods/
distr_PATH   = ../src/distributions/
example_PATH = ../examples/

# make qstart
qstart: qstart.texi

# qstart.texi
qstart.texi: scripts
	tex qstart.texi
	makeinfo --html qstart.texi
	makeinfo --no-headers -o qstart.txt qstart.texi
	makeinfo qstart.texi

# run perl scripts
scripts:
	./examples_to_texi.perl $(example_PATH)example[0-9]*.c
	./geninfo.perl $(method_PATH)x_gen.c
	./change_distr.perl $(method_PATH)distr.h
	./description.perl $(method_PATH)*.h
	./distr.perl $(distr_PATH)c_*.c
	./distr_discrete.perl $(distr_PATH)d_*.c
	./distrref.perl $(distr_PATH)*.h
	./funcref.perl $(method_PATH)*.h
	./internals.perl $(method_PATH)*.h

# make html files
html:
	texi2html -verbose -split_chapter $(info_TEXINFOS)

# clean backup files
clean-local:
	rm -f *~

# clean html files
maintainer-clean:
	rm -f *.html



