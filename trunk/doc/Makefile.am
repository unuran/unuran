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

# make html files
html:	version-ref.texi
	texi2html -verbose -split_chapter $(info_TEXINFOS)

# clean backup files
clean-local:
	rm -f *~

# clean html files
maintainer-clean:
	rm -f *.html

