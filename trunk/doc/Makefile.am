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
 fig-pdf-cauchy.png \
 fig-pdf-chisquare.png \
 fig-pdf-exponential.png \
 fig-pdf-gamma.png \
 fig-pdf-laplace.png \
 fig-pdf-lognormal.png \
 fig-pdf-normal.png \
 fig-pdf-pareto.png \
 fig-pdf-powerexponential.png \
 fig-pdf-rayleigh.png \
 fig-pdf-student.png

# make html files
html:	version-ref.texi
	texi2html -verbose -split_chapter $(info_TEXINFOS)

# clean backup files
clean-local:
	rm -f *~

# clean html files
maintainer-clean:
	rm -f *.html

