/*---------------------------------------------------------------------------*/

=NODE  Methods   Methods for generating non-uniform random variates
=UP TOP [40]
=EON

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

=NODE  Methods_for_CONT    Methods for continuous univariate distributions
=UP Methods [10]

=DESCRIPTION

@subheading Overview of methods

@cartouche
@noindent 
Methods for @b{continuous univariate distributions}@*
sample with @command{unur_sample_cont}
@sp 1
@multitable {method} {PDF} {dPDF} {mode} {area} {build-in standard distributionxx}
@item method @tab PDF @tab dPDF @tab mode @tab area @tab other
@item AROU   @tab  x  @tab  x   @tab  [x] @tab      @tab T-concave
@item CSTD   @tab     @tab      @tab      @tab      @tab build-in standard distribution
@item NINV   @tab [x] @tab      @tab      @tab      @tab CDF 
@item SROU   @tab  x  @tab      @tab   x  @tab   x  @tab T-concave
@item SSR    @tab  x  @tab      @tab   x  @tab   x  @tab T-concave
@item TABLE  @tab  x  @tab      @tab   x  @tab  [~] @tab all local extrema
@item TDR    @tab  x  @tab  x   @tab      @tab      @tab T-concave
@item UTDR   @tab  x  @tab      @tab   x  @tab   ~  @tab T-concave
@end multitable
@end cartouche

@subheading Example

@smallexample
@include ref_example_cont.texi
@end smallexample

@subheading Example (String API)

@smallexample
@include ref_example_cont_str.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Methods_for_CEMP    Methods for continuous empirical univariate distributions
=UP Methods [20]
=DESCRIPTION

@subheading Overview of methods

@cartouche
@noindent
Methods for @b{continuous empirical univariate distributions}@*
sample with @command{unur_sample_cont}
@sp 1
@noindent
EMPK:  Requires an observed sample.
@end cartouche

@subheading Example

@smallexample
@include ref_example_emp.texi
@end smallexample

@subheading Example (String API)

@smallexample
@include ref_example_emp_str.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Methods_for_CVEC    Methods for continuous multivariate distributions
=UP Methods [30]
=DESCRIPTION

@subheading Overview of methods

@cartouche
@noindent
Methods for @b{continuous multivariate distributions}@*
sample with @command{unur_sample_vec}
@sp 1
@noindent
VMT: Requires the mean vector and the covariance matrix.
@end cartouche

=EON

/*---------------------------------------------------------------------------*/

=NODE  Methods_for_CVEMP   Methods for continuous empirical multivariate distributions
=UP Methods [40]
=DESCRIPTION

@subheading Overview of methods

@cartouche
@noindent
Methods for @b{continuous empirical multivariate distributions}@*
sample with @command{unur_sample_vec}
@sp 1
@noindent
VEMPK: Requires an observed sample.
@end cartouche

@subheading Example

@smallexample
@include ref_example_vemp.texi
@end smallexample

@subheading Example (String API)

(not implemented)

=EON

/*---------------------------------------------------------------------------*/

=NODE  Methods_for_DISCR   Methods for discrete univariate distributions
=UP Methods [50]
=DESCRIPTION

@subheading Overview of methods

@cartouche
@noindent 
Methods for @b{discrete univariate distributions}@*
sample with @command{unur_sample_discr}
@sp 1
@multitable {method} {PMF} {PV} {mode} {area} {build-in standard distributionxx}
@item method@tab PMF @tab PV @tab mode @tab sum @tab other
@item DARI  @tab  x  @tab    @tab  x   @tab  ~   @tab T-concave
@item DAU   @tab [x] @tab  x @tab      @tab      @tab 
@item DGT   @tab [x] @tab  x @tab      @tab      @tab 
@item DSTD  @tab     @tab    @tab      @tab      @tab build-in standard distribution
@end multitable
@end cartouche

@subheading Example

@smallexample
@include ref_example_discr.texi
@end smallexample

@subheading Example (String API)

@smallexample
@include ref_example_discr_str.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Methods_for_UNID    Methods for uniform univariate distributions
=UP Methods [60]
=DESCRIPTION
=EON

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

