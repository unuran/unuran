
/*---------------------------------------------------------------------------*/

=NODE  Stddist   UNURAN Library of standard distributions
=UP TOP [55]

=DESCRIPTION
Although it is not its primary target, many
distributions are already implemented in UNURAN.
This section presents these available distributions
and their parameters.

The syntax to get a distribuion object for distributions
@code{<dname>} is:

@deftypefn -- {UNUR_DISTR*} unur_distr_<dname> (double* @var{params}, int @var{n_params})
@var{params} is an array of doubles of size 
@var{n_params} holding the parameters.
@end deftypefn

@noindent
E.g. to get an object for the gamma distribution (with shape parameter) use
@smallexample
unur_distr_gamma( params, 1 );
@end smallexample

Distributions may have default parameters with need not be given
explicitely.
E.g. The gamma distribution has three parameters: the
shape, scale and location parameter. Only the (first) shape parameter
is required. The others can be omitted and are then set by default
values.

@sp 1
@smallexample
@cc_start alpha = 5; default: beta = 1, gamma = 0 @cc_stop
double fpar[] = @{5.@};
unur_distr_gamma( fpar, 1 );

@cc_start alpha = 5, beta = 3; default: gamma = 0 @cc_stop
double fpar[] = @{5., 3.@};
unur_distr_gamma( fpar, 2 );

@cc_start alpha = 5, beta = 3, gamma = -2
double fpar[] = @{5., 3., -2.@};
unur_distr_gamma( fpar, 3 );
@end smallexample

@sp 1
@strong{Important:} Naturally the computational accuracy
limits the possible parameters. There shouldn't be problems
when the parameters of a distribution are in a ``reasonable'' range but
e.g. the normal distribution N(10^15,1) won't yield the desired results.
(In this case it would be better generating N(0,1) and @emph{then}
transform the results.)
@* Of course computational inaccuracy is not specific to UNURAN
and should always be kept in mind when working with computers.

@emph{Important:} The routines of the standard library are included
for non-uniform random variate generation and not to provide special
functions for statistical computations.

@subheading Remark

The following keywords are used in the tables:
@table @i
@item PDF
probability density function,
with variable @i{x}.

@item PMF
probability mass function,
with variable @i{k}.

@item constant
given PDF and PMF, resp., must be divided by @i{constant}.

@item domain
domain PDF and PMF, resp.

@item parameters  @var{n_std} (@var{n_total}): @r{list}
list of parameters for distribution, where @var{n_std} is the number
of parameters for the standard form of the distribution and
@var{n_total} the number for the non-standard form. 
@var{list} is the list of parameter in the order as they are stored in
the array of parameters. Optional parameter that can be omitted are
enclosed in square brackets @code{[@dots{}]}.

A detailed list of these parameters gives then the range of valid
parameters and defaults for optional parameters that are used when
these are omitted.

@item reference
gives reference for distribution
(@pxref{Bibliography}).

@end table

=EON

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

=NODE  Stddist_CONT    UNURAN Library of continuous univariate distributions
=UP Stddist [10]
=EON

/*---------------------------------------------------------------------------*/

=NODE  Stddist_CVEC    UNURAN Library of continuous multivariate distributions
=UP Stddist [20]
=EON

/*---------------------------------------------------------------------------*/

=NODE  Stddist_DISCR    UNURAN Library of discrete univariate distributions
=UP Stddist [30]
=EON

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

