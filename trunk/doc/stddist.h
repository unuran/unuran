
/*---------------------------------------------------------------------------*/

=NODE  Stddist   UNURAN Library of standard distributions
=UP TOP [55]

=DESCRIPTION
Although it is not its primary target, many
distributions are already implemented in UNURAN.
This section presents these available distributions
and the belonging parameters.

The syntax of using one of these distributions is always the same (see ).
@c @ref{Ex1} and @ref{Ex3}
If standard values for parameters of a specific distribution
are available, it is possible to pass as second argument
(the number of passed arguments) a number smaller than the
number of required parameters.
Of course you have to provide at least as much as parameters
as the second argument indicates.
The first paramters (as much as the value of the second argument)
are set to the provided values and the remaining are set to the
standard values in the order specified in this section.
@c (See @ref{Ex1} (0 parameters provided) and @ref{Ex3} (2 parameters provided).)

Some methods allow arbitrary domains with predefined densities
(e.g. NINV -- numerical inversion):
And if the domain is changed so might the mode. How to change the
domain and due to that the mode or mayby another parameter can be
looked up in 
@c @ref{Do_it_yourself}.


@noindent @strong{WARNING:} Naturally the computational accuracy
limits the possible parameters. There shouldn't be problems
when the parameters of a distribution are in a ``reasonable'' range but
e.g. the normal distribution N(10^50,1) won't yield the desired results.
(In this case it would be better generating N(0,1) and @emph{then}
transform the results.)
@* Of course computational inaccuracy is not specific to UNURAN
and should always be kept in mind when working with computers.

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

