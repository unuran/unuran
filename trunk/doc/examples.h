
=NODE  Examples  Examples
=UP TOP [20]

=DESCRIPTION

   The examples in this chapter should compile cleanly and can be
   found in the directory @file{examples} of the source tree of
   UNURAN. Assuming that UNURAN as well as the PRNG and Cephes
   libraries have been installed properly (@pxref{Installation}) each
   of these can be compiled (using the GCC in this example) with
   @example
      gcc -Wall -O2 -o example example.c -lunuran -lprng -lm -lmd
   @end example 
   @noindent
   @emph{Remark:} @code{-lprng} must be omitted when the PRNG library
   is not installed. @code{-lmd} must be omitted when the Cephes
   library is not installed. In both cases however some of the
   examples might not work.

   The library uses three objects:
   @code{UNUR_DISTR}, @code{UNUR_PAR} and @code{UNUR_GEN}.
   It is not important to understand the details of these objects but
   it is important not to changed the order of their creation.
   The distribution object can be destroyed @emph{after} the generator
   object has been made. (The parameter object is freed automatically
   by the unur_init() call.) It is also important to check the result
   of the unur_init() call. If it has failed the NULL pointer is
   returned and causes a segmentation fault when used for sampling.

=EON

/*---------------------------------------------------------------------------*/

=NODE  Example_1  As short as possible
=UP Examples [10]

=DESCRIPTION

@smallexample
@include ref_example1.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Example_2  Arbitrary distributions
=UP Examples [20]

=DESCRIPTION

If you want to sample from a non-standard distribution,
UNURAN might be exactly what you need. 
Depending on the information is available, a method
must be choosen for sampling, 
@pxref{Concepts} for an overview and 
@ref{Methods} for details.

@smallexample
@include ref_example2.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Example_3  Change parameters of the method
=UP Examples [30]

=DESCRIPTION

Each method for generating random numbers allows several
parameters to be modified. If you do not want to use default values,
it is possible to change them.
The following example illustrates how to change parameters.
For details @pxref{Methods}.

@smallexample
@include ref_example3.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Example_4  Change uniform random generator
=UP Examples [40]

=DESCRIPTION

All generator object use the same default uniform random number
generator by default. This can be changed to any generator of your
choice such that each generator object has its own random number
generator or can share it with some other objects.
It is also possible to change the default generator at any time.
See @ref{URNG,Using uniform random number generators,Using uniform random number generators}
for details.

The following example shows how the uniform random number generator
can be set or changed for a generator object. It requires the PRNG
library to be installed and used. Otherwise the example must be
modified accordingly.

@smallexample
@include ref_example4.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Example_More  More examples
=UP Examples [50]
=DESCRIPTION

@xref{Methods_for_CONT,Methods for continuous univariate distributions,Methods for continuous univariate distributions}.

@xref{Methods_for_CEMP,Methods for continuous empirical univariate distributions,Methods for continuous empirical univariate distributions}.

@xref{Methods_for_CVEMP,Methods for continuous empirical multivariate distributions,Methods for continuous empirical multivariate distributions}.

@xref{Methods_for_DISCR,Methods for discrete univariate distributions,Methods for discrete univariate distributions}.

=EON

/*---------------------------------------------------------------------------*/

