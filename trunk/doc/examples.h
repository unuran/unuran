
=NODE  Examples  Examples
=UP TOP [20]

=DESCRIPTION

Assuming UNURAN has been installed successfully, here are some examples
showing how to obtain random numbers. Beginning with the basics of using
UNURAN we evolve to examples using more advanced features.@*
Call the compiler with:
@smallexample
  gcc -Wall -O2 example.c -lunuran -lprng -lm -lmd -o outfile
@end smallexample 
@*
If the PRNG-library is not installed, @code{-lprng} must be
omitted and of course can't be used.
@*
The internals of UNURAN are hidden behind three basic structures:
The @emph{distribution object},
the @emph{parameter object}
and the @emph{generator object}.
It is not important to understand the details of these objects but
to observe the order of creating them which is the order given above.
The actual sampling is done with the @emph{generator object} while
the other two objects are mainly aids for creating it. 
(The @emph{paramter object} is destroyed during the creation of
the @emph{generator object} and can't be reused.)
It is also recommended not to change the @emph{distribution object}
between the creation of the other two objects.

=EON

/*---------------------------------------------------------------------------*/

=NODE  Ex1  Ex1: As short as possible
=UP Examples [10]

=DESCRIPTION

@smallexample
@include ref_example1.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Ex2  Ex2: Arbitrary distributions
=UP Examples [20]

=DESCRIPTION

If you want to sample from a non-standard distribution,
UNURAN might be exactly what you need. 
Depending on the information is available, a method
must be choosen for sampling.
--> REFERENCE to table!!!

@smallexample
@include ref_example2.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Ex3  Ex3: Changing Parameters of the method
=UP Examples [30]

=DESCRIPTION

Each method for generating random numbers allows several
parameters to be modified. If you don't wand UNURAN to choose
default values, it is possible to change them.
--> REFERENCE to Methods
The following example will illustrate how to change parameters.

@smallexample
@include ref_example3.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  Ex4  Ex4: Changing the uniform random generator
=UP Examples [40]

=DESCRIPTION

To generate special uniformly distributed random numbers,
the PRNG-package can be used (it must be installed seperately).
This can be done by changing two lines in the file
@file{unuran_config.h}. Find the following lines in the file
and modify them as shown below:
@smallexample
  /* set type of uniform generator             */
  /* #define UNUR_URNG_TYPE UNUR_URNG_POINTER  */
  #DEFINE UNUR_URNG_TYPE UNUR_URNG_PRGN
@end smallexample
@* (Just comment one line and uncomment the next) 

@smallexample
@include ref_example4.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  ExCont  ExCont: 
=UP Examples [50]

=DESCRIPTION

@smallexample
@include ref_example_cont.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  ExDiscr  ExDiscr:
=UP Examples [60]

=DESCRIPTION

@smallexample
@include ref_example_discr.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  ExEmp  ExEmp:
=UP Examples [70]

=DESCRIPTION

@smallexample
@include ref_example_emp.texi
@end smallexample

=EON

/*---------------------------------------------------------------------------*/

=NODE  ExVemp  ExVemp:
=UP Examples [80]

=DESCRIPTION

@smallexample
 ref_example_vemp.texi
@end smallexample


/*---------------------------------------------------------------------------*/

