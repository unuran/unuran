
=NODE Intro Introduction
=UP TOP [10]
=EON

/*---------------------------------------------------------------------------*/

=NODE  Overview  Overview
=UP Intro [10]
=EON

/*---------------------------------------------------------------------------*/

=NODE  Explain  What is UNURAN?
=UP  Overview [10]

=DESCRIPTION
   UNURAN (Universal Non-Uniform RAndom Number generator) is
   designed for applications using non-standard distributions
   where no seperate generator is yet available.
   UNURAN is also useful if the parameters of a standard distribution
   should be changed frequently.
   Additionaly there are lots of predefined distributions available
@c    (see @ref{Distributions}).

   UNURAN implements several methods for generating random numbers.
   The choice depents primary on what information about the
   distribution can be provided and---if the user is familar with
   the different methods---on the preferences of the user.
=EON

/*---------------------------------------------------------------------------*/

=NODE  Usage  Usage of this document
=UP  Overview [20]

=DESCRIPTION
   We tried do design this document in a way one can
   use UNURAN with reading as little as possible.
   On basis of the available information about the density,
   one should choose a suitable method.
   @c ------------ (see @ref{Requirements}).
   Then the Examples @ref{Ex1} and @ref{Ex2} in addition
   with a short look at section 
   @c ----------- @ref{Do_it_yourself}
   should be enough.
   
   The function for changing the parameters while sampling can be
   found in 
   @c ------------ @ref{Function Reference}.

   To use predefined distributions take a look on the list of
   availabe Distributions 
@c    @ref{Distributions}, 
   Examples @ref{Ex1}
   and @ref{Ex1} and if required also on 
   @c ------------- @ref{Readymade}.   

   The greatest deal of this document is for users with knowledge about
   the methods who want to change method-specific parameters and can
   be ignored by all others.
=EON

/*---------------------------------------------------------------------------*/

=NODE  Install Installation
=UP Intro [20]

=DESCRIPTION
   First unzip and untar the package. Then change to the directory
   @file{/unuran} and execute the following commands:
   @smallexample
   configure --prefix=HOMEDIRECTORY
   make
   make install
   @end smallexample
   @c
   @noindent This installs the following files:
   @c
   @smallexample
   $(prefix)/include/unuran_tests.h
   $(prefix)/include/unuran.h
   $(prefix)/include/unuran_config.h
   $(prefix)/lib/libunuran.a
   @end smallexample
   
   @noindent
   If you set the follwing systemvariables:
   @smallexample
   LIBRARY_PATH="HOMEDIRECTORY/lib"
   C_INCLURE_PATH="HOMEDIRECTORY/include"
   @end smallexample
   @noindent
   the relevant files are where the C-compiler can find them.
   If every user on your machine should be able to use UNURAN,
   you omit @file{--prefix=HOMEDIRECTORY} when calling @file{configure}
   (as usual the standard value of @file{$(prefix)} is @file{/usr/local}.
   Of course you can do this only if you have the right permissions.
=EON

/*---------------------------------------------------------------------------*/

=NODE  Contact Contact the authors
=UP  Intro [30]

=DESCRIPTION
   If you have any problems with UNURAN, suggestions how to improve
   the program or find a bug, please contact us via the UNURAN-homepage:
   @uref{http://statistik.wu-wien.ac.at}
   or via email:
   @email{leydold@@statistik.wu-wien.ac.at}
=EON

/*---------------------------------------------------------------------------*/


