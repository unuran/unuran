
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
   Additionally there are lots of predefined distributions available
   (@pxref{Stddist,Standard distributions,Standard distributions}).

   UNURAN implements several methods for generating random numbers.
   The choice depends primary on the information about the
   distribution can be provided and -- if the user is familar with
   the different methods -- on the preferences of the user.

=EON 

/*---------------------------------------------------------------------------*/

=NODE  Usage  Usage of this document
=UP  Overview [20]

=DESCRIPTION
   We designed this document in a way such that one can
   use UNURAN with reading as little as necessary.
   On basis of the available information about the distribution
   one should choose a suitable method.
   Thus we give examples that can be copied and modified
   (@pxref{Examples}).

   Predefined distributions can be found in
   @pxref{Stddist,Standard distributions,Standard distributions}.
   If these should be changed or if it is necessary to build a
   distribution object from scratch
   @pxref{Distribution_objects,Handling distribution objects,Handling distribution objects}.

   The various methods are described in @ref{Methods}.
   For each of possible distribution classes 
   (continuous, discrete, empirical, multivariate)
   there exists a short overview section that can be used to choose an
   appropriate method followed by sections that describe each of the
   particular methods in detail.
   These are merely for users with some knowledge about
   the methods who want to change method-specific parameters and can
   be ignored by others.

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


