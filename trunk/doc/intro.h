
=NODE Intro Introduction
=UP TOP [10]
=EON

/*---------------------------------------------------------------------------*/

=NODE  Explain  What is UNURAN?
=UP  Intro [10]

=DESCRIPTION
   UNURAN (Universal Non-Uniform RAndom Number generator) is
   especially designed for such situations where 

   @itemize @minus
   @item a non-standard distribution or a truncated distribution is
	 needed.

   @item experiments with different types of distributions are made.

   @item random variates for variance reduction techniques are used.

   @item fast generators of predictable quality are necessary.

   @end itemize

   Of course it is also well suited for standard distributions. 
   However due to its more sophisticated programming interface it
   might not be as easy to use if you only look for a generator for
   the standard normal distribution. (Although UNURAN provides
   generators that are superior in many aspects to those found in
   quite a number of other libraries.)

   UNURAN implements several methods for generating random numbers.
   The choice depends primary on the information about the
   distribution can be provided and -- if the user is familar with
   the different methods -- on the preferences of the user.

=EON 

/*---------------------------------------------------------------------------*/

=NODE  Usage  Usage of this document
=UP  Intro [20]

=DESCRIPTION
   We designed this document in a way such that one can
   use UNURAN with reading as little as necessary.
   Read @ref{Installation} for the instructions to
   install the library.
   @ref{Concepts,Concepts of UNURAN,Concepts of UNURAN}
   discribes the basics of UNURAN.
   It also has a short guideline for choosing an appropriate method.
   In @ref{Examples} gives examples that can be copied and modified.
   They also can be found in the directory @file{examples} in the
   source tree.

   Further information are given in consecutive chapters.
   @ref{Distribution_objects,Handling distribution objects,Handling distribution objects}
   describes how to create and manipulate distribution objects.
   @ref{Stddist,standard distributions,standard distributions}
   describes predefined distribution objects that are ready to use.
   @ref{Methods} describes the various methods in detail.
   For each of possible distribution classes 
   (continuous, discrete, empirical, multivariate)
   there exists a short overview section that can be used to choose an
   appropriate method followed by sections that describe each of the
   particular methods in detail.
   These are merely for users with some knowledge about
   the methods who want to change method-specific parameters and can
   be ignored by others.

   Abbreviations and explanation of some basic terms can be found in
   @ref{Glossary}.
=EON

/*---------------------------------------------------------------------------*/

=NODE  Installation Installation
=UP Intro [30]

=DESCRIPTION
   UNURAN was develloped on an Intel architecture under Linux with
   the gnu C compiler.


   @subsubheading Uniform random number generator

   It can be used with any uniform random number generator but (at the
   moment) some features work best with Otmar Lendl's @emph{prng}
   library (see @url{http://statistik.wu-wien.ac.at/prng/} for
   description and downloading).
   For more details on using uniform random number in UNURAN
   @pxref{URNG,Using uniform random number generators,Using uniform random number generators}.


   @subsubheading Special mathematical functions

   We used Stephen L. Moshier's @emph{Cephes} mathematical library for
   special function (eg. Gamma function). So we recommend to install
   this library too. It can be found at
   @url{http://www.netlib.org/cephes/}. 

   A local copy that we have used can be downloaded from
   @url{ftp://statistik.wu-wien.ac.at/src/cephes-math-28.tar.gz}.
   Unfortunately there is no comfortable installation script. So you
   have 

   @enumerate
   @item Unpack the tar archiv:
	 @smallexample
	    tar zxvf cephes-math-28.tar.gz
	 @end smallexample
   @item Compile the double library:
	 @smallexample
	    cd cephes/double
	    make
	 @end smallexample
   @item Copy the library to a place in the search path of your
	 compiler (eg. @file{/usr/local/lib}):
	 @smallexample
	    cp libmd.a /usr/local/lib
	    chmod 644 /usr/local/lib/libmd.a
	 @end smallexample
	 (You must be root when you want to install it into 
	 @file{/usr/local/lib}).
   @end enumerate

   If you do not install the Cephes library then some of the predefined
   distrubutions have some of the described features.
   In future releases we plan to include such special functions in to
   the library.

   
   @subsubheading UNURAN

   @enumerate 

   @item First unzip and untar the package and change to the directory:
	 @smallexample
	    tar zxvf unuran-@value{VERSION}.tar.gz
            cd unuran-@value{VERSION}
	 @end smallexample

   @item Edit the file @file{src/unuran_config.h}.
	 Set the appropriate source of uniform random numbers:
	 @code{#define UNUR_URNG_TYPE}
	 (@pxref{URNG} for details).

	 @emph{Important:} If the prng library is not installed you
	 must not use @code{UNUR_URNG_PRNG}.

	 @emph{Warning:} If @code{UNUR_URNG_POINTER} is used then the
	 buildin default uniform random number generators should be
	 used only for small sample sizes or for demonstration. They
	 are not state of the art any more.
	 
   @item Run a configuration script:
	 @smallexample
	    sh ./configure --prefix=<prefix>
	 @end smallexample
	 where @code{<prefix>} is the root of the installation tree.
	 When omitted @file{/usr/local} is used.

	 Use @code{configure --help} to get a list of other options.

	 @emph{Important:} You must install prng and Cephes
	 @emph{before} @code{configure} is executed.

   @item Compile and install the libray:
	 @smallexample
	    make
	    make install
	 @end smallexample
	 This installs the following files:
	 @smallexample
	    $(prefix)/include/unuran.h
	    $(prefix)/include/unuran_config.h
	    $(prefix)/include/unuran_tests.h
	    $(prefix)/lib/libunuran.a
	    $(prefix)/info/unuran.info
         @end smallexample
	 Obviously @code{$(prefix)/include} and @code{$(prefix)/lib}
	 must be in the search path of your compiler. You can use environment
	 variables to add these directories to the search path. If you
	 are using the bash type (or add to your profile):
	 @smallexample
	    export LIBRARY_PATH="HOMEDIRECTORY/lib"
	    export C_INCLURE_PATH="HOMEDIRECTORY/include"
	@end smallexample

   @item Documentation in various formats (PS, PDF, info, dvi, HTML,
	 plain text) can be found in the directory @file{doc}.

   @item You can run some tests my 
	 @smallexample
	    make check
	 @end smallexample
	 However this test suite requires the usage of prng.
	 Notice a number of these tests fail if you have not installed
	 the Cephes library. Moreover it might happen that some of the
	 tests might fail due to roundoff errors or the mysteries of
	 floating point arithmetic, since we have used some extreme
	 settings to test the library.

   @end enumerate

=EON



/*---------------------------------------------------------------------------*/

=NODE  Concepts Concepts of UNURAN
=UP Intro [40]

=DESCRIPTION
   UNURAN is a C library for generating non-uniformly distributed
   random variates. Its emphasis is on the generation of non-standard
   distribution and on streams of random variates of special purposes.
   It is designed to provide a consistent tool to 
   sample from distributions with various properties. 
   Since there is no universal method that fits for all situations,
   various methods for sampling are implemented.

   UNURAN solves this complex task by means of an object oriented
   programming interface. Three basic objects are used:

   @itemize @bullet
   @item distribution object @code{UNUR_DISTR}@*
      Hold all information about the random variates that should be
      generated.

   @item generator object @code{UNUR_GEN}@*
      Hold the generators for the given distributions.
      Two generator objects are completely independent of each other.
      They may share a common uniform random number generator or have
      their owns.

   @item parameter object @code{UNUR_PAR}@*
      Hold all information for creating a generator object. It is
      necessary due to various parameters and switches for each of
      these generation methods. 
   
   @end itemize

   The idea behind these structures is to make distributions, choosing
   a generation method and sampling to orthogonal (ie. independent)
   functions of the library. 
   The parameter object is only introduced due to the necessity to
   deal with various parameters and switches for
   each of these generation methods which are required to adjust the
   algorithms to unusual distributions with extreme properties but
   have defaut values that are suitable for most applications.
   These parameters and the data for distributions are set by various
   functions.

   Once a generator object has been created sampling (from the
   univariate continuous distribution) can be done by
   the following command:
   @example
      double x = unur_sample_cont(generator);
   @end example
   @noindent
   Analogous commands exist for discrete and multivariate
   distributions.
   For detailed examples that can be copied and modified
   @pxref{Examples}.


   @subheading Distribution objects
   
   All information about a distribution are stored in objects
   (structures) of type @code{UNUR_DISTR}.
   UNURAN has five different types of distribution objects:

   @table @code
   @item cont
      Continuous univariate distributions.
   @item cvec
      Continuous multivariate distributions.
   @item discr
      Discrete univariate distributions.
   @item cemp
      Continuous empirical univariate distribution, ie. given by sample.
   @item cvemp
      Continuous empirical multivariate distribution, ie. given by sample.

   @end table

   @noindent
   Distribution objects can be
   created from scratch by the following call
   @example
      distr = unur_sample_<type>_new();
   @end example
   @noindent
   where @code{<type>} is one of the five possible types from the
   above table.
   Notice that these commands only create an @emph{empty} object which
   still must be filled by means of calls for each type of
   distribution object
   (@pxref{Distribution_objects,Handling distribution objects,Handling distribution objects}).
   The naming scheme of these functions is designed to indicate the
   corresponding type of the distribution object and the task to be
   performed. It is demonstated on the following example.
   @example
     unur_distr_cont_set_pdf(distr, mypdf);
   @end example
   @noindent
   This command stores a PDF named @code{mypdf} in the distribution
   object @code{distr} which must have the type @code{cont}.

   Of course UNURAN provides an easier way to use standard distribution.
   Instead of using @command{unur_distr_<type>_new} calls and fuctions
   @command{unur_distr_<type>_set_<@dots{}>} for setting data
   objects for standard distribution can be created by a single call.
   Eg. to get an object for the normal distribution with mean 2 and
   standard deviation 5 use
   @example
   double parameter[2] = @{2.0 ,5.0@};
   UNUR_DISTR *distr = unur_distr_normal(parameter, 2);
   @end example
   @noindent
   For a list of standard distributions
   @pxref{Stddist,Standard distributions,Standard distributions}.


   @subheading Generation methods

   The information a distribution object must contain depends
   heavily on the method choosen for sampling random variates.

   @include unuran_method_requirements.texi

   @sp 2

   Because of tremendous variety of possible problems, UNURAN provides many
   methods. All information for creating an generator object have to
   collected in a parameter first.
   For example if the task is to sample from a continuous distribution
   the method AROU might be a good choice. Then the call
   @example
     UNUR_PAR *par = unur_arou_new(distribution);
   @end example
   @noindent
   creates an parameter object @code{par} with a pointer to the
   distribution object and default values for all necessary parameters
   for method AROU.
   Other methods can be used by replacing @code{arou} with the name
   of the desired methods (in lower case letters):
   @example
     UNUR_PAR *par = unur_<method>_new(distribution);
   @end example
   @noindent
   This sets the default values for all necessary parameters for the
   chosen methods. These are suitable for almost all
   applications. Nevertheless it is possible to control the behaviour
   of the method using corresponding @command{set} calls for each method.
   This might be necessary to adjust the algorithm for an unusual
   distribution with extreme properties, or just for fine tuning the
   perforence of the algorithm.
   The following example demonstrates how to change the maximum
   number of iterations for method NINV to the value 50:
   @example
    unur_ninv_set_max_iteration(par, 50);
   @end example
   All available methods are described in details in 
   @ref{Methods,Methods,Methods}.
  

   @subheading Creating a generator object

   Now it is possible to create a generator object:
   @example
     UNUR_GEN *generator = unur_init(par);
     if (generator == NULL) exit(EXIT_FAILURE);
   @end example

   @noindent
   @strong{Important:} You must always check whether unur_init() has
   been executed successfully. Otherwise the NULL pointer is returned
   which causes a segmentation fault when used for sampling.

   @noindent
   @strong{Important:}
   The call of unur_init() @strong{destroys} the parameter object!@*
   Moreover it is recommended to call unur_init() immediately after
   the parameter object @code{par} has created and modified.

   @sp 1
   An existing generator object is a rather static construct.
   Nevertheless some of the parameters can still be modified by
   @command{chg} calls, e.g.
   @example
     unur_ninv_chg_max_iteration(gen, 30);
   @end example
  
   @sp 1
   Notice that it is important @emph{when} pararameters are
   changed because different functions must be used:

   To change the parameters @emph{before} creating the generator object,
   the function name includes the term @command{set} and the first
   argument must be of type @code{UNUR_PAR}. 

   To change the parameters for an @emph{existing} generator object,
   the function name includes the term @command{chg} and the first
   argument must be of type @code{UNUR_GEN}.

   For details @pxref{Methods,Methods,Methods}.


   @subheading Sampling

   You can now use your generator object in any place of your program
   to sample from your distribution. You only have take about the type
   of number it computes: @code{double}, @code{int} or a vector (array
   of @code{double}s). 
   Notice that at this point it does not matter whether you are
   sampling from a gamma distribution, a truncated normal distribution
   or even an empirical distribution.


   @subheading Destroy

   When you do not need your generator object any more, you should
   destroy it:
   @example
     unur_free(generator);
   @end example

   
   @subheading Uniform random numbers

   Each generator object can have its own uniform random number
   generator or share one with others.
   When created a parameter object the pointer for the uniform random
   number generator is set to the default generator. However it can be
   changed at any time to any other generator:
   @example
     unur_set_urng(par, urng);
   @end example
   @noindent
   or
   @example
     unur_chg_urng(generator, urng);
   @end example
   @noindent
   respectively.
   See @ref{URNG,Using uniform random number generators,Using uniform random number generators}
   for details.

=EON

/*---------------------------------------------------------------------------*/

=NODE  Contact Contact the authors
=UP  Intro [50]

=DESCRIPTION
   If you have any problems with UNURAN, suggestions how to improve
   the library or find a bug, please contact us via email
   @email{unuran@@statistik.wu-wien.ac.at}.

   For news please visit out homepage at
   @uref{http://statistik.wu-wien.ac.at/unuran/}.

=EON

/*---------------------------------------------------------------------------*/


