
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
   UNURAN is a C library for generating nonuniformly distributed
   random numbers.
   It is designed to provide a simple tool to produce random numbers with
   various properties. Because of the divergent requirements of the
   random numbers various methods for sampling are needed.
   Nevertheless UNURAN solves this complex task in an easy and
   universal way.
   
   To achive this, an object oriented approach with three basic
   structures is introduced:
   @itemize @bullet
   @item distribution object
     (name of the corresponding structure: @code{UNUR_DISTR})
   @item generator object
     (name of the corresponding structure: @code{UNUR_GEN})
   @item parameter object
     (name of the corresponding structure: @code{UNUR_PAR})
   @end itemize

   The intention behind these structures is that the user
   should not change the values within these structures staightly
   but only via the provides functions. This enables anyone to use
   UNURAN without knowing anything about internal structures.

   When these structures are created it is possible to sample the
   desired random numbers with the following or a similar command:
   @example
    x = unur_sample_cont(generator);
   @end example
   The variable @code{x} is a double, @code{generator}
   is a structure of the type @code{UNUR_GEN} and 
   @code{unur_sample_cont()} is a double valued function.
   
   @sp 2

   Of course the user has to provide information about the
   internal properties the random numbers should meet. It is a C
   structure named @code{UNUR_DISTR} which holds this kind of
   information.

   Depending on the properties of the distribution UNURAN uses
   five different functions to create a structure of the kind
   @code{UNUR_DISTR}:
   @example
      distr = unur_sample_<type>_new();
   @end example
   The variable @code{distr} is a structure of type @code{UNUR_DISTR}
   and instead of <type> one of the type @code{cont},   
   @code{cvec},   @code{discr},   @code{cemp} or @code{cvemp}
   as shown below must be used. Depending on <type> you have to use
   a function @code{unur_sample_<cont|discr|vec>()}
   to sample random numbers.
   @itemize @bullet
   @item @code{unur_distr_cont_new()}:
     univariate continous distributions
   @item  @code{unur_distr_cvec_new()}:
     continuous multivariate distributions
   @item @code{unur_distr_discr_new()}:
     discrete univariate distributions
   @item @code{unur_distr_cemp_new()}:
     empirical continuous univariate distribuions
   @item  @code{unur_distr_cvemp_new()}:
     empirical continuous multivariate distributions
   @end itemize

   Notice that these commands only create an @emph{empty} structure which
   still must be filled with the help of a special set of
   fuctions---depending on the type
   @code{<type>} of the distribution object. 
   The function reference of these types of distribution objects
   are found in the corresponding sections
   (@pxref{Distribution_objects,Handling
   distribution objects,Handling distribution objects}). The naming scheme
   of these functions is designed to indicate the corresponding type of
   the distribution object and the task to be performed. To demonstrate
   this an example will be given:
   @example
     unur_distr_cont_set_pdf(distr, mypdf);
   @end example
   This command stores a pdf named @code{mypdf} in the distribution
   object @code{distr} which has the type @code{cont}.
   The function @code{mypdf} which returns a
   @code{double} must be provided by the user.

   The information the user must fill into the distribution object depends
   heavily on the method choosen for sampling random numbers. 
  
   The requirements of the methods are indicated in the following tables: 
    

     @include unuran_method_requirements.texi


   If we want to sample from a continuous distribution and
   the choosen method needs the mode, it can be set with
   the following line:
   @example
    unur_distr_cont_set_mode(distr, mode);  
   @end example
   The argument @code{distr} is of the type @code{UNUR_DISTR} and
   @code{mode} is a double variable holding the position of the
   mode.

  Of course UNURAN provides an easier way to use standard distribution.
  Instead of the usage of @code{unur_distr_<type>_new()} and fuctions of
  the kind @code{unur_distr_<type>_set_<...>()} only one function call
  is necessary, e.g.:
  @example
   UNUR_DISTR distr;
   double parameter[2] = @{2.0 ,5.0@};
   distr = unur_distr_normal(parameter, 2);
  @end example
  This defines a structure @code{distr} of the kind @code{UNUR_DISTR},
  a double array holding parameters of the desired distribution and
  creates @emph{and} initializes a distribution object to the
  normal distribution with mean 2 and variance 5.
  When using ready made distributions the user must know the
  corresponding type <type> of distribution to use the correct
  function for sampling.
  In combination with the previous example,
  the sampling routine `@code{unur_distr_sample_cont()}' must be used.
  All implemented standard distribution are described below.( @pxref{Stddist,standard distributions,standard distributions}).
 
  @sp 2

  Because of tremendous variety of possible problems, UNURAN provides many
  state of the art methods. 
  There exist parameter objects which contain all relevant informations
  about the choosen method similar to the distribution object.
  For example if the task is to sample from a continuous distribution
  the method AROU would be a good choice:
  @example
    par = unur_arou_new(distr);
  @end example
  The variable @code{par} is a structure of the type @code{UNUR_PAR} and
  the function @code{unur_arou_new()} creates @emph{and} initializes
  a parameter object for the method AROU. Of course @code{distr}
  is of the type @code{UNUR_DISTR} as described above.
  Other methods can be used by replacing `@code{arou}' with the name
  of the desired methods (in lower case letters):
  @example
    unur_<method>_new();
  @end example
  This sets all necessary parameters of the methods. Nevertheless
  it is possible to control the behaviour of the method using
  the corresponding functions---a set for each method---(@ref{Methods,Methods,Methods}).
  The following example demonstrates how to change the maximum
  number of iterations of the method NINV to the value 50:
  @example
    unur_ninv_set_max_iteration(par, 50);
  @end example
  The available methods are described later (@pxref{Methods,Methods,Methods}).
  
  Now it is possible to create a generator object:
  @example
    generator = unur_init(par);
  @end example
  The variable @code{generator} is a structure of the type @code{UNUR_GEN}
  and @code{par} is of the type @code{UNUR_PAR}.

  ATTENTION: The call of @code{unur_init()} @emph{destroys} the
  parameter object!  

  Nevertheless it is still possible to adjust some of the
  parameters of the distribution and the method within the
  generator object, e.g:
  @example
    unur_ninv_chg_max_iteration(gen, 30);
  @end example
  
  Please notice that it is important @emph{when} the pararameters are
  changed because different functions must be used:
  To change the parameters @emph{before} creating the generator object,
  the function name includes the term `@code{set}' and the first argument
  must be of type @code{UNUR_PAR}. To change the parameters
  @emph{afterwards}, the function name includes the term `@code{chg}'
  and the first argument must be of type @code{UNUR_GEN}.
  The main advantage of this concept is the possibility to change
  essential parameters during sampling.
  All functions corresponding to a specific method are
  explained later (@pxref{Methods,Methods,Methods}).

  Finally you can sample random number with:
  @example
    x = unur_sample_<type>(gen);
  @end example
  The variable @code{x} holding the random number is of type double,
  @code{gen} is of type @code{UNUR_GEN} and @code{unur_sample<_type>()}
  is one of the sampling fuctions described above.

  After sampling you should release the allocated memory
  using  the functions @code{unur_free(gen)} and @code{unur_distr_free(distr)}.
=EON

/*---------------------------------------------------------------------------*/

=NODE  Contact Contact the authors
=UP  Intro [50]

=DESCRIPTION
   If you have any problems with UNURAN, suggestions how to improve
   the program or find a bug, please contact us via the UNURAN-homepage:
   @uref{http://statistik.wu-wien.ac.at}
   or via email:
   @email{leydold@@statistik.wu-wien.ac.at}
=EON

/*---------------------------------------------------------------------------*/


