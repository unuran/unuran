
=NODE Intro Introduction
=UP TOP [10]
=EON

/*---------------------------------------------------------------------------*/

=NODE  Explain  What is UNURAN?
=UP  Intro [10]

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
=UP  Intro [20]

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
   UNURAN was develloped on an Intel architecture under Linux with
   the gnu C compiler.
   
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

=NODE  Concepts Concepts of UNURAN
=UP Intro [40]

=DESCRIPTION
   UNURAN is a C library for generating nonuniformly distributed
   random numbers.
   It is designed to provide a simple tool to produce random numbers with
   various properties. Because of the divergent requirements to the
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

   When these structures are created it is possible to sample the
   disired random numbers with the following or a similar command:
   @example
    x = unur_sample_cont(generator);
   @end example
   Where the variable @code{x} is a double, @code{generator}
   is a structure of the type @code{UNUR_GEN} and 
   @code{unur_sample_cont()} is a double valued function.


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
   and instead of <type> one of the type @code{cont}   
   @code{cvec}   @code{discr}   @code{cemp}   @code{vemp}
   as shown below must be used. Depending on <type> you have to use
   a function @code{unur_saple_<cont|discr|vec>()}
   to sample random numbers.
   @itemize @bullet
   @item @code{unur_distr_cont_new()}@*
        This function is used to create a distribution object
	if the cdf is continuous.@*
	Use for sampling: @code{unur_sample_cont()} 
   @item  @code{unur_distr_cvec_new()}@*
	This function is used to create a distribution object
	if you want to sample from a continuous multivariate
	distribution.@*
        Use for sampling: @code{unur_sample_vec()}
   @item @code{unur_distr_discr_new()}@*
	If your distribution is descrete, use this function to
	create the distribution object.@*
        Use for sampling: @code{unur_sample_discr()}
   @item @code{unur_distr_cemp_new()}@*
        If you want to create random numbers similar to a given
	set of numbers, use this function to get a 
	distribution object.@*
	Use for sampling: @code{unur_sample_cont()}
   @item  @code{unur_distr_vemp_new()}@*
        If you want to create a vector of random numbers similar
	to a given set of vectors, this function will produce
	the distribution object.@*
	Use for sampling: @code{unur_sample_vec()}
   @end itemize

   Notice that these commands only create an @emph{empty} structure which
   still must be filled with the help of a special set of
   fuctions---depending on the type
   @code{<type>} of the distribution object. 
   The function reference of these types of distribution objects
   are found in the corresponding sections XXXXX-XXXXX. The naming scheme
   of these functions is designed to indicate the corresponding type of
   the distribution object and the task to be performed. To demonstrate
   this an example will be given:
   @example
     unur_distr_cont_set_pdf(distr, mypdf);
   @end example
   This command stores a pdf named @code{mypdf} to the distribution
   object @code{distr} which has the type @code{cont}.
   The function @code{mypdf} which returns a
   @code{double} must be provided by the user.

   The information the user must fill into the distribution object depends
   heavily on the method choosen for sampling random numbers. 

   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 

   Tabelle!!!   

   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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
  All implemented standard distribution are described in section XXXXXX
 
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
  the functions---a set for each method---describes in section XXXXX.
  The following example demonstrates how to change the maximum
  number of iterations of the method NINV to the value 50:
  @example
    unur_ninv_set_max_iteration(par, 50);
  @end example
  The available methods are described in section XXXXX.
  
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
  All functions corresponding to a specific method are found in
  section XXXXX.

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


