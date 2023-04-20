
# UNU.RAN -- Universal Non-Uniform Random number generator

-------------------------------------------------------------------------------

__UNU.RAN__ is an ANSI C library licensed under GPL. 
It contains universal (also called automatic or black-box) algorithms
that can generate random numbers from large classes of continuous or
discrete distributions, and also from practically all standard
distributions.

The library and an extensive online documentation are available at:

<https://statmath.wu.ac.at/unuran/>

---------------------------------------------------------

## A short overview

To generate random numbers the user must supply some information about
the desired distribution, especially a C-function that computes the
density and - depending on the chosen methods - some additional
information (like the borders of the domain, the mode, the derivative
of the density ...). After a user has given this information an
init-program computes all tables and constants necessary for the
random variate generation. The sample program can then generate
variates from the desired distribution.

The main part of __UNU.RAN__ are different universal algorithms 
(called methods).
The choice of the method depends on the information available for
the distribution and on the desired characteristics of the algorithm
(fast initialization and slow sampling, slow initialization and
fast sampling).

A second important part of __UNU.RAN__ is the distribution module
containing all necessary functions for many continuous and discrete
univariate standard distributions. Thus __UNU.RAN__ can be used  without
extra coding to obtain very fast generators for the best known
standard distributions.

__UNU.RAN__ is coded in ANSI C but uses an object oriented programming
interface. There are three objects:

+ *distribution objects*  
  contain all information of the distribution.
		
+ *parameter objects*  
  contain all input parameters (and defaults) for the
  different methods.
		
+ *generator objects*  
  are produced by the initialization program and contain all
  table for sampling from the given distribution.everything

Of course a uniform random number generator is necessary to use
__UNU.RAN__. We provide an interface to use the PRNG uniform package from
the pLab team from the University of Salzburg (Austria),
available at <https://statmath.wu.ac.at/prng/>, 
to the GSL library and a library Rngstream
(<https://statmath.wu.ac.at/software/RngStreams/>) by Pierre L'Ecuyer.
It is also no problem to include any other uniform random number
generator.

---------------------------------------------------------

## Advantages of UNU.RAN


Why can it be worth the time to download __UNU.RAN__ and to understand the
concept of its interface? Isn't it much faster to implement a simple
standard method for the distribution I am interested in?

- The first and main advantage lies in the modeling flexibility you
  gain for your simulation. Once you have installed __UNU.RAN__ you can
  sample from practically all uni-modal (and other) distributions
  without coding more than the density functions. For a big number of
  standard distributions (and truncated versions of these standard
  distributions) you need not even code the densities as these are
  already included in __UNU.RAN__.

- It is possible to sample from non-standard distribution. In fact
  only a pointer to a function that returns e.g. the density at a
  given point *x* is required.

- Distributions can be exchanged easily. For example it is not
  difficult at all to start your simulation with the normal
  distribution, and switch to an empirical distribution later.

- The library contains reliable and fast generation algorithms. The
  characteristics of some these algorithms (like speed, expected
  number of uniforms required etc, ...) are only slightly influenced
  by the chosen distribution. (However numerical inversion is included
  as a (very slow) brute force algorithm for the rare cases where the
  more sophisticated methods do not work.)

- Correlation induction facilities are included. 

---------------------------------------------------------

## Packages

The __UNU.RAN__ is also accessible from other languages:

+ __[R](https://www.r-project.org/)__: 
  CRAN package
  *[Runuran](https://CRAN.R-project.org/package=Runuran)*.  
  The sources are available from <https://github.com/unuran/Runuran>.

+ __[Python](https://www.python.org/)__:
  __UNU.RAN__ is in included in *[SciPy](https://scipy.org/)*.
  The sources are available from <https://github.com/scipy/unuran>.
  
+ Packages for __[GNU Octave](https://octave.org/)__ and 
  __[MATLAB](https://de.mathworks.com/products/matlab.html)__ do exist
  but are not released yet.
  
---------------------------------------------------------

## Build and check


> `./autogen.sh` (or download the tar ball from <https://statmath.wu.ac.at/unuran/>)
> `./configure`
> `make`

Run `./configure --help` for building options.

---------------------------------------------------------

## References

* J. Leydold and W. Hörmann:
  [UNU.RAN User Manual](https://statmath.wu.ac.at/unuran/).

* W. Hörmann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg
