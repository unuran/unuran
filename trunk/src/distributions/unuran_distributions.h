/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_distributions.h                                              *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for PDF, CDF, etc.         *
 *         of distribtions.                                                  *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [1] N.L. Johnson, S. Kotz, and A.W. Kemp                                *
 *       Univariate Discrete Distributions,                                  *
 *       2nd edition                                                         *
 *       John Wiley & Sons, Inc., New York, 1992                             *
 *                                                                           *
 *   [2] N.L. Johnson, S. Kotz, and N. Balakrishnan                          *
 *       Continuous Univariate Distributions,                                *
 *       Volume 1, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1994                             *
 *                                                                           *
 *   [3] N.L. Johnson, S. Kotz, and N. Balakrishnan                          *
 *       Continuous Univariate Distributions,                                *
 *       Volume 2, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1995                             *
 *                                                                           *
 *   [4] N.L. Johnson, S. Kotz, and N. Balakrishnan                          *
 *       Discrete Multivariate Distributions,                                *
 *       John Wiley & Sons, Inc., New York, 1997                             *
 *                                                                           *
 *   [5] S. Kotz, N. Balakrishnan, and N.L. Johnson                          *
 *       Continuous Multivariate Distributions,                              *
 *       Volume 1: Models and Applications                                   *
 *       John Wiley & Sons, Inc., New York, 2000                             *
 *                                                                           *
 *   [0] S. Kotz and N.L. Johnson                                            *
 *       Encyclopedia of Statistical Sciences                                *
 *       Volumes 1-9                                                         *
 *       John Wiley & Sons, Inc., New York, 1982-1988                        *
 *                                                                           *
 *****************************************************************************
     $Id$
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifndef __UNURAN_DISTRIBUTIONS_H_SEEN
#define __UNURAN_DISTRIBUTIONS_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_stddistr.h>

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                    Continuous univariate distributions                    *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/*  Beta distribution  [3; ch.25, p.210]                                     */
/* 
   =DISTR    beta  Beta distribution
   =UP       Stddist_CONT
   =REF      [JKBc95]   ch.25, p.210
   =PDF      (x-a)^{p-1} * (b-x)^{q-1}
   =CONST    1 / (Beta(p,q) * (b-a)^{p+q-1})
   =DOMAIN   a < x < b
   =FPARAM    0  : p : > 0 :   : scale    :
              1  : q : > 0 :   : scale    :
             [2] : a :     : 0 : location, scale :
             [3] : b : > a : 1 : location, scale :
   =EON
*/
UNUR_DISTR *unur_distr_beta(double *params, int n_params);
/** TODO: STDGEN **/

/*---------------------------------------------------------------------------*/
/* Burr family of distributions  [2; ch.12, p.54]                            */
UNUR_DISTR *unur_distr_burr(double *params, int n_params);
/** under construction **/

/*---------------------------------------------------------------------------*/
/* Cauchy distribution  [2; ch.16, p.299]                                    */
/* 
   =DISTR    cauchy Cauchy distribution
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.16, p.299
   =PDF      \frac{1}{1 + ((x-theta)/lambda)^2}
   =CONST    \frac{1}{pi * lambda}
   =DOMAIN   -infinity < x < infinity 
   =FPARAM    [0]   : theta  :     : 0 : location :
             [[1]]  : lambda : > 0 : 1 : scale    :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_cauchy(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Chi distribution  [2; ch.18, p.417]                                       */
/* 
   =DISTR    chi Chi distribution
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.18, p.417
   =PDF      x^{nu-1} * exp( -x^2/2 )
   =CONST    1 / (2^{(nu/2)-1} * Gamma(nu/2))
   =DOMAIN   0 <= x < infinity 
   =FPARAM   0 : nu : > 0 : : shape :
   =STDGEN   DEF  Ratio of Uniforms with shift (only for @tex $nu >= 1$@end tex) [MJa87]
   =EON
*/
UNUR_DISTR *unur_distr_chi(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Chisquare distribution  [2; ch.18, p.416]                                 */
/* 
   =DISTR    chisquare Chisquare distribution
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.18, p.416
   =PDF      x^{(nu/2)-1} * exp( -x/2 )
   =CONST    1 / (2^{nu/2} * Gamma(nu/2))
   =DOMAIN   0 <= x < infinity 
   =FPARAM   0 : nu : > 0 : : shape (degrees of freedom) :
   =EON
*/
UNUR_DISTR *unur_distr_chisquare(double *params, int n_params);
/** TODO: STDGEN **/

/*---------------------------------------------------------------------------*/
/* Erlang distribution                                                       */
/* not implemented */

/*---------------------------------------------------------------------------*/
/*  Exponential distribution  [2; ch.19, p.494]                              */
/* 
   =DISTR    exponential Exponential distribution
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.19, p.494
   =PDF      exp( -\frac{x-theta}{sigma})
   =CONST    \frac{1}{sigma}
   =DOMAIN   theta <= x < infinity 
   =FPARAM    [0]  : sigma : > 0 : 1 : scale    :
	     [[1]] : theta :     : 0 : location :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_exponential(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Extreme value type I distribution  [3; ch.22, p.2]                       */
/* 
   =DISTR    extremeI  Extreme value type I (Gumbel-type) distribution
   =UP       Stddist_CONT
   =REF      [JKBc95]   ch.22, p.2
   =PDF      exp( -exp( -\frac{x-zeta}{theta} ) - \frac{x-zeta}{theta} )
   =CONST    \frac{1}{theta}
   =DOMAIN   -infinity < x <infinity
   =FPARAM    [0]  : zeta  :     : 0 : location :
	     [[1]] : theta : > 0 : 1 : scale    :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_extremeI(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Extreme value type II distribution  [3; ch.22, p.2]                      */
/* 
   =DISTR    extremeII  Extreme value type II (Frechet-type) distribution
   =UP       Stddist_CONT
   =REF      [JKBc95]   ch.22, p.2
   =PDF      exp( -(\frac{x-zeta}{theta})^{-k}) * (\frac{x-zeta}{theta})^{-k-1}
   =CONST    \frac{k}{theta}
   =DOMAIN   zeta < x <infinity
   =FPARAM     0   : k     : > 0 :   : shape    :
              [1]  : zeta  :     : 0 : location :
	     [[2]] : theta : > 0 : 1 : scale    :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_extremeII(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Gamma distribution  [2; ch.17, p.337]                                     */
/* 
   =DISTR    gamma  Gamma distribution
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.17, p.337
   =PDF      (\frac{x-gamma}{beta})^{alpha-1} * exp( -\frac{x-gamma}{beta} )
   =CONST    1 / (beta * Gamma(alpha))
   =DOMAIN   gamma < x < infinity 
   =FPARAM     0    : alpha : > 0 :   : shape    :
              [1]   : beta  : > 0 : 1 : scale    :
	     [[2]]  : gamma :     : 0 : location :
   =STDGEN   DEF  Acceptance Rejection combined with Acceptance Complement [ADa74] [ADa82]
             2    Rejection from log-logistic envelopes [CHa77]
   =EON
*/
UNUR_DISTR *unur_distr_gamma(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Generalized inverse Gaussian distribution  [2; ch.15, p.284]              */
UNUR_DISTR *unur_distr_gig(double *params, int n_params);
/** TODO: STDGEN **/

/*---------------------------------------------------------------------------*/
/*  Laplace distribution  [3; ch.24, p.164]                                  */
/* 
   =DISTR    laplace  Laplace distribution
   =UP       Stddist_CONT
   =REF      [JKBc95]   ch.24, p.164
   =PDF      exp( -\frac{|x-theta|}{phi} )
   =CONST    \frac{1}{2 * phi}
   =DOMAIN   -infinity < x <infinity
   =FPARAM    [0]  : theta :     : 0 : location :
	     [[1]] : phi   : > 0 : 1 : scale    :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_laplace(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Logistic distribution  [3; ch.23, p.115]                                  */
/* 
   =DISTR    logistic  Logistic distribution
   =UP       Stddist_CONT
   =REF      [JKBc95]   ch.23, p.115
   =PDF      exp(-\frac{x-alpha}{beta} * (1 + exp(-\frac{x-alpha}{beta}))^{-2}
   =CONST    \frac{1}{beta}
   =DOMAIN   -infinity < x <infinity
   =FPARAM    [0]  : alpha :     : 0 : location :
	     [[1]] : beta  : > 0 : 1 : scale    :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_logistic(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Lognormal distribution  [2; ch.14, p.208]                                */
UNUR_DISTR *unur_distr_lognormal(double *params, int n_params);
/** TODO: STDGEN **/

/*---------------------------------------------------------------------------*/
/*  Lomax distribution (Pareto distr. of second kind)  [2; ch.20, p.575]     */
/* 
   =DISTR    lomax  Lomax distribution (Pareto distribution of second kind)
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.20, p.575
   =PDF      (x+C)^{-(a+1)}
   =CONST    a * C^a
   =DOMAIN   0 <= x < infinity 
   =FPARAM    0  : a : > 0 :   : shape :
             [1] : C : > 0 : 1 : scale :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_lomax(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Normal distribution  [2; ch.13, p.80]                                     */
/* 
   =DISTR    normal  Normal distribution
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.13, p.80
   =PDF      exp( -\frac{1}{2} * (\frac{x-mu}{sigma})^2 )
   =CONST    \frac{1}{sigma * sqrt{2 pi}}
   =DOMAIN   -infinity < x < infinity 
   =FPARAM    [0]   : mu    :     : 0 : location :
             [[1]]  : sigma : > 0 : 1 : scale    :
   =STDGEN   DEF  ACR method (Acceptance-Complement Ratio) [HDa90]
             1    Box-Muller method [BMa58]
	     2    Polar method with rejection [MGa62]
	     3    Kindermann-Ramage method [KRa76]
             INV  Inversion method (slow)
   =EON
*/
UNUR_DISTR *unur_distr_normal( double *params, int n_params );

/*---------------------------------------------------------------------------*/
/* Pareto distribution (of first kind)  [2; ch.20, p.574]                    */
/* 
   =DISTR    pareto  Pareto distribution (of first kind)
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.20, p.574
   =PDF      x^{-(a+1)}
   =CONST    a * k^a
   =DOMAIN   k < x < infinity 
   =FPARAM   0 : k : > 0 :  : shape, location :
             1 : a : > 0 :  : shape           :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_pareto( double *params, int n_params );

/*---------------------------------------------------------------------------*/
/* Pearson VI distribution                                                   */
/* not implemented */

/*---------------------------------------------------------------------------*/
/* Perks distribution                                                        */
/* not implemented */

/*---------------------------------------------------------------------------*/
/* Planck distribution                                                       */
/* not implemented **/

/*---------------------------------------------------------------------------*/
/*  Power-exponential (Subbotin) distribution  [3; ch.24, p.195]             */
/* 
   =DISTR    powerexponential  Powerexponential (Subbotin) distribution
   =UP       Stddist_CONT
   =REF      [JKBc95]   ch.24, p.195
   =PDF      exp( -|x|^tau )
   =CONST    1 / (2 * Gamma(1+1/tau))
   =DOMAIN   -infinity < x < infinity
   =FPARAM   0 : tau : > 0 : : shape :
   =STDGEN   DEF  Transformed density rejection (only for @tex $tau >= 1$@end tex) [DLa86]
   =EON
*/
UNUR_DISTR *unur_distr_powerexponential(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Rayleigh distribution  [2; ch.18, p.456]                                  */
/* 
   =DISTR    rayleigh  Rayleigh distribution
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.18, p.456
   =PDF      x * exp( -1/2 * (\frac{x}{sigma})^2 )
   =CONST    \frac{1}{sigma^2}
   =DOMAIN   0 <= x < infinity
   =FPARAM   0 : sigma : > 0 : : scale :
   =EON
*/
UNUR_DISTR *unur_distr_rayleigh(double *params, int n_params);
/** TODO: STDGEN **/

/*---------------------------------------------------------------------------*/
/* Snedecor's F distribution                                                 */
/* not implemented */


/*---------------------------------------------------------------------------*/
/* Student's t distribution  [3; ch. 28; p. 362]                             */
UNUR_DISTR *unur_distr_student(double *params, int n_params);
/** CDF not implemented !!!!!! */
/** TODO: STDGEN **/

/*---------------------------------------------------------------------------*/
/* Slash distribution  [2; ch.12, p.63]                                      */
UNUR_DISTR *unur_distr_slash(double *params, int n_params);
/** TODO: STDGEN **/

/*---------------------------------------------------------------------------*/
/*  Triangular distribution  [3; ch.26, p.297]                               */
/* 
   =DISTR    triangular  Triangular distribution
   =UP       Stddist_CONT
   =REF      [JKBc95]   ch.26, p.297
   =PDF      2*x / H,          \hbox{ for } 0 <= x <= H \hfill\break
             2*(1-x) / (1-H),  \hbox{ for } H <= x <= 1 
   =CONST    1
   =DOMAIN   0 <= x <= 1
   =FPARAM   [0] : H : 0 <= H <= 1 : 1/2 : shape :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_triangular(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Uniform distribution  [3; ch.26, p.276]                                   */
/* 
   =DISTR    uniform  Uniform distribution
   =UP       Stddist_CONT
   =REF      [JKBc95]   ch.26, p.276
   =PDF      \frac{1}{b-a}
   =CONST    1
   =DOMAIN   a < x < b
   =FPARAM   [0] : a :     : 0 : location :
             [1] : b : > a : 1 : location :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_uniform(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Weibull distribution  [2; ch.21, p.628]                                   */
/* 
   =DISTR    weibull  Weibull distribution
   =UP       Stddist_CONT
   =REF      [JKBb94]   ch.21, p.628
   =PDF      (\frac{x-zeta}{alpha})^{c-1} * exp( -(\frac{x-zeta}{alpha})^c )
   =CONST    \frac{c}{alpha}
   =DOMAIN   zeta < x < infinity 
   =FPARAM     0    : c     : > 0 :   : shape    :
              [1]   : alpha : > 0 : 1 : scale    :
	     [[2]]  : zeta  :     : 0 : location :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_weibull(double *params, int n_params);

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                   Continuous multivariate distributions                   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Multinormal distribution  [5; ch.45, p.107]                               */
UNUR_DISTR *unur_distr_multinormal(int dim, double *mean, double *covar);
/* 
   Creates a distribution object for the multinormal distribution with
   @var{dim} components. @var{mean} is an array of size @var{dim}
   A NULL pointer for @var{mean} is interpreted as the zero
   vector (0,@dots{},0).
   @var{covar} is an array of size @var{dim}x@var{dim} and holds the
   covariance matrix, where the rows of the matrix are stored
   consecutively in this array. The NULL pointer can be used
   instead the identity matrix.

   For standard form of the distribution use NULL for @var{mean} and 
   @var{covar}.
*/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                     Discrete univariate distributions                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Binomial distribution  [1; ch.3, p.105]                                   */
/*
   =DISTR    binomial  Binomial distribution
   =UP       Stddist_DISCR
   =REF      [JKKa92]   ch.3, p.105
   =PMF      {n \choose k} * p^k * (1-p)^{n-k}
   =CONST    1
   =DOMAIN   0 <= k <= n
   =FPARAM   0 : n : >= 1      :  : no. of elements  :
             1 : p : 0 < p < 1 :  : shape            :
   =STDGEN   DEF  Ratio of Uniforms/Inversion [STa89]
   =EON
*/
UNUR_DISTR *unur_distr_binomial(double *params, int n_params);
/** No CDF !!! **/

/*---------------------------------------------------------------------------*/
/* Geometric distribution  [1; ch.5.2, p.201]                                */
/*
   =DISTR    geometric  Geometric distribution
   =UP       Stddist_DISCR
   =REF      [JKKa92]   ch.5.2, p.201
   =PMF      p * (1-p)^k
   =CONST    1
   =DOMAIN   0 <= k < infinity
   =FPARAM   0 : p : 0 < p < 1 :  : shape :
   =STDGEN   INV  Inversion method
   =EON
*/
UNUR_DISTR *unur_distr_geometric(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Hypergeometric distribution  [1; ch.6, p.237]                             */
/*
   =DISTR    hypergeometric  Hypergeometric distribution
   =UP       Stddist_DISCR
   =REF      [JKKa92]   ch.6, p.237
   =PMF      {M \choose k} * {N-M \choose n-k} / {N \choose n}
   =CONST    1
   =DOMAIN   max(0,n-N+M) <= k <= min(n,M)
   =FPARAM   0 : N : >= 1        :  : no. of elements  :
             1 : M : 1 <= M <= N :  : shape            :
             2 : n : 1 <= n <= N :  : shape            :
   =STDGEN   DEF  Ratio of Uniforms/Inversion [STa89]
   =EON
*/
UNUR_DISTR *unur_distr_hypergeometric(double *params, int n_params);
/** No CDF !!! **/

/*---------------------------------------------------------------------------*/
/* Logarithmic distribution  [1; ch.7, p.285]                                */
/*
   =DISTR    logarithmic  Logarithmic distribution
   =UP       Stddist_DISCR
   =REF      [JKKa92]   ch.7, p.285
   =PMF      theta^k / k
   =CONST    - log( 1.-theta);
   =DOMAIN   1 <= k < infinity
   =FPARAM   0 : theta : 0 < theta < 1 :  : shape :
   =STDGEN   DEF   Inversion/Transformation [KAa81]
   =EON
*/
UNUR_DISTR *unur_distr_logarithmic(double *params, int n_params);
/** No CDF !!! **/

/*---------------------------------------------------------------------------*/
/* Negative Binomial distribution  [1; ch.5.1, p.200]                        */
/*
   =DISTR    negativebinomial  Negative Binomial distribution
   =UP       Stddist_DISCR
   =REF      [JKKa92]   ch.5.1, p.200
   =PMF      {k+r-1 \choose r-1} * p^r * (1-p)^k
   =CONST    1
   =DOMAIN   0 <= k < infinity
   =FPARAM   0 : p : 0 < p < 1 :  : shape :
             1 : r : > 0       :  : shape :
   =EON
*/
UNUR_DISTR *unur_distr_negativebinomial(double *params, int n_params);
/** No CDF !!! **/
/** TODO: STDGEN **/

/*---------------------------------------------------------------------------*/
/* Poisson distribution  [1; ch.4, p.151]                                    */
/*
   =DISTR    poisson  Poisson distribution
   =UP       Stddist_DISCR
   =REF      [JKKa92]   ch.4, p.151
   =PMF      theta^k / k!
   =CONST    exp(theta)
   =DOMAIN   0 <= k < infinity
   =FPARAM   0 : theta : > 0 :  : shape :
   =STDGEN   DEF  Tabulated Inversion combined with Acceptance Complement [ADb82]
             2    Tabulated Inversion combined with Patchwork Rejection [ZHa94]
   =EON
*/
UNUR_DISTR *unur_distr_poisson(double *params, int n_params);
/** No CDF !!! **/

/*---------------------------------------------------------------------------*/
/* Zipf (or Zeta) distribution  [1; ch.11.20, p.465]                         */
UNUR_DISTR *unur_distr_zipf(double *params, int n_params);
/** No CDF !!! **/
/** TODO: STDGEN **/

/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_DISTRIBUTIONS_H_SEEN */
/*---------------------------------------------------------------------------*/

