##################################################################
#
# Configuration file for ANURAN tests.
#
#################################################################

#################################################################
# Constants

$sample_size = 2000;
$accuracy = 1.0e-7;

#################################################################
# Methods

@method_list =
    ( 'method = tdr; c = -0.5',
      'method = tdr; c = 0' );

#################################################################
# Distributions

@distr_list = 
    ( 
      'beta (1 .. 10, 1 .. 10)',
      'beta (10 .. 1000, 10 .. 1000)',
      'beta (0.1 .. 1, 0.1 .. 1)',
      'beta (1 .. 10, 1 .. 10, -100 .. -10, -9 .. 100)',
      'beta (1 .. 10, 1 .. 10); domain=(0,0.5)',
      'beta (1 .. 10, 1 .. 10); domain=(0.5,1.5)',

      'cauchy ()',
      'cauchy (-10 .. 10)',
      'cauchy (1, 0.1 .. 10)',
      'cauchy (); domain=(0,inf)',
      'cauchy (); domain=(-inf,0)',

      'chi (0.1 .. 1)',
      'chi (1 .. 100)',

      'chisquare (0.1 .. 1)',
      'chisquare (1 .. 100)',

      'exponential ()',
      'exponential (0.1 .. 10)',
      'exponential (0.1 .. 10, -10 .. 100)',
      
      'extremeI ()',
      'extremeI (-10 .. 100)',
      'extremeI (-10 .. 100, 0.1 .. 10)',
      
      'extremeII (0.1 .. 10)',
      'extremeII (0.1 .. 10, -10 .. 100)',
      'extremeII (0.1 .. 10, -10 .. 100, 0.1 .. 10)',
      
      'gamma (0.1 .. 1)',
      'gamma (1 .. 10)',
      'gamma (1 .. 10, 0.1 .. 10)',
      'gamma (1 .. 10, 0.1 .. 10, -10 .. 10)',
      
      'laplace ()',
      'laplace (-10 .. 100)',
      'laplace (-10 .. 100, 0.1 .. 10)',
      
      'logistic ()',
      'logistic (-10 .. 100)',
      'logistic (-10 .. 100, 0.1 .. 10)',

      'lomax (0.1 .. 10)',
      'lomax (0.1 .. 10, 0.001 .. 1000)',

      'normal ()',
      'normal (2.)',
      'normal (-5 .. 5, 1 .. 10)',
      'normal (); domain=(-3,0.23)',
      'normal (); domain=(-inf,inf)',
      'normal (); domain=(0,inf)',

      'pareto (0.001 .. 1000, 0.1 .. 3)',
      'pareto (0.001 .. 1000, 0.1 .. 3)',

      'powerexponential (0.1 .. 1)',
      'powerexponential (1 .. 5)',

      'rayleigh (0.1 .. 10)',

      'student (0.1 .. 0.9)',
      'student (1 .. 100)',
      
      'triangular (0 .. 1)',
      
      'uniform ()',
      'uniform (-100 .. -10, -9 .. 10)',

      'weibull (0.1 .. 10)',
      'weibull (0.1 .. 10, 0.1 .. 10)',
      'weibull (0.1 .. 10, 0.1 .. 10, -10 .. 10)',

      );

# End
1;
