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
    ( '-c -0.5',
      '-c 0' );

#################################################################
# Distributions

@distr_list = 
    ( 
      '-d beta -p "1 .. 10, 1 .. 10"',
      '-d beta -p "10 .. 1000, 10 .. 1000"',
      '-d beta -p "0.1 .. 1, 0.1 .. 1"',
      '-d beta -p "1 .. 10, 1 .. 10, -100 .. -10, -9 .. 100"',
      '-d beta -p "1 .. 10, 1 .. 10" -D "0 0.5"',
      '-d beta -p "1 .. 10, 1 .. 10" -D "0.5 1.5"',

      '-d cauchy',
      '-d cauchy -p "-10 .. 10"',
      '-d cauchy -p "1, 0.1 .. 10"',
      '-d cauchy -D "0 inf"',
      '-d cauchy -D "-inf 0"',

      '-d chi -p "0.1 .. 1"',
      '-d chi -p "1 .. 100"',

      '-d chisquare -p "0.1 .. 1"',
      '-d chisquare -p "1 .. 100"',

      '-d exponential',
      '-d exponential -p "0.1 .. 10"',
      '-d exponential -p "0.1 .. 10, -10 .. 100"',
      
      '-d extremeI',
      '-d extremeI -p "-10 .. 100"',
      '-d extremeI -p "-10 .. 100, 0.1 .. 10"',
      
      '-d extremeII -p "0.1 .. 10"',
      '-d extremeII -p "0.1 .. 10, -10 .. 100"',
      '-d extremeII -p "0.1 .. 10, -10 .. 100, 0.1 .. 10"',
      
      '-d gamma -p "0.1 .. 1"',
      '-d gamma -p "1 .. 10"',
      '-d gamma -p "1 .. 10, 0.1 .. 10"',
      '-d gamma -p "1 .. 10, 0.1 .. 10, -10 .. 10"',
      
      '-d laplace',
      '-d laplace -p "-10 .. 100"',
      '-d laplace -p "-10 .. 100, 0.1 .. 10"',
      
      '-d logistic',
      '-d logistic -p "-10 .. 100"',
      '-d logistic -p "-10 .. 100, 0.1 .. 10"',

      '-d lomax -p "0.1 .. 10"',
      '-d lomax -p "0.1 .. 10, 0.001 .. 1000"',

      '-d normal',
      '-d normal -p "2."',
      '-d normal -p "-5 .. 5, 1 .. 10"',
      '-d normal -D "-3 0.23"',
      '-d normal -D "-inf inf"',
      '-d normal -D "0 inf"',

      '-d pareto -p "0.001 .. 1000, 0.1 .. 3"',
      '-d pareto -p "0.001 .. 1000, 0.1 .. 3"',

      '-d powerexponential -p "0.1 .. 1"',
      '-d powerexponential -p "1 .. 5"',

      '-d rayleigh -p "0.1 .. 10"',
      
      '-d triangular -p "0 .. 1"',
      
      '-d uniform',
      '-d uniform -p "-100 .. -10, -9 .. 10"',

      '-d weibull -p "0.1 .. 10"',
      '-d weibull -p "0.1 .. 10, 0.1 .. 10"',
      '-d weibull -p "0.1 .. 10, 0.1 .. 10, -10 .. 10"',

      );

# End
1;
