#!/usr/bin/perl
####################################################

# Compiler
$GCC = "gcc -Wall";
$G77 = "g77 -Wall";

# Files
$C_src = "./uniftest_c.c";
$C_exec = "./uniftest_c";

$FORTRAN_src = "./uniftest_f.f";
$FORTRAN_exec = "./uniftest_f";

$JAVA_src = "./uniftest_java.java";
$JAVA_exec = "./Uniftest";

####################################################

$seed = int(rand 12345678) + 1;
$sample_size = 10000;
$accuracy = 1.0e-15;

####################################################

# number of different results
$n_diffs = 0;

# Make source files
make_C_src();
make_FORTRAN_src();
make_JAVA_src();

# Compile sources
system "$GCC -o $C_exec $C_src";
system "$G77 -o $FORTRAN_exec $FORTRAN_src";
system "javac $JAVA_src";

# Print Test data
print "seed = $seed\n";
print "sample size = $sample_size\n";
print "accuracy = $accuracy\n";
print "languages = C, FORTRAN, JAVA\n";

# Start generators
open C, "$C_exec |" or die "cannot run $C_exec"; 
open FORTRAN, "$FORTRAN_exec |" or die "cannot run $FORTRAN_exec"; 
open JAVA, "java $JAVA_exec |" or die "cannot run $JAVA_exec"; 

# Run generatores and compare output
$FORTRAN_n_diffs = 0;
$JAVA_n_diffs = 0;
while ($C_out = <C>) {
    $FORTRAN_out = <FORTRAN>;
    $JAVA_out = <JAVA>;
    chomp $C_out;
    chomp $FORTRAN_out;
    chomp $JAVA_out;

    $FORTRAN_diff = abs($C_out - $FORTRAN_out);
    $JAVA_diff = abs($C_out - $JAVA_out);

    if ($FORTRAN_diff > $accuracy) {
	++$FORTRAN_n_diffs;
	print "C = $C_out\tFORTRAN = $FORTRAN_out\tdifference = $FORTRAN_diff\n";
    }
    if ($JAVA_diff > $accuracy) {
	++$JAVA_n_diffs;
	print "C = $C_out\tJAVA = $JAVA_out\tdifference = $JAVA_diff\n";
    }
}

# End
close C;
close FORTRAN;
close JAVA;

####################################################

$exitcode = 0;

if ($FORTRAN_n_diffs > 0) {
    print "FORTRAN Test FAILED\n";
    ++$exitcode;
}
if ($JAVA_n_diffs > 0) {
    print "JAVA Test FAILED\n";
    ++$exitcode;
}
if ($exitcode == 0) {
    print "Test PASSED\n";
}

####################################################

exit ($exitcode);


####################################################
# Make C version of generator

sub make_C_src
{
    open C_src, ">$C_src" or die "cannot open $C_src for writing";

    print C_src <<EOX;
    
/* ---------------------------------------------------------------- */
/* C version                                                        */
/* ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- */
/* LCG (Linear Congruential Generator) by G. Marsaglia (1972).      */
/*   x_(n+1) = 69069 * x_n + 1 mod 2^32                             */
/* ---------------------------------------------------------------- */

static unsigned int xn = $seed;   /* seed  */

static double uniform (void)
{
  
  xn = 69069 * xn + 1;
  
  if (sizeof(unsigned int) > 4) 
    /* no overflow. 2^32 = 4294967296 */
    xn %= 4294967296;
  
  return (xn * 2.32830643653869628906e-10 + 1.16415321826934814453e-10);
}

static void useed(unsigned int seed)
{
  xn = seed;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>

int main ()
{
  int i;

  useed($seed);

  for (i=0; i<$sample_size; i++)
    printf("%.17e\\n",uniform());

  exit (0);
}

/* ---------------------------------------------------------------- */

EOX

    close C_src;

} # end of make_C_src() 

####################################################
# Make FORTRAN version of generator

sub make_FORTRAN_src
{
    open FORTRAN_src, ">$FORTRAN_src" or die "cannot open $FORTRAN_src for writing";

    print FORTRAN_src <<EOX;

* ------------------------------------------------------------------ *
* FORTRAN version                                                    *
* ------------------------------------------------------------------ *

* ------------------------------------------------------------------ *
* LCG (Linear Congruential Generator) by G. Marsaglia (1972).        *
*   x_(n+1) = 69069 * x_n + 1 mod 2^32                               *
* ------------------------------------------------------------------ *

      DOUBLE PRECISION FUNCTION unif()

      INTEGER xn
      DOUBLE PRECISION f1, f2
      PARAMETER (f1=2.d0**(-32))
      PARAMETER (f2=2.d0**(-33))

C     state variable
      COMMON /state/xn
      DATA xn/$seed/
      SAVE /state/

      xn = 69069 * xn + 1

      unif = xn * f1 + f2
      IF (xn .LT. 0) unif = unif + 1.d0

      END


      SUBROUTINE useed(seed)

      INTEGER seed, xn

C     state variable
      COMMON /state/xn
      SAVE /state/

C     seed generator
      xn = seed

      RETURN
      END

* ------------------------------------------------------------------ *
* ------------------------------------------------------------------ *

      PROGRAM MAIN

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      CALL useed($seed)

      DO 1 i=1,$sample_size
         u = unif()
         WRITE (*,'(d24.18)') u
 1    CONTINUE
      END

* ------------------------------------------------------------------ *

EOX

    close FORTRAN_src;

} # end of make_FORTRAN_src() 

####################################################
# Make JAVA version of generator

sub make_JAVA_src
{
    open JAVA_src, ">$JAVA_src" or die "cannot open $JAVA_src for writing";

    print JAVA_src <<EOX;

/* ---------------------------------------------------------------- */
/* JAVA version                                                     */
/* ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- */
/* LCG (Linear Congruential Generator) by G. Marsaglia (1972).      */
/*   x_(n+1) = 69069 * x_n + 1 mod 2^32                             */
/* ---------------------------------------------------------------- */

public class Urand {

    /* class variables */

    private static int xn = $seed;

    /* member functions */

    public static double myrandom() 
    {

        private double retval;

        xn = 69069 * xn + 1;
        retval = 1.0 * xn;

        retval = 2.32830643653869628906e-10 * retval + 1.16415321826934814453e-10;

	if ( retval  < 0 ) {
           retval += 1.0;
        }

        return (retval);
    }


    public static void useed(int seed) 
    {
	xn = seed;
    }

} /* end of class Urand */

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

public class Uniftest {
   public static void main(String[] args) throws Exception {

      /* set new seed */
      Urand.useed($seed);

      for (private int i = 0; i<$sample_size;i++) {
        System.out.println( Urand.myrandom() );
      }
   }

}  /* end of class Test */

/* ---------------------------------------------------------------- */

EOX

    close JAVA_src;

} # end of make_JAVA_src() 

####################################################
