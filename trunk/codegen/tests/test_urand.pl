#!/usr/bin/perl
# ****************************************************************

# ----------------------------------------------------------------
# Compiler

$GCC = "gcc -Wall -ansi -pedantic";
$G77 = "g77 -Wall";
$JAVAC = "javac -w1";

# ----------------------------------------------------------------
# Constants

$seed = int(rand 12345678) + 1;
$sample_size = 10000;
$accuracy = 1.0e-15;

# ----------------------------------------------------------------
# Prefix for file names

$file_prefix = "./run_test_urand";

# ----------------------------------------------------------------
# Global log file

open LOG, ">$file_prefix.log" or die "Cannot open log file $file_prefix.log";

# ----------------------------------------------------------------
# Files

$PRNG_exec = "$file_prefix\_PRNG";
$PRNG_src = "$PRNG_exec.c";

$C_exec = "$file_prefix\_C";
$C_src = "$C_exec.c";

$FORTRAN_exec = "$file_prefix\_FORTRAN";
$FORTRAN_src = "$FORTRAN_exec.f";

$JAVA_urand_src = "./Urand.java";
$JAVA_src = "$file_prefix\_JAVA.java";
$JAVA_exec = "$file_prefix\_JAVA";

# ----------------------------------------------------------------
# number of different results

$n_diffs = 0;

# ----------------------------------------------------------------
# Make source files

make_PRNG_src();
make_C_src();
make_FORTRAN_src();
make_JAVA_src();

# ----------------------------------------------------------------
# Compile sources

system "$GCC -o $PRNG_exec $PRNG_src -lprng -lm";
system "$GCC -o $C_exec $C_src";
system "$G77 -o $FORTRAN_exec $FORTRAN_src";
system "$JAVAC $JAVA_src $JAVA_urand_src";

# ----------------------------------------------------------------
# Print Test data

print_log("seed = $seed\n");
print_log("sample size = $sample_size\n");
print_log("accuracy = $accuracy\n");
print_log("languages = C, FORTRAN, JAVA\n");

# ----------------------------------------------------------------
# Start generators

open PRNG, "$PRNG_exec |" or die "cannot run $PRNG_exec"; 
open C, "$C_exec |" or die "cannot run $C_exec"; 
open FORTRAN, "$FORTRAN_exec |" or die "cannot run $FORTRAN_exec"; 
open JAVA, "java $JAVA_exec |" or die "cannot run $JAVA_exec"; 

# ----------------------------------------------------------------
# Run generatores and compare output

$C_n_diffs = 0;
$FORTRAN_n_diffs = 0;
$JAVA_n_diffs = 0;

while ($PRNG_out = <PRNG>) {
    $C_out = <C>;
    $FORTRAN_out = <FORTRAN>;
    $JAVA_out = <JAVA>;

    chomp $PRNG_out;
    chomp $C_out;
    chomp $FORTRAN_out;
    chomp $JAVA_out;

    $C_diff = abs($PRNG_out - $C_out);
    $FORTRAN_diff = abs($PRNG_out - $FORTRAN_out);
    $JAVA_diff = abs($PRNG_out - $JAVA_out);

    if ($C_diff > $accuracy) {
	++$C_n_diffs;
	print LOG "PRNG = $PRNG_out\tC = $C_out\tdifference = $C_diff\n";
    }
    if ($FORTRAN_diff > $accuracy) {
	++$FORTRAN_n_diffs;
	print LOG "PRNG = $PRNG_out\tFORTRAN = $FORTRAN_out\tdifference = $FORTRAN_diff\n";
    }
    if ($JAVA_diff > $accuracy) {
	++$JAVA_n_diffs;
	print LOG "PRNG = $PRNG_out\tJAVA = $JAVA_out\tdifference = $JAVA_diff\n";
    }
}

# ----------------------------------------------------------------
# End

close PRNG;
close C;
close FORTRAN;
close JAVA;

# ----------------------------------------------------------------
# Results

$exitcode = 0;

if ($C_n_diffs > 0) {
    print_log("C Test FAILED\n");
    ++$exitcode;
}
if ($FORTRAN_n_diffs > 0) {
    print_log("FORTRAN Test FAILED\n");
    ++$exitcode;
}
if ($JAVA_n_diffs > 0) {
    print_log("JAVA Test FAILED\n");
    ++$exitcode;
}
if ($exitcode == 0) {
    print_log("Test PASSED\n");
}

# ----------------------------------------------------------------

close LOG;

if ($exitcode) {
    $exitcode = 1;
}

exit $exitcode;

# ****************************************************************
#
# Subroutines
#
# ****************************************************************

# ----------------------------------------------------------------
# Print on screen and log file

sub print_log
{
    my $msg = $_[0];
    print $msg;
    print LOG $msg;
} # end of print_log()

# ****************************************************************
#
# Make PRNG version of generator
#
# ****************************************************************

sub make_PRNG_src
{
    open PRNG_src, ">$PRNG_src" or die "cannot open $PRNG_src for writing";

    print PRNG_src <<EOX;
    
/* ---------------------------------------------------------------- */
/* PRNG version                                                     */
/* ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- */
/* LCG (Linear Congruential Generator) by Park & Miller (1988).     */
/*   x_(n+1) = 16807 * x_n mod 2^31 - 1    (Minimal Standard)       */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <prng.h>

int main ()
{
  int i;
  struct prng *urng;

  urng = prng_new(\"LCG(2147483647,16807,0,$seed)\");
  prng_get_next(urng);   /* synchronize with other generators */

  for (i=0; i<$sample_size; i++)
    printf("%.17e\\n",prng_get_next(urng));

  prng_free(urng);
  exit (0);
}

/* ---------------------------------------------------------------- */

EOX

    close PRNG_src;

} # end of make_PRNG_src() 

# ****************************************************************
#
# Make C version of generator
#
# ****************************************************************

sub make_C_src
{
    open C_src, ">$C_src" or die "cannot open $C_src for writing";

    print C_src <<EOX;
    
/* ---------------------------------------------------------------- */
/* C version                                                        */
/* ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- */
/* LCG (Linear Congruential Generator) by Park & Miller (1988).     */
/*   x_(n+1) = 16807 * x_n mod 2^31 - 1    (Minimal Standard)       */
/* ---------------------------------------------------------------- */

static unsigned int xn = $seed;   /* seed  */

static double urand (void)
{
#define a 16807
#define m 2147483647
#define q 127773      /* m / a */
#define r 2836        /* m % a */

  int hi, lo, test;

  hi = xn / q;
  lo = xn % q;
  test = a * lo - r * hi;
  
  xn = (test > 0 ) ? test : test + m;

  return (xn * 4.656612875245796924105750827e-10);

#undef a
#undef m
#undef q
#undef r
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
    printf("%.17e\\n",urand());

  exit (0);
}

/* ---------------------------------------------------------------- */

EOX

    close C_src;

} # end of make_C_src() 

# ****************************************************************
#
# Make FORTRAN version of generator
#
# ****************************************************************

sub make_FORTRAN_src
{
    open FORTRAN_src, ">$FORTRAN_src" or die "cannot open $FORTRAN_src for writing";

    print FORTRAN_src <<EOX;

* ------------------------------------------------------------------ *
* FORTRAN version                                                    *
* ------------------------------------------------------------------ *

* ------------------------------------------------------------------ *
* LCG (Linear Congruential Generator) by Park & Miller (1988).       *
*   x_(n+1) = 16807 * x_n mod 2^31 - 1    (Minimal Standard)         *
* ------------------------------------------------------------------ *

      DOUBLE PRECISION FUNCTION urand()

      INTEGER a, m, q, r, xn, hi, lo, test
      PARAMETER (a = 16807)
      PARAMETER (m = 2147483647)
      PARAMETER (q = 127773)
      PARAMETER (r = 2836)

C     state variable
      COMMON /state/xn
      DATA xn/$seed/
      SAVE /state/

      hi = xn / q
      lo = MOD(xn,q)

      test = a * lo - r * hi
      IF (test .gt. 0) THEN
	  xn = test
      ELSE
          xn = test + m
      END IF

      urand = xn * 4.656612875245796924105750827D-10

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
         u = urand()
         WRITE (*,'(d24.18)') u
 1    CONTINUE
      END

* ------------------------------------------------------------------ *

EOX

    close FORTRAN_src;

} # end of make_FORTRAN_src() 

# ****************************************************************
#
# Make JAVA version of generator
#
# ****************************************************************

sub make_JAVA_src
{
    open JAVA_urand_src, ">$JAVA_urand_src" or die "cannot open $JAVA_urand_src for writing";

    print JAVA_urand_src <<EOS;

/* ---------------------------------------------------------------- */
/* JAVA version                                                     */
/* ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- */
/* LCG (Linear Congruential Generator) by Park & Miller (1988).     */
/*   x_(n+1) = 16807 * x_n mod 2^31 - 1    (Minimal Standard)       */
/* ---------------------------------------------------------------- */

public class Urand {

    /* class variables */

    private static int xn = $seed;

    /* member functions */

    public static double random() 
    {

        static final int a = 16807;
        static final int m = 2147483647;
        static final int q = 127773;      /* m / a */
        static final int r = 2836;        /* m % a */

	private int hi, lo, test;

        hi = xn / q;
        lo = xn % q;
        test = a * lo - r * hi;
  
        xn = (test > 0 ) ? test : test + m;

        return (xn * 4.656612875245796924105750827e-10);
    }


    public static void useed(int seed) 
    {
	xn = seed;
    }

} /* end of class Urand */

/* ---------------------------------------------------------------- */

EOS
    close JAVA_urand_src;

    open JAVA_src, ">$JAVA_src" or die "cannot open $JAVA_src for writing";

    print JAVA_src <<EOS;

/* ---------------------------------------------------------------- */

public class run_test_urand_JAVA {
   public static void main(String[] args) throws Exception {

      /* set new seed */
      Urand.useed($seed);

      for (private int i = 0; i<$sample_size;i++) {
        System.out.println( Urand.random() );
      }
   }

}  /* end of class Test */

/* ---------------------------------------------------------------- */

EOS

    close JAVA_src;

} # end of make_JAVA_src() 

# ----------------------------------------------------------------
