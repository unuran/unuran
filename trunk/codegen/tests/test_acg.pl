#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

$DEBUG = 0;

# ----------------------------------------------------------------

require "../read_PDF.pl";

# ----------------------------------------------------------------

my $ACG = "../acg";

# ----------------------------------------------------------------
# Constants

$sample_size = 10;
$accuracy = 1.0e-7;

# ----------------------------------------------------------------
# Prefix for file names

$file_prefix = "./run_test_acg";

# ----------------------------------------------------------------
# Global log file

open LOG, ">$file_prefix.log" or die "Cannot open log file $file_prefix.log";

# ----------------------------------------------------------------
# Compiler

$GCC = "gcc -Wall -ansi -pedantic -I../../src -L../../src";
$G77 = "g77 -Wall";
$JAVAC = "javac -w1";

# ----------------------------------------------------------------

# Read configuration file name for tests from argument list ...
my $test_conf_file = shift
    or die "no argument given";

# C file for making code generator tests
my $make_test_codegen = "make_test_codegen.c";

# C file for tests
my $test_codegen = "test_codegen.c";

# ----------------------------------------------------------------
# List of distributions
my $DISTR = read_PDFdata('../..');

# For description of data fields in this list see file `read_PDF.pl'.

# ................................................................

# Print result on screen
    if ($DEBUG) {
	print_data($DISTR);
    }

# ................................................................

# Print Test data
print_log("sample size = $sample_size\n");
print_log("accuracy = $accuracy\n");
print_log("languages = C, FORTRAN\n\n");

# ----------------------------------------------------------------
# Get list of distributions 

require $test_conf_file; 

# ----------------------------------------------------------------
# Check for missing CONTinuous distributions

my %distr_names;
foreach my $distr (@distr_list) {
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    $distr_names{$1} = 1;
}

foreach my $d (sort keys %{$DISTR}) {
    next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
    unless ($distr_names{$d}) {
	print_log("test missing for distribution \"$d\"\n\n");
    }
}

# ----------------------------------------------------------------
# Fill in parameters for PDF at random

foreach my $distr (@distr_list) {
    while ($distr =~ /([\d\+\-\.]+)\s*\.\.\s*([\d\+\-\.]+)/) {
	my $param = ($1>0) ? exp(log($1)+rand()*(log($2)-log($1))) : $1+rand()*($2-$1);
	$distr =~ s/([\d\+\-\.]+)\s*\.\.\s*([\d\+\-\.]+)/$param/;
    }
}

# ----------------------------------------------------------------
# Process test for each distribution

my $test_nr = 0;
my $errorcode = 0;

foreach my $distr (@distr_list) {
    foreach my $gen (@gen_list) {
	$errorcode += run_test(++$test_nr,"$distr $gen");
    }
}

# ----------------------------------------------------------------
# End

if ($errorcode) {
    $errorcode = 1;
    print_log("\n\tTEST(S) FAILED\n");
}
else {
    print_log("\n\tALL TESTS PASSED\n");
}
close LOG;

exit $errorcode;

# ****************************************************************

# ****************************************************************
#
# Run test
#
# ****************************************************************

# ----------------------------------------------------------------
# Run test for a particular distribution

sub run_test
{
    $test = $_[0];
    $distr = $_[1];

    # Print on screen
    my $test_key = sprintf "%03d", $test;
    print_log("[$test_key] $distr");

    # Remove commas
    $distr =~ s/\,/ /g;

    # Seed for uniform rng
    my $seed = int(rand 12345678) + 1;

    # Files
    my $file_name = "$file_prefix\_$test_key";

    my $UNURAN_exec = "$file_name\_UNURAN";
    my $UNURAN_src = "$UNURAN_exec.c";
    my $UNURAN_log = "$UNURAN_exec.log";

    my $C_exec = "$file_name\_C";
    my $C_src = "$C_exec.c";
    my $C_log = "$C_exec.log";

    my $FORTRAN_exec = "$file_name\_FORTRAN";
    my $FORTRAN_src = "$FORTRAN_exec.f";
    my $FORTRAN_log = "$FORTRAN_exec.log";

    my $JAVA_exec = "$file_name\_JAVA";
    my $JAVA_src = "$JAVA_exec.java";
    my $JAVA_log = "$JAVA_exec.log";

    # Get random variate generators

    # UNURAN version
    my $UNURAN_code = make_UNURAN_code($UNURAN_log,$distr,$seed);
    unless ($UNURAN_code) {
	print_log("  .........  cannot create generator.\n");
	next;
    }
    print_log("\n");
    make_UNURAN_exec($UNURAN_code,$UNURAN_src,$UNURAN_exec);

    # C version
    my $C_code = make_C_code($C_log,$distr,$seed);
    make_C_exec($C_code,$C_src,$C_exec);

    # FORTRAN version
    my $FORTRAN_code = make_FORTRAN_code($FORTRAN_log,$distr,$seed);
    make_FORTRAN_exec($FORTRAN_code,$FORTRAN_src,$FORTRAN_exec);

    # JAVA version
    make_JAVA_gen($JAVA_log,$distr,$seed);
    make_JAVA_exec($JAVA_exec,$distr,$seed);

    # Start generators
    open UNURAN, "$UNURAN_exec |" or die "cannot run $UNURAN_exec"; 
    open C, "$C_exec |" or die "cannot run $C_exec"; 
    open FORTRAN, "$FORTRAN_exec |" or die "cannot run $FORTRAN_exec"; 

    # Run generatores and compare output
    my $C_n_diffs = 0;
    my $FORTRAN_n_diffs = 0;
    my $n_sample = 0;

    while (my $UNURAN_out = <UNURAN>) {
	my $C_out = <C>;
	my $FORTRAN_out = <FORTRAN>;
	
	chomp $UNURAN_out;
	chomp $C_out;
	chomp $FORTRAN_out;
	
	++$n_sample;
	
	(my $UNURAN_x,  my $UNURAN_pdfx)  = split /\s+/, $UNURAN_out, 2;
	(my $C_x,       my $C_pdfx)       = split /\s+/, $C_out, 2;
	(my $FORTRAN_x, my $FORTRAN_pdfx) = split /\s+/, $FORTRAN_out, 2;
	
	my $C_x_diff          = $UNURAN_x    - $C_x;
	my $C_pdfx_diff       = $UNURAN_pdfx - $C_pdfx;
	my $FORTRAN_x_diff    = $UNURAN_x    - $FORTRAN_x;
	my $FORTRAN_pdfx_diff = $UNURAN_pdfx - $FORTRAN_pdfx;
	
	if ( !FP_equal($C_x,$UNURAN_x) or !FP_equal($C_pdfx,$UNURAN_pdfx) ) {
	    ++$C_n_diffs;
	    print LOG "  C: x    = $C_x\tdifference = $C_x_diff\n";
	    print LOG "  C: pdfx = $C_pdfx\tdifference = $C_pdfx_diff\n";
	}
	if ( !FP_equal($FORTRAN_x,$UNURAN_x) or !FP_equal($FORTRAN_pdfx,$UNURAN_pdfx) ) {
	    ++$FORTRAN_n_diffs;
	    print LOG "  FORTRAN: x    = $FORTRAN_x\tdifference = $FORTRAN_x_diff\n";
	    print LOG "  FORTRAN: pdfx = $FORTRAN_pdfx\tdifference = $FORTRAN_pdfx_diff\n";
	}
	
    }

    # End
    close UNURAN;
    close C;
    
    my $errorcode = $n_sample ? 0 : 1;
    
    if ($C_n_diffs > 0) {
	print_log("\t...  C Test FAILED\n");
	++$errorcode;
    }
    if ($FORTRAN_n_diffs > 0) {
	print_log("\t...  FORTRAN Test FAILED\n");
	++$errorcode;
    }
    
    if ($errorcode == 0) {
	print_log("\t...  PASSED\n");
    }

    return $errorcode;

} # run_test()

# ----------------------------------------------------------------
# When two floats are equal

sub FP_equal
{
    my $a = $_[0];
    my $b = $_[1];

    if ($a==$b || abs($a-$b) <= ((abs($a)<abs($b)) ? abs($a) : abs($b)) * $accuracy) {
	return 1;
    }
    else {
	return 0;
    }
} # end of FP_equal()

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
# UNURAN version
#
# ****************************************************************

# ----------------------------------------------------------------
# Make generator code for test file (UNURAN version)

sub make_UNURAN_code
{
    my $logfile = $_[0];
    my $distr = $_[1];
    my $seed = $_[2];

    my $acg_query = "$ACG $distr -l UNURAN -L $logfile";

    my $generator = `$acg_query`;

    return "" if $?;

    # Get name of distribution
    my $distr_name;
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    $distr_name = $1;

    my $urng = make_UNURAN_urng($seed);
    my $main = make_UNURAN_main($distr_name,$seed);

    return $urng.$generator.$main;

} # end of make_UNURAN_code()

# ----------------------------------------------------------------
# uniform rng (UNURAN/PRNG version)

sub make_UNURAN_urng
{
    my $seed = $_[0];

    my $code = <<EOS;
    
/* ---------------------------------------------------------------- */
/* UNURAN / PRNG version                                            */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <unuran.h>

/* ---------------------------------------------------------------- */
/* LCG (Linear Congruential Generator) by Park & Miller (1988).     */
/*   x_(n+1) = 16807 * x_n mod 2^31 - 1    (Minimal Standard)       */
/* ---------------------------------------------------------------- */

void useed(unsigned seed)
{
    static struct prng *urng = NULL;

    if (urng == NULL)
	urng = prng_new("LCG(2147483647,16807,0,1)");

    prng_seed(urng,seed);

    /* synchronize with other generators */
    prng_get_next(urng);
	    
    unur_set_default_urng(urng);
}

/* ---------------------------------------------------------------- */

EOS

    return $code;

} # end of make_UNURAN_urng() 

# ----------------------------------------------------------------
# Make main for test file (UNURAN version)

sub make_UNURAN_main
{
    my $distr = $_[0];
    my $seed = $_[1];

    my $code = <<EOS;

#include <stdio.h>
#include <stdlib.h>

int main()
{
    int i;
    double x, fx;

    useed($seed);

    for (i=0; i<$sample_size; i++) {
	x = rand\_$distr();
	fx = pdf\_$distr(x);
	printf("%.17e\\t%.17e\\n",x,fx);
    }

    exit (0);
}

EOS

    return $code;

} # end of make_UNURAN_main()

# ----------------------------------------------------------------
# Make executable from test file (UNURAN version)

sub make_UNURAN_exec
{
    my $code = $_[0];
    my $src = $_[1];
    my $exec = $_[2];

    # make source file
    open SRC, ">$src" or die "cannot open $src for writing";
    print SRC $code;
    close SRC;

    # compile
    system "$GCC -o $exec $src -lunuran -lprng -lm";

} # end of make_UNURAN_exec()

# ----------------------------------------------------------------

# ****************************************************************
#
# C version
#
# ****************************************************************

# ----------------------------------------------------------------
# Make generator code for test file (C version)

sub make_C_code
{
    my $logfile = $_[0];
    my $distr = $_[1];
    my $seed = $_[2];

    my $acg_query = "$ACG $distr -l C -L $logfile";

    my $generator = `$acg_query`;

    return "" if $?;

    # Get name of distribution
    my $distr_name;
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    $distr_name = $1;

    my $urng = make_C_urng($seed);
    my $main = make_C_main($distr_name,$seed);

    return $urng.$generator.$main;

} # end of make_C_code()

# ----------------------------------------------------------------
# uniform rng (C version)

sub make_C_urng
{
    my $seed = $_[0];

    my $code = <<EOS;
    
/* ---------------------------------------------------------------- */
/* C version                                                        */
/* ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- */
/* LCG (Linear Congruential Generator) by Park & Miller (1988).     */
/*   x_(n+1) = 16807 * x_n mod 2^31 - 1    (Minimal Standard)       */
/* ---------------------------------------------------------------- */

#define uniform() urand()

static unsigned int xn = $seed;   /* seed  */

double urand (void)
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

void useed(unsigned int seed)
{
  xn = seed;
}

EOS

    return $code;

} # end of make_C_urng() 

# ----------------------------------------------------------------
# Make main for test file (C version)

sub make_C_main
{
    my $distr = $_[0];
    my $seed = $_[1];

    my $code = <<EOS;

#include <stdio.h>
#include <stdlib.h>

int main()
{
    int i;
    double x, fx;

    useed($seed);

    for (i=0; i<$sample_size; i++) {
	x = rand\_$distr();
	fx = pdf\_$distr(x);
	printf("%.17e\\t%.17e\\n",x,fx);
    }

    exit (0);
}

EOS

    return $code;

} # end of make_C_main()

# ----------------------------------------------------------------
# Make executable from test file (C version)

sub make_C_exec
{
    my $code = $_[0];
    my $src = $_[1];
    my $exec = $_[2];

    # make source file
    open SRC, ">$src" or die "cannot open $src for writing";
    print SRC $code;
    close SRC;

    # compile
    system "$GCC -o $exec $src -lm";

} # end of make_C_exec()

# ----------------------------------------------------------------

# ****************************************************************
#
# FORTRAN version
#
# ****************************************************************

# ----------------------------------------------------------------
# Make generator code for test file (FORTRAN version)

sub make_FORTRAN_code
{
    my $logfile = $_[0];
    my $distr = $_[1];
    my $seed = $_[2];

    my $acg_query = "$ACG $distr -l FORTRAN -L $logfile";

    my $generator = `$acg_query`;

    return "" if $?;

    # remove built-in uniform rng
    $urng_pattern = "\\s+DOUBLE PRECISION FUNCTION urand\\s*\\(\\s*\\).*?END\\s*";
    $generator =~ s/($urng_pattern)/\n/s;

    # Get name of distribution
    my $distr_name;
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    $distr_name = $1;

    my $urng = make_FORTRAN_urng($seed);
    my $main = make_FORTRAN_main($distr_name,$seed);

    return $urng.$generator.$main;

} # end of make_FORTRAN_code()

# ----------------------------------------------------------------
# uniform rng (FORTRAN version)

sub make_FORTRAN_urng
{
    my $seed = $_[0];

    my $code = <<EOS;
    
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

EOS

    return $code;

} # end of make_FORTRAN_urng() 

# ----------------------------------------------------------------
# Make main for test file (FORTRAN version)

sub make_FORTRAN_main
{
    my $distr = $_[0];
    my $seed = $_[1];

    my $rand_name = substr "r$distr", 0, 6;
    my $pdf_name = substr "f$distr", 0, 6;

    my $code = <<EOS;

* ------------------------------------------------------------------ *

      PROGRAM MAIN

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      CALL useed($seed)

      DO 1 i=1,$sample_size
         x = $rand_name()
	 fx = $pdf_name(x)
         WRITE (*,'(d24.18, 1x, d24.18)') x, fx
 1    CONTINUE
      END

* ------------------------------------------------------------------ *

EOS

    return $code;

} # end of make_FORTRAN_main()

# ----------------------------------------------------------------
# Make executable from test file (FORTRAN version)

sub make_FORTRAN_exec
{
    my $code = $_[0];
    my $src = $_[1];
    my $exec = $_[2];

    # make source file
    open SRC, ">$src" or die "cannot open $src for writing";
    print SRC $code;
    close SRC;

    # compile
    system "$G77 -o $exec $src";

} # end of make_FORTRAN_exec()

# ----------------------------------------------------------------

# ****************************************************************
#
# Java version
#
# ****************************************************************

# ----------------------------------------------------------------
# uniform rng (Java version)

#
# run test_urng.pl to create generator
#

# ----------------------------------------------------------------
# Make executable for Generator (Java version)

sub make_JAVA_gen
{
    my $logfile = $_[0];
    my $distr = $_[1];
    my $seed = $_[2];

    my $acg_query = "$ACG $distr -l JAVA -L $logfile";

    my $generator = `$acg_query`;

    return "" if $?;

    # Replace Java URNG
    $generator =~ s/random\s*\(\s*\)\s*\;/Urand.random\(\)\;/g;

    # Get name of Generator class
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    my $distr_name = $1;
    my $gen_src = "Generator_$distr_name.java";

    # make source file
    open SRC, ">$gen_src" or die "cannot open $gen_src for writing";
    print SRC $generator;
    close SRC;

    # compile
    system "$JAVAC $gen_src";

} # end of make_JAVA_gen()

# ----------------------------------------------------------------
# Make executable for Generator test (Java version)

sub make_JAVA_exec
{
    my $test = $_[0];
    my $distr = $_[1];
    my $seed = $_[2];

    my $src = "$test.java";

    # remove leading "./"
    $test =~ s/\.\///;

    # Get name of Generator class
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    my $distr_name = $1;
    my $gen = "Generator_$distr_name";

    # Make code
    my $code = <<EOS;

/* ---------------------------------------------------------------- */

public class $test {
   public static void main(String[] args) throws Exception {

      private double x;
      private double pdfx;

      /* set seed */
      Urand.useed($seed);

      for (private int i = 0; i<$sample_size;i++) {
	x = $gen.sample();
	pdfx = $gen.pdf(x);
        System.out.println( x );
      }
   }

}  /* end of class Test */

EOS

    # make source file
    open SRC, ">$src" or die "cannot open $src for writing";
    print SRC $code;
    close SRC;

    # compile
    system "$JAVAC $src";

} # end of make_JAVA_exec()

# ----------------------------------------------------------------

