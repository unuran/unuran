#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

$DEBUG = 0;

# ----------------------------------------------------------------

require "../read_PDF.pl";
require "read_test_conf.pl";

# ----------------------------------------------------------------

my $ACG = "../acg";

# ----------------------------------------------------------------
# Constants

$sample_size = 10000;
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
print "sample size = $sample_size\n";
print LOG "sample size = $sample_size\n";

print "accuracy = $accuracy\n";
print LOG "accuracy = $accuracy\n";

print "languages = C, FORTRAN\n\n";
print LOG "languages = C, FORTRAN\n\n";

# ----------------------------------------------------------------
# Get list of distributions 

my $list_distr = get_test_distributions( $test_conf_file, $DISTR );

#.................................................................
# Check for missing CONTinuous distributions

# list of names of distributions without distribution number
my $list_distr_short;
foreach my $d (sort keys %{$list_distr}) {
    my $distr_short = $d;
    $distr_short =~ s/\_[^\_]+$//;
    $list_distr_short->{$distr_short} = $d;
}

foreach my $d (sort keys %{$DISTR}) {
    next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
    unless ($list_distr_short->{$d}) {
	print STDERR "test missing for distribution \"$d\"\n";
    }
}

# ----------------------------------------------------------------
# Process test for each distribution

my $test_nr = 0;

foreach my $d (sort keys %{$list_distr}) {
    # number of test
    ++$test_nr;
    my $test_key = sprintf "%03d", $test_nr;

    # get name of distribution
    my $distr_key = $d;
    my $distr = $d;
    $distr =~ s/\_(.*)$//;

    # get parameter list
    my $fparam = $list_distr->{$distr_key};

    # print on screen
    print "[$test_key] $distr($fparam)";
    print LOG "[$test_key] $distr($fparam)";

    # seed for uniform rng
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
    my $UNURAN_code = make_UNURAN_code($UNURAN_log,$distr,$fparam,$seed);
    unless ($UNURAN_code) {
	print "  .........  cannot create generator.\n";
	print LOG "  .........  cannot create generator.\n";
	next;
    }
    make_UNURAN_exec($UNURAN_code,$UNURAN_src,$UNURAN_exec);

    # C version
    my $C_code = make_C_code($C_log,$distr,$fparam,$seed);
    make_C_exec($C_code,$C_src,$C_exec);

    # FORTRAN version
    my $FORTRAN_code = make_FORTRAN_code($FORTRAN_log,$distr,$fparam,$seed);
    make_FORTRAN_exec($FORTRAN_code,$FORTRAN_src,$FORTRAN_exec);

    # JAVA version
    my $JAVA_code = make_JAVA_code($JAVA_log,$distr,$fparam,$seed);
    make_JAVA_exec($JAVA_code,$JAVA_src,$JAVA_exec);

    # Start generators
    open UNURAN, "$UNURAN_exec |" or die "cannot run $UNURAN_exec"; 
    open C, "$C_exec |" or die "cannot run $C_exec"; 
    open FORTRAN, "$FORTRAN_exec |" or die "cannot run $FORTRAN_exec"; 

    # Run generatores and compare output
    $C_n_diffs = 0;
    $FORTRAN_n_diffs = 0;
    $n_sample = 0;

    while ($UNURAN_out = <UNURAN>) {
	$C_out = <C>;
	$FORTRAN_out = <FORTRAN>;
	
	chomp $UNURAN_out;
	chomp $C_out;
	chomp $FORTRAN_out;
	
	++$n_sample;
	
	($UNURAN_x, $UNURAN_pdfx) = split /\s+/, $UNURAN_out, 2;
	($C_x, $C_pdfx) = split /\s+/, $C_out, 2;
	($FORTRAN_x, $FORTRAN_pdfx) = split /\s+/, $FORTRAN_out, 2;
	
	$C_x_diff = abs($UNURAN_x - $C_x);
	$C_pdfx_diff = abs($UNURAN_pdfx - $C_pdfx);
	$FORTRAN_x_diff = abs($UNURAN_x - $FORTRAN_x);
	$FORTRAN_pdfx_diff = abs($UNURAN_pdfx - $FORTRAN_pdfx);
	
	if ( !FP_equal($C_x,$UNURAN_x) or !FP_equal($C_pdfx,$UNURAN_pdfx) ) {
	    ++$C_n_diffs;
	    print LOG "\n  C: x    = $C_x\tdifference = $C_x_diff\n";
	    print LOG "  C: pdfx = $C_pdfx\tdifference = $C_pdfx_diff";
	}
	if ( !FP_equal($FORTRAN_x,$UNURAN_x) or !FP_equal($FORTRAN_pdfx,$UNURAN_pdfx) ) {
	    ++$FORTRAN_n_diffs;
	    print LOG "\n  FORTRAN: x    = $FORTRAN_x\tdifference = $FORTRAN_x_diff\n";
	    print LOG "  FORTRAN: pdfx = $FORTRAN_pdfx\tdifference = $FORTRAN_pdfx_diff";
	}
	
    }
    
    # End
    close UNURAN;
    close C;
    
    $errorcode = $n_sample ? 0 : 1;
    
    if ($C_n_diffs > 0) {
	print "\n\t ...  C Test FAILED\n";
	print LOG "\n\t ...  C Test FAILED\n";
	++$errorcode;
    }
    if ($FORTRAN_n_diffs > 0) {
	print "\n\t ...  FORTRAN Test FAILED\n";
	print LOG "\n\t ...  FORTRAN Test FAILED\n";
	++$errorcode;
    }
    
    if ($errorcode == 0) {
	print "  ...  PASSED\n";
	print LOG "  ...  PASSED\n";
    }
}

# ----------------------------------------------------------------
# End

close LOG;

exit 0;

# ****************************************************************

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
    my $fparam = $_[2];
    my $seed = $_[3];

    my $acg_query = "$ACG -l UNURAN -d $distr -L $logfile";
    $acg_query .= " -p \"$fparam\"" if $fparam; 

    my $generator = `$acg_query`;

    return "" if $?;

    my $urng = make_UNURAN_urng($seed);
    my $main = make_UNURAN_main($distr,$seed);

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

    my $code = <<EOS

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
    my $fparam = $_[2];
    my $seed = $_[3];

    my $acg_query = "$ACG -l C -d $distr -L $logfile";
    $acg_query .= " -p \"$fparam\"" if $fparam; 

    my $generator = `$acg_query`;

    return "" if $?;

    my $urng = make_C_urng($seed);
    my $main = make_C_main($distr,$seed);

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

    my $code = <<EOS

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
    my $fparam = $_[2];
    my $seed = $_[3];

    my $acg_query = "$ACG -l FORTRAN -d $distr -L $logfile";
    $acg_query .= " -p \"$fparam\"" if $fparam; 

    my $generator = `$acg_query`;

    return "" if $?;

    # remove built-in uniform rng
    $urng_pattern = "\\s+DOUBLE PRECISION FUNCTION unif\\s*\\(\\s*\\).*?END\\s*";
    $generator =~ s/($urng_pattern)/\n/s;

    my $urng = make_FORTRAN_urng($seed);
    my $main = make_FORTRAN_main($distr,$seed);

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

      DOUBLE PRECISION FUNCTION unif()

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

      unif = xn * 4.656612875245796924105750827D-10

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

    my $code = <<EOS

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
# Make generator code for test file (C version)

sub make_JAVA_code
{
    my $logfile = $_[0];
    my $distr = $_[1];
    my $fparam = $_[2];
    my $seed = $_[3];

    my $acg_query = "$ACG -l JAVA -d $distr -L $logfile";
    $acg_query .= " -p \"$fparam\"" if $fparam; 

    my $generator = `$acg_query`;

    return "" if $?;

    my $main = "";
##    my $main = make_JAVA_main($distr,$seed);

    return $generator.$main;

} # end of make_JAVA_code()

# ----------------------------------------------------------------
# uniform rng (Java version)

#
# run test_urng.pl to create generator
#

# ----------------------------------------------------------------
# Make main for test file (Java version)

sub make_JAVA_main
{
    my $distr = $_[0];
    my $seed = $_[1];

    my $code = <<EOS

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

} # end of make_JAVA_main()

# ----------------------------------------------------------------
# Make executable from test file (Java version)

sub make_JAVA_exec
{
    my $code = $_[0];
    my $src = $_[1];
    my $exec = $_[2];

    # make source file
    open SRC, ">$src" or die "cannot open $src for writing";
    print SRC $code;
    close SRC;

    # compile
    system "$JAVAC $src";

} # end of make_JAVA_exec()

# ----------------------------------------------------------------

