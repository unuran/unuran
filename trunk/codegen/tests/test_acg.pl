#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

$DEBUG = 0;
$| = 1;

# ----------------------------------------------------------------

require "../read_PDF.pl";

# ----------------------------------------------------------------

my $ACG = "../acg";

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
$JAVAC = "javac";

# ----------------------------------------------------------------

# Read configuration file name for tests from argument list ...
my $test_conf_file = shift
    or die "no argument given";

# C file for making code generator tests
my $make_test_codegen = "make_test_codegen.c";

# C file for tests
my $test_codegen = "test_codegen.c";

# ----------------------------------------------------------------
# Read data from test config file

print_log("Read data ...\n\n");
require $test_conf_file; 

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
print_log("languages = C, FORTRAN, JAVA\n\n");

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
# File names

print_log("Make sources for tests (".(($#distr_list+1)*($#method_list+1)).") ...\n\n");

my $UNURAN_exec = "$file_prefix\_UNURAN";
my $UNURAN_src = "$UNURAN_exec.c";

my $C_exec = "$file_prefix\_C";
my $C_src = "$C_exec.c";

my $FORTRAN_exec = "$file_prefix\_FORTRAN";
my $FORTRAN_src = "$FORTRAN_exec.f";

my $JAVA_urand_src = "./Urand.java";
my $JAVA_exec = "$file_prefix\_JAVA";
my $JAVA_src = "$JAVA_exec.java";
my $JAVA_class = $JAVA_exec;
$JAVA_class =~ s/^\.\///;

# ----------------------------------------------------------------
# Get code for Uniform random number generator

my $UNURAN_urng  = make_UNURAN_urng();
my $C_urng       = make_C_urng();
my $FORTRAN_urng = make_FORTRAN_urng();
make_JAVA_urng();

# ----------------------------------------------------------------
# Make Generators for all tests

my $UNURAN_gen;
my $UNURAN_main;
my $C_gen;
my $C_main;

my $test_nr = 0;
$test_runs = 0;

foreach my $distr (@distr_list) {
    foreach my $method (@method_list) {
	
	# Increment counter
	++$test_nr;
	
	# Make a 3 digit string from test number
	my $test_key = sprintf "%03d", $test_nr;
	print_log("[$test_key]");
	
	# Seed for uniform rng
	my $seed = int(rand 12345678) + 1;
	
	# Remove commas from ./acg command line
	$distr =~ s/\,/ /g;
	
	# Get random variate generators
	
	# UNURAN version
	my $test_name = make_UNURAN_gen(\$UNURAN_gen,"$distr $method",$test_key,$seed);
	$UNURAN_main .= "\t$test_name ();\n";
	
	# C version
	$test_name = make_C_gen(\$C_gen,"$distr $method",$test_key,$seed);
	$C_main .= "\t$test_name ();\n";
	
	# FORTRAN version
	$test_name = make_FORTRAN_gen(\$FORTRAN_gen,"$distr $method",$test_key,$seed);
	$FORTRAN_main .= "      CALL $test_name\n";
	
	# JAVA version
	$JAVA_main .= make_JAVA_gen(\$JAVA_gen,"$distr $method",$test_key,$seed);
    }
}

# ----------------------------------------------------------------
# Print number of tests

my $test_runs_rel = 100.0 * $test_runs / $test_nr;
print_log("\n\nFound $test_nr tests.  Build $test_runs tests ($test_runs_rel \%) ...\n\n");

# ----------------------------------------------------------------
# Make source files

open UNURAN, ">$UNURAN_src" or die "Cannot open file $UNURAN_src";
print UNURAN $UNURAN_urng;
print UNURAN $UNURAN_gen;
print UNURAN make_UNURAN_main($UNURAN_main);
close UNURAN;

open C, ">$C_src" or die "Cannot open file $C_src";
print C $C_urng;
print C $C_gen;
print C make_C_main($C_main);
close C;

open FORTRAN, ">$FORTRAN_src" or die "Cannot open file $FORTRAN_src";
print FORTRAN $FORTRAN_urng;
print FORTRAN $FORTRAN_gen;
print FORTRAN make_FORTRAN_main($FORTRAN_main);
close FORTRAN;

open JAVA, ">$JAVA_src" or die "Cannot open file $JAVA_src";
print JAVA make_JAVA_main($JAVA_main);
close JAVA;

# ----------------------------------------------------------------
# Compile sources

print_log("Compiling sources ...\n\n");

system "$GCC -o $UNURAN_exec $UNURAN_src -lunuran -lprng -lm";
system "$GCC -o $C_exec $C_src -lm";
system "$G77 -o $FORTRAN_exec $FORTRAN_src -lm";
system "$JAVAC $JAVA_src";

# ----------------------------------------------------------------
# Run tests

print_log("Run tests ...\n\n");

# Start generators
open UNURAN, "$UNURAN_exec |";
open C, "$C_exec |";
open FORTRAN, "$FORTRAN_exec |";

open JAVA, "java $JAVA_class |";
$HAVE_JAVA = ($?) ? 0 : 1;
unless ($HAVE_JAVA) {
    print_log("Cannot run JAVA tests!\n\n");
}

# Run generatores and compare output
my $data_mode = 0;

my $C_errors = 0;
my $FORTRAN_errors = 0;
my $JAVA_errors = 0;

my $n_sample = 0;
my $C_n_diffs = 0;
my $FORTRAN_n_diffs = 0;
my $JAVA_n_diffs = 0;

while (my $UNURAN_out = <UNURAN>) {
    my $C_out = <C>;
    my $FORTRAN_out = <FORTRAN>;
    my $JAVA_out = <JAVA>;
    
    chomp $UNURAN_out;
    chomp $C_out;
    chomp $FORTRAN_out;
    chomp $JAVA_out;

    # print section headers
    if (not $data_mode) {
	if ($UNURAN_out eq "start") {
	    $n_sample = 0;
	    $C_n_diffs = 0;
 	    $FORTRAN_n_diffs = 0;
	    $JAVA_n_diffs = 0;
	    
	    $data_mode = 1;
	    $errors = 0;
	}
	else {
	    print_log("$UNURAN_out\n");
	}
	next;
    }
    
    # no more data, compute results
    if ($UNURAN_out eq "stop") {
	$data_mode = 0;
	my $errors = 0;
	
	if ($C_n_diffs > 0) {
	    my $quote = 100.0 * $C_n_diffs / $n_sample;
	    print_log("\t...  C Test FAILED  ($quote \%)\n");
	    ++$C_errors;
	    ++$errors;
	}
	if ($FORTRAN_n_diffs > 0) {
	    my $quote = 100.0 * $FORTRAN_n_diffs / $n_sample;
	    print_log("\t...  FORTRAN Test FAILED  ($quote \%)\n");
	    ++$FORTRAN_errors;
	    ++$errors;
	}
	if ($HAVE_JAVA and $JAVA_n_diffs > 0) {
	    my $quote = 100.0 * $JAVA_n_diffs / $n_sample;
	    print_log("\t...  JAVA Test FAILED  ($quote \%)\n");
	    ++$JAVA_errors;
	    ++$errors;
	}
	
	unless ($errors) {
	    print_log("\t...  PASSED\n") if $n_sample;
	}

	next;
    }
    
    # reading data
    ++$n_sample;
    
    (my $UNURAN_x,  my $UNURAN_pdfx)  = split /\s+/, $UNURAN_out, 2;
    (my $C_x,       my $C_pdfx)       = split /\s+/, $C_out, 2;
    (my $FORTRAN_x, my $FORTRAN_pdfx) = split /\s+/, $FORTRAN_out, 2;
    (my $JAVA_x,    my $JAVA_pdfx)    = split /\s+/, $JAVA_out, 2;
    
    my $C_x_diff          = $UNURAN_x    - $C_x;
    my $C_pdfx_diff       = $UNURAN_pdfx - $C_pdfx;
    my $FORTRAN_x_diff    = $UNURAN_x    - $FORTRAN_x;
    my $FORTRAN_pdfx_diff = $UNURAN_pdfx - $FORTRAN_pdfx;
    my $JAVA_x_diff       = $UNURAN_x    - $JAVA_x;
    my $JAVA_pdfx_diff    = $UNURAN_pdfx - $JAVA_pdfx;
    
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
    if ( $HAVE_JAVA and (!FP_equal($JAVA_x,$UNURAN_x) or !FP_equal($JAVA_pdfx,$UNURAN_pdfx)) ) {
	++$JAVA_n_diffs;
	print LOG "  JAVA: x    = $JAVA_x\tdifference = $JAVA_x_diff\n";
	print LOG "  JAVA: pdfx = $JAVA_pdfx\tdifference = $JAVA_pdfx_diff\n";
    }
    
}

# End
close UNURAN;
close C;
close FORTRAN;
close JAVA;

# ----------------------------------------------------------------
# End

my $errors = 0;
if ($C_errors) {
    print_log("\n\tC TEST(S) FAILED ($C_errors)\n");
    ++$errors;
}
if ($FORTRAN_errors) {
    print_log("\n\tFORTRAN TEST(S) FAILED ($FORTRAN_errors)\n");
    ++$errors;
}
if ($JAVA_errors) {
    print_log("\n\tJAVA TEST(S) FAILED ($JAVA_errors)\n");
    ++$errors;
}

unless ($HAVE_JAVA) {
    print_log("\nCannot run JAVA tests!\n\n");
}

unless ($errors) {
    print_log("\n\tALL TESTS PASSED\n");
}

exit ($errors ? 1 : 0);

# ----------------------------------------------------------------
# End


# ****************************************************************

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
# Uniform random number generators
#
# ****************************************************************

# ----------------------------------------------------------------
# uniform rng (UNURAN/PRNG version)

sub make_UNURAN_urng
{
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
# uniform rng (C version)

sub make_C_urng
{
    my $code = <<EOS;
    
/* ---------------------------------------------------------------- */
/* C version                                                        */
/* ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- */
/* LCG (Linear Congruential Generator) by Park & Miller (1988).     */
/*   x_(n+1) = 16807 * x_n mod 2^31 - 1    (Minimal Standard)       */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>

#define uniform() urand()

static unsigned int xn = 1;   /* seed  */

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

/* ---------------------------------------------------------------- */

EOS

    return $code;

} # end of make_C_urng() 

# ----------------------------------------------------------------
# uniform rng (FORTRAN version)

sub make_FORTRAN_urng
{
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
      DATA xn/1/
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

EOS

    return $code;

} # end of make_FORTRAN_urng() 

# ----------------------------------------------------------------
# uniform rng (JAVA version)

sub make_JAVA_urng
{
    # Source
    my $code = <<EOS;
    
/* ---------------------------------------------------------------- */
/* JAVA version                                                     */
/* ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- */
/* LCG (Linear Congruential Generator) by Park & Miller (1988).     */
/*   x_(n+1) = 16807 * x_n mod 2^31 - 1    (Minimal Standard)       */
/* ---------------------------------------------------------------- */

public class Urand {

    /* class variables */

    private static int xn = 1;  /* seed */
    private static final int a = 16807;
    private static final int m = 2147483647;
    private static final int q = 127773;      /* m / a */
    private static final int r = 2836;        /* m % a */

    /* member functions */
    public static double random() 
    {
	int hi, lo, test;

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

    # Make source file
    open JAVA_urng, ">$JAVA_urand_src" or die "Cannot open file $JAVA_urand_src";
    print JAVA_urng $code;
    close JAVA_urng;

    # Compile
    system "$JAVAC $JAVA_urand_src";

} # end of make_JAVA_urng() 


# ****************************************************************
#
# Make generator code
#
# ****************************************************************

# ----------------------------------------------------------------
# Make generator code for test file (UNURAN version)

sub make_UNURAN_gen
{
    my $UNURAN_code = $_[0];
    my $distr = $_[1];
    my $test_key = $_[2];
    my $seed = $_[3];
    
    # Programming language
    my $l = "UNURAN";

    # Log file
    my $Log = "$file_prefix.log.$test_key.$l";

    # Get name of distribution
    my $distr_name;
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    $distr_name = $1."\_$test_key";

    # Get code for generator
    my $gencode = `$ACG $distr -N $distr_name -l $l -L $Log`;
    my $failed = $?;

    # Name of test routine
    my $test_name = "test\_$distr_name";

    # Name of generator and pdf
    my $gen_name = "rand\_$distr_name";
    my $pdf_name = "pdf\_$distr_name";

    # We have to mask "
    $distr =~ s/\"/\\\"/g;

    # Make code for test routine
    my $testcode = "int $test_name (void)\n{\n";
    
    if ($failed) {
	$testcode .= "\tprintf(\"[$test_key] $distr  .........  cannot create generator.\\n\");\n";
	$testcode .= "\tprintf(\"start\\nstop\\n\");\n";
	$testcode .= "\treturn 0;\n}\n\n";
    }
    else {
	++$test_runs;
	$testcode .= <<EOS;
\tint i;
\tdouble x, fx;

\tuseed($seed);

\tprintf("[$test_key] $distr\\n");
\tprintf("start\\n");

\tfor (i=0; i<$sample_size; i++) {
\t\tx = $gen_name ();
\t\tfx = $pdf_name (x);
\t\tprintf("%.17e  %.17e\\n",x,fx);
\t}

\tprintf("stop\\n");
\treturn 1;
\}

/* ---------------------------------------------------------------- */

EOS
}
    # Store code
    $$UNURAN_code .= $gencode;
    $$UNURAN_code .= $testcode;

    # return name of test routine
    return $test_name;

} # end of make_UNURAN_gen()


# ----------------------------------------------------------------
# Make generator code for test file (C version)

sub make_C_gen
{
    my $C_code = $_[0];
    my $distr = $_[1];
    my $test_key = $_[2];
    my $seed = $_[3];
    
    # Programming language
    my $l = "C";

    # Log file
    my $Log = "$file_prefix.log.$test_key.$l";

    # Get name of distribution
    my $distr_name;
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    $distr_name = $1."\_$test_key";

    # Get code for generator
    my $gencode = `$ACG $distr -N $distr_name -l $l -L $Log`;
    my $failed = $?;

    # Name of test routine
    my $test_name = "test\_$distr_name";

    # Name of generator and pdf
    my $gen_name = "rand\_$distr_name";
    my $pdf_name = "pdf\_$distr_name";

    # We have to mask "
    $distr =~ s/\"/\\\"/g;

    # Make code for test routine
    my $testcode = "int $test_name (void)\n{\n";
    
    if ($failed) {
	$testcode .= "\tprintf(\"[$test_key] $distr  .........  cannot create generator.\\n\");\n";
	$testcode .= "\tprintf(\"start\\nstop\\n\");\n";
	$testcode .= "\treturn 0;\n}\n\n";
    }
    else {
	$testcode .= <<EOS;
\tint i;
\tdouble x, fx;

\tuseed($seed);

\tprintf("[$test_key] $distr\\n");
\tprintf("start\\n");

\tfor (i=0; i<$sample_size; i++) {
\t\tx = $gen_name ();
\t\tfx = $pdf_name (x);
\t\tprintf("%.17e  %.17e\\n",x,fx);
\t}

\tprintf("stop\\n");
\treturn 1;
\}

/* ---------------------------------------------------------------- */

EOS

}

    # Store code
    $$C_code .= $gencode;
    $$C_code .= $testcode;

    # return name of test routine
    return $test_name;

} # end of make_C_gen()


# ----------------------------------------------------------------
# Make generator code for test file (FORTRAN version)

sub make_FORTRAN_gen
{
    my $FORTRAN_code = $_[0];
    my $distr = $_[1];
    my $test_key = $_[2];
    my $seed = $_[3];
    
    # Programming language
    my $l = "FORTRAN";
    
    # Log file
    my $Log = "$file_prefix.log.$test_key.$l";
    
    # Get name of distribution
    my $distr_name;
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    $distr_name = substr($1,0,1)."$test_key";
    
    # Get code for generator
    my $gencode = `$ACG $distr -N $distr_name -l $l -L $Log`;
    my $failed = $?;
    
    # Remove built-in uniform rng
    my $urng_pattern = "\\s+DOUBLE PRECISION FUNCTION urand\\s*\\(\\s*\\).*?END.*?END\\s*";
    $gencode =~ s/($urng_pattern)/\n/s;
    
    # Name of test routine
    my $test_name = "t$test_key";

    # Name of generator and pdf
    my $gen_name = "r$distr_name";
    my $pdf_name = "f$distr_name";
    
    # Make code for test routine
    my $testcode = "      SUBROUTINE $test_name\n\n";
    
    if ($failed) {
	$testcode .= "      WRITE (*,'(''[$test_key] '',\n";
	foreach $p (split (/\s+/, $distr)) {
	    $testcode .= "     #   ''$p '',\n";
        }
	$testcode .= "     #   '' .........  cannot create generator.'')')\n";
	$testcode .= "      WRITE (*,'(''start'')')\n";
	$testcode .= "      WRITE (*,'(''stop'')')\n";
	$testcode .= "      RETURN\n";
	$testcode .= "      END\n";
    }
    else {
	$testcode .= "      IMPLICIT DOUBLE PRECISION (a-h,o-z)\n\n";
	$testcode .= "      WRITE (*,'(''[$test_key] '',\n";
	foreach $p (split /\s+/, $distr) {
	    $testcode .= "     #   ''$p '',\n";
	}
	$testcode .= "     #   '''')')\n";
	$testcode .= "      WRITE (*,'(''start'')')\n";

        $testcode .= "      CALL useed($seed) \n";

	$testcode .= "      DO 1 i=1,$sample_size \n";
	$testcode .= "         x = $gen_name() \n";
	$testcode .= "         fx = $pdf_name(x) \n";
	$testcode .= "         WRITE (*,'(d24.18, 1x, d24.18)') x, fx \n";
	$testcode .= "1     CONTINUE \n\n";

	$testcode .= "      WRITE (*,'(''stop'')')\n";
	$testcode .= "      RETURN\n";
	$testcode .= "      END\n\n";
	$testcode .= "* ------------------------------------------------------------------ *\n";
    }

    # Store code
    $$FORTRAN_code .= $gencode;
    $$FORTRAN_code .= $testcode;

    # return name of test routine
    return $test_name;

} # end of make_FORTRAN_gen()


# ----------------------------------------------------------------
# Make generator code for test file (JAVA version)

sub make_JAVA_gen
{
    my $JAVA_code = $_[0];
    my $distr = $_[1];
    my $test_key = $_[2];
    my $seed = $_[3];
    
    # Programming language
    my $l = "JAVA";

    # Log file
    my $Log = "$file_prefix.log.$test_key.$l";

    # Get name of distribution
    my $distr_name;
    die unless "$distr " =~ /-d\s+(\w+)\s+/;
    $distr_name = $1."\_$test_key";

    # Get code for generator
    my $gencode = `$ACG $distr -N $distr_name -l $l -L $Log`;
    my $failed = $?;

    # We have to mask "
    $distr =~ s/\"/\\\"/g;

    # We have nothing to do when ACG did not work
    if ($failed) {
	return <<EOS;
\t\tSystem.out.println("[$test_key] $distr  .........  cannot create generator.");
\t\tSystem.out.println("start");
\t\tSystem.out.println("stop");

EOS
    }

    # Replace Java URNG
    $gencode =~ s/random\s*\(\s*\)\s*\;/Urand.random\(\)\;/g;

    # Make class for generator
    my $gen_src = "Generator_$distr_name.java";
    open SRC, ">$gen_src" or die "cannot open $gen_src for writing";
    print SRC $gencode;
    close SRC;

    # Compile
##    system "$JAVAC $gen_src";

    # Name of generator and pdf
    my $gen_name = "Generator\_$distr_name";

    # Now make tests
    return <<EOS;
\t\tSystem.out.println("[$test_key] $distr");
\t\tUrand.useed($seed);
\t\tSystem.out.println("start");
\t\tfor (int i = 0; i<$sample_size;i++) {
\t\t\tx = $gen_name.sample();
\t\t\tpdfx = $gen_name.pdf(x);
\t\t\tSystem.out.println( x +" "+ pdfx);
\t\t}
\t\tSystem.out.println("stop");

/* ---------------------------------------------------------------- */

EOS

} # end of make_JAVA_gen()


# ****************************************************************
#
# Make Main
#
# ****************************************************************

# ----------------------------------------------------------------
# Make code for main of test file (UNURAN version)

sub make_UNURAN_main
{
    my $body = $_[0];

    return "int main()\n\{\n$body\n\texit (0);\n\}\n";

} # end of make_UNURAN_main()

# ----------------------------------------------------------------
# Make code for main of test file (C version)

sub make_C_main
{
    my $body = $_[0];

    return "int main()\n\{\n$body\n\texit (0);\n\}\n";

} # end of make_C_main()

# ----------------------------------------------------------------
# Make code for main of test file (FORTRAN version)

sub make_FORTRAN_main
{
    my $body = $_[0];

    return <<EOS;

* ------------------------------------------------------------------ *

      PROGRAM MAIN

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

$body

      END

* ------------------------------------------------------------------ *

EOS

} # end of make_FORTRAN_main()


# ----------------------------------------------------------------
# Make code for main of test file (JAVA version)

sub make_JAVA_main
{
    my $body = $_[0];

    my $class = $JAVA_exec;
    $class =~ s/^\.\///;

    # Make code
    return <<EOS;

/* ---------------------------------------------------------------- */

public class $class {
\tpublic static void main(String[] args) throws Exception {

\t\tdouble x;
\t\tdouble pdfx;

$body

\t}

}  /* end of class Test */

EOS

} # end of make_JAVA_main()

# ----------------------------------------------------------------

