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

$sample_size = 10;
$accuracy = 1.0e-7;

# ----------------------------------------------------------------
# Compiler

$GCC = "gcc -Wall -ansi -pedantic -I../../src -L../../src";
$G77 = "g77 -Wall";

# ----------------------------------------------------------------

# Read configuration file name for tests from argument list ...
my $test_conf_file = shift
    or die "no argument given";

# C file for making code generator tests
my $make_test_codegen = "make_test_codegen.c";

# C file for tests
my $test_codegen = "test_codegen.c";

# Sample size for test
my $SAMPLE_SIZE = 100000;

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
print "accuracy = $accuracy\n";
print "languages = C\n\n";

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

foreach my $d (sort keys %{$list_distr}) {

    # get name of distribution
    my $distr_key = $d;
    my $distr = $d;
    $distr =~ s/\_(.*)$//;
    my $distr_nr = $1;

    # get parameter list
    my $fparam = $list_distr->{$distr_key};

    # print on screen
    print "[$distr_nr] $distr($fparam)";

    # seed for uniform rng
    $seed = int(rand 12345678) + 1;

    # Files
    $file_name = "./run_test\_acg\_$distr_nr";

    $UNURAN_exec = "$file_name\_UNURAN";
    $UNURAN_src = "$UNURAN_exec.c";

    $C_exec = "$file_name\_C";
    $C_src = "$C_exec.c";

    # Get random variate generators

    # UNURAN version
    my $UNURAN_code = make_UNURAN_code($distr,$fparam,$seed);
    unless ($UNURAN_code) {
	print "  .........  cannot create generator.\n";
	next;
    }
    print "\n";
    make_UNURAN_exec($UNURAN_code,$UNURAN_src,$UNURAN_exec);

    # C version
    my $C_code = make_C_code($distr,$fparam,$seed);
    make_C_exec($C_code,$C_src,$C_exec);

    # Start generators
    open UNURAN, "$UNURAN_exec |" or die "cannot run $UNURAN_exec"; 
    open C, "$C_exec |" or die "cannot run $C_exec"; 

    # Run generatores and compare output
    $C_n_diffs = 0;

    while ($UNURAN_out = <UNURAN>) {
	$C_out = <C>;

	chomp $UNURAN_out;
	chomp $C_out;
	
	($UNURAN_x, $UNURAN_pdfx) = split /\s+/, $UNURAN_out, 2;
	($C_x, $C_pdfx) = split /\s+/, $C_out, 2;

	$C_x_diff = abs($UNURAN_x - $C_x);
	$C_pdfx_diff = abs($UNURAN_pdfx - $C_pdfx);

	if ( !FP_equal($C_x,$UNURAN_x) or !FP_equal($C_pdfx,$UNURAN_pdfx) ) {
	    ++$C_n_diffs;
	    print "x    = $UNURAN_x\tdifference = $C_x_diff\n";
	    print "pdfx = $UNURAN_pdfx\tdifference = $C_pdfx_diff\n";
	}

    }

    # End
    close UNURAN;
    close C;

}

# ----------------------------------------------------------------
# End

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
    my $distr = $_[0];
    my $fparam = $_[1];
    my $seed = $_[2];

    my $acg_query = "$ACG -l UNURAN -d $distr";
    $acg_query .= " -p \"$fparam\"" if $fparam; 

    my $generator = `$acg_query`;

    return "" if $?;

    # read code for generating distribution object from code
    die unless $generator =~ /(double\s+\*?fpar[^;]*;[^;]+;)/;
    my $distr_code = $1;
    
    my $urng = make_UNURAN_urng($seed);
    my $main = make_UNURAN_main($distr, $distr_code);

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

void init_urand(void)
{
    static struct prng *urng = NULL;

    if (urng == NULL)
	urng = prng_new("LCG(2147483647,16807,0,$seed)");
    else
	prng_seed(urng,$seed);

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
    my $distr_code = $_[1];

    my $code = <<EOS

int main()
{
    int i;
    double x, fx;
    UNUR_DISTR *distr;

    {
	$distr_code
    }

    init_urand();

    for (i=0; i<$sample_size; i++) {
	x = rand\_$distr();
	fx = unur_distr_cont_eval_pdf (x,distr);
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
    my $distr = $_[0];
    my $fparam = $_[1];
    my $seed = $_[2];

    my $acg_query = "$ACG -l C -d $distr";
    $acg_query .= " -p \"$fparam\"" if $fparam; 

    my $generator = `$acg_query`;

    return "" if $?;

    my $urng = make_C_urng($seed);
    my $main = make_C_main($distr);

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

EOS

    return $code;

} # end of make_C_urng() 

# ----------------------------------------------------------------
# Make main for test file (C version)

sub make_C_main
{
    my $distr = $_[0];

    my $code = <<EOS

#include <stdio.h>
#include <stdlib.h>

int main()
{
    int i;
    double x, fx;

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
# Make executable from test file (UNURAN version)

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

