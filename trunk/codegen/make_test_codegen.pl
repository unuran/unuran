#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

$DEBUG = 0;

# ----------------------------------------------------------------

require "readPDF.pl";

# ----------------------------------------------------------------

# Read configuration file name for tests from argument list ...
my $test_conf_file = shift
    or die "no argument given";

# C file for making code generator tests
my $make_test_PDFgen = "make_test_codegen.c";

# C file for tests
my $test_PDFgen = "test_codegen.c";

# Sample size for test
my $SAMPLE_SIZE = 10000;

# ----------------------------------------------------------------
# List of distributions
my $DISTR = read_PDFdata();

# For description of data fields in this list see file `readPDF.pl'.

# ................................................................

# Print result on screen
    if ($DEBUG) {
	print_data($DISTR);
    }

# ----------------------------------------------------------------

# Make test files for code generator
make_codegen_tests();

# ----------------------------------------------------------------
# End

exit 0;

# ----------------------------------------------------------------
#
# Subroutines
#
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Make C file for testing code generator

sub make_codegen_tests
{
    my $test_file;
    my $test_list;
    my $test_number = 0;

    # Mark distributions for which tests exist
    my %test_distr;

    # Read the test.conf file
    open CONF, $test_conf_file or die "cannot open file $test_conf_file\n";
    my $conf_content;
    while (<CONF>) {
	next if /^\#/;
	$conf_content .= $_;
    }
    close CONF;
    
    # Get tests
    my @tests = split /\n\s*\n/, $conf_content;

    # rule for C files 
    my $hrule = "/* ----------------------------------------------------------------- */";

#.................................................................
# C file for tests

    # The C file header 
    my $test_header = <<EOX;
\#include <string.h>
\#include <float.h>
\#include <unuran.h>
\#include <config.h>

\#ifdef WITH_DMALLOC
\#  include <dmalloc.h>
\#endif

\#if UNUR_URNG_TYPE != UNUR_URNG_PRNG
\#  error UNUR_URNG_TYPE must be set to UNUR_URNG_PRNG in unuran_config.h
\#endif

\#define FP_equal(a,b)  ((a)==(b) ||  fabs((a)-(b)) <= ((fabs(a)<fabs(b))?fabs(a):fabs(b))*100*DBL_EPSILON)

/* We use global variables for the uniform random number generators */
static struct prng *urng1 = NULL;
static struct prng *urng2 = NULL;

/* Define uniform(); use the second generator */
\#define uniform()  (prng_get_next(urng2))

$hrule

EOX

#.................................................................
# Process each test

    foreach my $t (@tests) {
	# There might be empty entries in @tests
	next unless $t;

	# Get name and parameters of distribution
	$t =~ /DISTR:\s*(\w+)\s*\(([^\)]*)\)/
	    or die "cannot find valid distribution tag for test";
	
	# Store data
	my $distr = $1;
	my $params = $2;

	# Mark distribution
	$test_distr{$distr} = 1;

	# Check for existing distribution
	die "Unknown distribution: $distr" unless $DISTR->{$distr};

	# Name of test routine
	++$test_number;
	my $testroutine = "test\_$distr\_".$test_number;

	# Name of distribution for tests
	my $distr_name = "$distr\_".$test_number;

	# Process parameters
	my @param_list = split /\,/, $params;
	my $n_params = $#param_list + 1;
	my $fpm;
	if ($n_params > 0) {
	    $fpm = "double fpm[] = { ";
	    foreach my $p (@param_list) {
		if ($p =~ /(.+)\.\.(.+)/) {
		    $fpm .= ($1>0) ? exp(log($1)+rand()*(log($2)-log($1))) : $1+rand()*($2-$1);
		    $fpm .= ", ";
		}
		else {
		    $fpm .= "$p, ";
		}
	    }
	    $fpm =~ s/,\s*$/ \}/;
	}
	else {
	    $fpm = "double *fpm = NULL";
	}

	# Make distribution object
	my $distribution = 
	    "\t{\n".
	    "\t\t$fpm;\n".
	    "\t\tdistr = unur_distr_$distr(fpm,$n_params);\n".
	    "\t}\n";    

	# Make Generator object
	my $generator = 
	    "\tpar = unur_tdr_new (distr);\n".
	    "\tunur_tdr_set_max_sqhratio (par, 0.);\n".
            "\tgen = unur_init( par );\n";

# The test routine
	# Begin of test routines 
	my $test_test_routine = "int $testroutine (void)\n\{\n";
	
	# Declarations for test routines
	my $test_test_decl = 
	    "\tUNUR_DISTR *distr;\n".
            "\tUNUR_PAR *par;\n".
            "\tUNUR_GEN *gen;\n".
	    "\tint i;\n".
	    "\tdouble x1, x2;\n".
	    "\tint n_failed = 0;\n".
	    "\n";

	# Body of test routine
	my $test_test_body =
	    "$distribution\n".                          # the distribution object
	    "$generator\n";                             # the generator object

	# Set uniform random number generator
	$test_test_body .= "\tunur_chg_urng(gen,urng1);\n\n";

	# Print info on screen
	$test_test_body .= "\tprintf(\\\"$distr \\\");\n";
	$test_test_body .= "\tfflush(stdout);\n\n";

	# Compare generator output
	$test_test_body .= <<EOX;
\tfor (i=0; i<$SAMPLE_SIZE; i++) {
\t\tx1 = unur_sample_cont(gen);
\t\tx2 = rand_$distr_name();
\t\tif (!FP_equal(x1,x2)) {
\t\t\tfprintf(stderr,\\\"error! %%g, %%g, diff = %%g\\\\n\\\",x1,x2,x1-x2);
\t\t\t++n_failed;
\t\t}
\t}

EOX

        # End of test routine
	$test_test_body .= "\tunur_distr_free(distr);\n";
	$test_test_body .= "\tunur_free(gen);\n";
	$test_test_body .= "\treturn n_failed;\n";
	$test_test_body .= "\}\n\n";

	# The test routine
	my $test_test = 
	    $test_test_routine.
	    $test_test_decl.
	    $test_test_body;

# The make test routine
	# Begin of test routines 
	my $make_test_routine = "void $testroutine (FILE *out)\n\{\n";

	# Declarations for make test routine
	my $make_test_decl = 
	    "\tUNUR_DISTR *distr;\n".
            "\tUNUR_PAR *par;\n".
            "\tUNUR_GEN *gen;\n";

	# Body of make test routine
	my $make_test_body = 
	    "$distribution\n".                          # the distribution object
	    "$generator\n".                             # the generator object
	    "\tunurgen(gen,out,\"$distr_name\");\n\n".  # make code 
	    "\tunur_distr_free(distr);\n".              # free distribution object
	    "\tunur_free(gen);\n\n";                    # free generator object

	# Header for test routine for distribution
	foreach $l (split /\n/, $test_test) { 
	    $l =~ s/\t/\\t/g;
	    $make_test_body .= "\tfprintf(out,\"$l\\n\");\n";
	}

	# End of make test routine
	$make_test_body .= "}\n\n";

	# The make test routine
	my $make_test = 
	    $make_test_routine.
	    $make_test_decl.
	    $make_test_body;

	# Add test routine to output
	$test_file .= $make_test; 

# Add test to list
	$test_list .= "$testroutine\n";
    }

#.................................................................
# Check for missing CONTinuous distributions

    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	unless ($test_distr{$d}) {
	    print STDERR "test missing for distribution \"$d\"\n";
	}
    }

#.................................................................
# C file that makes test file

    # The C file header 
    my $make_header = <<EOX;
\#include <string.h>
\#include <unuran.h>
\#include <unurgen.h>
\#ifdef WITH_DMALLOC
\#  include <dmalloc.h>
\#endif

$hrule

EOX

# The test main()

    # Begin of main()
    my $test_main_routine = "int main(void)\n\{\n";

    # Declarations for main()
    my $test_main_decl = 
	"\tFILE *LOG;\n".
	"\tint n_failed = 0;\n";

    # Body of make main()
    my $test_main_body = '';

    # Log files 
    $test_main_body =
	"\tLOG = fopen( \\\"test_codegen.log\\\",\\\"w\\\" );\n".
	"\tunur_set_stream( LOG );\n".
        "unur_set_default_debug(UNUR_DEBUG_ALL);\n\n";

    # Initialize URNG
    my $seed = int(rand 123456) + 1;
    $test_main_body .= "\turng1 = prng_new(\\\"mt19937($seed)\\\");\n";
    $test_main_body .= "\turng2 = prng_new(\\\"mt19937($seed)\\\");\n";
    $test_main_body .= "\n";

    # Execute the tests
    foreach my $t (split /\n/, $test_list) {
	$test_main_body .= "\tn_failed += $t();\n";
    }
    $test_main_body .= "\n";

    # End of main()
    $test_main_body .= "\tprintf(\\\"\\\\n\\\");\n\n";
    $test_main_body .= "\tfclose(LOG);\n\n";
    $test_main_body .= "\texit ((n_failed) ? 1 : 0);\n";
    $test_main_body .= "}\n\n";

    # test main()
    my $test_main = 
	$test_main_routine.
	    $test_main_decl."\n".
	$test_main_body;

# The make test main() 

    # Begin of main()
    my $make_main_routine = "int main(void)\n\{\n";

    # Declarations for make main()
    my $make_main_decl = 
	"\tFILE *out;\n".
	"\tFILE *LOG;\n";

    # Log files 
    my $make_main_body =
	"\tLOG = fopen( \"make_test_codegen.log\",\"w\" );\n".
	"\tunur_set_stream( LOG );\n".
        "unur_set_default_debug(UNUR_DEBUG_ALL);\n\n";

    # Body of make main()
    $make_main_body .= 
	"\tout = fopen(\"$test_PDFgen\",\"w\");\n";
    $make_main_body .= "\n";

    # Write the header for the test file
    foreach my $l (split /\n/, $test_header) {
	$l =~ s/\t/\\t/g;
	$make_main_body .= "\tfprintf(out,\"$l\\n\");\n";
    }
    $make_main_body .= "\n";

    # Write the calls for the test routines
    foreach my $t (split /\n/, $test_list) {
	$make_main_body .= "\t$t(out);\n";
    }
    $make_main_body .= "\n";

    # Write main for test file
    foreach my $l (split /\n/, $test_main) {
	$l =~ s/\t/\\t/g;
	$make_main_body .= "\tfprintf(out,\"$l\\n\");\n";
    }
    $make_main_body .= "\n";

    # End of make main()
    $make_main_body .= "\tfclose(out);\n";
    $make_main_body .= "\tfclose(LOG);\n\n";
    $make_main_body .= "\texit (EXIT_SUCCESS);\n";
    $make_main_body .= "}\n\n";

    # make main()
    my $make_main =
	$make_main_routine.
	$make_main_decl."\n".
	$make_main_body;

#.................................................................
# Make C file for creating tests

    open TEST, ">$make_test_PDFgen" or die "cannot open file $make_test_PDFgen\n";
    print TEST $make_header;
    print TEST $test_file;
    print TEST $make_main;
    close TEST;

} # end of make_codegen_tests()

# ----------------------------------------------------------------
