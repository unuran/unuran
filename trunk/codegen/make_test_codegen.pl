#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

$DEBUG = 0;

# ----------------------------------------------------------------

require "read_PDF.pl";
require "read_test_conf.pl";

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
my $DISTR = read_PDFdata();

# For description of data fields in this list see file `read_PDF.pl'.

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

    # rule for C files 
    my $hrule = "/* ----------------------------------------------------------------- */";

#.................................................................
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

#.................................................................
# Get list of generation methods

    my $list_gen = get_test_methods( $test_conf_file );

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

\#define FP_equal(a,b)  ((a)==(b) ||  fabs((a)-(b)) <= ((fabs(a)<fabs(b))?fabs(a):fabs(b))*1000*DBL_EPSILON)

/* We use global variables for the uniform random number generators */
static struct prng *urng1 = NULL;
static struct prng *urng2 = NULL;

/* Define uniform(); use the second generator */
\#define uniform()  (prng_get_next(urng2))

$hrule

EOX

#.................................................................
# Process each test

    # short name for distribution
    my $last_distr_short;

    # Each distribution 
    foreach my $d (sort keys %{$list_distr}) {

	# use local copy
	my $distr = $d;

	# short name for distribution
	my $distr_short = $distr;
	$distr_short =~ s/\_[^\_]+$//;
	my $distr_print_name;
	if ($last_distr_short ne $distr_short) {
	    $distr_print_name = "\tprintf(\" $distr_short \");\n";
	    $last_distr_short = $distr_short;
	}
	else {
	    $distr_print_name = "";
	}

	# Each generation method
	foreach my $g (keys %{$list_gen}) {

	    # use local copy
	    my $gen = $g;

	    # Name of test routine
	    my $testroutine = "test\_$distr\_$gen";

	    # Distribution object
	    my $distribution = $list_distr->{$distr};
	
	    # Generator object
	    my $generator;
	    foreach my $l (split /\;/, $list_gen->{$gen}) {
		$l =~ s/^\s+//;
		$l =~ s/\s+$//;
		$generator .= "\t$l;\n";
	    }
	    $generator .=  "\tgen = unur_init( par );\n";

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
		"$distribution\n".         # the distribution object
		"$generator\n";            # the generator object

	    # Init of generator failed
	    $test_test_body .=
		"\tif (gen == NULL) return 1;\n\n";

	    # Set uniform random number generator
	    $test_test_body .= "\tunur_chg_urng(gen,urng1);\n\n";

	    # Print info on screen
	    $test_test_body .= $distr_print_name; 
	    $test_test_body .= "\tfflush(stdout);\n\n";

	    # Compare generator output
	    $test_test_body .= <<EOX;
\tfor (i=0; i<$SAMPLE_SIZE; i++) {
\t\tx1 = unur_sample_cont(gen);
\t\tx2 = rand_$distr();
\t\tif (!FP_equal(x1,x2)) {
\t\t\tfprintf(stderr,\"error! %%g, %%g, diff = %%g\\n\",x1,x2,x1-x2);
\t\t\t++n_failed;
\t\t}
\t}

EOX

            # End of test routine
            $test_test_body .= 
		"\tunur_distr_free(distr);\n".
		"\tunur_free(gen);\n".
		"\t(n_failed > 0) ? printf(\"!\") : printf(\"+\");\n".
                "\tfflush(stdout);\n\n".
		"\treturn n_failed;\n".
	        "\}\n\n";
	    
	    # The test routine
	    my $test_test = 
		$test_test_routine.
		$test_test_decl.
		$test_test_body;

	    # Test routine when init of generator object failed
	    my $test_test_failed = 
		$test_test_routine.
		$distr_print_name. 
		"\tprintf(\".\");\n".
                "\tfflush(stdout);\n\n".
		"\treturn 0;\n".
		"\}\n";

# The make test routine
	    # Begin of make test routines 
	    my $make_test_routine = "void $testroutine (FILE *out)\n\{\n";

	    # Declarations for make test routine
	    my $make_test_decl = 
		"\tUNUR_DISTR *distr;\n".
                "\tUNUR_PAR *par;\n".
                "\tUNUR_GEN *gen;\n";

	    # Body of make test routine
	    my $make_test_body = 
		"$distribution\n".                       # the distribution object
		"$generator\n";                          # the generator object

	    # Init of generator failed
	    $make_test_body .=
		"\tif (gen == NULL) {\n".
		"\t\tunur_distr_free(distr);\n";
	    foreach my $l (split /\n/, $test_test_failed) { 
		$l =~ s/\t/\\t/g;
		$l =~ s/\\n/\\\\n/g;
		$l =~ s/\"/\\\"/g;
		$make_test_body .= "\t\tfprintf(out,\"$l\\n\");\n";
	    }
	    $make_test_body .=
		"\t\treturn;\n".
		"\t}\n\n";

	    # Make code
	    $make_test_body .=
		"\tunur_acg(gen,out,\"$distr\");\n\n".   # make code 
		"\tunur_distr_free(distr);\n".           # free distribution object
		"\tunur_free(gen);\n\n";                 # free generator object

	    # write test file
	    foreach my $l (split /\n/, $test_test) { 
		$l =~ s/\t/\\t/g;
		$l =~ s/\\n/\\\\n/g;
		$l =~ s/\"/\\\"/g;
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
    }

#.................................................................
# C file that makes test file

    # The C file header 
    my $make_header = <<EOX;
\#include <string.h>
\#include <unuran.h>
\#include <unuran_acg.h>
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

    # Body of main()
    my $test_main_body = '';

    # Log files 
    $test_main_body =
	"\tLOG = fopen( \"test_codegen.log\",\"w\" );\n".
	"\tunur_set_stream( LOG );\n".
        "\tunur_set_default_debug(UNUR_DEBUG_ALL);\n\n";

    # Initialize URNG
    my $seed = int(rand 123456) + 1;
    $test_main_body .= "\turng1 = prng_new(\"mt19937($seed)\");\n";
    $test_main_body .= "\turng2 = prng_new(\"mt19937($seed)\");\n";
    $test_main_body .= "\n";

    # Execute the tests
    foreach my $t (split /\n/, $test_list) {
	$test_main_body .= "\tn_failed += $t();\n";
    }
    $test_main_body .= "\n";

    # End of main()
    $test_main_body .= "\tprintf(\"\\n\");\n\n";
    $test_main_body .= "\tfclose(LOG);\n\n";
    $test_main_body .= "\texit ((n_failed) ? EXIT_FAILURE : EXIT_SUCCESS);\n";
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
        "\tunur_set_default_debug(UNUR_DEBUG_ALL);\n\n";

    # Body of make main()
    $make_main_body .= 
	"\tout = fopen(\"$test_codegen\",\"w\");\n\n";

    # Write the header for the test file
    foreach my $l (split /\n/, $test_header) {
	$l =~ s/\t/\\t/g;
	$l =~ s/\\n/\\\\n/g;
	$l =~ s/\"/\\\"/g;
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
	$l =~ s/\\n/\\\\n/g;
	$l =~ s/\"/\\\"/g;
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

    open TEST, ">$make_test_codegen" or die "cannot open file $make_test_codegen\n";
    print TEST $make_header;
    print TEST $test_file;
    print TEST $make_main;
    close TEST;

} # end of make_codegen_tests()

# ----------------------------------------------------------------


