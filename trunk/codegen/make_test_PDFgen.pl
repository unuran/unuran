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
my $make_test_PDFgen = "make_test_PDFgen.c";

# C file for tests
my $test_PDFgen = "test_PDFgen.c";

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
make_PDFgen_tests();

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

sub make_PDFgen_tests
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
# C file for tests

    # The C file header 
    my $test_header = <<EOX;
\#include <string.h>
\#include <float.h>
\#include <unuran.h>
\#include <config.h>
\#include \"PDFgen_source.h\"
\#ifdef WITH_DMALLOC
\#  include <dmalloc.h>
\#endif

\#if UNUR_URNG_TYPE != UNUR_URNG_PRNG
\#  error UNUR_URNG_TYPE must be set to UNUR_URNG_PRNG in unuran_config.h
\#endif

\#define FP_equal(a,b)  ((a)==(b) ||  fabs((a)-(b)) <= ((fabs(a)<fabs(b))?fabs(a):fabs(b))*100*DBL_EPSILON)

$hrule

EOX

#.................................................................
# Process each test

    # Each distribution 
    foreach my $d (keys %{$list_distr}) {

	# use local copy
	my $distr = $d;

	# Name of test routine
	my $testroutine = "test\_$distr";

	# Name of PDF function
	my $PDFroutine = "pdf\_$distr";

	# Distribution object
	my $distribution = $list_distr->{$distr};
	
# The test routine

	# Begin of test routines 
	my $test_test_routine = "int $testroutine (void)\n\{\n";

	# Declarations for test routines
	my $test_test_decl = 
	    "\tUNUR_DISTR *distr;\n".
	    "\tUNUR_GEN *gen;\n".
            "\tint i;\n".
	    "\tdouble x, f1, f2;\n".
	    "\tint n_failed = 0;\n".
	    "\n";

	# Body of test routine
	my $test_test_body =
	    "$distribution\n";         # the distribution object

	# Print info on screen
	$test_test_body .= "\tprintf(\"$distr \");\n";
	$test_test_body .= "\tfflush(stdout);\n\n";

	# We need a generator for importance sampling
	$test_test_body .= "\tgen = unur_init( unur_tdr_new(distr) );\n\n";

	# Compare PDFs
	$test_test_body .= <<EOX;
\tfor (i=0; i<$SAMPLE_SIZE; i++) \{
\t\tx  = unur_sample_cont(gen);
\t\tf1 = unur_distr_cont_eval_pdf (x,distr);
\t\tf2 = $PDFroutine (x);
\t\tif (!FP_equal(f1,f2)) \{
\t\t\tfprintf(stderr,\"error! %%g, %%g, diff = %%g\\n\",f1,f2,f1-f2);
\t\t\t++n_failed;
\t\t\}
\t\}

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
	# Begin of make test routines 
	my $make_test_routine = "void $testroutine (FILE *out)\n\{\n";

	# Declarations for make test routine
	my $make_test_decl = 
	    "\tUNUR_DISTR *distr;\n";

	# Body of make test routine
	my $make_test_body = 
	    "$distribution\n".                       # the distribution object
	    "\t_unurgen_C_PDF(distr,out,\"$PDFroutine\");\n".  # make code 
	    "\tunur_distr_free(distr);\n\n";         # free distribution object

	# write test file
	foreach $l (split /\n/, $test_test) { 
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

#.................................................................
# C file that makes test file

    # The C file header 
    my $make_header = <<EOX;
\#include <string.h>
\#include <unuran.h>
\#include \"PDFgen_source.h\"
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
    $test_main_body .=
	"\tLOG = fopen( \"test_PDFgen.log\",\"w\" );\n".
	"\tunur_set_stream( LOG );\n".
        "\tunur_set_default_debug(UNUR_DEBUG_ALL);\n\n";

    # Execute the tests
    foreach my $t (split /\n/, $test_list) {
	$test_main_body .= "\tn_failed += $t();\n";
    }
    $test_main_body .= "\n";

    # End of main()
    $test_main_body .= "\tprintf(\"\\n\");\n\n";
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

    # Declarations for main()
    my $make_main_decl = 
	"\tFILE *out;\n".
	"\tFILE *LOG;\n";

    # Body of make main()
    my $make_main_body =
	"\tLOG = fopen( \"make_test_PDFgen.log\",\"w\" );\n".
	"\tunur_set_stream( LOG );\n".
        "\tunur_set_default_debug(UNUR_DEBUG_ALL);\n\n";

    # open output stream
    $make_main_body .=
	"\tout = fopen(\"$test_PDFgen\",\"w\");\n\n";

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

} # end of make_PDFgen_tests()

# ----------------------------------------------------------------
