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

# Configuration file for tests
my $test_conf_file = "test.conf";

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
	$test_file .= "void $testroutine (FILE *out)\n\{\n";

	# Name of PDF function
	my $PDFroutine = "pdf\_$distr\_".$test_number;

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
	    $fpm =~ s/,\s*$/ \}\;/;
	}
	else {
	    $fpm = "double *fpm = NULL;";
	}

	# Make distribution object
	my $distribution = "\tUNUR_DISTR *distr;\n";
	$distribution .= "\t$fpm\n";
	$distribution .= "\tdistr = unur\_distr\_$distr(fpm,$n_params);\n";

	# Add Distribution object to test file
	$test_file .= $distribution;
	$test_file .= "\t_unurgen_C_PDF(distr,out,\"$PDFroutine\");\n";
	$test_file .= "\tunur_distr_free(distr);\n\n";
	$test_file .= "\tfprintf(out,\"\\n\");\n";

	# Header for test routine for distribution
	$test_file .= "\tfprintf(out,\"int $testroutine (void)\\n\{\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\tint i;\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\tdouble x, f1, f2;\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\tint n_failed = 0;\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\tUNUR_GEN *gen;\\n\");\n";

	# Write distribution object into test routine
	foreach $l (split /\n/, $distribution) {
	    $l =~ s/\t/\\t/g;
	    $test_file .= "\tfprintf(out,\"$l\\n\");\n";
	}

	# Generator for importance sampling
	$test_file .= "\tfprintf(out,\"\\tgen = unur_init( unur_tdr_new(distr) );\\n\");\n\n";

	# Print info on screen
	$test_file .= "\tfprintf(out,\"\\tprintf(\\\"$distr \\\");\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\tfflush(stdout);\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\n\");\n";

	# Compare PDFs
	$test_file .= "\tfprintf(out,\"\\tfor (i=0; i<$SAMPLE_SIZE; i++) \{\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\t\\tx  = unur_sample_cont(gen);\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\t\\tf1 = unur_distr_cont_eval_pdf (x,distr);\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\t\\tf2 = $PDFroutine (x);\\n\");\n";

	$test_file .= "\tfprintf(out,\"\\t\\tif (!FP_equal(f1,f2)) {\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\t\\t\\tfprintf(stderr,\\\"error! %%g, %%g, diff = %%g\\\\n\\\",f1,f2,f1-f2);\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\t\\t\\t++n_failed;\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\t\\t\}\\n\");\n";

	$test_file .= "\tfprintf(out,\"\\t\}\\n\");\n\n";
	$test_file .= "\tfprintf(out,\"\\n\");\n";

	# End of test routine
	$test_file .= "\tfprintf(out,\"\\tunur_distr_free(distr);\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\tunur_free(gen);\\n\");\n";
	$test_file .= "\tfprintf(out,\"\\treturn n_failed;\\n\");\n";
	$test_file .= "\tfprintf(out,\"\}\\n\\n\");\n";
	$test_file .= "\n";

	$test_file .= "}\n\n";

	# Add test to list
	$test_list .= "$testroutine\n";
    }

    # Check for missing CONTinuous distributions
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	unless ($test_distr{$d}) {
	    print STDERR "test missing for distribution \"$d\"\n";
	}
    }

    # Make header for C file that creates tests
    my $test_file_header = 
	 "\#include <stdio.h>\n"
	."\#include <stdlib.h>\n"
	."\#include <string.h>\n"
        ."\#include <unuran.h>\n"
        ."\#include \"PDFgen_source.h\"\n\n";

    # Make main()
    my $test_file_main = "int main(void)\n\{\n";

    # The header for the resulting test file
    $test_file_main .= "\tFILE *out;\n";
    $test_file_main .= "\n";
    $test_file_main .= "\tout = fopen(\"$test_PDFgen\",\"w\");\n";
    $test_file_main .= "\n";

    # The header for the test file
    $test_file_main .= "\tfprintf(out,\"#include <stdio.h>\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#include <stdlib.h>\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#include <string.h>\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#include <math.h>\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#include <float.h>\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#include <prng.h>\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#include <unuran.h>\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#include \\\"PDFgen_source.h\\\"\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#include <config.h>\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#ifdef WITH_DMALLOC\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#  include <dmalloc.h>\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#endif\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#if UNUR_URNG_TYPE != UNUR_URNG_PRNG\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#  error UNUR_URNG_TYPE must be set to UNUR_URNG_PRNG in unuran_config.h\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"#endif\\n\");\n";

    $test_file_main .= "\tfprintf(out,\"#define FP_equal(a,b) \");\n";
    $test_file_main .= "\tfprintf(out,\" ((a)==(b) || \");\n";
    $test_file_main .= "\tfprintf(out,\" fabs((a)-(b)) <= ((fabs(a)<fabs(b))?fabs(a):fabs(b))*100*DBL_EPSILON)\\n\\n\");\n\n";


    # Call the test routines
    foreach my $t (split /\n/, $test_list) {
	$test_file_main .= "\t$t(out);\n";
    }
    $test_file_main .= "\n";

    # Make the main file for the test file
    $test_file_main .= "\tfprintf(out,\"int main(void)\\n\{\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"\\tint n_failed = 0;\\n\\n\");\n";
    foreach my $t (split /\n/, $test_list) {
	$test_file_main .= "\tfprintf(out,\"\\tn_failed += $t();\\n\");\n";
    }

    $test_file_main .= "\tfprintf(out,\"\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"\\tprintf(\\\"\\\\n\\\");\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"\\n\");\n";

    $test_file_main .= "\tfprintf(out,\"\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"\\texit ((n_failed) ? 1 : 0);\\n\");\n";
    $test_file_main .= "\tfprintf(out,\"}\\n\");\n";

    # End if C file
    $test_file_main .= "\n";
    $test_file_main .= "\tfclose(out);\n";
    $test_file_main .= "\n";
    $test_file_main .= "\texit (0);\n";
    $test_file_main .= "\}\n\n";

    # Make C file for creating tests
    open TEST, ">$make_test_PDFgen" or die "cannot open file $make_test_PDFgen\n";
    print TEST $test_file_header;
    print TEST $test_file;
    print TEST $test_file_main;
    close TEST;

} # end of make_PDFgen_tests()

# ----------------------------------------------------------------
