#!/usr/bin/perl

############################################################
# $Id$

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;

    print STDERR <<EOM;
usage: $progname <test.conf> 
      
Scans <test.conf> and generates a C file with test routines.
File must have suffix "conf".

EOM

    exit;
}

############################################################

use English;

############################################################

$unuran_h_file = "../src/unuran.h";

#############################################################

############################################################
#                                                          #
#  main                                                    #
#                                                          #
############################################################

# get file names

# program name ...
$name_program = $0;
$name_program =~ s#^.*/##g;

# read in input file name from argument list ...
$file_in = shift;
(usage and die) unless $file_in;

#check suffix
(usage and die) unless ($file_in =~ /\.conf$/);

# compose name of output file
$file_out = $file_in;
$file_out =~ s#^.*/##g;
$file_out =~ s/\.conf$/\.c/;
$file_out = "t_$file_out";

# name of log files
$file_testlog = $file_out;
$file_testlog =~ s/\.c$/_test\.log/;
$file_unuranlog = $file_out;
$file_unuranlog =~ s/\.c$/_unuran\.log/;

#get name of file
$file_name = $file_in;
$file_name =~ s#^.*/##g;
$file_name =~ s/\.conf$//;

#open files ...
open (IN,"$file_in")    or die "Cannot open file $file_in for reading";
open (OUT,">$file_out") or die "Cannot open file $file_out for writing";

# write out put header ...
print OUT "/*\n\tfile automatically generated by $name_program\n\t";
print OUT scalar localtime;
print OUT "\n*/\n\n";

############################################################
#                                                          #
#  scan main section                                       #
#                                                          #
############################################################

# data we want ...
undef $method;
undef $gen_type;
undef $distr_type;
undef $urng;
undef $C_header_aux;

my $section = "main";
undef my $next_section;

# search for begin of [main] section ...
$_ = "";
$_ = <IN> until /^\[main/;

# skip section marker
$_ = <IN> if /^\[main\]/;

# scan section ...
while (1) {

    # search for begin of next (sub) section ...
    $_ = <IN> until /^\[/;

    # check (sub) section name ...
    unless (/^\[($section)/) {
	die "wrong subsection" if /-/;
	# next section ...
	/\[(.*)\]/;
	$next_section = $1;
	last;
    }

    # get subsection name ...
    die "wrong section" unless /^\[$section\s+\-\s+(.*):/;
    my $subsection = $1;

    # search for closing ] ...
    $_ = <IN> until /\]/;
    $_ = <IN>;

    # scan subsection ...
    while (1) {
	if ( $subsection =~ /data/ ) {
	    #read data ...
	    if (/^\s*method\s*:\s*(\w+)/)        { $method = $1; }
	    if (/^\s*urng\s*:\s*(.+)[\s\n]/)     { $urng = $1; }
	}
	elsif ( $subsection =~ /header/ ) {
	    # add verbatim to C header (except conf comments) ...
	    $C_header_aux .= $_ unless /^\#/; 
	}
	else {
	    die "Unknown subsection in [main]";
	}
	# next line ...
	$_ = <IN>;
	last if /^\[/;       # start of next (sub) section
    }
}

# transform \# --> #
$C_header_aux =~ s/\\#/#/g;

# check data ...
die "Data missing" unless (defined $method and
			   defined $urng);

# name of method 
$method =~ tr/[A-Z]/[a-z]/;
$METHOD = $method;
$METHOD =~ tr/[a-z]/[A-Z]/;


############################################################
#                                                          #
#  C header and main()                                     #
#                                                          #
############################################################

# write ...
print_C_prototypes();

print OUT <<EOM;

/*---------------------------------------------------------------------------*/

int main()
{ 
	/* open log file for unuran and set output stream for unuran messages */
	UNURANLOG = fopen( "$file_unuranlog","w" );
	abort_if_NULL( stderr,-1, UNURANLOG );
	unur_set_stream( UNURANLOG );

	/* open log file for testing */
	TESTLOG = fopen( "$file_testlog","w" );
	abort_if_NULL( stderr,-1, TESTLOG );

	/* write header into log file */
  	{
		time_t started;  
		fprintf(TESTLOG,"\\nUNURAN - Universal Non-Uniform RANdom number generator\\n\\n");
		if (time( \&started ) != -1)
			fprintf(TESTLOG,"%s",ctime(\&started));
		fprintf(TESTLOG,"\\n====================================================\\n\\n");
	}

	/* set uniform random number generator */
	urng = prng_new("$urng");
	unur_set_default_urng(urng);

	/* set default debugging flag */
	unur_set_default_debug(UNUR_DEBUG_ALL);

	/* start test */
	printf("$method: ");

	/* run tests */
	test_new();
	test_set();
	test_get();
	test_chg();
	test_init();
	test_sample();
	test_reinit();
	test_validate();

	/* test finished */
	printf("\\n");  fflush(stdout);

	/* close log files */
	fprintf(TESTLOG,"\\n====================================================\\n\\n");
	if (test_ok)
		fprintf(TESTLOG,"All tests PASSED.\\n");
	else
		fprintf(TESTLOG,"Test(s) FAILED.\\n");

	fclose(UNURANLOG);
	fclose(TESTLOG);

	/* exit */
	exit( (test_ok) ? 0 : -1 );

} /* end of main */

EOM

############################################################
#                                                          #
#  scan all section until [validate]                       #
#                                                          #
############################################################

    $section = $next_section;

until ($section =~ /validate/) {
    $next_section = scan_section($section);
    $section = $next_section;
}

############################################################
#                                                          #
#  validate                                                #
#                                                          #
############################################################

scan_validate();

############################################################
#                                                          #
#  verbatim                                                #
#                                                          #
############################################################

print OUT "/*---------------------------------------------------------------------------*/\n\n";

# scan section ...
while (<IN>) {
    next if /^\#/;  # comment line
    s/\\#/#/g;
    print OUT;
}

print_C_routines();

add_unur_set_verify_routine();

############################################################
#                                                          #
#  end                                                     #
#                                                          #
############################################################

print OUT <<EOM;
/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
int main() { exit(77); } /* ignore test */
/*---------------------------------------------------------------------------*/
#endif  /* T_$METHOD */
/*---------------------------------------------------------------------------*/
EOM

############################################################

# close files ...
close(IN);
close(OUT);

############################################################

exit 0;

############################################################

############################################################
############################################################
############################################################
############################################################

############################################################
#                                                          #
#  subroutine: scan section                                #
#                                                          #
############################################################

sub scan_section {
    my $section = $_[0];   # name of section
    my $next_section;

    print OUT <<EOM;
/*---------------------------------------------------------------------------*/

void test_$section (void)
{
	/* start test */
	printf("[$section "); fflush(stdout);
	fprintf(TESTLOG,"\\n[$section]\\n");

	/* reset counter */
	n_tests_failed = 0;
  
EOM

    # search for begin of next (sub) section ...
    $_ = <IN> until /^\[/;
	
    # skip section marker
    $_ = <IN> if /^\[$section\]/;

    # scan section ...
    while (1) {

	# search for begin of next (sub) section ...
	$_ = <IN> until /^\[/;

	# check (sub) section name ...
	unless (/^\[($section)/) {
	    die "wrong subsection" if /-/;
	    # next section ...
	    /\[(.*)\]/;
	    $next_section = $1;
	    last;
	}

	# skip section marker
	next if /^\[$section\]/;

	# get subsection name ...
	die "wrong section" unless /^\[$section\s+\-\s+(.*):/;
	print OUT "/* $1 */\n";

	# clear variables ...
	$line = "";
	$subsection_closing = "";

	unless ( /\]\s*$/ ) {
	    while (<IN>) { 
		next if /^\#/;  # comment line
		# read until next empty line
		$line .= $_;
		last unless /[^\s]+/;
	    }

            # transform \# --> #
	    $line =~ s/\\#/#/g;

	    # there should be a closing ] ...
	    die "closing ] missing" unless $line =~ /\]\s*$/;
	
	    # remove this bracket ...
	    $line =~ s/\]\s*$/\n/;
	
	    # add declarations for variables distr, par and gen ...
	    # (and free these at the end of the subsection)
	    if ($line =~ /\s+gen\s+=/) { 
		$line = "UNUR_GEN   *gen = NULL\;\n".$line;
		$subsection_closing .= "unur_free(gen)\;\n";
	    }
	    if ($line =~ /\s+par\s+=/) {
		$line = "UNUR_PAR   *par = NULL\;\n".$line; 
	    }
	    if ($line =~ /\s+distr\s+=/) {
		$line = "UNUR_DISTR *distr = NULL\;\n".$line; 
		$subsection_closing .= "unur_distr_free(distr)\;\n";
	    }
	
	    # lines indicated with "<-- ! NULL" must not produce a NULL pointer 
	    $line =~ s/^(.*)=(.*)<--\s+!\s*NULL\s*\n/$1=$2\nabort_if_NULL\(TESTLOG, $INPUT_LINE_NUMBER, $1\)\;\n/mg;

	}
	else {
	    $_ = <IN>;
	}

	# print subsection header ...
	print OUT "{\n$line";

	# scan subsection body ...
	while (1) {
	    # search for non-empty line ...
	    $_ = <IN> until /\w+/;

	    # scan till next empty line of opening [ ...
	    my $line = "";
	    while (1) {
		if ( /^\#/ ) {
		    # comment line
		    $_ = <IN>;		# next line ...
		    next;
		}
		last unless /\w+/;   # empty line
		last if /^\[/;       # start of next (sub) section
		
		# convert line starting with ~ ...
		s/~/unur_$method\_$section/ if /^~/;

		# append line ...
		$line .= $_;

		# next line ...
		$_ = <IN>;
	    }
	    last if /^\[/;       # start of next (sub) section

            # transform \# --> #
	    $line =~ s/\\#/#/g;

	    # analyze string ...
	    my ($code,$test_command,$errno) = split /-->/, $line;

	    # is there any test ?
	    unless ($test_command) {
		print OUT $code if defined $code;
		next;
	    }
	    
	    # split into lines again ...
	    my @lines = split /\n/, $code;

	    # get last C command ...
	    undef my $last_C_line;
	    while ( $last_C_line = pop @lines ) {
		last if $last_C_line =~ /\w+/;
		last if !defined $last_C_line;
	    }
	    $last_C_line =~ s/;\s*$// if defined $last_C_line;

	    # print ...
	    print OUT "\nunur_errno = 0\;\n";

	    # C lines ...
	    foreach (@lines) {
		print OUT $_,"\n";
	    }

	    # test ...
	    print_test_command( $test_command, $last_C_line );

	    # error code ...
	    $errno =~ s/[\s\n]+//g if defined $errno;
	    print OUT "n_tests_failed += check_errorcode( TESTLOG, $INPUT_LINE_NUMBER, $errno )\;\n" if $errno;

	}

	# close subsection ...
	print OUT $subsection_closing,"}\n\n";
    }

  # print section closing ...
  print OUT <<EOM;

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_$section() */

EOM

    # return name of next section ...
    return $next_section;
} # end of scan_section()

############################################################
#                                                          #
#  subroutine: scan validate section                       #
#                                                          #
############################################################

sub scan_validate {
    my $section = "validate";   # name of section
    my $next_section;

    my @generators;
    my @distributions;
    my $chi2;
    my $timing;

    # search for begin of next (sub) section ...
    $_ = <IN> until /^\[/;
	
    # skip section marker
    $_ = <IN> if /^\[$section\]/;

    # scan section ...
    while (1) {
	
	# search for begin of next (sub) section ...
	$_ = <IN> until /^\[/;
	
	# check (sub) section name ...
	unless (/^\[($section)/) {
	    die "wrong subsection" if /-/;
	    # next section ...
	    /\[(.*)\]/;
	    $next_section = $1;
	    last;
	}
	
	# skip section marker
	next if /^\[$section\]/;
	
	# get subsection name ...
	die "wrong section" unless /^\[$section\s+\-\s+(.*):/;
	$subsection = $1;
	
	# there should be a closing ] ...
	die "closing ] missing" unless /\]\s*$/;
	
	# go to next line
	$_ = <IN>;
	
	# scan subsection body ...
	while (1) {
	    # search for non-empty line ...
	    $_ = <IN> until /\w+/;
	    
	    # scan till next empty line of opening [ ...
	    my $line = "";
	    while (1) {
		if ( /^\#/ ) {
		    # comment line
		    $_ = <IN>;		# next line ...
		    next;
		}
		last unless /\w+/;   # empty line
		last if /^\[/;       # start of next (sub) section
		
		# convert line starting with ~ ...
		s/~/unur_$method\_$section/ if /^~/;

		# append line ...
		$line .= $_;
		
		# next line ...
		$_ = <IN>;
	    }
	    last if /^\[/;       # start of next (sub) section

	    # skip over empty lines between comment lines
	    next if $line =~ /^[\s\n]*$/;

            # transform \# --> #
	    $line =~ s/\\#/#/g;
	    
	    # subsection generators
	    if ($subsection eq "generators") {
		push @generators, $line;
		next;
	    }
	    
	    if ($subsection eq "distributions") {
		push @distributions, $line;
		next;
	    }

	    if ($subsection eq "test chi2") {
		$chi2 .= $line;
		next;
	    }

	    if ($subsection eq "timing") {
		$timing .= $line;
		next;
	    }

	}
    }

    # number of given distributions
    $n_distributions = $#distributions + 1;

    # number of generators
    $n_generators = $#generators + 1;

    # print out ...

    ## header ##

    print OUT "/*---------------------------------------------------------------------------*/\n\n";
    print OUT "void test_validate (void)\n{\n";

    print OUT "\tUNUR_DISTR *distr[$n_distributions];\n";
    print OUT "\tUNUR_PAR *par;\n";
    print OUT "\tUNUR_GEN *gen;\n";
    print OUT "\tint n_tests_failed = 0;\n";

    print OUT "\tdouble *darray;\n";
    print OUT "\tdouble fpm[10];\n";

    print OUT "\n\t/* start test */\n";
    print OUT "\tprintf(\"[validate \"); fflush(stdout);\n";
    print OUT "\tfprintf(TESTLOG,\"\\n[validate]\\n\");\n\n";
    print OUT "\t/* reset counter */\n\tn_tests_failed = 0;\n\n";

    ## list of distributions ##

    print OUT "\n/* distributions: $n_distributions */\n";
    foreach (@distributions) {
	print OUT "$_\n\n";
    }

    ## chi^2 test ##

    if ($chi2) {

	# analyse chi^2 tests
	my @chi2tests = split /\n/, $chi2; 
	die "wrong number of chi2 tests" unless ($#chi2tests <= $#distributions);

	print OUT "\tprintf(\"\\n(chi^2) \"); fflush(stdout);\n";
	print OUT "\n/* chi^2 tests: ".($#generators+1)*$n_distributions." */\n\n";
	print OUT "\tunur_set_default_debug(~UNUR_DEBUG_SAMPLE);\n";
	print OUT "\tfprintf( TESTLOG,\"\\nChi^2 Test:\\n\");";
	
	foreach $test (@chi2tests) {
	    die "invalide test line" unless ($test =~ /<(\d+)>/);
	    my $n_distr = $1;
	    print OUT "/* distribution [$n_distr] */\n\n";
	    $test =~ s/\#.+$//;    # remove comments
	    $test =~ s/^\s+//;
	    $test =~ s/\s+$//;
	    my @gentest = split /\s+/, $test;
	    shift @gentest;
	    die "invalide number of test indicators" unless ($#gentest == $#generators);
	    foreach (@generators) {
		# get entry for generator
		my $genline = $_;
		
		# remove [..] from par[..]
		# (it is just to number generators for convience)
		$genline =~ s/par\[(\d+)\]/par/g;
		
		# insert distribution object
		$genline =~ s/\@distr\@/distr\[$n_distr\]/g;
		
		# read what we have to test
		$todo = shift @gentest;
		
		# split into lines again
		my @lines = split /\n/, $genline;
		
		# print lines 
		print OUT "\tunur_errno = 0;\n";
		
		my $have_gen_lines = 0;
		foreach $l (@lines) {
		    if ($l =~ /gen/ and !$have_gen_lines) {
			$have_gen_lines = 1;
			if ( $todo eq '.' ) {
			    # nothing to do
			    print OUT "\tgen = NULL; if (0) {"; 
			    last;
			}
			else {
			    print OUT "\tgen = unur_init(par);\n\tif (gen) {\n";
			}
		    }
		    
		    print OUT "$l\n";
		}
		
		if ($have_gen_lines) {
		    print OUT "\t;}\n";
		}
		else {
		    if ( $todo eq '.' ) {
			# nothing to do
			print OUT "\tgen = NULL;\n"; }
		    else {
			print OUT "\tgen = unur_init(par);\n"; }
		}
		print OUT "\tn_tests_failed += run_validate_chi2( TESTLOG, 0, gen, '$todo' );\n";
		print OUT "\tunur_free(gen);\n\n";
	    }	    
	}
    }

    ## timing ##

    if ($timing) {

	# analyse chi^2 tests
	my @timingtests = split /\n/, $timing; 
	die "wrong number of timing tests" unless ($#timingtests <= $#distributions);


	print OUT "\tprintf(\"\\n(timing) \"); fflush(stdout);\n";
	print OUT "\n/* timing tests: ".($#generators+1)*$n_distributions." */\n\n";
	print OUT "\tunur_set_default_debug(~1u);\n";

	print OUT "{\n\tdouble time_setup, time_sample;\n";
	print OUT "\tdouble timing_result[$n_generators];\n";
	print OUT "\tint i;\n\n";

	print OUT "\tfprintf( TESTLOG,\"\\nTimings (marginal generation times in micro seconds):\\n\");\n";
	print OUT "\tfor (i=0; i<$n_generators; i++)\n";
	print OUT "\t\tfprintf( TESTLOG, \"  [%2d]\", i);";
	print OUT "\tfprintf( TESTLOG,\"\\n\");\n";
	
	foreach $test (@timingtests) {
	    die "invalide test line" unless ($test =~ /<(\d+)>/);
	    my $n_distr = $1;

	    print OUT "/* distribution [$n_distr] */\n\n";
	    print OUT "\tfor (i=0; i<$n_generators; i++) timing_result[i] = -1.;\n";

	    $test =~ s/\#.+$//;    # remove comments
	    $test =~ s/^\s+//;
	    $test =~ s/\s+$//;
	    my @gentest = split /\s+/, $test;
	    shift @gentest;
	    die "invalide number of test indicators" unless ($#gentest == $#generators);
	    my $n_gen = 0;
	    foreach (@generators) {
		# get entry for generator
		my $genline = $_;
		
		# remove [..] from par[..]
		# (it is just to number generators for convience)
		$genline =~ s/par\[(\d+)\]/par/g;
		
		# insert distribution object
		$genline =~ s/\@distr\@/distr\[$n_distr\]/g;
		
		# read what we have to test
		$todo = shift @gentest;
		
		# split into lines again
		my @lines = split /\n/, $genline;
		
		# print lines 
		print OUT "\tunur_errno = 0;\n";
		
		my $have_gen_lines = 0;
		foreach $l (@lines) {
		    if ($l =~ /gen/ and !$have_gen_lines) {
			$have_gen_lines = 1;
			if ( $todo eq '.' ) {
			    # nothing to do
			    print OUT "\tgen = NULL; if (0) {"; 
			    last;
			}
			else {
			    print OUT "\tgen = unur_test_timing(par,5,&time_setup,&time_sample,0);\n"; 
			    print OUT "\tif (gen) {\n";
			}
		    }
		    
		    print OUT "$l\n";
		}
		
		if ($have_gen_lines) {
		    print OUT "\t;}\n";
		}
		else {
		    if ( $todo eq '.' ) {
			# nothing to do
			print OUT "\tgen = NULL;\n"; }
		    else {
			print OUT "\tgen = unur_test_timing(par,5,&time_setup,&time_sample,0);\n"; 
			print OUT "\ttiming_result[$n_gen] = time_sample;\n";
		    }
		}
		print OUT "\tif (gen) timing_result[$n_gen] = time_sample;\n";
		print OUT "\tunur_free(gen);\n\n";

		# increment counter for generator
		++$n_gen;
	    }
	    # print result of timings 
	    print OUT "\tprint_timing_results( TESTLOG, 0, distr\[$n_distr\], timing_result, $n_generators );\n\n";
	}
	
	print OUT "}\n";

    }

    ## end ##

    print OUT "\n\t/* test finished */\n";
    print OUT "\tif (n_tests_failed>0) n_tests_failed--;  /* we accept one failure */\n";
    print OUT "\ttest_ok &= (n_tests_failed) ? 0 : 1;\n";
    print OUT "\t(n_tests_failed) ? printf(\" --> failed] \") : printf(\" --> ok] \");\n";
    print OUT "\n} /* end of test_validate */\n\n";

    # return name of next section ...
    return $next_section;
} # end of scan_validate()

############################################################
#                                                          #
#  subroutine: print test command                          #
#                                                          #
############################################################

sub print_test_command {
    my $test_command = $_[0];
    my $last_C_line = $_[1];

  SWITCH: {
      if ($test_command =~ /^\s*none\s*/ ) {
	  $test_command =~ s/\s+//g;
	  print OUT "$last_C_line\;\n";
	  last SWITCH;
      }
      if ($test_command =~ /^\s*expected_NULL\s*/ or 
	  $test_command =~ /^\s*expected_setfailed\s*/ or 
	  $test_command =~ /^\s*expected_INFINITY\s*/ or 
	  $test_command =~ /^\s*expected_reinit\s*/) {
	  $test_command =~ s/\s+//g;
	  print OUT "n_tests_failed += check_$test_command\( TESTLOG, $INPUT_LINE_NUMBER, ($last_C_line) )\;\n";
	  last SWITCH;
      }
      if ($test_command =~ /^\s*compare_double_sequence_par\s*$/ or
	  $test_command =~ /^\s*compare_double_sequence_par_start\s*$/ or
	  $test_command =~ /^\s*compare_int_sequence_par\s*$/ or
	  $test_command =~ /^\s*compare_int_sequence_par_start\s*$/ ) {
	  $test_command =~ s/\s+//g;
	  print OUT "$last_C_line\;\n";
	  print OUT "n_tests_failed += $test_command\( TESTLOG, $INPUT_LINE_NUMBER, urng, par, COMPARE_SAMPLE_SIZE );\n";
	  last SWITCH;
      }
      if ($test_command =~ /^\s*compare_double_sequence_urng_start\s*$/ ) {
	  $test_command =~ s/\s+//g;
	  print OUT "$last_C_line\;\n";
	  print OUT "$test_command\( TESTLOG, $INPUT_LINE_NUMBER, urng, COMPARE_SAMPLE_SIZE );\n";
	  last SWITCH;
      }
      if ($test_command =~ /^\s*run_verify_generator\s*$/) {
	  print OUT "$last_C_line\;\n";
	  print OUT "$test_command( TESTLOG,$INPUT_LINE_NUMBER, par );\n";
	  last SWITCH;
      }

      # otherwise
      die "Unknown test command: \"$test_command\"";
  }

} # end print_test_command() 
		

############################################################
#                                                          #
#  subroutine: print common C functions                    #
#                                                          #
############################################################

sub print_C_prototypes {

    print OUT <<EOM;
/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************/
    
/**
 ** Tests for $METHOD
 **/
    
/*---------------------------------------------------------------------------*/
#include "testunuran.h"
/*---------------------------------------------------------------------------*/
#ifdef T_$METHOD
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

static struct prng *urng;           /* uniform random number generator       */

static FILE *TESTLOG;               /* test log file                         */
static FILE *UNURANLOG;             /* unuran log file                       */

static int test_ok = TRUE;          /* all tests ok (boolean)                */
static int n_tests_failed;          /* number of failed tests                */

/*---------------------------------------------------------------------------*/

void test_new( void );
void test_set( void );
void test_get( void );
void test_chg( void );
void test_init( void );
void test_reinit( void );
void test_sample( void );
void test_validate( void );

void run_verify_generator( FILE *LOG, int line, UNUR_PAR *par );

int unur_$method\_set_verify( UNUR_PAR *par, int );

/*---------------------------------------------------------------------------*/

$C_header_aux

EOM

} # end of print_C_prototypes()

############################################################

sub print_C_routines {

    print OUT <<EOM;


/*---------------------------------------------------------------------------*/
/* run generator in verifying mode */

void run_verify_generator( FILE *LOG, int line, UNUR_PAR *par )
{
	UNUR_GEN *gen;
	int i;

	/* switch to verifying mode */
	unur_$method\_set_verify(par,1);

	/* initialize generator */
	gen = unur_init( par ); abort_if_NULL(LOG, line, gen);

	/* run generator */
	for (i=0; i<VIOLATE_SAMPLE_SIZE; i++)
		unur_sample_cont(gen);

	/* destroy generator */
	unur_free(gen); 

} /* end of run_verify_generator() */

EOM

} # end of print_C_routines()

############################################################

sub add_unur_set_verify_routine {
    my $verify = "unur_$method\_set_verify";
    undef my $have_found;

    open (H,"$unuran_h_file")    or die "Cannot open file $unuran_h_file for reading";
    while (<H>) {
	$have_found = 1 if /$verify/;
    }

    unless ($have_found) {
	print OUT "int $verify(UNUR_PAR *par, int verify) {return 0;}\n";
    }
}

############################################################
