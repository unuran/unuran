#!/usr/bin/perl
# ----------------------------------------------------------------

use File::Find;

# ----------------------------------------------------------------

$DEBUG = 0;

# ----------------------------------------------------------------

# Start the search for the files in these directories
my @Startdirs = ('..');

# C file for code generator  
my $PDFgen_C_file = "PDFgen.c";
my $PDFgen_H_file = "PDFgen_source.h";

# Configuration file for tests
my $test_conf_file = "test.conf";

# C file for making code generator tests
my $make_test_PDFgen = "make_test_PDFgen.c";

# C file for tests
my $test_PDFgen = "test_PDFgen.c";

# Sample size for test
my $SAMPLE_SIZE = 10000;

# Header file for UNURAN standard distributions
my $h_stddistr = "unuran_distributions.h";

# List of distribution types
%distr_types =
    ( "CONT"  => { "file_prefix" => "c",
		   "PDF_prefix"  => "_unur_pdf_",
		   "PDF_type"    => "double" },

      "DISCR" => { "file_prefix" => "d",
		   "PDF_prefix"  => "_unur_pmf_",
		   "PDF_type"    => "double" } );

# ----------------------------------------------------------------
# List of distributions
my $DISTR;

# Description of data fields (with beta distribution as example)
#
#   $DISTR->{"beta"}                           ... entries for distribution "beta"
#
#   $DISTR->{"beta"}->{"=NAME"}                ... name of distribution
#   $DISTR->{"beta"}->{"=TYPE"}                ... type of distribution (CONT|DISCR)
#   $DISTR->{"beta"}->{"=FILE"}                ... file name + path for C file
#   $DISTR->{"beta"}->{"=ID"}                  ... id of distribution
#
#   $DISTR->{"beta"}->{"=DOC"}                 ... documentation for distribution
#   $DISTR->{"beta"}->{"=DOC"}->{"=PDF"}       ... formula for PDF
#   $DISTR->{"beta"}->{"=DOC"}->{"=CONST"}     ... normalization constant for PDF
#   $DISTR->{"beta"}->{"=DOC"}->{"=DOMAIN"}    ... domain for PDF
#   $DISTR->{"beta"}->{"=DOC"}->{"=FPARAM"}    ... list of parameters with constraints
#      Remark: There exist other fields which are not relevant here
#              (see src/distributions/unuran_distributions.h).
#
#   $DISTR->{"beta"}->{"=PDF"}                 ... PDF of distribution
#   $DISTR->{"beta"}->{"=PDF"}->{"=NAME"}      ... name of PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=RTYPE"}     ... return type of PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=ARGS"}      ... list of arguments for PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=N_PARAMS"}  ... number of parameters for PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=PARAMS"}[2] ... parameter #2 for PDF (starting at 0)
#   $DISTR->{"beta"}->{"=PDF"}->{"=BODY"}      ... function body of PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=CONST"}     ... macro expansion for (LOG)NORMCONSTANT
#   $DISTR->{"beta"}->{"=PDF"}->{"=DISTR"}     ... macro expansion for DISTR
#

# ----------------------------------------------------------------
# List of files
my %file_list;
find (\&find_files, @Startdirs);

# ----------------------------------------------------------------

# Scan header file for UNURAN standard distributions
# and get a list of valid distributions
scan_h_stddistr( $file_list{$h_stddistr} );

# ----------------------------------------------------------------

# Print result on screen
if ($DEBUG) {
    print_data();
}

# ----------------------------------------------------------------

# Make C file with code generator
make_PDFgen_files();

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
# Finds recursively all C- or H-files containing the relevant
# informations. 

sub find_files 
{
    if ($File::Find::name =~ /(.*\/)?(.+\.[ch])$/ ) {
	my $file_path =  $File::Find::name;
	my $file_name =  $2;
	$file_list{$file_name} = $file_path;
    }
} # end of find_files()


# ----------------------------------------------------------------
# Scan header file for UNURAN standard distributions

sub scan_h_stddistr
{
    my $file = $_[0];
    
    # Read header file
    open HFILE, $file or die "cannot open file $file\n";
    my $file_content;
    while (<HFILE>) {
	chomp;
	$_ =~ s/^\s*(.*)\s*$/$1/;
	$file_content .= "$_\n";
    }
    close HFILE;

    # Split into sections
    my @sections = split /=EON/, $file_content;

    # Scan sections
    foreach $s (@sections) {
	# Remove trailing part of sections
	(my $dummy, $s) = split /=DISTR/, $s, 2;
	next unless $s;

	# Get name of distribution
	$s =~ s/(\w+)\s+(.+)\n//;
	my $distr = $1;
	$DISTR->{$distr}->{"=NAME"} = $2;

	# Scan distribution section
	scan_distr($distr,$s);

	# Get type of distribution and path of distribution file
	get_distr_file($distr);

	# Read distribution file
	read_distr_file($distr);

    }
    
} # end of scan_h_stddistr()

# ----------------------------------------------------------------
# Scan distribution

sub scan_distr
{
    my $distr = $_[0];
    my $distr_text = $_[1];

    # add =END TAG to text (for convenience)
    $distr_text .= "\n=END";

    # split into lines
    my @lines = split /\n/, $distr_text;

    # scan all TAGs (node sections)
    my $this_TAG = "=END";  # add =END tag to node TAG

    foreach my $l (@lines) {
	# next TAG ?
	if ($l =~ /^\s*(\/\*)?\s*(=[A-Z]+)\s*(.*)$/) {
	    # store next TAG
	    $this_TAG = $2;
	    # save rest of line
	    $l = $3;
	}

	# append to stored lines
	# (except for =END TAG)
	unless ($this_TAG eq "=END") {
	    $DISTR->{$distr}->{"=DOC"}->{$this_TAG} .= $l."\n";
	}
    }

} # end of scan_distr() 

# ----------------------------------------------------------------
# Get type of distribution and path of distribution file

sub get_distr_file
{
    my $distr = $_[0];

    my $found = 0;

    foreach my $type (keys %distr_types) {
	next unless $type;
	my $file_name = $distr_types{$type}{"file_prefix"}."\_$distr\.c";
	if ($file_list{$file_name}) {
	    $found = 1;
	    $DISTR->{$distr}->{"=FILE"} = $file_list{$file_name};
	    $DISTR->{$distr}->{"=TYPE"} = $type;
	    last;
	}
    }

    die "Cannot find file for $distr" unless $found;

} # end of get_distr_file()

# ----------------------------------------------------------------
# Read distribution file

sub read_distr_file
{
    my $distr = $_[0];

    # Read file
    my $file = $DISTR->{$distr}->{"=FILE"};
    open CFILE, $file or die "cannot open file $file\n";
    my $file_content;
    while (<CFILE>) {
	$file_content .= $_;
    }
    close CFILE;

    # Check distribution name
    $file_content =~ /static\s+const\s+char\s+distr\_name\s*\[\s*\]\s*=\s*\"$distr\"/
	or die "$distr: distr_name inconsistent";

    # Type of distribution
    $type = $DISTR->{$distr}->{"=TYPE"};

    # Get PDF source
    my $PDF_name = $distr_types{$type}{"PDF_prefix"}.$distr;
    my $PDF_pattern = 
	"(int|double)\\s+"                   # $1: return type  
	.$PDF_name                           #     name of function
        ."\\s*\\(([^\\)]*)\\)\\s*"           # $2: arguments of function
	."([^;])"                            # $3: first character (to distinguish from prototype) 
	."(.*)"                              # $4: function body
        ."\\/\\*\\s+end\\s+of\\s+$PDF_name"; # end of function marker
    $file_content =~ /$PDF_pattern/s
	or die "cannot find PDF for $distr";

    # Store data
    $DISTR->{$distr}->{"=PDF"}->{"=NAME"}  = $PDF_name;  # name of PDF function
    $DISTR->{$distr}->{"=PDF"}->{"=RTYPE"} = $1;         # return type 
    $DISTR->{$distr}->{"=PDF"}->{"=ARGS"}  = $2;         # arguments for function
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  = $3.$4;      # function body

    # Modify function arguments:
    #   remove DISTR from argument list
    $DISTR->{$distr}->{"=PDF"}->{"=ARGS"}  =~
	s /\,\s*(UNUR_DISTR|struct unur_distr)\s*\*\s*distr\s*//;

    # Modify function body:
    #   remove comments
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s {/\*.*?\*/} []gsx;
    #   remove enclosing brackets
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s /^\s*\{(.*)\}\s*$/$1/s;
    #   remove empty lines
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s /\n\s*\n/\n/gx;
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s /^\s*\n//;
    #   remove declaration of params
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ 
	s /\s*(register)?\s+double\s*\*\s*params\W.*\n//;
    #   remove all `DISTR.' from body
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s /DISTR\.//g;

    # Get parameters for PDF
    $file_content =~ s {/\*.*?\*/} []gsx;    # remove all comments
    my @lines = split /\n/, $file_content;
    my $n_params = -1;
    foreach my $l (@lines) {
	next unless $l =~ /\#define\s+(\w+)\s+(.*)$/;
	my $macro_name = $1;
	my $macro_body = $2;

	if ($macro_body =~ /params\s*\[(\d)+\]/) {
	    $n_params = ($n_params < $1) ? $1 : $n_params;
	    $DISTR->{$distr}->{"=PDF"}->{"=PARAMS"}[$1] = $macro_name;
	    next;
	}
	
	if ($macro_name =~ /(.*NORMCONSTANT)/) {
	    $DISTR->{$distr}->{"=PDF"}->{"=CONST"} = $macro_body;
	    next;
	}

	if ($macro_name =~ /(DISTR)/) {
	    $DISTR->{$distr}->{"=PDF"}->{"=DISTR"} = $macro_body;
	    next;
	}
    }

    # Number of parameters for PDF
    $DISTR->{$distr}->{"=PDF"}->{"=N_PARAMS"} = $n_params+1;

    # Id of distribution
    $file_content =~ /distr\-\>id\s*=\s*(\w+)/
	or die "cannot find ID for $distr";
    $DISTR->{$distr}->{"=ID"} = $1;
	
} # end of read_distr_file()

# ----------------------------------------------------------------
# Print data on screen (for debugging)

sub print_data
{
    foreach my $d (keys %{$DISTR}) {
	print "-------------------------------------\n";
	print "distribution = \"$d\"\n";

	print "NAME: ".$DISTR->{$d}->{"=NAME"}."\n";
	print "TYPE: ".$DISTR->{$d}->{"=TYPE"}."\n";
	print "ID  : ".$DISTR->{$d}->{"=ID"}."\n";
	print "FILE: ".$DISTR->{$d}->{"=FILE"}."\n\n";

	print "DOC: PDF    : ".$DISTR->{$d}->{"=DOC"}->{"=PDF"};
	print "DOC: CONST  : ".$DISTR->{$d}->{"=DOC"}->{"=CONST"};
	print "DOC: DOMAIN : ".$DISTR->{$d}->{"=DOC"}->{"=DOMAIN"};
	print "DOC: FPARAM :\n".$DISTR->{$d}->{"=DOC"}->{"=FPARAM"}."\n";

	print "PDF: NAME  : ".$DISTR->{$d}->{"=PDF"}->{"=NAME"}."\n";
	print "PDF: RTYPE : ".$DISTR->{$d}->{"=PDF"}->{"=RTYPE"}."\n";
	print "PDF: ARGS  : ".$DISTR->{$d}->{"=PDF"}->{"=ARGS"}."\n";

	print "PDF: PARAM : ".$DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"}."\n";
	foreach my $i (0 .. $DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"} - 1) {
	    print "\t[$i]: ".$DISTR->{$d}->{"=PDF"}->{"=PARAMS"}[$i]."\n"
	}

	print "PDF: CONST : ".$DISTR->{$d}->{"=PDF"}->{"=CONST"}."\n";
	print "PDF: DISTR : ".$DISTR->{$d}->{"=PDF"}->{"=DISTR"}."\n";

	print "PDF: BODY :\n".$DISTR->{$d}->{"=PDF"}->{"=BODY"}."\n";
    }

} # end of print_data()

# ----------------------------------------------------------------
# Make C file with code generator

sub make_PDFgen_files
{
    my $PDFgen;
    my $PDFgen_static_prototypes;
    my $PDFgen_global_prototypes;
    

    my $empty_line = "\tfprintf (out,\"\\n\");\n\n";

# ................................................................
# Main

    # Mark begin of Main
    $PDFgen .= "/* ------------------------------------- */\n";
    $PDFgen .= "/* Main */\n\n";
    
    # Function name
    my $function = "int \_unurgen\_C\_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf)";
    
    # Function prototype
    $PDFgen_global_prototypes .= "$function;\n";
    
    # Function header
    $PDFgen .= "$function\n{\n";

    # Check for invalid NULL pointer
    $PDFgen .= "\t_unur_check_NULL(\"unurgen\", distr, 0 );\n";

    # Switch
    $PDFgen .= "\tswitch (distr->id) {\n";
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	$PDFgen .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";
	$PDFgen .= "\t\treturn \_unurgen\_C\_PDF_$d (distr,out,pdf);\n";
    }
    $PDFgen .= "\tdefault:\n";
    $PDFgen .= "\t\treturn 0;\n";
    $PDFgen .= "\t}\n";
    
    # End of function
    $PDFgen .= "}\n\n";

    # Mark end of Main
    $PDFgen .= "/* end of Main */\n";
    $PDFgen .= "/* ------------------------------------- */\n";

# ................................................................
# Continuous distributions
    
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	# Mark begin of distribution
	$PDFgen .= "/* ------------------------------------- */\n";
	$PDFgen .= "/* $d: ".$DISTR->{$d}->{"=NAME"}." */\n\n";

	# Function name
	my $function = "int \_unurgen\_C\_PDF\_$d (UNUR_DISTR *distr, FILE *out, const char *pdf)";

	# Function prototype
	$PDFgen_static_prototypes .= "static $function;\n";

	# Function header
	$PDFgen .= "$function\n{\n";

	# Macro definitions for parameters
	#   (total) number of parameters
	$PDFgen .= "\tfprintf (out,\"\#define  n_params  %d\\n\","
	    .$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".n_params);\n";
	#   parameter list
	foreach my $i (0 .. $DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"} - 1) {
	    $PDFgen .= "\tfprintf (out,\"\#define  "
		.$DISTR->{$d}->{"=PDF"}->{"=PARAMS"}[$i]
		."  (%.20g)\\n\","
	        .$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".params[$i]);\n";
	}
	#   normalization constant
	if ($DISTR->{$d}->{"=PDF"}->{"=CONST"}) {
	    $PDFgen .= "\tfprintf (out,\"\#define  NORMCONSTANT  (%.20g)\\n\","
		.$DISTR->{$d}->{"=PDF"}->{"=CONST"}.");\n";
	    $PDFgen .= "\tfprintf (out,\"\#define  LOGNORMCONSTANT  NORMCONSTANT\\n\");\n";
	}
	$PDFgen .= $empty_line;

	# compose PDF name
	my $PDFname = "\tfprintf (out,\"static ";
	$PDFname .= $DISTR->{$d}->{"=PDF"}->{"=RTYPE"};
	$PDFname .= " %s (";
	$PDFname .= $DISTR->{$d}->{"=PDF"}->{"=ARGS"};
	$PDFname .= ")\\n{\\n\",";
	$PDFname .= " ((pdf) ? pdf : \"pdf_$d\") ";
	$PDFname .= ");\n";

	# Write PDF
	$PDFgen .= $PDFname; 
	foreach my $l (split /\n/, $DISTR->{$d}->{"=PDF"}->{"=BODY"}) {
	    $PDFgen .= "\tfprintf (out,\"$l\\n\");\n";
	}

	$PDFgen .= "\tfprintf (out,\"}\\n\");\n";
	$PDFgen .= $empty_line;

	# Undefine macros
	#   (total) number of parameters
	$PDFgen .= "\tfprintf (out,\"\#undef  n_params\\n\");\n";
	#   parameter list
	foreach my $i (0 .. $DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"} - 1) {
	    $PDFgen .= "\tfprintf (out,\"\#undef  "
		.$DISTR->{$d}->{"=PDF"}->{"=PARAMS"}[$i]."\\n\");\n";
	}
	#   normalization constant
	if ($DISTR->{$d}->{"=PDF"}->{"=CONST"}) {
	    $PDFgen .= "\tfprintf (out,\"\#undef  NORMCONSTANT\\n\");\n";
	    $PDFgen .= "\tfprintf (out,\"\#undef  LOGNORMCONSTANT\\n\");\n";
	}

	# End of function
	$PDFgen .= "\n";
	$PDFgen .= "\treturn 1;\n";
	$PDFgen .= "}\n\n";

	# Mark end of distribution
	$PDFgen .= "/* end of $d */\n";
	$PDFgen .= "/* ------------------------------------- */\n";
    }

# ................................................................
# Print C code int files

    # Header file
    open HFILE, ">$PDFgen_H_file" or die "cannot open file $PDFgen_H_file\n";
    print HFILE $PDFgen_global_prototypes;
    close HFILE;

    # C file
    open CFILE, ">$PDFgen_C_file" or die "cannot open file $PDFgen_C_file\n";
    print CFILE "\#include <source_unuran.h>\n";
    print CFILE "\#include \"PDFgen_source.h\"\n\n";
    print CFILE $PDFgen_static_prototypes;
    print CFILE $PDFgen;
    close CFILE;

} # end of make_code_generator()

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
	$test_file .= "\tfprintf(out,\"\\n\\n\");\n\n";

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
