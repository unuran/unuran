#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

$DEBUG = 0;

# ----------------------------------------------------------------

require "read_PDF.pl";

# ----------------------------------------------------------------

# C file for code generator  
my $headergen_file = "headergen.c";

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

# List of prototypes
my $headergen_prototypes;

# Generator codes
my $headergen;

# ................................................................

# Make C version of code generator
$headergen .= make_headergen_C($DISTR);

# ................................................................

# Print code into file
open FILE, ">$headergen_file" or die "cannot open file $headergen_file\n";

print FILE make_bar("include header file");
print FILE "\#include <codegen_source.h>\n";

print FILE $headergen;
close FILE;

# ----------------------------------------------------------------
# End

exit 0;


# ----------------------------------------------------------------
#
# Subroutines
#
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Make routines for header code generator (C version)

sub make_headergen_C
{
    my $DISTR = $_[0];   # data for distributions

    my $gencode;
    
    # C version
    $gencode .= make_bar("C version");

    # Mark begin of Main
    $gencode .= make_bar("Make header");

    # Function header
    $gencode .= "int _unur_acg_C_header (UNUR_DISTR *distr, FILE *out, const char *pdf)\n{\n";

    # Define horizontal rule
    $gencode .= "\tconst char hrule[] = \"/* ---------------------------------------------------------------- */\\n\";\n";

    # Format for printing string
    $gencode .= "\tconst char sformat[] = \"/* %-64.64s */\\n\";\n\n";

    # Check for invalid NULL pointer
    $gencode .= "\t_unur_check_NULL(\"ACG\", distr, 0 );\n\n";

    # Start of header
    $gencode .= "\tfprintf(out,hrule);\n\n";

    # Switch (distribution)
    $gencode .= "\tswitch (distr->id) {\n";

    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	$gencode .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";
	$gencode .= "\t\tfprintf(out,sformat,\"Generator for $d distribution\");\n";
	$gencode .= "\t\tfprintf(out,sformat,\"   with parameters\");\n";

##    #   List of parameters
##    $gencode .= make_PDF_params_C($DISTR,$d);
	
	$gencode .= "\t\tbreak;\n";
    }

    $gencode .= "\tdefault:\n";
    $gencode .= "\t\t_unur_error(distr->name,UNUR_ERR_GEN_DATA,\"Unknown distribution\");\n";
    $gencode .= "\t\treturn 0;\n";
    $gencode .= "\t}\n\n";

    # End of header
    $gencode .= "\tfprintf(out,hrule);\n\n";
    
    # End of function
    $gencode .= "\treturn 1;\n";
    $gencode .= "}\n";

    # End
    return $gencode;

} # end of make_headergen_C()




# ----------------------------------------------------------------
# Process parameter list (C version)

sub make_PDF_params
{
    my $DISTR = $_[0];    # data for distributions
    my $d = $_[1];        # name of distribution

    # whether we have variable number of parameters
    my $have_n_params = ($DISTR->{$d}->{"=PDF"}->{"=BODY"} =~ /n_params/);

    # parameter list
    my $params;

    my $n_in_params = $DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"};
    my $in_params = $DISTR->{$d}->{"=PDF"}->{"=PARAMS"};

    foreach my $i (0 .. $n_in_params - 1) {

	if ($have_n_params) {
	    $params .=
		"\tif (".$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".n_params > $i)\n".
		"\t\tfprintf (out,\"\\tconst double ".
		$in_params->[$i].
		" = %.20e;\\n\",".
		$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".params[$i]);\n";
	}

	else {
	    $params .= 
		"\tfprintf (out,\"\\tconst double ".
	        $in_params->[$i].
	        " = %.20e;\\n\",".
	        $DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".params[$i]);\n";
	}
    }

    return $params;

} # end of make_PDF_params_C()


# ----------------------------------------------------------------
# Process normalization constants (C version)

sub make_PDF_normconstant_C
{
    my $DISTR = $_[0];    # data for distributions
    my $distr = $_[1];    # name of distribution

    if ($DISTR->{$distr}->{"=PDF"}->{"=CONST"}) {
	if ($DISTR->{$distr}->{"=PDF"}->{"=BODY"} =~ /((LOG)?NORMCONSTANT)/) {
	    return 
		"\tfprintf (out,\"\\tconst double $1 = %.20e;\\n\\n\",".
		$DISTR->{$distr}->{"=PDF"}->{"=CONST"}.");\n";
	}
    }

    # else
    return "";

} # end of make_PDF_normconstant_C()


# ----------------------------------------------------------------
# Make bar in headergen file

sub make_bar
{
    my $string = $_[0];

    my $hline = "/* ";
    for (1..64) { $hline .= "-"; }
    $hline .= " */\n";

    for (1..64) { $string .= " "; }
    $string = substr $string, 0, 64;

    my $bar = "\n";
    $bar .= $hline;
    $bar .= "/* $string */\n";
    $bar .= $hline;
    $bar .= "\n";

    return $bar;

} # end of make_bar()

# ----------------------------------------------------------------

