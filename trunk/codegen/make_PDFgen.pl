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
my $PDFgen_C_file = "PDFgen.c";

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

# Make C file with code generator
make_PDFgen_C();

# ----------------------------------------------------------------
# End

exit 0;

# ----------------------------------------------------------------
#
# Subroutines
#
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Make C file with code generator

sub make_PDFgen_C
{
    my $PDFgen;
    my $PDFgen_prototypes;
    
    my $empty_line = "\tfprintf (out,\"\\n\");\n";

# ................................................................
# Main

    # Mark begin of Main
    $PDFgen .= make_bar_C("PDF main");

    # Function name
    my $function = "int \_unur_acg\_C\_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf)";
    
    # Function header
    $PDFgen .= "$function\n{\n";

    # Check for invalid NULL pointer
    $PDFgen .= "\t_unur_check_NULL(\"ACG\", distr, 0 );\n\n";

    # Switch
    $PDFgen .= "\tswitch (distr->id) {\n";
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	$PDFgen .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";
	$PDFgen .= "\t\treturn \_unur_acg\_C\_PDF_$d (distr,out,pdf);\n";
    }
    $PDFgen .= "\tdefault:\n";
    $PDFgen .= "\t\t_unur_error(distr->name,UNUR_ERR_GEN_DATA,\"Cannot make PDF\");\n";
    $PDFgen .= "\t\treturn 0;\n";
    $PDFgen .= "\t}\n";
    
    # End of function
    $PDFgen .= "}\n";

# ................................................................
# Continuous distributions
    
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	# Mark begin of distribution
	$PDFgen .= make_bar_C("PDF $d: ".$DISTR->{$d}->{"=NAME"});

	# Function name
	my $function = "int \_unur_acg\_C\_PDF\_$d (UNUR_DISTR *distr, FILE *out, const char *pdf)";

	# Function prototype
	$PDFgen_prototypes .= "static $function;\n";

	# Function header
	$PDFgen .= "$function\n{\n";

	# compose PDF name
	my $PDFname = 
	    "\tfprintf (out,\"static ".
	    $DISTR->{$d}->{"=PDF"}->{"=RTYPE"}.
	    " %s (".
	    $DISTR->{$d}->{"=PDF"}->{"=ARGS"}.
	    ")\\n{\\n\",".
	    " ((pdf) ? pdf : \"pdf_$d\") ".
	    ");\n";

	# Constants (parameters)
	my $PDFconst = "\tfprintf (out,\"\\t/* parameters for PDF */\\n\");\n";

	# List of parameters
	$PDFconst .= make_PDF_params($DISTR,$d);

	# Normalization constant
	$PDFconst .= make_PDF_normconstant($DISTR,$d);

	# Make body of PDF
	my $PDFbody = make_PDF_body($DISTR,$d);

	# Write PDF
	$PDFgen .= 
	    "\t_unur_acg_C_print_sectionheader(out, 1, \"PDF for $d distribution.\");\n\n".
	    $PDFname.
	    $PDFconst.
	    $PDFbody;

	$PDFgen .= "\tfprintf (out,\"}\\n\");\n";
	$PDFgen .= $empty_line;

	# End of function
	$PDFgen .= "\n\treturn 1;\n";
	$PDFgen .= "}\n";
    }

# ................................................................
# Print C code into files

    open CFILE, ">$PDFgen_C_file" or die "cannot open file $PDFgen_C_file\n";

    print CFILE make_bar_C("include header file");
    print CFILE "\#include <codegen_source.h>\n";

    print CFILE make_bar_C("local prototypes");
    print CFILE $PDFgen_prototypes;

    print CFILE $PDFgen;
    close CFILE;

} # end of make_PDFgen_C()

# ----------------------------------------------------------------
# Process parameter list

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
} # end of make_PDF_params()

# ----------------------------------------------------------------
# Process normalization constants

sub make_PDF_normconstant
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
} # end of make_PDF_normconstant()

# ----------------------------------------------------------------
# Process PDF body

sub make_PDF_body
{
    my $DISTR = $_[0];    # data for distributions
    my $d = $_[1];        # name of distribution

    # code for PDF as extracted from UNURAN library
    my $in_body = $DISTR->{$d}->{"=PDF"}->{"=BODY"};

    # at the moment we do DO NOT handle strings the following regular expression
    if ($in_body =~ /if\s*\(\s*n_params\s*[<>=!]+\s*\d+\s*\)\s*\{/) {
	die "cannot handle PDF body";
    } 

    # string for make PDF function
    my $body = "\tfprintf (out,\"\\t/* compute PDF */\\n\");\n";

    foreach my $l (split /\n/, $in_body) {
	if ($l =~ /if\s*\(\s*n_params\s*[<>=!]+\s*\d+\s*\)/) {
	    $l =~ s/n_params/$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.n_params/;
	    $l =~ s/  /\t/g;
	    $body .= "$l\n\t";
	    next;
	}
	$l =~ s/  /\\t/g;
	$body .= "\tfprintf (out,\"$l\\n\");\n";
    }

    return $body;
} # end of make_PDF_body()

# ----------------------------------------------------------------
# Make bar in PDFgen file

sub make_bar_C
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
} # end of make_bar_C()

# ----------------------------------------------------------------

