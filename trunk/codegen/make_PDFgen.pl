#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

use strict;
my $DEBUG = 0;

# ----------------------------------------------------------------
# Languages

my @languages = ('C', 'FORTRAN', 'JAVA');

# ----------------------------------------------------------------
# Directory with sources

my $top_srcdir = $ENV{'srcdir'} ? $ENV{'srcdir'} : '.';
$top_srcdir .= "/..";

# ----------------------------------------------------------------

require "$top_srcdir/scripts/read_PDF.pl";

# ----------------------------------------------------------------

# C file for code generator  
#my $PDFgen_file = "PDFgen.c";

# ----------------------------------------------------------------
# List of distributions
my $DISTR = read_PDFdata( $top_srcdir );

# For description of data fields in this list see file `readPDF.pl'.

# ................................................................

# Print result on screen
    if ($DEBUG) {
	print_data($DISTR);
    }

# ----------------------------------------------------------------
# Read template C file from STDIN and insert C code for string interpreter 
#
while ( <STDIN> ){

    unless (/^s*=INPUT\s*(\w+)\s+(\w*)/){
	print $_;
    }

    else {
	my $key   = $1;    # INPUT key value
	my $value = $2;

	if ( $key eq "prototypes" ){
	    print make_PDFgen_prototypes($DISTR);
	}
	elsif ( $key eq "distributions" ) {
	    print make_list_of_distributions($DISTR,$value);
	}
	elsif ( $key eq "PDF" ) {
	    print make_PDFgen($DISTR,$value);
	}
	else {
	    die "Error: unknown qualifier after =INPUT: $key\n";
	}
    }
}

# ----------------------------------------------------------------
# End

exit 0;


# ----------------------------------------------------------------
#
# Subroutines
#
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# Make list of prototypes

sub make_PDFgen_prototypes
{
    my $DISTR = $_[0];   # data for distributions

    my %prototypes;

    # Continuous distributions
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	# Function prototype
	foreach my $l (@languages) {
	    $prototypes{$l} .= "static ".get_PDFgen_funct($d,$l).";\n";
	}
    }

    my $gencode;
    foreach my $l (@languages) {
	$gencode .= $prototypes{$l}."\n";
    }

    # End
    return $gencode;

} # end of make_PDFgen_prototypes()


# ----------------------------------------------------------------
# Make switch list for distributions

sub make_list_of_distributions
{
    my $DISTR    = $_[0];   # data for distributions
    my $language = $_[1];   # chosen programming language

    my $gencode;            # code for switch list

    # Switch (distribution)
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$gencode .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";
	$gencode .= "\t\treturn _unur_acg_$language\_PDF_$d (distr,out,pdf);\n";
    }

    return $gencode;

} # end of make_list_of_distributions()


# ----------------------------------------------------------------
# Make routines for PDF code generator (C version)

sub make_PDFgen
{
    my $DISTR    = $_[0];   # data for distributions
    my $language = $_[1];   # chosen programming language

    my $gencode;            # code for switch list

    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

        # Function for distribution
        $gencode .= eval("make_PDF_distr_$language(\$DISTR,\$d)");
    }

    return $gencode;

} # end of make_PDF_main_C()


# ----------------------------------------------------------------
# Get name of generation routine

sub get_PDFgen_funct
{
    my $d = $_[0];       # distribution
    my $l = $_[1];       # chosen programming language

    return "int _unur_acg_$l\_PDF_$d (UNUR_DISTR *distr, FILE *out, const char *pdf)";
} # end of get_PDFgen_funct()


# ----------------------------------------------------------------
# Print PDF code (C version)

sub make_PDF_distr_C
{
    my $DISTR = $_[0];   # data for distributions
    my $d = $_[1];       # distribution

    my $gencode;         # code for creating PDFs

    # Mark begin of distribution
    $gencode .= make_bar("PDF $d: ".$DISTR->{$d}->{"=NAME"});

    # Function header
    $gencode .= get_PDFgen_funct($d,'C')."\n{\n";

    # Print a short description of PDF
    $gencode .= 
	"\t_unur_acg_C_print_section_title(out,\"PDF for $d distribution.\");\n\n";

    # include header files
    $gencode .=
	"\tfprintf (out,\"#include <math.h>\\n\\n\");\n";	

    # compose PDF name
    $gencode .= 
	"\tfprintf (out,\"static double %s (double x)\\n{\\n\",".
        " ((pdf) ? pdf : \"pdf_$d\") );\n";

    # Constants (parameters)
    $gencode .= 
	"\tfprintf (out,\"\\t/* parameters for PDF */\\n\");\n";

    #   List of parameters
    $gencode .= make_PDF_params_C($DISTR,$d);

    #   Normalization constant
    $gencode .= make_PDF_normconstant_C($DISTR,$d);

    # Body of PDF
    $gencode .= make_PDF_body_C($DISTR,$d);

    # End of function
    $gencode .= "\tfprintf (out,\"}\\n\");\n";

    $gencode .= "\n\treturn 1;\n";
    $gencode .= "}\n";

    return $gencode;
    
} # make_PDF_distr_C()


# ----------------------------------------------------------------
# Print PDF code (FORTRAN version)

sub make_PDF_distr_FORTRAN
{
    my $DISTR = $_[0];   # data for distributions
    my $d = $_[1];       # distribution

    my $gencode;         # code for creating PDFs

    # Mark begin of distribution
    $gencode .= make_bar("PDF $d: ".$DISTR->{$d}->{"=NAME"});

    # Function header
    $gencode .= get_PDFgen_funct($d,'FORTRAN')."\n{\n";

    # compose PDF name
    $gencode .= "\tchar pdf_name[7];\n";
    $gencode .= "\tsprintf (pdf_name,\"%.6s\",((pdf) ? pdf : \"f$d\") );\n\n";

    # Print a short description of PDF
    $gencode .= 
	"\t_unur_acg_FORTRAN_print_section_title(out, \"PDF for $d distribution.\");\n\n";

    # print PDF function name
    $gencode .= 
	"\tfprintf (out,\"      DOUBLE PRECISION FUNCTION %.6s(xin)\\n\\n\", pdf_name);\n";
	
    # Declarations
    $gencode .= 
	"\tfprintf (out,\"      IMPLICIT DOUBLE PRECISION (A-Z)\\n\");\n";

    # Define pow()
    $gencode .= 
	"\tfprintf (out,\"      pow(base,exponent) = base**exponent\\n\");\n";

    # Constants (parameters)
    $gencode .= 
	"\tfprintf (out,\"C\\n\");\n".
	"\tfprintf (out,\"C     parameters for PDF\\n\");\n".
	"\tfprintf (out,\"C\\n\");\n";

    #   List of parameters
    $gencode .= make_PDF_params_FORTRAN($DISTR,$d);

    #   Normalization constant
    $gencode .= make_PDF_normconstant_FORTRAN($DISTR,$d);

    # Body of PDF
    $gencode .= make_PDF_body_FORTRAN($DISTR,$d);

    # End of function
    $gencode .= "\tfprintf (out,\"\\n\");\n";
    $gencode .= "\tfprintf (out,\"      END\\n\");\n";
    $gencode .= "\tfprintf (out,\"\\n\");\n";

    $gencode .= "\n\treturn 1;\n";
    $gencode .= "}\n";

    return $gencode;
    
} # make_PDF_distr_FORTRAN()


# ----------------------------------------------------------------
# Print PDF code (JAVA version)

sub make_PDF_distr_JAVA
{
    my $DISTR = $_[0];   # data for distributions
    my $d = $_[1];       # distribution

    my $gencode;         # code for creating PDFs

    # Mark begin of distribution
    $gencode .= make_bar("PDF $d: ".$DISTR->{$d}->{"=NAME"});

    # Function header
    $gencode .= get_PDFgen_funct($d,'JAVA')."\n{\n";

    # Print a short description of PDF
    $gencode .= 
	"\t_unur_acg_JAVA_print_section_title(out, \"PDF for $d distribution.\");\n\n";

    # Constants (parameters)
    $gencode .= 
	"\tfprintf (out,\"\\t/* parameters for PDF */\\n\");\n";

    #   List of parameters
    $gencode .= make_PDF_params_JAVA($DISTR,$d);

    #   Normalization constant
    $gencode .= make_PDF_normconstant_JAVA($DISTR,$d);

    # compose PDF name
    $gencode .= 
	"\tfprintf (out,\"\\tstatic double %s (double x)\\n\\t{\\n\",".
        " ((pdf) ? pdf : \"pdf_$d\") );\n";

    # Body of PDF
    $gencode .= make_PDF_body_JAVA($DISTR,$d);

    # End of function
    $gencode .= "\tfprintf (out,\"\\t}\\n\");\n";

    $gencode .= "\n\treturn 1;\n";
    $gencode .= "}\n";

    return $gencode;
    
} # make_PDF_distr_JAVA()


# ----------------------------------------------------------------
# Process parameter list (C version)

sub make_PDF_params_C
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
# Process parameter list (FORTRAN version)

sub make_PDF_params_FORTRAN
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
		"\tif (".$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".n_params > $i) {\n".
		"\t\tfprintf (out,\"      PARAMETER (".
	        $in_params->[$i]." = \");\n".
		"\t\t_unur_acg_FORTRAN_print_double(out,".
	        $DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".params[$i]);\n".
	        "\t\tfprintf (out,\")\\n\");\n\t}\n";
	}

	else {
	    $params .= 
		"\tfprintf (out,\"      PARAMETER (".
	        $in_params->[$i]." = \");\n".
		"\t_unur_acg_FORTRAN_print_double(out,".
	        $DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".params[$i]);\n".
	        "\tfprintf (out,\")\\n\");\n";
	}
    }

    return $params;

} # end of make_PDF_params_FORTRAN()



# ----------------------------------------------------------------
# Process parameter list (JAVA version)

sub make_PDF_params_JAVA
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
		"\t\tfprintf (out,\"\\tprivate static final double ".
		$in_params->[$i].
		" = %.20e;\\n\",".
		$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".params[$i]);\n";
	}

	else {
	    $params .= 
		"\tfprintf (out,\"\\t\\tstatic final double ".
	        $in_params->[$i].
	        " = %.20e;\\n\",".
	        $DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".params[$i]);\n";
	}
    }

    return $params;

} # end of make_PDF_params_JAVA()


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
# Process normalization constants (FORTRAN version)

sub make_PDF_normconstant_FORTRAN
{
    my $DISTR = $_[0];    # data for distributions
    my $distr = $_[1];    # name of distribution

    if ($DISTR->{$distr}->{"=PDF"}->{"=CONST"}) {
	if ($DISTR->{$distr}->{"=PDF"}->{"=BODY"} =~ /((LOG)?NORMCONSTANT)/) {
	    my $norm = ($1 =~ /LOG/) ? "lncnst" : "const";
	    return 
		"\tfprintf (out,\"      PARAMETER ($norm = \");\n".
		"\t_unur_acg_FORTRAN_print_double(out,".$DISTR->{$distr}->{"=PDF"}->{"=CONST"}.");\n".
	        "\tfprintf (out,\")\\n\");\n";
	}
    }

    # else
    return "";

} # end of make_PDF_normconstant_FORTRAN()


# ----------------------------------------------------------------
# Process normalization constants (JAVA version)

sub make_PDF_normconstant_JAVA
{
    my $DISTR = $_[0];    # data for distributions
    my $distr = $_[1];    # name of distribution

    if ($DISTR->{$distr}->{"=PDF"}->{"=CONST"}) {
	if ($DISTR->{$distr}->{"=PDF"}->{"=BODY"} =~ /((LOG)?NORMCONSTANT)/) {
	    return 
		"\tfprintf (out,\"\\tprivate static final double $1 = %.20e;\\n\\n\",".
		$DISTR->{$distr}->{"=PDF"}->{"=CONST"}.");\n";
	}
    }

    # else
    return "";

} # end of make_PDF_normconstant_JAVA()


# ----------------------------------------------------------------
# Process PDF body (C version)

sub make_PDF_body_C
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

} # end of make_PDF_body_C()


# ----------------------------------------------------------------
# Process PDF body (FORTRAN version)

sub make_PDF_body_FORTRAN
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
    my $body = 
	"\tfprintf (out,\"C\\n\");\n".
	"\tfprintf (out,\"C     compute PDF\\n\");\n".
	"\tfprintf (out,\"C\\n\");\n";

    # Make a copy of the argument
    $body .=
	"\tfprintf (out,\"      x = xin\\n\\n\");\n";

    # translate C -> FORTRAN
    my $is_if_block = 0;
    my $next_is_if_block = 0;

    foreach my $l (split /\n/, $in_body) {

	# if (n_params ...) statement
	if ($l =~ /if\s*\(\s*n_params\s*[<>=!]+\s*\d+\s*\)/) {
	    $l =~ s/n_params/$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.n_params/;
	    $l =~ s/  /\t/g;
	    $body .= "$l\n\t";
	    next;
	}

	# remove variable declarations 
	next if $l =~ /^\s*(register)?\s*(double|int)/;

	my $indent = ($is_if_block) ? "        " : "     ";
	my $return = "";

	# rename normalization constant
	$l =~ s/LOGNORMCONSTANT/lncnst/g;
	$l =~ s/NORMCONSTANT/const/g;

	# double constants
	$l =~ s/([\d\.]+)/$1D0 /g;

	# FORTRAN commands
	$l =~ s/<=/.LE./g;
	$l =~ s/>=/.GE./g;
	$l =~ s/</.LT./g;
	$l =~ s/>/.GT./g;

	$l =~ s/==/.EQ./g;
	$l =~ s/!=/.NE./g;

	$l =~ s/\&\&/.AND./g;
	$l =~ s/\|\|/.OR./g;

	$l =~ s/fabs\(([^\)]+)\)/ABS($1)/g;

	# return statement
	if ($l =~ /return\s+(.*);/) {
	    my $val = $1;
	    $val =~ s/^\s*/$indent %.6s = /g;
	    if (length $val > 71) {
		substr($val,71,0) = "\\n     #  ";
	    }
	    $body .= "\tfprintf (out,\"$val\\n\", pdf_name);\n";
	    $return = "\tfprintf (out,\"$indent RETURN\\n\");\n";
	}

	else {
	    # remove `;' at and of line
	    $l =~ s/;\s*$//;

	    # `if' statement ?
	    if ($l =~ /^\s*if\W/) {
		$l =~ s/^\s*if/IF/;
		$l .= " THEN";
		$next_is_if_block = 1;
	    }

	    # change beginning of line
	    $l =~ s/^\s*/$indent /g;

	    # print
	    if (length $l > 71) {
		substr($l,71,0) = "\\n     #  ";
	    }
	    $body .= "\tfprintf (out,\"$l\\n\");\n";
	}

	# return statement (if there is one)
	$body .= $return;

	# inside if block ?
	#       we assume that each if block consists of one line
	if ($is_if_block) {
	    $body .= "\tfprintf (out,\"      ENDIF\\n\");\n";
	}
	$is_if_block = $next_is_if_block;
	$next_is_if_block = 0;

    }

    return $body;

} # end of make_PDF_body_FORTRAN()


# ----------------------------------------------------------------
# Process PDF body (JAVA version)

sub make_PDF_body_JAVA
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
    my $body = "\tfprintf (out,\"\\t\\t/* compute PDF */\\n\");\n";

    foreach my $l (split /\n/, $in_body) {
	if ($l =~ /if\s*\(\s*n_params\s*[<>=!]+\s*\d+\s*\)/) {
	    $l =~ s/n_params/$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.n_params/;
	    $l =~ s/  /\t/g;
	    $body .= "$l\n\t";
	    next;
	}
	$l =~ s/  /\\t/g;
       
        #remove declaratien "register"
        $l =~ s/register//;

	# modify mathematical functions
	$l =~ s/(\W)(exp|log|pow)(\W)/$1Math\.$2$3/g;
	$l =~ s/(\W)fabs(\W)/$1Math\.abs$2$3/g;

	$body .= "\tfprintf (out,\"\\t$l\\n\");\n";
    }

    return $body;

} # end of make_PDF_body_JAVA()


# ----------------------------------------------------------------
# Make bar in PDFgen file

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
