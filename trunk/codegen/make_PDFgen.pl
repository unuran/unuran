#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

$DEBUG = 0;

# ----------------------------------------------------------------
# Directory with sources

my $top_srcdir = $ENV{'srcdir'} ? $ENV{'srcdir'} : '.';
$top_srcdir .= "/..";

# ----------------------------------------------------------------

require "$top_srcdir/scripts/read_PDF.pl";

# ----------------------------------------------------------------

# C file for code generator  
my $PDFgen_file = "PDFgen.c";

# ----------------------------------------------------------------
# List of distributions
my $DISTR = read_PDFdata('..');

# For description of data fields in this list see file `readPDF.pl'.

# ................................................................

# Print result on screen
    if ($DEBUG) {
	print_data($DISTR);
    }

# ----------------------------------------------------------------

# List of prototypes
my $PDFgen_prototypes;

# Generator codes
my $PDFgen;

# ................................................................

# Make C version of code generator
$PDFgen .= make_PDFgen_C($DISTR);

# Make FORTRAN version of code generator
$PDFgen .= make_PDFgen_FORTRAN($DISTR);

# Make JAVA version of code generator
$PDFgen .= make_PDFgen_JAVA($DISTR);

# ................................................................

# Print code into file
open FILE, ">$PDFgen_file" or die "cannot open file $PDFgen_file\n";

print FILE make_bar("include header file");
print FILE "\#include <codegen_source.h>\n";

print FILE make_bar("local prototypes");
print FILE $PDFgen_prototypes;

print FILE $PDFgen;
close FILE;

# ----------------------------------------------------------------
# End

exit 0;


# ----------------------------------------------------------------
#
# Subroutines
#
# ----------------------------------------------------------------

# print a blank line
my $empty_line = "\tfprintf (out,\"\\n\");\n";


# ----------------------------------------------------------------
# Make routines for PDF code generator (C version)

sub make_PDFgen_C
{
    my $DISTR = $_[0];   # data for distributions

    my $PDFgen;
    
    # C version
    $PDFgen .= make_bar("C version");

    # Main (C)
    $PDFgen .= make_PDF_main_C( $DISTR );

    # Continuous distributions
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	# Function prototype
	$PDFgen_prototypes .= "static ".get_PDFgen_funct_C($d).";\n";

	# Function for distribution
	$PDFgen .= make_PDF_distr_C($DISTR,$d);
    }

    # End
    return $PDFgen;

} # end of make_PDFgen_C()


# ----------------------------------------------------------------
# Make routines for PDF code generator (FORTRAN version)

sub make_PDFgen_FORTRAN
{
    my $DISTR = $_[0];   # data for distributions

    my $PDFgen;
    
    # C version
    $PDFgen .= make_bar("FORTRAN version");

    # Main (C)
    $PDFgen .= make_PDF_main_FORTRAN( $DISTR );

    # Continuous distributions
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	# Function prototype
	$PDFgen_prototypes .= "static ".get_PDFgen_funct_FORTRAN($d).";\n";

	# Function for distribution
	$PDFgen .= make_PDF_distr_FORTRAN($DISTR,$d);
    }

    # End
    return $PDFgen;

} # end of make_PDFgen_FORTRAN()


# ----------------------------------------------------------------
# Make routines for PDF code generator (JAVA version)

sub make_PDFgen_JAVA
{
    my $DISTR = $_[0];   # data for distributions

    my $PDFgen;
    
    # JAVA version
    $PDFgen .= make_bar("JAVA version");

    # Main (JAVA)
    $PDFgen .= make_PDF_main_JAVA( $DISTR );

    # Continuous distributions
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	# Function prototype
	$PDFgen_prototypes .= get_PDFgen_funct_JAVA($d).";\n";

	# Function for distribution
	$PDFgen .= make_PDF_distr_JAVA($DISTR,$d);
    }

    # End
    return $PDFgen;

} # end of make_PDFgen_JAVA()



# ----------------------------------------------------------------
# Main (C version)

sub make_PDF_main_C
{
    my $DISTR = $_[0];   # data for distributions
    
    my $gencode;         # code for creating PDFs

    # Mark begin of Main
    $gencode .= make_bar("PDF main");

    # Function header
    $gencode .= "int _unur_acg_C_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf)\n{\n";

    # Check for invalid NULL pointer
    $gencode .= "\t_unur_check_NULL(\"ACG\", distr, 0 );\n\n";

    # Switch (distribution)
    $gencode .= "\tswitch (distr->id) {\n";
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	$gencode .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";
	$gencode .= "\t\treturn _unur_acg_C_PDF_$d (distr,out,pdf);\n";
    }
    $gencode .= "\tdefault:\n";
    $gencode .= "\t\t_unur_error(distr->name,UNUR_ERR_GEN_DATA,\"Cannot make PDF\");\n";
    $gencode .= "\t\treturn 0;\n";
    $gencode .= "\t}\n";
    
    # End of function
    $gencode .= "}\n";

    return $gencode;

} # end of make_PDF_main_C()


# ----------------------------------------------------------------
# Main (FORTRAN version)

sub make_PDF_main_FORTRAN
{
    my $DISTR = $_[0];   # data for distributions
    
    my $gencode;         # code for creating PDFs

    # Mark begin of Main
    $gencode .= make_bar("PDF main");

    # Function header
    $gencode .= "int _unur_acg_FORTRAN_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf)\n{\n";

    # Check for invalid NULL pointer
    $gencode .= "\t_unur_check_NULL(\"ACG\", distr, 0 );\n\n";

    # Switch (distribution)
    $gencode .= "\tswitch (distr->id) {\n";
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	$gencode .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";
	$gencode .= "\t\treturn _unur_acg_FORTRAN_PDF_$d (distr,out,pdf);\n";
    }
    $gencode .= "\tdefault:\n";
    $gencode .= "\t\t_unur_error(distr->name,UNUR_ERR_GEN_DATA,\"Cannot make PDF\");\n";
    $gencode .= "\t\treturn 0;\n";
    $gencode .= "\t}\n";
    
    # End of function
    $gencode .= "}\n";

    return $gencode;

} # end of make_PDF_main_FORTRAN()


# ----------------------------------------------------------------
# Main (JAVA version)

sub make_PDF_main_JAVA
{
    my $DISTR = $_[0];   # data for distributions
    
    my $gencode;         # code for creating PDFs

    # Mark begin of Main
    $gencode .= make_bar("PDF main");

    # Function header
    $gencode .= "int _unur_acg_JAVA_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf)\n{\n";

    # Check for invalid NULL pointer
    $gencode .= "\t_unur_check_NULL(\"ACG\", distr, 0 );\n\n";

    # Switch (distribution)
    $gencode .= "\tswitch (distr->id) {\n";
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	$gencode .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";
	$gencode .= "\t\treturn _unur_acg_JAVA_PDF_$d (distr,out,pdf);\n";
    }
    $gencode .= "\tdefault:\n";
    $gencode .= "\t\t_unur_error(distr->name,UNUR_ERR_GEN_DATA,\"Cannot make PDF\");\n";
    $gencode .= "\t\treturn 0;\n";
    $gencode .= "\t}\n";
    
    # End of function
    $gencode .= "}\n";

    return $gencode;

} # end of make_PDF_main_JAVA()


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
    $gencode .= get_PDFgen_funct_C($d)."\n{\n";

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
    $gencode .= $empty_line;

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
    $gencode .= get_PDFgen_funct_FORTRAN($d)."\n{\n";

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
    $gencode .= get_PDFgen_funct_JAVA($d)."\n{\n";

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
    $gencode .= $empty_line;

    $gencode .= "\n\treturn 1;\n";
    $gencode .= "}\n";

    return $gencode;
    
} # make_PDF_distr_JAVA()


# ----------------------------------------------------------------
# Get name of generation routine (C version)

sub get_PDFgen_funct_C
{
    my $d = $_[0];       # distribution

    return "int _unur_acg_C_PDF_$d (UNUR_DISTR *distr, FILE *out, const char *pdf)";

} # end of get_PDFgen_funct_C()


# ----------------------------------------------------------------
# Get name of generation routine (FORTRAN version)

sub get_PDFgen_funct_FORTRAN
{
    my $d = $_[0];       # distribution

    return "int _unur_acg_FORTRAN_PDF_$d (UNUR_DISTR *distr, FILE *out, const char *pdf)";

} # end of get_PDFgen_funct_FORTRAN()


# ----------------------------------------------------------------
# Get name of generation routine (JAVA version)

sub get_PDFgen_funct_JAVA
{
    my $d = $_[0];       # distribution

    return "int _unur_acg_JAVA_PDF_$d (UNUR_DISTR *distr, FILE *out, const char *pdf)";

} # end of get_PDFgen_funct_JAVA()



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













