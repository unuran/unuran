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

# Make FORTRAN version of code generator
## $headergen .= make_headergen_FORTRAN($DISTR);

# Make JAVA version of code generator
## $headergen .= make_headergen_JAVA($DISTR);

# ................................................................

# Print code into file
open FILE, ">$headergen_file" or die "cannot open file $headergen_file\n";

print FILE make_bar("include header file");
print FILE "\#include <codegen_source.h>\n";

print FILE make_bar("local prototypes");
print FILE $headergen_prototypes;

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

# print a blank line
my $empty_line = "\tfprintf (out,\"\\n\");\n";


# ----------------------------------------------------------------
# Make routines for ACG header code generator (C version)

sub make_headergen_C
{
    my $DISTR = $_[0];   # data for distributions

    my $headergen;
    
    # C version
    $headergen .= make_bar("C version");

    # Constants 
    $headergen .= make_header_constants();

    # Main (C)
    $headergen .= make_header_main_C( $DISTR );

    # Print date when generated
    $headergen .= make_header_date_C();

    # Copyright notice
    $headergen .= make_header_copyright_C();

    # Function for printing parameter values
    $headergen .= print_param();

    # Function for printing of PDF
    $headergen .= make_header_PDFdomain_C();

    # Continuous distributions
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	# Function prototype
	$headergen_prototypes .= "static ".get_headergen_funct_C($d).";\n";

	# Function for distribution
	$headergen .= make_header_distr_C($DISTR,$d);
    }

    # End
    return $headergen;

} # end of make_headergen_C()


# ----------------------------------------------------------------
# Main (C version)

sub make_header_main_C
{
    my $DISTR = $_[0];   # data for distributions
    
    my $gencode;         # code for creating ACG header

    # Mark begin of Main
    $gencode .= make_bar("ACG header main");

    # Function header
    $gencode .= "int _unur_acg_C_header (UNUR_DISTR *distr, FILE *out, const char *rand)\n{\n";

    # Check for invalid NULL pointer
    $gencode .= "\t_unur_check_NULL(\"ACG\", distr, 0 );\n\n";

    # Switch (distribution)
    $gencode .= "\tswitch (distr->id) {\n";
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	$gencode .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";
	$gencode .= "\t\treturn _unur_acg_C_header_$d (distr,out,rand);\n";
    }
    $gencode .= "\tdefault:\n";
    $gencode .= "\t\t_unur_error(distr->name,UNUR_ERR_GEN_DATA,\"Cannot make ACG header\");\n";
    $gencode .= "\t\treturn 0;\n";
    $gencode .= "\t}\n";
    
    # End of function
    $gencode .= "}\n";

    return $gencode;

} # end of make_header_main_C()


# ----------------------------------------------------------------
# Print header code for distribution (C version)

sub make_header_distr_C
{
    my $DISTR = $_[0];   # data for distributions
    my $d = $_[1];       # distribution

    my $gencode;         # code for creating ACG header

    # Mark begin of distribution
    $gencode .= make_bar("PDF $d: ".$DISTR->{$d}->{"=NAME"});

    # Function header
    $gencode .= get_headergen_funct_C($d)."\n{\n";

    # Name of distribution
    $gencode .= 
	 "\tfprintf (out,hrule);\n"
	."\tif (distr->set \& UNUR_DISTR_SET_STDDOMAIN) {\n"
	."\t\tfprintf (out,sformat,\"Generator for $d distribution.\");\n"
	."\t} else {\n"
	."\t\tfprintf (out,sformat,\"Generator for truncated $d distribution.\");\n"
	."\t}\n";

    # PDF
    my $PDF = scan_PDF($DISTR,$d);
    $gencode .= 
	 "\tfprintf (out,hrule);\n"
	."\tfprintf (out,sformat,\"  PDF(x) = $PDF\");\n";

    # List of parameters for PDF
    $gencode .= 
	make_header_PDFparams_C($DISTR,$d);

    # Domain of PDF
    $gencode .= 
	"\t_unur_acg_C_header_domain (out, distr);\n";

    # Synopsis
    $gencode .= 
	 "\tfprintf (out,hrule);\n"
	."\tfprintf (out,sformat,\"Synopsis:\");\n"
	."\tfprintf (out,\"\\ndouble %s (void);\\n\\n\",rand);\n";

    # Copyright notice
    $gencode .= 
	"\t_unur_acg_C_header_copyright (out);\n";

    $gencode .= $empty_line;

    $gencode .= "\n\treturn 1;\n";
    $gencode .= "}\n";

    return $gencode;
    
} # make_header_distr_C()


# ----------------------------------------------------------------
# Get name of generation routine (C version)

sub get_headergen_funct_C
{
    my $d = $_[0];       # distribution

    return "int _unur_acg_C_header_$d (UNUR_DISTR *distr, FILE *out, const char *rand)";

} # end of get_headergen_funct_C()


# ----------------------------------------------------------------
# Process parameter list

sub make_header_PDFparams_C
{
    my $DISTR = $_[0];    # data for distributions
    my $d = $_[1];        # name of distribution

    my $gencode;         # code for creating ACG header

    # parameters for PDF
    my $n_in_params = $DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"};
    my $in_params = $DISTR->{$d}->{"=PDF"}->{"=PARAMS"};

    foreach my $i (0 .. $n_in_params - 1) {
	my $param = $in_params->[$i];
	$gencode .= 
	    "\t_unur_acg_C_header_param (out,\"$param\",distr->data.cont.params[$i]);\n";
    }

    # end
    return $gencode;

} # end of make_header_PDFparams_C()


# ----------------------------------------------------------------
# Print domain

sub make_header_PDFdomain_C
{
    my $gencode;         # code for creating ACG header

    # Mark begin of Copyright notice
    $gencode .= make_bar("Print domain for PDF");

    # Function prototype
    $headergen_prototypes .= "static int _unur_acg_C_header_domain (FILE *out, UNUR_DISTR *distr);\n";

    # Print copyright notice
    $gencode .= <<EOX;
int _unur_acg_C_header_domain (FILE *out, UNUR_DISTR *distr)
{
\tchar buf[1024];
\tsprintf (buf,"      domain = [ %g, %g ]", distr->data.cont.domain[0], distr->data.cont.domain[1]);
\tfprintf (out,sformat,buf);
\treturn 1;
}
EOX

    # end
    return $gencode;

} # end of make_header_PDFdomain_C()


# ----------------------------------------------------------------
# scan PDF for distribution
#

sub scan_PDF {
    my $DISTR = $_[0];    # data for distributions
    my $d = $_[1];        # name of distribution

    # content of PDF description
    my $PDF = $DISTR->{$d}->{"=DOC"}->{"=PDF"};

    # empty ? 
    return unless $PDF;

    # trim heading blanks
    $PDF =~ s/^\s*//;

    # chop off trailing blanks
    $PDF =~ s/\s+$//;

    # remove newlines
    $PDF =~ s/\n+/ /g;

    # format other output
    $PDF =~ s/\\over\s+/\//g;

    $PDF =~ s/\\hbox{\s*(\w+)\s*}/ $1 /g;
    $PDF =~ s/\\hfil+\\break/\\n/g;

    $PDF =~ s/\\frac\{([^\}]+[\s\+\-]+[^\}]+)\}\{([^\}]+[\s\+\-]+[^\}]+)\}/\($1\)\/\($2\)/g;
    $PDF =~ s/\\frac\{([^\}]+)\}\{([^\}]+[\s\+\-]+[^\}]+)\}/$1\/\($2\)/g;
    $PDF =~ s/\\frac\{([^\}]+[\s\+\-]+[^\}]+)\}\{([^\}]+)\}/\($1\)\/$2/g;
    $PDF =~ s/\\frac\{([^\}]+)\}\{([^\}]+)\}/$1\/$2/g;

    $PDF =~ s/\{/\(/g;
    $PDF =~ s/\}/\)/g;

    # return result
    return $PDF;

} # end of scan_PDF()


# ----------------------------------------------------------------
# Constants for headergen file

sub make_header_constants
{
    my $gencode;         # code for creating ACG header

    # Mark begin of Copyright notice
    $gencode .= make_bar("constants");

    # Print copyright notice
    $gencode .= <<EOX;
/* C version */
const static char hrule_C[] = "/* ---------------------------------------------------------------- */\\n";
const static char sformat_C[] = "/* %-64.64s */\\n";

/* FORTRAN version */
const static char hrule_FORTRAN[] = "* ------------------------------------------------------------------ *\\n";
const static char sformat_FORTRAN[] = "* %-66.66s *\\n";

/* JAVA version */
const static char *hrule_JAVA = hrule_C;
const static char *sformat_JAVA = sformat_C;

/* switch */
const static char *hrule = hrule_C;
const static char *sformat = sformat_C;
EOX

    return $gencode;
} # end of make_header_constants_C()

# ----------------------------------------------------------------
# Print copyright notice (C version)

sub make_header_copyright_C
{
    my $gencode;         # code for creating ACG header

    # Mark begin of Copyright notice
    $gencode .= make_bar("ACG copyright notice");

    # Function prototype
    $headergen_prototypes .= "static int _unur_acg_C_header_copyright (FILE *out);\n";

    # Print copyright notice
    $gencode .= <<EOX;
int _unur_acg_C_header_copyright (FILE *out)
{
\tfprintf (out,hrule);
\tfprintf (out,sformat,"See NOTICE for uniform random number generator below!");
\tfprintf (out,hrule);
\t_unur_acg_C_header_date (out);
\tfprintf (out,sformat,"");
\tfprintf (out,sformat,"Copyright (c) 2001   Wolfgang Hoermann and Josef Leydold");
\tfprintf (out,sformat,"Dept. for Statistics, University of Economics, Vienna, Austria");
\tfprintf (out,sformat,"All rights reserved.");
\tfprintf (out,sformat,"");
\tfprintf (out,sformat,"Redistribution and use in source and binary forms, with or");
\tfprintf (out,sformat,"without modification, are permitted without a fee.");
\tfprintf (out,sformat,"Redistributions of source code must retain this copyright notice");
\tfprintf (out,sformat,"and the following disclaimer.");
\tfprintf (out,sformat,"");
\tfprintf (out,sformat,"THIS PROGRAM IS PROVIDED \\\"AS IS\\\" AND WITHOUT ANY EXPRESS OR");
\tfprintf (out,sformat,"IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED");
\tfprintf (out,sformat,"WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR");
\tfprintf (out,sformat,"PURPOSE.");
\tfprintf (out,sformat,"");
\tfprintf (out,sformat,"IN NO EVENT SHALL JOSEF LEYDOLD OR WOLFGANG HOERMANN BE LIABLE");
\tfprintf (out,sformat,"FOR ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES");
\tfprintf (out,sformat,"OF ANY KIND, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF");
\tfprintf (out,sformat,"USE, DATA OR PROFITS, WHETHER OR NOT ADVISED OF THE POSSIBILITY");
\tfprintf (out,sformat,"OF DAMAGE, AND ON ANY THEORY OF LIABILITY, ARISING OUT OF OR IN");
\tfprintf (out,sformat,"CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.");
\tfprintf (out,hrule);

\treturn 1;
}
EOX

    return $gencode;

} # end of make_header_copyright_C()

# ----------------------------------------------------------------
# Print date notice (C version)

sub make_header_date_C
{
    my $gencode;         # code for creating ACG header

    # Mark begin of Copyright notice
    $gencode .= make_bar("ACG copyright notice");

    # Function prototype
    $headergen_prototypes .= "static int _unur_acg_C_header_date (FILE *out);\n";

    # Print copyright notice
    $gencode .= <<EOX;
int _unur_acg_C_header_date (FILE *out)
{
\tchar buf[1024];
\tchar date[64];
\ttime_t timer;
\ttimer = time(&timer);
\tstrftime(date,60,"%B %d %Y %H:%M:%S",localtime(&timer));
\tsprintf (buf,"Generated by ANURAN (v.%s).   %s", VERSION, date);
\tfprintf (out,sformat,buf);
\treturn 1;
}
EOX

    return $gencode;

} # end of  make_header_date_C()

# ----------------------------------------------------------------
# Print parameter

sub print_param
{

    my $gencode;         # code for creating ACG header

    # Mark begin of Copyright notice
    $gencode .= make_bar("Print parameter for PDF");

    # Function prototype
    $headergen_prototypes .= "static int _unur_acg_C_header_param (FILE *out, char *par, double val);\n";

    # Print copyright notice
    $gencode .= <<EOX;
int _unur_acg_C_header_param (FILE *out, char *par, double val)
{
\tchar buf[1024];
\tsprintf (buf,"      %s = %g", par, val);
\tfprintf (out,sformat,buf);
\treturn 1;
}
EOX

    return $gencode;

} # end of print_param()

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


