#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

use strict;
my $DEBUG = 0;

# ----------------------------------------------------------------
# Directory with sources

my $top_srcdir = $ENV{'srcdir'} ? $ENV{'srcdir'} : '.';
$top_srcdir .= "/..";

# ----------------------------------------------------------------

require "$top_srcdir/scripts/read_PDF.pl";

# ----------------------------------------------------------------

# C file for code generator  
#my $headergen_file = "headergen.c";

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

    unless (/^s*=INPUT\s*(\w*)/){
	print $_;
    }

    else {
	my $type = $1;    # INPUT type 

	if ( $type eq "list_of_distributions" ){
	    print make_list_of_distributions($DISTR);
	}
	else{
	    die "Error: unknown qualifier after =INPUT: $type\n";
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
# Make list of distributions 

sub make_list_of_distributions
{
    my $DISTR = $_[0];   # data for distributions

    my $gencode;         # code for creating ACG header

    # Switch list for distributions
    foreach my $d (sort keys %{$DISTR}) {
	# at the moment we include continuous distributions only
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	# add distribution
	$gencode .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";

	# PDF
	my $PDF = scan_PDF($DISTR,$d);
	$gencode .= "\t\tfprintf (out,sformat,\"  PDF(x) = $PDF\");\n";

	# List of parameters for PDF
	$gencode .= make_header_PDFparams($DISTR,$d);

	$gencode .= "\t\tbreak;\n";

    }

    return $gencode;

} # end of make_list_of_distributions()


# ----------------------------------------------------------------
# Process parameter list

sub make_header_PDFparams
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
	    "\t\t_unur_acg_header_param (out, \"$param\", distr->data.cont.params[$i], hrule, sformat);\n";
    }

    # end
    return $gencode;

} # end of make_header_PDFparams()


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
    $PDF =~ s/\\hfil+\\break/;/g;

    $PDF =~ s/\\frac\{([^\}]+[\s\+\-]+[^\}]+)\}\{([^\}]+[\s\+\-]+[^\}]+)\}/\($1\)\/\($2\)/g;
    $PDF =~ s/\\frac\{([^\}]+)\}\{([^\}]+[\s\+\-]+[^\}]+)\}/$1\/\($2\)/g;
    $PDF =~ s/\\frac\{([^\}]+[\s\+\-]+[^\}]+)\}\{([^\}]+)\}/\($1\)\/$2/g;
    $PDF =~ s/\\frac\{([^\}]+)\}\{([^\}]+)\}/$1\/$2/g;

    $PDF =~ s/\{/\(/g;
    $PDF =~ s/\}/\)/g;

    # replace tabs and multiple blanks by single blanks
    $PDF =~ s/\t+/ /g;
    $PDF =~ s/ +/ /g;

    # return result
    return $PDF;

} # end of scan_PDF()

# ----------------------------------------------------------------

