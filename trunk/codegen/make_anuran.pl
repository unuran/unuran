#!/usr/bin/perl
# ----------------------------------------------------------------
# make CGI script for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

use strict;

# ----------------------------------------------------------------
# Directory with sources

my $top_srcdir = $ENV{'srcdir'} ? $ENV{'srcdir'} : '.';
$top_srcdir .= "/..";

# ----------------------------------------------------------------
# List of standard distributions

require "$top_srcdir/scripts/read_PDF.pl";

my $DISTR = read_PDFdata( $top_srcdir );

# For description of data fields in this list see file `readPDF.pl'.

# we do not need the uniform distribution
undef $DISTR->{'uniform'};

# ----------------------------------------------------------------
# Read template C file from STDIN and insert C code for string interpreter 

while ( <STDIN> ){

    unless (/^s*=INPUT\s*(\w+)\s+(\w*)/){
	print $_;
    }

    else {
	my $key   = $1;    # INPUT key value
	my $value = $2;

	if ( $key eq "labels_menue_cstd_distributions" ){
	    # Popup menue for list of distributions
	    print make_label_menue_cstd_distribution($DISTR);
	}
	elsif ( $key eq "data_distributions" ) {
	    # Data for distributions
	    print make_data_distributions($DISTR,$value);
	}
	else {
	    die "Error: unknown qualifier after =INPUT: $key\n";
	}
    }
}

# ----------------------------------------------------------------
# End
# ----------------------------------------------------------------
exit 0;
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# make entries for popup menue for list of distributions
# ----------------------------------------------------------------

sub make_label_menue_cstd_distribution
{
    my $DISTR = $_[0];

    # associate array for CGI script
    my $label;

    # list of distributions
    foreach my $d (sort keys %{$DISTR}) { 
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$label .= "\t'$d' => \"".$DISTR->{$d}->{"=NAME"}."\",\n\t";
    }

    # end
    return $label;

} # end of make_label_menue_distribution()

# ----------------------------------------------------------------
# make table with data of distributions
# ----------------------------------------------------------------

sub make_data_distributions
{
    my $DISTR = $_[0];

    # associate array for CGI script
    my $data;

    # generic distributions
    $data .= "\t'cont' => {\n\t\t";
    $data .= "'=NAME' => 'my distribution',\n\t\t";
    $data .= "'=N_REQ' => '0',\n\t\t";
    $data .= "'=N_TOT' => '0',\n\t\t";
    $data .= "'=LEFT' => '-infinity',\n\t\t";
    $data .= "'=RIGHT' => 'infinity',\n\t\t";
    $data .= "},\n";

    # list of standard distributions
    foreach my $d (sort keys %{$DISTR}) { 
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$data .= "\t'$d' => {\n\t\t";

	# name
	$data .= "'=NAME' => \"".$DISTR->{$d}->{"=NAME"}."\",\n\t\t";

	# PDF
	$data .= "'=PDF' => '".format_PDF($DISTR->{$d}->{"=DOC"}->{"=PDF"})."',\n\t\t";

	# parameters for distribution
	$data .= scan_FPARAM($DISTR->{$d}->{"=DOC"}->{"=FPARAM"});

	# domain
	my $domain = $DISTR->{$d}->{"=DOC"}->{"=DOMAIN"};
	$domain =~ s/\s//g;
	(my $left, my $right) = split /[<>=]+x[<>=]+/, $domain, 2;
	$data .= "'=LEFT' => '$left',\n\t\t";
	$data .= "'=RIGHT' => '$right',\n\t\t";

	$data .= "},\n";
    }

    # end
    return $data;

} # end of make_data_distributions()

# ----------------------------------------------------------------
# Format PDF for distribution
# ----------------------------------------------------------------

sub format_PDF 
{
    my $PDF = $_[0];

    # empty ? 
    return '' unless $PDF;

    # trim heading blanks
    $PDF =~ s/^\s*//;

    # chop off trailing blanks
    $PDF =~ s/\s+$//;

    # remove newlines
    $PDF =~ s/\n+/ /g;

    # format output
    $PDF =~ s/\\over\s+/\//g;

    $PDF =~ s/\\hbox{\s*(\w+)\s*}/ $1 /g;
    $PDF =~ s/\\hfil+\\break/<BR>/g;

    $PDF =~ s/\\frac\{([^\}]+[\s\+\-]+[^\}]+)\}\{([^\}]+[\s\+\-]+[^\}]+)\}/\($1\)\/\($2\)/g;
    $PDF =~ s/\\frac\{([^\}]+)\}\{([^\}]+[\s\+\-]+[^\}]+)\}/$1\/\($2\)/g;
    $PDF =~ s/\\frac\{([^\}]+[\s\+\-]+[^\}]+)\}\{([^\}]+)\}/\($1\)\/$2/g;
    $PDF =~ s/\\frac\{([^\}]+)\}\{([^\}]+)\}/$1\/$2/g;

    $PDF =~ s/\^\s*\{([^\}]+)\}/<SUP>($1)<\/SUP>/g;
    $PDF =~ s/\^\s*(\w+)/<SUP>$1<\/SUP>/g;

    $PDF =~ s/\{/\(/g;
    $PDF =~ s/\}/\)/g;

    # return result
    return $PDF;

} # end of format_PDF()

# ----------------------------------------------------------------
# scan list of parameters for distribution
# ----------------------------------------------------------------

sub scan_FPARAM {
    my $DOC_params = $_[0];

    # empty ? 
    return "" unless $DOC_params;

    # trim heading blanks
    $DOC_params =~ s/^\s*//;

    # chop off trailing blanks
    $DOC_params =~ s/\s+$//;

    # remove blanks around `:'
    $DOC_params =~ s/[ \t]*\:[ \t]*/\:/g;

    # split into lines
    my @lines = split /\n+/, $DOC_params;

    # process lines
    my $params = "[\n\t\t\t";
    my $n_total = 0;
    my $n_optional = 0;

    foreach my $l (@lines) {
	# each line must start with a [\d+]
	next unless $l =~ /\[*\d+/;

	# split into columns
	my @cols = split /\:/, $l;
	die "\nwrong number of columns for =FPARAM: $#cols" if $#cols != 4;

	# get entries
	my $number  = $cols[0];
	$number     =~ s/(.*)(\d+).*/$2/;
	my $name    = $cols[1];
	my $cond    = $cols[2];
	my $default = $cols[3];
	my $type    = $cols[4];

	# append to list of parameters
	$params .= "{";
	$params .= "'=NAME' => '$name', ";
	$params .= "'=DEF'  => '".((length $default) ? $default : "*")."', ";
	$params .= "'=COND' => '$cond', ";

	# upper and lower bounds for parameter
	$cond =~ s/\s//g;
	my $upper = '';
	my $lower = '';
      SWITCH: {
	  if ($cond =~ /^>([^<>=]+)$/) {
	      $lower = $1;
	      last SWITCH;
	  }
	  if ($cond =~ /^([^<>=]+)[<=]+([^<>=]+)[<=]+([^<>=]+)$/) {
	      die "wrong variable in condition" if $2 ne $name; 
	      $lower = $1;
	      $upper = $3;
	      last SWITCH;
	  }
      }
	
	if (length $upper > 0) {
	    $params .= "'=UPPER' => '$upper', ";
	}
	if (length $lower > 0) {
	    $params .= "'=LOWER' => '$lower', ";
	}

	# close
	$params .= "},\n\t\t\t";


	# increment counter for parameters
	++$n_total;
	++$n_optional if length $default;
    }

    # close array of parameters
    $params .= "],\n\t\t";

    # number of parameters for standard form
    my $n_required = $n_total - $n_optional;

    # Return result
    return 
	"'=N_REQ' => '$n_required',\n\t\t".
        "'=N_TOT' => '$n_total',\n\t\t".
	"'=FPARAMS' => $params,\n\t\t";
	    

} # end of scan_FPARAM()

# ----------------------------------------------------------------
