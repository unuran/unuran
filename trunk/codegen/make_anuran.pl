#!/usr/bin/perl
# ----------------------------------------------------------------
# make CGI script for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

# ----------------------------------------------------------------

require "read_PDF.pl";

# ----------------------------------------------------------------
# Gray Color for disabled text

my $Gray = "\#C0C0C0";

# ----------------------------------------------------------------
# List of distributions
my $DISTR = read_PDFdata();

# For description of data fields in this list see file `readPDF.pl'.

# ----------------------------------------------------------------
# Popup menue for list of distributions

my $label_menue_distributions = make_label_menue_distribution($DISTR);

# ----------------------------------------------------------------
# Data for distributions

my $data_distr = make_data_distributions($DISTR);

# ----------------------------------------------------------------
# Print main part of CGI script
# ----------------------------------------------------------------

print <<EOX;
#!/usr/bin/perl
# ----------------------------------------------------------------
# CGI script for code generator
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Load packages

use CGI qw/:standard *blockquote *font/;
use CGI::Pretty;

# ----------------------------------------------------------------
# Create CGI object

my \$q = new CGI;

# ----------------------------------------------------------------
# Data for distributions

\%labels_menue_distributions = $label_menue_distributions;

\@menue_distributions = sort keys \%labels_menue_distributions;

\%data_distr = $data_distr;

# ----------------------------------------------------------------
# Get some data from query

\$step = \$q->param('step');
Delete('step');

\$distr = \$q->param('distribution');
\$Stdform = \$q->param('Stdform');

\$n_req = \$data_distr{\$distr}{'=N_REQ'};
\$n_tot = \$data_distr{\$distr}{'=N_TOT'};

# ----------------------------------------------------------------

anuran_start();
anuran_select_distribution();
anuran_params_distribution();
anuran_domain_distribution();
anuran_properties_generator();
anuran_language();
anuran_code();
anuran_end();

# ----------------------------------------------------------------
exit 0;
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# Start of HTML page
# ----------------------------------------------------------------

sub anuran_start
{
    print
	
	# HTTP header
	header(),

	# HTML header and title
	start_html(-title=>'ANURAN Test Page',
		   -author=>'unuran\@statistik.wu-wien.ac.at',
		   -meta=>{'keywords'=>'nonuniform random number generator',
			   'copyright'=>'copyright 2001 Institut fuer Statistik, WU Wien'},
		   -BGCOLOR=>'#FFFFF5'),

	# Heading
	h1('Automatic Code Generator'),
	b('Experimental version'),

	# ruler
	p().hr().p();

} # end of anuran_start()

# ----------------------------------------------------------------
# Step 1: Select a distribution
# ----------------------------------------------------------------

sub anuran_select_distribution
{
    my \$warning;

    # distribution properly selected ?
    if (\$step > 1 && \$distr eq '--none--') {
	# no distribution selected at step 1
	\$step = 1;
	\$warning = p().b('You must select a distribution first!').'\n';
    }

# ................................................................
# CASE: do first step
# ................................................................

    if (\$step <=1) { 
	print
	    
	    # heading
	    'Step 1: '.b('Distribution').br(),
	    
	    # form
	    blockquote( \$warning,
			startform(-method=>'GET',
				  -action=>\$q->url()),
			popup_menu(-name=>'distribution',
				   -values=>\\\@menue_distributions,
				   -labels=>\\\%labels_menue_distributions),
			'&nbsp;&nbsp;&nbsp;\n',
			hidden('step','2'),'\n',
			submit('Continue'),'\n',
			endform() ),
	    
	    # ruler
	    p().hr().p();
    }

# ................................................................
# CASE: first step done
# ................................................................

    if (\$step > 1) {
	print

	    # selected distribution
	    'Step 1: '.b(\$data_distr{\$distr}{'=NAME'}).br(),

	    # PDF of distribution
	    blockquote( 'PDF(x) = [const] * ',
			\$data_distr{\$distr}{'=PDF'} ),

	    # ruler
	    p().hr().p();
    }

} # end of anuran_select_distribution()

# ----------------------------------------------------------------
# Step 2: Parameters of distribution
# ----------------------------------------------------------------

sub anuran_params_distribution
{

# ................................................................
# CASE: no distribution selected till now    
# ................................................................

    if (\$step < 2) {
	print

	    # heading
	    font({-color=>"$Gray"},
		 'Step 2: '.b('Parameters for distribution')),
	    br(),

	    # ruler
	    p().hr().p();
    }

# ................................................................
# Distribution selected
# ................................................................
    
    if (\$step >= 2) {
	print 
	    'Step 2: ',
	    b('Parameters for '.\$data_distr{\$distr}{'=NAME'}).br();
    }

# ................................................................
# CASE: Distribution selected and Parameters given
# ................................................................

    if (\$step > 2) {

	# Read parameters for distribution
	if (\$Stdform eq 'yes') {
	    # Standardform

	    \$n_param = \$n_req;
	    for (my \$i=0; \$i < \$n_tot; \$i++) {
		\$param[\$i] = \$q->param("Std_param_\$i");
	    }
	}
	else {
	    # Non-Standardform

	    \$n_param = \$n_tot;
	    for (my \$i=0; \$i < \$n_tot; \$i++) {
		\$param[\$i] = \$q->param("param_\$i");
	    }
	}

	# String to be printed on screen
	for (my \$i=0; \$i < \$n_tot; \$i++) {
	    \$print_param[\$i] = \$param[\$i];
	}
	
	# Check parameters
	for (my \$i=0; \$i < \$n_tot; \$i++) {
	    unless (\$param[\$i] =~ /\\d+/) {
		\$step = 2;
		\$print_param[\$i] = '<FONT COLOR="red"><BLINK>missing</BLINK></FONT>';
	    }
	}

	# Print parameters into web page
	print start_blockquote();

	for (my \$i=0; \$i < \$n_param; \$i++) {
	    print 
		\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'},
		' = ',
		\$print_param[\$i], '&nbsp;&nbsp;&nbsp;';
	}

	if (\$n_tot > \$n_param) {
	    print 
		start_font({-color=>'$Gray'}),
		'(&nbsp;&nbsp;&nbsp;';
	    for (my \$i=\$n_param; \$i < \$n_tot; \$i++) {
		print 
		    \$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'},
		    ' = ',
		    \$print_param[\$i], '&nbsp;&nbsp;&nbsp;';
	    }
	    print 
		')',
		end_font();
	}

	print end_blockquote();

	# ruler
	if (\$step > 2) {
	    print p().hr().p();
	}
    }

# ................................................................
# CASE: Distribution selected, input Parameters
# ................................................................

    if (\$step == 2) {

	# start form
	print 
	    start_blockquote(),
	    startform(-method=>'GET',
		      -action=>\$q->url());

	if (\$n_req < \$n_tot) {
	    # Standardform exists

	    if (\$Stdform ne 'no') {
		print '\t<INPUT TYPE="radio" NAME="Stdform" VALUE="yes" CHECKED> Standardform\n';
		print start_blockquote();
		for (my \$i = 0; \$i < \$n_req; \$i++) {
		    print 
			\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'},
			' = ',
			textfield(-name=>"Std_param_\$i",
				  -size=>10,
				  -maxlength=>10),
			( (\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=COND'})
			  ? "    &nbsp;&nbsp;&nbsp;(\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=COND'})"
			  : ''),
			br();
		}
		print 
		    p(), start_font({-color=>'$Gray'}),
		    '(&nbsp;&nbsp;&nbsp;';
		for (my \$i = \$n_req; \$i < \$n_tot; \$i++) {
		    print 
			hidden(-name=>"Std_param_\$i",
			       -value=>\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=DEF'}),
			\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'},
			' = ',
			\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=DEF'},
			'&nbsp;&nbsp;&nbsp;';
		}
		print ')', end_font();
		print end_blockquote();
	    }

	    if (\$Stdform ne 'yes') {
		my \$checked = (\$Stdform eq 'no') ? 'CHECKED' : ''; 
		print '\t<INPUT TYPE="radio" NAME="Stdform" VALUE="no" '.\$checked.'> Non-Standardform\n';
		print start_blockquote();
		for (my \$i = 0; \$i < \$n_tot; \$i++) {
		    print 
			\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'},
			' = ',
			textfield(-name=>"param_\$i",
				  -default=>((\$i<\$n_req)
					     ? ''
					     : \$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=DEF'}),
				  -size=>10,
				  -maxlength=>10),
			( (\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=COND'})
			  ? "    &nbsp;&nbsp;&nbsp;(\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=COND'})"
			  : ''),
			br();
		}
		print end_blockquote();
	    }
	}
	
	else {
	    # No special standard form
	    print hidden('Stdform','no'),'\n';
	    for (my \$i = 0; \$i < \$n_tot; \$i++) {
		print 
		    \$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'},
		    ' = ',
		    textfield(-name=>"Std_param_\$i",
			      -size=>10,
			      -maxlength=>10),
		    ( (\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=COND'})
		      ? "    &nbsp;&nbsp;&nbsp;(\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=COND'})"
		      : ''),
		    br();
	    }
	}

	# make a check box for truncated distribution
	print 
	    p(),
	    checkbox(-name=>'truncated',
		     -value=>'yes',
		     -label=>' Truncated domain');

	# Continue
	print
	    p(),
	    hidden('step','3'),'\n',
	    hidden('distribution'),
	    submit('Continue'),'\n';

	# end of form
	print 
	    end_form(),
	    end_blockquote();

	# ruler
	print p().hr().p();
    }


} # end of anuran_params_distribution()

# ----------------------------------------------------------------
# Step 3: Domain of distribution
# ----------------------------------------------------------------

sub anuran_domain_distribution
{

# ................................................................
# CASE: no distribution or no parameters selected till now    
# ................................................................

    if (\$step < 3) {
	print
	    
	    # heading
	    font({-color=>"$Gray"},
		 'Step 3: '.b('Domain for distribution')),
	    br(),
	    
	    # ruler
	    p().hr().p();

	return;
    }

# ................................................................
# CASE: Distribution selected and Parameters given
# ................................................................

    print 
	'Step 3: ',
	b('Domain for '.\$data_distr{\$distr}{'=NAME'}).br();


    # compute standard domain
    \$std_left = '-infinity';
    \$std_right = 'infinity';

    for (my \$i=0; \$i < \$n_tot; \$i++) {
	if (\$data_distr{\$distr}{'=LEFT'} eq \$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'}) {
	    \$std_left = \$param[\$i];
	}
	if (\$data_distr{\$distr}{'=RIGHT'} eq \$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'}) {
	    \$std_right = \$param[\$i];
	}
    }

    # Get inserted domain
    \$left = \$std_left;
    if (defined \$q->param('left') and \$std_left > \$q->param('left')) {
	\$left = \$q->param('left');
    }
    \$right = \$std_right;
    if (defined \$q->param('right') and \$std_right > \$q->param('right')) {
	\$right = \$q->param('right');
    }




    # ruler
    print p().hr().p();


} # end of anuran_domain_distribution()

# ----------------------------------------------------------------
# Step 4: Properties of generator
# ----------------------------------------------------------------

sub anuran_properties_generator
{

    print

	# heading
	font({-color=>"$Gray"},
	     'Step 4: '.b('Properties of generator')),
	br(),

	# ruler
	p().hr().p();

} # end of anuran_properties_generator()

# ----------------------------------------------------------------
# Step 5: Select a programming language
# ----------------------------------------------------------------

sub anuran_language
{

    print

	# heading
	font({-color=>"$Gray"},
	     'Step 5: '.b('Programming language')),
	br(),

	# ruler
	p().hr().p();

} # end of anuran_language()

# ----------------------------------------------------------------
# Step 6: Finally the code
# ----------------------------------------------------------------

sub anuran_code
{

    print

	# heading
	font({-color=>"$Gray"},
	     'Step 6: '.b('Generator code for distribution')),
	br(),

	# ruler
	p().hr().p();

} # end of anuran_code()

# ----------------------------------------------------------------
# End of HTML page
# ----------------------------------------------------------------

sub anuran_end
{
    print 
	
	# HTML page closing
	end_html();

} # end of anuran_end()

# ----------------------------------------------------------------
EOX
# ----------------------------------------------------------------
exit 0;
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# make entries for popup menue for list of distributions
# ----------------------------------------------------------------

sub make_label_menue_distribution
{
    my $DISTR = $_[0];

    # associate array for CGI script
    my $label = 
	"(\n\t".
	"'--none--' => '-- Select a distribution --',\n\t";

    # list of distributions
    foreach my $d (sort keys %{$DISTR}) { 
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$label .= "'$d' => '".$DISTR->{$d}->{"=NAME"}."',\n\t";
    }

    # close array
    $label .= ")";

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
    my $data = "(\n\t";

    # list of distributions
    foreach my $d (sort keys %{$DISTR}) { 
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$data .= "'$d' => {\n\t\t";

	# name 
	$data .= "'=NAME' => '".$DISTR->{$d}->{"=NAME"}."',\n\t\t";

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

	$data .= "},\n\t";
    }

    # close array
    $data .= ")";

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
    my $params;

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

    foreach $l (@lines) {
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
