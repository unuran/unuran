#!/usr/bin/perl
# ----------------------------------------------------------------
# make CGI script for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

use strict;

# ----------------------------------------------------------------
# File names and pathes

# Log file
my $log_dir = "/home/staff/leydold/public_html";
my $anuran_log_file = "$log_dir/anuran.log";

# Automatic code generator
my $ACG = './acg';

# Frogs
my $frog_dir = "http://statistik.wu-wien.ac.at/anuran/frogs";
my $n_frogs = 12;

# ----------------------------------------------------------------
# Directory with sources

my $top_srcdir = $ENV{'srcdir'} ? $ENV{'srcdir'} : '.';
$top_srcdir .= "/..";

# ----------------------------------------------------------------

require "$top_srcdir/scripts/read_PDF.pl";

# ----------------------------------------------------------------
# Gray Color for disabled text

my $Gray = "\#C0C0C0";

# ----------------------------------------------------------------
# C compiler

my $C_inc_dir = '-I/home/staff/leydold/include';
my $C_lib_dir = '-L/home/staff/leydold/lib';
my $C_libs    = '-lunuracg -lunuran -lm -lprng';
my $CC        = 'gcc';

# ----------------------------------------------------------------
# List of distributions
my $DISTR = read_PDFdata( $top_srcdir );

# For description of data fields in this list see file `readPDF.pl'.

# ----------------------------------------------------------------
# Popup menue for list of distributions

my $label_menue_distributions = make_label_menue_distribution($DISTR);

# ----------------------------------------------------------------
# Data for distributions

my $data_distr = make_data_distributions($DISTR);

# ----------------------------------------------------------------
# Bar with links
# ----------------------------------------------------------------

my $logo = <<EOX;
"<TABLE BORDER='0' CELLSPACING='1'>
<TR>
    <TD ALIGN=LEFT><IMG SRC='../arvag/images/logo_arvag_40.png' ALT='' ALIGN=LEFT>
    <TD ALIGN=LEFT>&nbsp;&nbsp;&nbsp;
    <FONT FACE='Arial, Helvetica, sans-serif' SIZE='4' COLOR='#1E426E'>
    ANURAN - Automatic Non-Uniform RANdom number generators
    </FONT></TD>
</TR></TABLE>"
EOX

my $links = <<EOX;
"<TABLE WIDTH='100%' BORDER='0' CELLSPACING='2'>
  <TR> 
    <TD WIDTH='17%' BGCOLOR='#1E426E'> 
      <A HREF='/anuran/index.html'>
      <DIV ALIGN='CENTER'>
      <FONT FACE='Arial, Helvetica, sans-serif' SIZE='2' COLOR='#FFFFFF'>
      <B>ANURAN</B></FONT></A></DIV>
    </TD>
    <TD WIDTH='17%' BGCOLOR='#BCD2EE'> 
      <DIV ALIGN='CENTER'>
      <FONT FACE='Arial, Helvetica, sans-serif' SIZE='2' COLOR='#1E426E'>
      <B>Code Generator</B></FONT></DIV>
    </TD>
    <TD WIDTH='17%' BGCOLOR='#1E426E'> 
      <DIV ALIGN='CENTER'>
      <A HREF='/anuran/project.html'>
      <FONT FACE='Arial, Helvetica, sans-serif' SIZE='2' COLOR='#FFFFFF'>
      <B>Project</B></FONT></A></DIV>
    </TD>
    <TD WIDTH='17%' BGCOLOR='#1E426E'> 
      <DIV ALIGN='CENTER'>
      <A HREF='/arvag/software.html'>
      <FONT FACE='Arial, Helvetica, sans-serif' SIZE='2' COLOR='#FFFFFF'>
      <B>Other Software</B></FONT></A></DIV>
    </TD>
    <TD WIDTH='17%' BGCOLOR='#1E426E'> 
      <DIV ALIGN='CENTER'>
      <A HREF='/arvag/index.html'>
      <FONT FACE='Arial, Helvetica, sans-serif' SIZE='2' COLOR='#FFFFFF'>
      <B>ARVAG</B></FONT></A></DIV>
    </TD>
    <TD WIDTH='17%' BGCOLOR='#1E426E'> 
      <DIV ALIGN='CENTER'>
      <A HREF='mailto:anuran\\\@statistik.wu-wien.ac.at'>
      <FONT FACE='Arial, Helvetica, sans-serif' SIZE='2' COLOR='#FFFFFF'>
      <B>Feedback</B></FONT></A></DIV>
    </TD>
    <TD WIDTH='1%'>&nbsp;</TD>
  </TR>
</TABLE>"
EOX

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

use CGI qw/:standard -nosticky *blockquote *font/;
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

\$step = param('step');
Delete('step');

\$distr = param('distribution');
\$Stdform = param('Stdform');

\$n_req = \$data_distr{\$distr}{'=N_REQ'};
\$n_tot = \$data_distr{\$distr}{'=N_TOT'};

# ----------------------------------------------------------------
# We need a frog

if (param('frog')) {
    \$frog = param('frog');
}
else {
    \$frog = sprintf("%02d",int(rand $n_frogs) + 1);
}

\$img_frog = 
    "&nbsp;&nbsp;&nbsp;".
    "<img SRC=\\"$frog_dir/frog\$frog.jpg\\"".
    "ALT=\\"[Miran: 'A Frog']\\"".
    "ALIGN = top>";

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
	start_html(-title  =>'ANURAN Test Page',
		   -author =>'anuran\@statistik.wu-wien.ac.at',
		   -meta   =>{'keywords'=>'nonuniform random number generator',
			      'copyright'=>'copyright 2001 Josef.Leydold\@statistik.wu-wien.ac.at'},
		   -BGCOLOR=>'#FFFFF5', 
		   -TEXT   =>'#000000',
		   -LINK   =>'#000000',
		   -VLINK  =>'#000000',
		   -ALINK  =>'#000000',
		   -MARGINWIDTH =>'5',
		   -MARGINHEIGHT=>'5',
		   -LEFTMARGIN  =>'5',
		   -TOPMARGIN   =>'5'),

        # Make bar with links
	$logo,
	$links,

	# Heading
	h1('Automatic Code Generator'),
	h2('Experimental version'),
	b('Resulting code is reliable but should used with care!'),br
	b('This code generator cannot produce generators for all distributions / parameters.'),br,

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
	\$warning = p().b('You must select a distribution first!');
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
				  -action=>url()),
			popup_menu(-name=>'distribution',
				   -values=>\\\@menue_distributions,
				   -labels=>\\\%labels_menue_distributions),
			'&nbsp;&nbsp;&nbsp;',
			hidden('step','2'),
			hidden('frog',\$frog),
			submit('Continue'),
			\$img_frog,
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

	return;
    }

# ................................................................
# Distribution selected
# ................................................................
    
    print 
	'Step 2: ',
	b('Parameters for '.\$data_distr{\$distr}{'=NAME'}).br();


    # Store id of parameter
    for (my \$i=0; \$i < \$n_tot; \$i++) {
	\$param_by_name{\$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'}} = \$i;
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
		# first transform string to number (if possible)
		if (param("Std_param_\$i") =~ /\\d+/) {
		    param("Std_param_\$i",param("Std_param_\$i")+0);
		}
		# copy parameter into temporary array
		\$param[\$i] = param("Std_param_\$i");
	    }
	}
	else {
	    # Non-Standardform
	    \$n_param = \$n_tot;
	    for (my \$i=0; \$i < \$n_tot; \$i++) {
		# first transform string to number (if possible)
		if (param("param_\$i") =~ /\\d+/) {
		    param("param_\$i",param("param_\$i")+0);
		}
		# copy parameter into temporary array
		\$param[\$i] = param("param_\$i");
	    }
	}

	# String to be printed on screen
	for (my \$i=0; \$i < \$n_tot; \$i++) {
	    \$print_param[\$i] = \$param[\$i];
	}
	
	# Check input for empty entries and transform string to numbers
	for (my \$i=0; \$i < \$n_tot; \$i++) {
	    unless (\$param[\$i] =~ /\\d+/) {
		\$step = 2;
		\$print_param[\$i] = '<FONT COLOR="red"><BLINK>missing</BLINK></FONT> ';
	    }
	}

	# Check validity of parameters
	for (my \$i=0; \$i < \$n_tot; \$i++) {
	    my \$lower = \$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=LOWER'};
	    if (defined \$lower) {
		if ( (\$lower =~ /\\d+/ and \$lower > \$param[\$i]) or
		     (\$lower !~ /\\d+/ and \$param[\$param_by_name{\$lower}] > \$param[\$i]) ) {
		    \$step = 2;
		    \$print_param[\$i] = ' <FONT COLOR="red">'
			.\$print_param[\$i]
			    .' <BLINK>invalid</BLINK></FONT> ';
		}
	    }
	    my \$upper = \$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=UPPER'};
	    if (defined \$upper) {
		if ( (\$upper =~ /\\d+/ and \$upper < \$param[\$i]) or
		     (\$upper !~ /\\d+/ and \$param[\$param_by_name{\$upper}] < \$param[\$i]) ) {
		    \$step = 2;
		    \$print_param[\$i] = ' <FONT COLOR="red">'
			.\$print_param[\$i]
			    .' <BLINK>invalid</BLINK></FONT> ';
		}
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
	    startform(-method => 'GET',
		      -action => url());

	if (\$n_req < \$n_tot) {
	    # Standardform exists

	    if (\$Stdform ne 'no') {
		print '\t<INPUT TYPE="radio" NAME="Stdform" VALUE="yes" CHECKED> Standardform'."\\n";
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
		print '\t<INPUT TYPE="radio" NAME="Stdform" VALUE="no" '.\$checked.'> Non-Standardform'."\\n";
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
	    print hidden('Stdform','yes'),"\\n";
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
	    hidden('step','3'),"\\n",
	    hidden('frog'),
	    hidden('distribution'),
	    submit('Continue'),
	    \$img_frog,"\\n";


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


    # Compute standard domain
    \$std_left = \$data_distr{\$distr}{'=LEFT'};
    if ( \$std_left !~ /\\d+/ and \$std_left ne '-infinity') {
	\$std_left = \$param[\$param_by_name{\$data_distr{\$distr}{'=LEFT'}}];
    }
    \$std_right = \$data_distr{\$distr}{'=RIGHT'};
    if ( \$std_right !~ /\\d+/ and \$std_right ne 'infinity') {
	\$std_right = \$param[\$param_by_name{\$data_distr{\$distr}{'=RIGHT'}}];
    }

    # Get domain
    if (!defined(param('left')) or param('left') < \$std_left) {
	param('left',\$std_left);
    }
    if (!defined(param('right')) or param('right') > \$std_right) {
	param('right',\$std_right);
    }

    # Transform string into number
    param('left',param('left')+0);
    param('right',param('right')+0);

    # Truncated domain ?
    if (param('truncated') ne 'yes' && \$step == 3) {
	\$step = 5; # nothing to do --> skip
    }
    if (param('truncated' ne 'yes')) {
	param('truncated','no');
    }

# ................................................................
# CASE: Distribution selected, Parameters and Domain given
# ................................................................

    if (\$step > 3) {

	# check domain
	if (param('left') > param('right')) {
	    \$step = 3;
	    print
		blockquote('<FONT COLOR="red"><BLINK>domain invalid</BLINK> (left > right)</FONT>');
	}

	else {
	    # print domain
	    print 
		blockquote( ' domain = [ ',
			    param('left'), 
			    ',',
			    param('right'),
			    ']',
			    '&nbsp;&nbsp;&nbsp;',
			    ((param('truncated')eq'yes') ? '(truncated)' : '') );
	}
    }

# ................................................................
# CASE: Distribution selected and Parameters given, input domain
# ................................................................

    if (\$step == 3) {

	print
	    start_blockquote(),
	    startform(-method => 'GET',
		      -action => url()),
	    
	    'left = ',
	    textfield(-name=>'left',
		      -size=>10,
		      -maxlength=>10),
	    '&nbsp;&nbsp;&nbsp;',
	    'right = ',
	    textfield(-name=>'right',
		      -size=>10,
		      -maxlength=>10),
	    p(),
	    hidden('step','5'),"\\n",
	    hidden('frog'),
	    hidden('distribution'),
	    hidden('truncated'),
	    hidden('Stdform');

	if (\$Stdform ne 'no') {  # i.e. 'yes'
	    for (my \$i = 0; \$i < \$n_tot; \$i++) {
		print hidden("Std_param_\$i");
	    }
	}
	if (\$Stdform ne 'yes') {  # i.e. 'no'
	    for (my \$i = 0; \$i < \$n_tot; \$i++) {
		print hidden("param_\$i");
	    }
	}
	
	print
	    submit('Continue'),
	    \$img_frog,"\\n",
	    end_form(),
	    end_blockquote();
    }

# ................................................................

    # ruler
    print p().hr().p();


} # end of anuran_domain_distribution()


# ----------------------------------------------------------------
# Step 4: Properties of generator
# ----------------------------------------------------------------

sub anuran_properties_generator
{

# ................................................................
# 
# ................................................................

    if (\$step < 4) {
	print

	    # heading
	    font({-color=>"$Gray"},
		 'Step 4: '.b('Properties of generator')),
	    br(),

	    # ruler
	    p().hr().p();

	return;
    }

# ................................................................
# 
# ................................................................

    print 
	'Step 4: ',
	b('Properties of generator for '.\$data_distr{\$distr}{'=NAME'}).br();

# ................................................................
# 
# ................................................................

    if (\$step > 4) {
	;
    }

# ................................................................
# 
# ................................................................

    if (\$step == 4) {
	;
    }

# ................................................................

    # ruler
    print p().hr().p();

} # end of anuran_properties_generator()


# ----------------------------------------------------------------
# Step 5: Select a programming language
# ----------------------------------------------------------------

sub anuran_language
{

# ................................................................
# 
# ................................................................

    if (\$step < 5) {
	print
	    
	    # heading
	    font({-color=>"$Gray"},
		 'Step 5: '.b('Programming language')),
	    br(),
	    
	    # ruler
	    p().hr().p();

	return;
    }

# ................................................................
# 
# ................................................................

    if (\$step > 5) {
	print 
	    'Step 5: ',
	    b('Programming language: '.param('language')).br();
    }

# ................................................................
# 
# ................................................................

    if (\$step == 5) {
	print 
	    startform(-method => 'GET',
		      -action => url()),
	    'Step 5: ',
	    b('Programming language: '),
	    popup_menu(-name=>'language',
		       -values=>['C','FORTRAN','JAVA']),
	    hidden('step','6'),"\\n",
	    hidden('frog'),
	    hidden('distribution'),
	    hidden('left'),
	    hidden('right'),
	    hidden('truncated'),
	    hidden('Stdform');

	if (\$Stdform ne 'no') {  # i.e. 'yes'
	    for (my \$i = 0; \$i < \$n_tot; \$i++) {
		print hidden("Std_param_\$i");
	    }
	}
	if (\$Stdform ne 'yes') {  # i.e. 'no'
	    for (my \$i = 0; \$i < \$n_tot; \$i++) {
		print hidden("param_\$i");
	    }
	}
	
	print 
	    submit('Continue'),
	    \$img_frog,
	    br(),"\\n";

	print start_blockquote();
	print radio_group(-name => 'codetype',
			  -values=>['generator','demo'],
			  -default=>'generator',
			  -rows=>2,
			  -columns=>1,
			  -labels=>{'generator'=>' Generator only',
				    'demo'=>' Make complete demo version'});
	print end_blockquote();

	print endform();
    }

# ................................................................

    # ruler
    print p().hr().p();

} # end of anuran_language()


# ----------------------------------------------------------------
# Step 6: Finally the code
# ----------------------------------------------------------------

sub anuran_code
{

# ................................................................
# 
# ................................................................

    if (\$step < 6) {
	print

	    # heading
	    font({-color=>"$Gray"},
		 'Step 6: '.b('Generator code for distribution')),
	    br(),

	    # ruler
	    p().hr().p();

	return;
    }

# ................................................................
# 
# ................................................................

    print 
	'Step 6: ',
	b('Generator for '.\$data_distr{\$distr}{'=NAME'});

    print '(Complete demo version)' if param('codetype') eq 'demo';
    print \$img_frog,br(),"\\n";

# ................................................................
# 
# ................................................................

    if (\$step >= 6) {

# ................................................................
# Make acg command
# ................................................................

	# check also for command line
	my \$command_ok = 1;

	# Command
	my \$acg_query = "$ACG";

	# Distribution
	\$acg_query .= " -d \$distr";
	\$command_ok = 0 if \$distr =~ /[^a-zA-Z]/;

	# List of parameters
	if (\$n_param) {
	    \$acg_query .= " -p \\\"\$param[0]";
	    for (my \$i=1; \$i < \$n_param; \$i++) {
		\$acg_query .= " \$param[\$i]";
		\$command_ok = 0 if \$param[\$i] =~ /[^\\d\\.\\-\\+eE]/;
	    }
	    \$acg_query .= "\\\"";
	}
	
	# Domain
	if (param('truncated')) {
	    \$acg_query .= ' -D "'.param('left').' '.param('right').'"';
	    \$command_ok = 0 if param('left') =~ /[^\\d\\.\\-\\+eEiInNfF]/;
	    \$command_ok = 0 if param('right') =~ /[^\\d\\.\\-\\+eEiInNfF]/;
	}

	# Programming language
	\$acg_query .= ' -l '.param('language');
	\$command_ok = 0 if param('language') =~ /[^a-zA-Z]/;
	
	# Make complete demo version ?
	\$acg_query .= ' -M' if param('codetype') eq 'demo';

# ................................................................
# Print result
# ................................................................

	if (\$command_ok) {
	    print '<PRE>';
	    print `\$acg_query`;
	    print '</PRE>';
	}

# ................................................................
# Exit status
# ................................................................

	my \$status;

	if (\$command_ok) {
	    if (\$?) {
		# error
		\$status = "failed"; 
		print p(),
		font({-color=>"red"},
		     b("Sorry. Cannot make generator for given distribution / parameters.")),
		p();
	    }

	    else {
		# successfull
		\$status = "ok"; 
	    }
	}

	else { # the composed command line might execute unsecure code
	    \$status = "dangerous command, not executed";
		print p(),
		font({-color=>"red"},
		     b("Sorry. Cannot make generator. Internal error.")),
		p();
	}

# ................................................................
# Make Entry into log file
# ................................................................

	my \$log;

	# date
	\$log .= "date = ".(scalar localtime)."\\n";

	# client
        \$log .= "client = ".remote_host()."  (".remote_addr().")\\n";

	# name of distribution
	\$log .= "distribution = \$distr\\n";

	# list of parameters (defaults settings start with '('
	for (my \$i=0; \$i < \$n_tot; \$i++) {
	    \$log .= ((\$i >= \$n_param) ? '(' : '')
		  .  "param[\$i] = "
		  .  \$data_distr{\$distr}{'=FPARAMS'}[\$i]{'=NAME'}
	          .  " = "
		  .  \$print_param[\$i]."\\n";
	}
	
	# Domain
	\$log .= 'domain = [ '.param('left').', '.param('right')." ]\\n";

	# Programming language
	\$log .= 'language = '.param('language')."\\n";

	# ACG code
	\$log .= "command = \$acg_query\\n";

	# Exit code
	\$log .= "status = \$status\\n";

	# Use blank line as separator
	\$log .= "\\n";

	# Write into log file
	# (This is not save !!)
	open LOG, ">>$anuran_log_file";
	print LOG \$log;
	close LOG;

    }

# ................................................................

    # ruler
    print p().hr().p();

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
