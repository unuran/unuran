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
my $PHP_file = "../../public_html/anuran.php";

# Gray Color for disabled text
my $Gray = "\#D0D0D0";

# ----------------------------------------------------------------
# List of distributions
my $DISTR = read_PDFdata();

# For description of data fields in this list see file `readPDF.pl'.

# ----------------------------------------------------------------

# Make C file with code generator
make_web();

# ----------------------------------------------------------------
# End

exit 0;

# ----------------------------------------------------------------
#
# Subroutines
#
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Make php file 

sub make_web
{
    my $PHP;

# ................................................................
# HTML header

    $PHP .= <<EOX;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<HTML>
<HEAD>
<TITLE>ANURAN Test Page</TITLE>
</HEAD>
<BODY BGCOLOR="#FFFFFF">

EOX

# ................................................................
# Header

    $PHP .= <<EOX;
<H1>Automatic Code Generator</H1>
<B>(Experimental version)</B>

<P><HR><P>

EOX

# ................................................................
# Select List for continuous distributions (required for step 1)
    
    my $HTML_distr_select;

    $HTML_distr_select .=
	"\t echo \"    <SELECT name=\\\"distribution\\\" size=1>\\n\";\n".
	"\t echo \"      <OPTION value=\\\"--none--\\\">-- Select a distribution --</OPTION>\\n\";\n";

    foreach my $d (sort keys %{$DISTR}) { 
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$HTML_distr_select .=
	    "\t echo \"      <OPTION value=".
	    "\\\"$d\\\">".
	    $DISTR->{$d}->{"=NAME"}.
	    "</OPTION>\\n\";\n";
    }

    $HTML_distr_select .=
	"\t echo \"    </SELECT>\\n\";";

# ................................................................
# Get name of continuous distribution (required for step 1)
    
    my $get_distr_name;

    $get_distr_name .=
	"\t # Get name and PDF of distribution\n".
	"\t switch (\$distribution) {\n";

    foreach my $d (sort keys %{$DISTR}) { 
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$get_distr_name .= 
	    "\t   case \"$d\":\n". 
	    "\t\t \$distr_name = \"".$DISTR->{$d}->{"=NAME"}."\";\n".
	    "\t\t \$distr_PDF = \"".format_PDF($DISTR->{$d}->{"=DOC"}->{"=PDF"})."\";\n".
	    "\t\t break;\n";
    }

    $get_distr_name .=
	"\t }\n";

# ................................................................
# Step 1: Select a distribution

    $PHP .= <<EOX;
<?php
 # Step 1

 if (\$step <= 0 || \$distribution == "--none--") {
\t # no distribution selected till now

\t echo "Step 1: <STRONG>Distribution</STRONG><BR>\\n";
\t echo "<BLOCKQUOTE>\\n";

\t if (\$step >= 1) {
\t\t # no distribution selected at step 1
\t\t \$step = 0;
\t\t echo "  <P><STRONG>You must select a distribution first!</STRONG></P>\\n";
\t }

\t echo "  <FORM method=\\"GET\\" action=\\"anuran.php\\">\\n";
$HTML_distr_select
\t echo "    &nbsp;&nbsp;&nbsp;\\n";
\t echo "    <INPUT type=\\"hidden\\" name=\\"step\\" value=\\"1\\">\\n";
\t echo "    <INPUT type=\\"submit\\" value=\\"Continue\\">\\n";
\t echo "  </FORM>\\n";

\t echo "</BLOCKQUOTE>\\n";
 }

 else {
\t #distribution selected

$get_distr_name
\t echo "Step 1: <STRONG>";
\t echo "\$distr_name";
\t echo "</STRONG><BR>\\n";
\t echo "<BLOCKQUOTE>\\n";
\t echo "PDF(x) = [const] * \$distr_PDF\\n";
\t echo "</BLOCKQUOTE>\\n";
 }
?>

<P><HR><P>

EOX

# ................................................................
# Get data about parameters for distribution (required for step 2)

    my $get_distr_PDFparams;

    $get_distr_PDFparams .=
	"\t switch (\$distribution) {\n";

    foreach my $d (sort keys %{$DISTR}) { 
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$get_distr_PDFparams .= 
	    "\t   case \"$d\":\n".
	    scan_FPARAM($DISTR->{$d}->{"=DOC"}->{"=FPARAM"}).
	    "\t\t break;\n";
    }

    $get_distr_PDFparams .= 
	"\t }\n";

# ................................................................
# Step 2: Insert parameters for PDF

    $PHP .= <<EOX;
<?php
 # Step 2

 if (\$step < 1) {
\t # no distribution selected till now    
\t echo "<FONT COLOR=\\"$Gray\\">\\n";
\t echo "Step 2: <STRONG>Parameters for distribution</STRONG><BR>\\n";
\t echo "</FONT>\\n";
 }

 else {

$get_distr_PDFparams

\t if (\$step == 1) {
\t\t # Insert parameters

\t\t echo "Step 2: <STRONG>Parameters for \$distr_name</STRONG><BR>\\n";
\t\t echo "<BLOCKQUOTE>\\n";
\t\t echo "  <FORM method=\\"GET\\" action=\\"anuran.php\\">\\n";

\t\t if (\$distr_n_requir < \$distr_n_total) {
\t\t\t # Standardform exists

\t\t\t echo "    <INPUT TYPE=\\"radio\\" NAME=\\"Stdform\\" VALUE=\\"yes\\" checked> Standardform\\n";
\t\t\t echo "    <BLOCKQUOTE>\\n";
\t\t\t for (\$i = 0; \$i < \$distr_n_requir; \$i++) {
\t\t\t\t echo "    \$distr_param_name[\$i] = \\n";
\t\t\t\t echo "    <INPUT TYPE=\\"text\\" NAME=\\"Std_param_\$i\\" SIZE=10 MAXLENGTH=10>\\n";
\t\t\t\t echo (\$distr_param_cond[\$i] ? "    &nbsp;&nbsp;&nbsp;(\$distr_param_cond[\$i])" : "");
\t\t\t\t echo "    <BR>\\n";
\t\t\t }
\t\t\t if (\$distr_n_total > \$distr_n_requir) {
\t\t\t\t echo "    <P><FONT COLOR=\\"$Gray\\">(&nbsp;&nbsp;&nbsp;\\n";
\t\t\t\t for (\$i = \$distr_n_requir; \$i < \$distr_n_total; \$i++) {
\t\t\t\t\t echo "    \$distr_param_name[\$i] = \$distr_param_default[\$i]\\n"; 
\t\t\t\t\t echo "&nbsp;&nbsp;&nbsp;";
\t\t\t\t }
\t\t\t\t echo ")</FONT></P>\\n";
\t\t\t }
\t\t\t echo "    </BLOCKQUOTE>\\n";

\t\t\t echo "    <INPUT TYPE=\\"radio\\" NAME=\\"Stdform\\" VALUE=\\"no\\"> Non-Standardform\\n";
\t\t\t echo "    <BLOCKQUOTE>\\n";
\t\t\t for (\$i = 0; \$i < \$distr_n_total; \$i++) {
\t\t\t\t echo "    \$distr_param_name[\$i] = \\n";
\t\t\t\t echo "    <INPUT TYPE=\\"text\\" NAME=\\"param_\$i\\" ";
\t\t\t\t echo ((\$i<\$distr_n_requir) ? "" : "VALUE=\\"\$distr_param_default[\$i]\\" ");
\t\t\t\t echo "SIZE=10 MAXLENGTH=10>\\n";
\t\t\t\t echo (\$distr_param_cond[\$i] ? "    &nbsp;&nbsp;&nbsp;(\$distr_param_cond[\$i])" : "");
\t\t\t\t echo "    <BR>\\n";
\t\t\t }
\t\t\t echo "    </BLOCKQUOTE>\\n";
\t\t }

\t\t else {
\t\t\t # No special standard form
\t\t\t echo "    <INPUT TYPE=\\"hidden\\" NAME=\\"Stdform\\" VALUE=\\"no\\">\\n";
\t\t\t for (\$i = 0; \$i < \$distr_n_total; \$i++) {
\t\t\t\t echo "    \$distr_param_name[\$i] = \\n";
\t\t\t\t echo "    <INPUT TYPE=\\"text\\" NAME=\\"param_\$i\\" SIZE=10 MAXLENGTH=10>\\n";
\t\t\t\t echo (\$distr_param_cond[\$i] ? "    &nbsp;&nbsp;&nbsp;(\$distr_param_cond[\$i])" : "");
\t\t\t\t echo "    <BR>\\n";
\t\t\t }
\t\t }

\t\t echo "    <P>\\n";
\t\t echo "    <INPUT type=\\"hidden\\" name=\\"step\\" value=\\"2\\">\\n";
\t\t echo "    <INPUT type=\\"hidden\\" name=\\"distribution\\" value=\\"\$distribution\\">\\n";
\t\t echo "    <INPUT type=\\"submit\\" value=\\"Continue\\">\\n";
\t\t echo "  </FORM>\\n";
\t\t echo "</BLOCKQUOTE>\\n";
\t }

\t else {
\t\t # Parameters given

\t\t \$distr_param = array();
\t\t if (\$Stdform == "yes") {
\t\t\t \$distr_n_params = \$distr_n_requir;
\t\t\t while (list (\$key, \$val) = each(\$HTTP_GET_VARS)) {
\t\t\t\t if ( preg_match ("/^Std_param_(\\d+)/", \$key, \$matches ) ) {
\t\t\t\t\t \$distr_param[\$matches[1]] = \$val;
\t\t\t\t }
\t\t\t }
\t\t\t for (\$i = \$distr_n_requir; \$i < \$distr_n_total; \$i++) {
\t\t\t\t \$distr_param[\$i] = \$distr_param_default[\$i];
\t\t\t }
\t\t }
\t\t else {
\t\t\t \$distr_n_params = \$distr_n_total;
\t\t\t while (list (\$key, \$val) = each(\$HTTP_GET_VARS)) {
\t\t\t\t if ( preg_match ("/^param_(\\d+)/", \$key, \$matches ) ) {
\t\t\t\t\t \$distr_param[\$matches[1]] = \$val;
\t\t\t\t }
\t\t\t }
\t\t }

\t\t echo "Step 2: <STRONG>Parameters for \$distr_name</STRONG><BR>\\n";
\t\t echo "<BLOCKQUOTE>\\n";
\t\t\t for (\$i = 0; \$i < \$distr_n_params; \$i++) {
\t\t\t\t echo "    \$distr_param_name[\$i] = \$distr_param[\$i]&nbsp;&nbsp;&nbsp;\\n";
\t\t\t }
\t\t\t if (\$distr_n_total > \$distr_n_params) {
\t\t\t\t echo "    <FONT COLOR=\\"$Gray\\">(&nbsp;&nbsp;&nbsp;\\n";
\t\t\t\t for (\$i = \$distr_n_requir; \$i < \$distr_n_total; \$i++) {
\t\t\t\t\t echo "    \$distr_param_name[\$i] = \$distr_param[\$i]&nbsp;&nbsp;&nbsp;\\n";
\t\t\t\t }
\t\t\t\t echo ")</FONT>\\n";
\t\t\t }
\t\t echo "</BLOCKQUOTE>\\n";

\t }
 }


?>

<P><HR><P>

EOX

# ................................................................
# End of HTML file
    $PHP .= <<EOX;
</BODY>
</HTML>

EOX

# ................................................................
# Print HTML file

    open FILE, ">$PHP_file" or die "cannot open file $PHP_file\n";
    print FILE $PHP;
    close FILE;

} # end of make_www()

# ----------------------------------------------------------------
# Format PDF for distribution
#

sub format_PDF {
    my $PDF = $_[0];

    # empty ? 
    return "" unless $PDF;

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
#

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
    my $plist_name = "array(";
    my $plist_def = "array(";
    my $plist_cond = "array(";
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
	$plist_name .= "\"$name\", ";
	$plist_def  .= ((length $default) ? $default : "\"*\"").", ";
	$plist_cond .= "\"$cond\", ";

	# increment counter for parameters
	++$n_total;
	++$n_optional if length $default;
    }

    # number of parameters for standard form
    my $n_required = $n_total - $n_optional;


    # store number of parameters
    $params .= "\t\t \$distr_n_requir = $n_required;\n";
    $params .= "\t\t \$distr_n_total = $n_total;\n";

    # close array definition
    $plist_name =~ s/,\s*$/\)/;
    $plist_def  =~ s/,\s*$/\)/;
    $plist_cond =~ s/,\s*$/\)/;

    $params .= "\t\t \$distr_param_name = $plist_name;\n";
    $params .= "\t\t \$distr_param_default = $plist_def;\n";
    $params .= "\t\t \$distr_param_cond = $plist_cond;\n";

    # Return result
    return $params;

} # end of scan_FPARAM()

# ----------------------------------------------------------------

