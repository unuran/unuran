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
# Continuous distributions (required for step 1)
    
    my $HTML_distr_select;
    my $HTML_distr_print;

    foreach my $d (sort keys %{$DISTR}) { 
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$HTML_distr_select .=
	    "    echo \"      <OPTION value=".
	    "\\\"$d\\\">".
	    $DISTR->{$d}->{"=NAME"}.
	    "</OPTION>\\n\";\n";

	$HTML_distr_print .= 
	    "        case (\"$d\"): echo \"".
	    $DISTR->{$d}->{"=NAME"}.
	    "\\n\"; break;\n";
    }

# ................................................................
# Step 1: Select a distribution


    $PHP .= <<EOX;
<?php

  echo "Step 1: <STRONG>Distribution</STRONG><BR>\\n";
  echo "<BLOCKQUOTE>\\n";

  if (\$step <= 0 || \$distribution == "--none--") {
    # no distribution selected till now

    if (\$step >= 1) {
      # no distribution selected at step 1
      echo "<P><STRONG>You must select a distribution first!</P>\\n";
    }

    echo "  <FORM method=\\"GET\\" action=\\"anuran.php\\">\\n";

    echo "    <SELECT name=\\"distribution\\" size=1>\\n";
    echo "      <OPTION value=\\"--none--\\">-- Select a distribution --</OPTION>\\n";

$HTML_distr_select
    echo "    </SELECT>\\n";
    echo "    &nbsp;&nbsp;&nbsp;\\n";
    echo "    <INPUT type=\\"hidden\\" name=\\"step\\" value=\\"1\\">\\n";
    echo "    <INPUT type=\\"submit\\" value=\\"Continue\\">\\n";
    echo "  </FORM>\\n";
  }

  else {
      # distribution selected
      switch (\$distribution) {
$HTML_distr_print
      }
  }

  echo "</BLOCKQUOTE>\\n";

?>

EOX


# ................................................................
# Continuous distributions
    
#    foreach my $d (sort keys %{$DISTR}) {
#	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

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
