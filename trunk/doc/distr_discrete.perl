#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl script generating a list of preset distributions
# available 
# 
# Call:  ./distr.perl ../src/distributions/*.c
# Input:  c-files with code of distributions 
# Output: description of distributions in texi-format
# 
# E.JANKA  September 2000
# $Id$
#
# This script scans c-files for description of distributions
# and produces texi-output
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# within a definition of a distribution
$DISTR = 0; 
# within a list of parameters
$PARAM = 0;

$distrcount = 0;

# file for  output:
open(OUTFILE, ">qstart_d_distributions.texi");


while($_ = <>)
{ 
    chomp;

    # first key word is "distr:"
    if ( ($distribution, $reference) =
         $_ =~ /^\s*\*\s*distr:\s*([^\[]*)\s+\[?([^\]]*)\]?\s*\*\s*$/ ){
        $DISTR = 1;
        $distrcount++;
	print OUTFILE "\@subsection ", $distribution, "\n";
  
        # print "distr:" explanation
	$_ = <>;
	chomp;
	($cut) = $_ =~ /^\s*\*(.*)\*\s*$/;
        while ($cut !~ /^\s*$/ && $cut !~ /^\s*\w+:/){
	    print OUTFILE "\@\* ", $cut, "\n";
	    $_ = <>;
	    chomp;
	    ($cut) = $_ =~ /^\s*\*(.*)\*\s*$/;
	}
    }
  
    # key words after "distr:"
    while ( $cut !~ /^\s*parameters:/ && $cut !~ /\*{10,}/ 
	    && $cut !~ /\.{10,}/){

	if ($DISTR == 1 && $cut =~ /^\s*(\w+):(.*?)\s*$/){
	    print OUTFILE "\@\*", $1, ": ", $2, "\n";

	    $_ = <>;
	    chomp;
	    ($cut) = $_ =~ /^\s*\*(.*)\*\s*$/;

            while ($cut !~ /^\s*$/ && $cut !~ /^\s*\w+:/){
		print OUTFILE "\@\* ", $cut, "\n";
	        $_ = <>;
	        chomp;
	        ($cut) = $_ =~ /^\s*\*(.*)\*\s*$/;
            }	    
	}
        else {
	    $_ = <>;
	    chomp;
	    ($cut) = $_ =~ /^\s*\*(.*)\*\s*$/;
        }

    }
 
    # lines with 10 or more starts end region with explanation
    if ( $_ =~ /\*{10,}/ || $_ =~ /\.{10,}/ ){
        # end multitable
	if ($PARAM == 1){
	    print OUTFILE "\@end multitable\n";
	}
	$DISTR = 0;
	$PARAM = 0;
    }


    if ($DISTR == 1 && $_ =~ /parameters:\s*(.*?)\s*\*\s*/){
	$DISTR = 0;
	$PARAM = 1;
        # begin multitable
	print OUTFILE "\@\*Parameters: $1";
	print OUTFILE "\@multitable \@columnfractions .04 .3 .2 .3\n";
	print OUTFILE
              "\@item Nr.\@tab Name \@tab Standard value \@tab Type\n";
   }   


    # parameters:
    if ( $PARAM == 1
	 && ( ($nr, $def) = $_ =~ /^\s*\*\s*([0-9]+):\s*(.*?)\s+\*\s*$/) ){

        ($def, $descr) = split  /\.{3}\s/, $def;
        ($def, $leer,$strd)  = $def =~ /([^\(]*)\s*(\(([^\)]+)\))?/; 

	print OUTFILE "\@item ", $nr, " \@tab ", $def, 
                      " \@tab ", $strd, " \@tab ", $descr, "\n";
    }


}




