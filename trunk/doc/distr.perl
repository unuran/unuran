#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl script for generating a list of preset distributions   
# 
# Call:  ./distr.perl ../src/distributions/*.c
# Input:  c-files with definition of the distributions 
# Output: list of distributions
# 
# E.JANKA  September 2000
# $Id$
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# file for  output:
open(OUTFILE, ">qstart_distributions.texi");

print OUTFILE "\@itemize \@bullet\n";

while($_ = <>)
{ 
    chomp;

       if ( ($distr) = ($_ =~/^unur_distr_(\w*)/) ){
	   print OUTFILE "\@item ", $distr, "\n";
          }

}
    print OUTFILE "\@end itemize\n";









