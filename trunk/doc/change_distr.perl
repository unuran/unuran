#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl script generating a list of changable properties of
# distributions   
# 
# Call:  ./distr.perl ../src/methods/*.c
# Input:  c-files with code to change pdf, mode, ... 
# Output: list fo properties
# 
# E.JANKA  September 2000
# $Id$
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# file for  output:
open(OUTFILE, ">qstart_change_distr.texi");

print OUTFILE "\@itemize \@bullet\n";

while($_ = <>)
{ 
    chomp;

       if ( ($distr) = ($_ =~/^unur_distr_cont_set_(\w*)/) ){
	   print OUTFILE "\@item ", $distr, "\n";
          }

}
    print OUTFILE "\@end itemize\n";









