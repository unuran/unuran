#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl script for generating a list of preset distributions   
# 
# Call:  ./distr.perl ../src/distributions/*.h
# Input:  c-files with definition of the distributions 
# Output: list of distributions
# 
# E.JANKA  September 2000
# $Id$
#
#  This script scans h-files for "type unur_distr_xxx()"
#  and produces texi output
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# known data types 
@TYPES = ("UNUR_PAR", "UNUR_GEN", "UNUR_DISTR", "struct", "void", "int", "double", "float", "long", "char", "short", "unsigned", "signed");

# file for  output:
open(OUTFILE, ">qstart_distributions.texi");

print OUTFILE "\@itemize \@minus\n";

while($_ = <>)
{ 
    chomp;

      # if ( ($distr) = ($_ =~/\s+\*unur_distr_(\w*)/) ){
	#   print OUTFILE "\@item ", $distr, "\n";
         # }

    foreach $type (@TYPES){

       if ( $_ =~/^\s*($type\s+\*unur_distr_\w+)\s*\((.*\))\s*;/ ){
           $CENABLE = 1;
           $DECL = $1;   # string before the braces
           $FUNC = $2;   # string between the braces 
           $DECL  =~ /(.*(\s+?|\*))(unur_distr_(\w+))/;
           $DECL1 = $1;
	   $DECL2 = $4;
           $DECL3 = $3;

           print OUTFILE  "\n\@item \@strong{", $DECL2 , "}\@*\n";
 
           print OUTFILE "\@code{", $DECL1 , "\@b{", $DECL3,"}("; 
           while ($FUNC =~ /(.*?)(\w*?)\s*?(,|\))/g){
	      print OUTFILE $1,"\@var{", $2,"}", $3;
	   }
           print OUTFILE "}\n"; 
       }

   }



}
    print OUTFILE "\@end itemize\n";













