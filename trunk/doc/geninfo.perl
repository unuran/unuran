#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl script generating a list of changable properties of
# distributions   
# 
# Call:  ./distr.perl ../src/methods/x_gen.h
# Input:  h-files with code relavant to generator 
# Output: list fo properties
# 
# E.JANKA  November 2000
# $Id$
#
# This script scans h-files for "type unur_xxx()"
# and produces texi-output
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# known data types 
@TYPES = ("UNUR_PAR", "UNUR_GEN", "UNUR_DISTR", "struct", "void", "int", "double", "float", "long", "char", "short", "unsigned", "signed");

# file for  output:
open(OUTFILE, ">qstart_geninfo.texi");

print OUTFILE "\@itemize \@minus\n";

while($_ = <>)
{ 
    chomp;
 
    foreach $type (@TYPES){

	if ( $_ =~/^\s*($type\s+\*?unur_[A-Za-z_]+)\s*\((.*\))\s*;/ ){
           $CENABLE = 1;
           $DECL = $1;   # string before the braces
           $FUNC = $2;   # string between the braces 
	   $DECL  =~ /(.*(\s+?|\*))(unur_([A-Za-z_]+))/;
           $DECL1 = $1;
	   $DECL2 = $4;
           $DECL3 = $3;

           print OUTFILE  "\n\@item \@strong{", $DECL2 , "}\@*\n";
 
           print OUTFILE "\@code{", $DECL1 , "\@b{", $DECL3,"}("; 
           while ($FUNC =~ /(.*?)(\w*?)\s*?(,|\))/g){
	      print OUTFILE $1,"\@var{", $2,"}", $3;
	   }
           print OUTFILE "}\n\@*";

           # print comments
           $_ = <>;
	   chomp;
	   while ($_ !~ /^\s*$/){
	       $print = join '', split /\s*\*\// , ( join '', split /\/\*/, $_);
               $print = join '', split /\(=>\)/ , $print;
	       print OUTFILE $print, "\n";
	       $_ = <>;
	       chomp;
	   }
       }

    }

}
    print OUTFILE "\@end itemize\n";










