#!/usr/bin/perl

# $ON = 1, wenn der "Dokumentations-Bereich" gefunden ist
$ON = 0;

# struct herausgenommen
@TYPES = ("void", "int", "double", "float", "long", "char", "short", "unsigned", "signed");


# file fuer  output:
open(OUTFILE, ">zzz");

while($_ = <>)
{
 
    # bestimmen des dokumentationsbereiches
    if ($_ =~/\/\* Routines for user interface/){
	$ON = 1;}
    if ($ON == 1 && $_ =~/\/\*-------------------/){
	$ON = 0;}
# Suche Function und DEfinitionszeilen

   if ( $ON == 1 && $_ =~/^\s*struct\s*unur.*(unur_[a-zA-z_]*)\(/ ){
       print OUTFILE $1, "\n";
       print OUTFILE $_, "\n";
   }


  foreach $type (@TYPES){
     # suche Funktionszeile
      
      if ( $ON == 1 && $_ =~/^\s*($type)\s*(unur_[a-zA-Z_]*)/){
       print OUTFILE $2 , "\n";
       print OUTFILE $_ , "\n";
      }

  }



}
