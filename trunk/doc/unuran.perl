#!/usr/bin/perl

# $ON = 1, wenn der "Dokumentations-Bereich" gefunden ist
$ON = 0;

@TYPES = ("void", "int", "double", "struct", "float", "long", "char", "short", "unsigned", "signed");


# file fuer  output:
open(OUTFILE, ">zzz");

while($_ = <>)
{
 
    # bestimmen des dokumentationsbereiches
    if ($_ =~/Routines.for.user.interface/){
	$ON = 1;}
    if ($ON == 1 && $_ =~/-------------------/){
	$ON = 0;}
# Suche Function und DEfinitionszeilen
  foreach $type (@TYPES){
 
      if ( $ON == 1 && $_=~/^($type)/){    # suche Funktionszeile
        print OUTFILE;
      }

  }



}

  
  print OUTFILE " Funktioniert doch \n";

