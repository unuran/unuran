#!/usr/bin/perl


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl spript zur automatischen Erstellung der   
# Funktions- Referenz von UNURAN im TEX-Info Format
# 
# Aufruf:  ./unuran.perl ../src/methods/*.h
# Input:  
# Output:
# 
# E.JANKA und G.TIRLER  August 2000
# $Id$
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# $ON = 1, wenn der "Dokumentations-Bereich" gefunden ist
$ON = 0;

# Schalter, um Kommentarzeilen hinter Funktionsdefinitionen zu erkennen 
$KOMMENT = 0;

# struct herausgenommen
@TYPES = ("void", "int", "double", "float", "long", "char", "short", "unsigned", "signed");


# file fuer  output:
open(OUTFILE, ">qstart_function_reference.texi");

while($_ = <>)
{ 
    # Kommentare enden mit einere Leerzeile
#    if !($_ =

    # bestimmen des dokumentationsbereiches
    if ($_ =~/\/\* Routines for user interface/){
	$ON = 1;}
    if ($ON == 1 && $_ =~/\/\*-------------------/){
	$ON = 0;}
# Suche Function und Definitionszeilen

   if ( $ON == 1 && $_ =~/^\s*struct\s*unur.*(unur_[a-zA-z_]*)\(/ ){
       $KOMMENT = 1;
       print   OUTFILE "\@unnumberedsubsubsec ", $1, "\n\n";
       print OUTFILE  "\@code\{" ,$_ , "\}\n";
      
   }


  foreach $type (@TYPES){
     # suche Funktionszeile
      
      if ( $ON == 1 && $_ =~/^\s*($type)\s*(unur_[a-zA-Z_]*)/){
       $KOMMENT = 1;
       print OUTFILE  "\@unnumberedsubsubsec ", $2 , "\n\n";
       print OUTFILE "\@code\{" ,$_ , "\}\n";
      }

  }



}
