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

#Name der aktuellen Methode (von unuran)
$MName = "a";

# struct herausgenommen, wird extra abgehandelt
@TYPES = ("void", "int", "double", "float", "long", "char", "short", "unsigned", "signed");



# file fuer  output:
open(OUTFILE, ">qstart_function_reference.texi");

while($_ = <>)
{ 
    # Kommentare enden mit einere Leerzeile
    if ($_ =~ /^\s$/ ){
	$KOMMENT = 0;
    }

    # bestimmen des dokumentationsbereiches
    if ($_ =~/\/\* Routines for user interface/){
	$ON = 1;}
    if ($ON == 1 && $_ =~/\/\*-------------------/){
	$ON = 0;}

    chomp;



# Suche Funktion und Definitionszeilen
   if ( $ON == 1 && $_ =~/^\s*(struct\s*unur_[a-zA-Z_]*)\s*(\*unur_[a-zA-z_]*)(.*)/ ){
       #Kommentar zur funktion moeglich       
       $KOMMENT = 1;

       # neues Texinfo-Kapitel, wenn sich methode ($MName) aendert
       /\*unur_(\w*)_/; 
       if ($MName ne $1){
           $MName = $1;
	   print OUTFILE "\n\n\@subsection ", $MName, "\n\n";
       }
       /^\s*(struct\s*unur_[a-zA-Z_]*)\s*(\*unur_[a-zA-z_]*)(.*)/;

       print   OUTFILE "\@unnumberedsubsubsec ", $2, "\n\n";
       print OUTFILE "\@code\{", "\@i\{", $1,"\} ", "\@b\{",$2, "\}", $3, "\}\n";      
   }

   foreach $type (@TYPES){
     # suche Funktionszeile
        if ( $ON == 1 && $_ =~/^\s*($type)\s*(unur_[a-zA-Z_]*)(.*)/){
       #Kommentar zur funktion moeglich
       $KOMMENT = 1;
       print OUTFILE  "\n\@unnumberedsubsubsec ", $2 , "\n\n";
       print OUTFILE "\@code\{", "\@i\{", $1,"\} ", "\@b\{",$2, "\}", $3, "\}\n";
      }
  }


# Suche zugehoerige Kommentare
      if($KOMMENT == 1 && $_ =~/^\/\*(.*)(\*\/)$/){
	  print OUTFILE $1, "\n\n";
      }



}










