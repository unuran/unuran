#!/usr/bin/perl


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl spript zur automatischen Erstellung der   
# (Methoden-)Funktions- Referenz von UNURAN im TEX-Info Format
# 
# 
# Aufruf:  ./unuran.perl ../src/methods/*.h
# Input:  
# Output:
# 
# E.JANKA und G.TIRLER  August 2000
# $Id$
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 
# $ON = 1, wenn der "Dokumentations-Bereich" gefunden ist
$ON = 0;

# Schalter, um Kommentarzeilen hinter Funktionsdefinitionen zu erkennen 
$KOMMENT = 0;

#Name der aktuellen Methode (von unuran)
$MName = "a";  # mit nicht vorkommendem Namen initialisiert

$ERST = 1;  # garantiert,dass keine itemize-umgebung beendet
            # wird bevor es ueberhaupt eine gibt
$ENDE = 0;  # ueberprueft, ob itemize-umgebung gesetzt wird
            # und garantiert, dass letztes itemize beendet wird

# 
@TYPES = ("struct", "void", "int", "double", "float", "long", "char", "short", "unsigned", "signed");



#  output file:
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
   foreach $type (@TYPES){

       # Wenn neuer MName vorkommt, dann neues Kapitel
       if ( $ON == 1 && $_ =~/^\s*($type.*)\s*\((.*\))\s*;/){
          # neues Texinfo-Kapitel, wenn sich methode ($MName) aendert
          /\*unur_(\w*)_/; 
          if ($type eq struct && $MName ne $1){
             $MName = $1;
             if ($ERST == 0){
		print OUTFILE "\@end itemize\n\n";
             }
             else{
		 $ERST = 0;
                 $ENDE = 1;
	     }
	     print OUTFILE "\n\n\@subsection ", $MName, "\n\n";
             print OUTFILE "\@itemize \@minus{} \n";
          }
       }

       # suche Funktionszeile
       if ( $ON == 1 && $_ =~/^\s*($type.*)\s*\((.*\))\s*;/){
           #Kommentar zur funktion moeglich wenn $KOMMENT = 1
           $KOMMENT = 1;
           $DECL = $1;   # vor der oeffnenden Funktionsklammer
           $FUNC = $2;   # zwischen den Funktionsklammern 
           $DECL  =~ /(.*(\s+?|\*))(\w+)/;
           $DECL1 = $1;
	   $DECL2 = $3;

           print OUTFILE  "\n\@item \@strong{", $DECL2 , "}\@*\n";
 
           print OUTFILE "\@code{", $DECL1 , "\@b{", $DECL2,"}("; 
           while ($FUNC =~ /(.*?)(\w*?)\s*?(,|\))/g){
	      print OUTFILE $1,"\@var{", $2,"}", $3;
	   }
           print OUTFILE "}\@\*\n"; 

       }
  }


# Suche zugehoerige Kommentare
  if($KOMMENT == 1 && $_ =~/^\/\*(.*)(\*\/)$/){
	  print OUTFILE $1, "\@\*\n";
  }



}

if ($ENDE == 1){
  print OUTFILE "\@end itemize";
}















