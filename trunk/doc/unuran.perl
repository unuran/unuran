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
# Das Perl-skrip durchsucht nach folgenden Schluesselworten      
# und und erzeugt Output wiefolgt:
# Folgende Schluesselwoerter beginnen mit '=...'
# Schluesselwort darf nur von Leerzeichen oder '/*'
# eingeleitet werden
#
#
# =METHODS name [LAngtext]
#          allg. Beschreibung der Methode in Header file
#
#          name      ... Name der Methode, Nur ein Wort!
#          [Langtext]... optional, Lange Beschreibung
#          Beispiel:
#             =METHODS NINV Numerical inversion 
#                     
#          folgende Kommentarzeilen/Block (bis zur ersten
#          Leerzeile) wird in TEXInfoformat ausgegeben
#          
# =ROUTINES
#          sucht C function zwischen =ROUTINES und =END 
#          beginnend mit 'unur_' und mit '()' abgeschlossen
#
#          Der Nachfolgende Kommentarblock/Zeilen werden mit
#          einer Leerzeile abgeschlossen.
#          Kommentare direkt vor jeder Funktion gelten als
#          interne Kommentare und werden nicht ausgegeben.
# =END     schliesst =ROUTINES ab (notwendig)
#
# (=>)     diese Zeichenfolge in Kommentarzeile wird in TEXINFO 
#          nicht ausgegegben (dient zur spaeteren Verwendung)
#
# =[A-Z,0-9] Unbekannte '=...' Zeichenfolgen wereden mit FEhler 
#            erkannt 
#
# =ERRORCODE derzeit ohne Funktion
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
@TYPES = ("UNUR_PAR", "UNUR_GEN", "struct", "void", "int", "double", "float", "long", "char", "short", "unsigned", "signed");



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
          if ($type eq UNUR_PAR && $MName ne $1){
             $MName = $1;
             if ($ERST == 0){
		print OUTFILE "\@end itemize\n\n";
             }
             else{
		 $ERST = 0;
                 $ENDE = 1;
	     }
	     print OUTFILE "\n\n\@node ", $MName, "\n";
	     print OUTFILE "\@subsection ", $MName, "\n\n";
             print OUTFILE "\@itemize \@minus{} \n";
          }
       }

       # suche Funktionszeile, funktion darf nicht
       # _unur enthalten (underscore!!)
       if ($_ =~ /_unur/ ){
	   $KOMMENT = 0;
       }
       elsif ( $ON == 1 && $_ =~/^\s*($type.*)\s*\((.*\))\s*;/){
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















