#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl spript zur automatischen Erstellung der   
# (Methoden-)Funktions- Referenz von UNURAN im TEX-Info Format
# 
# 
# Aufruf:  ./unuran.perl ../src/methods/*.h
#      (unuran.perl liegt im Verzeichnis /unuran/doc/
# Input:  
# Output:
# 
# E.JANKA und G.TIRLER  August 2000
# $Id$
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Das Perl-skrip durchsucht nach folgenden Schluesselworten      
# und und erzeugt Output wiefolgt:
# Folgende Schluesselwoerter beginnen mit '=...'
# Schluesselwort darf nur von Leerzeichen oder '/*'
# eingeleitet werden
#
#
# =METHOD  name [Langtext]
#          allgem. Beschreibung der Methode in Header file
#
#          name       ... Name der Methode, Nur ein Wort!
#          [Langtext] ... optional, Lange Beschreibung
#          Beispiel:
#             =METHODS NINV Numerical inversion 
#                     
#          folgende Kommentarzeilen/Block (bis zur ersten
#          Leerzeile) wird in TEXInfoformat ausgegeben
#          
# =ROUTINES
#          sucht C function zwischen =ROUTINES und =END 
#          beginnend mit 'unur_' und endend mit '()'
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
# =[A-Z,0-9] Unbekannte '=...' Zeichenfolgen wereden mit Fehler 
#            erkannt 
#
# =ERRORCODE derzeit ohne Funktion
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Kommentare sind nun erlaubt 
$CENABLE = 0;
# $_ ist Kommentarzeile (bei $CENABLE = 1) 
$KOMMENT = 0;
# Kennzeichnet Kommentarbloecke (zur Unterscheidung von Kommentarzeilen)
$BLOCK   = 0;
# Das soll ausgedruckt werden
$PRINT = "";


# Deklarations/Beschreibungsbereich der Routinen gefunden
$ROUTINES = 0;

# erlaubte Datentypen 
@TYPES = ("UNUR_PAR", "UNUR_GEN", "struct", "void", "int", "double", "float", "long", "char", "short", "unsigned", "signed");

# erlaubte Befehle
@COMMAND =("=METHOD", "=ROUTINES", "=END");

# output file:
open(OUTFILE, ">qstart_function_reference.texi");


while($_ = <>)
{ 
    chomp; 
   
    # Beschreibung der Routinen beginnt mit =ROUTINES
    if ($_ =~/(^\s*\/\*\s*|^\s*)=ROUTINES/){
     $ROUTINES = 1;  
    }

    # Beschreibung der Routinen endet mit =END
    if ($_ =~/(^\s*\/\*\s*|^\s*)=END/){
       # =END kann nur bei $ROUTINES = 1 erfolgen
       if ($ROUTINES = 0){
           print "\n\nERROR: =END before =ROUTINES\n\n";
       }
       $ROUTINES = 0;
       # letzte itemize-Umgebung muss beendet werden
       print OUTFILE "\@end itemize\n";
    }


    # Kommentare(Bloecke) beginnen mit /* (bei $CENABLE=1)
    if ($CENABLE == 1 && $_ =~/^\s*\/\*/){    
	$KOMMENT = 1;
    }

    # Kommentare enden mit Leerzeile 
    if ($BLOCK == 0 && $_ =~ /^\s*$/){
        $KOMMENT = 0;
        $CENABLE = 0;
    }
  
    # ueberpruefung ob falsches command (z.B.tippfehler?)
    if ($_ =~/^(\s*\/\*\s*|\s*)(=.*)\s/){
	foreach $command (@COMMAND){
	    if ($2 != $command){
		print "ERROR: unknown command: ", $2;
	    }
	}
    }
   
    # Suche Beschreibung der Methode (=METHOD)
    if ($_ =~/^(\s*\/\*\s*|\s*)=METHOD\s*(\w+)\s*(.*)/){
        # folgende Kommentare werden gedruckt 
        $CENABLE = 1;
        # sind bereits mitten im Kommentar
	$KOMMENT = 1;
        print OUTFILE "\n\n\@node ", $2, "\n";
        print OUTFILE "\@subsection ", $2, " ", $3, "\n\n";
        print OUTFILE "\@itemize \@minus\{\}\n";
	$BLOCK=1;
    }

   # Suche Funktion und Definitionszeilen
   foreach $type (@TYPES){

      # Suche Zeile mit erlauben Funktionsdeklarationen
      # (beginnen mit unur_  (_unur_ nicht erlaubt (-> intern) )
      if ($_ =~ /_unur/ ){
          $CENABLE = 0;
          $KOMMENT = 0; # eigentlich enden Kommentare mit Leerzeile
      }
      elsif ( $ROUTINES == 1 && $_ =~/^\s*($type.*)\s*\((.*\))\s*;/){
           $CENABLE = 1;
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


  # Ausgabe der Kommentare    
  if ($KOMMENT == 1){

     # einzelne Kommentarzeile  (beginnnen mit "/*" und enden mit "*/")
     if ( $_ =~ /^(\s*\/\*\s*)(.*)?(\s*\*\/)/ ){
        $BLOCK = 0;
	$PRINT = join  '',$2,"\n";
      }
     # Zeile beginnt mit Kommentarblock-Ende
     elsif ($BLOCK == 1 && $_ =~ /\s*\*\// ){
	 $BLOCK = 0;
         $PRINT = "";
     }
     # Kommentarzeile in einem Kommentarblock (nicht erste Zeile)
     elsif( $BLOCK == 1 && $_ =~ /^(\s*)(.*)\s*(\*\/)?/ ){
	 $PRINT = join '', $2, " ";
         # Kommentarblocke enten mit "*/"
	 if ($_ =~ /\*\//){
	     $BLOCK = 0;
	 }
     }
     # Beginn eines Kommentarblocks (Zeile beginnt mit "/*"
     # und darf "*/" nicht enthalten)
     elsif($_ =~ /^(\s*\/\*\s*)(.*)/ ){
       $PRINT = join '', $2, " ";
       $BLOCK = 1;
     }
     
     #  Ausdruck der Kommentare, wobei
     #  "(=>)" und zeilen die mit einem =BEFEHL beginnen
     #  nicht gedruckt werden sollen
     if( $PRINT !~ /(\s*\/\*\s*|\s*)=\w+/ ) { 
        print OUTFILE join ' ', split /\(=>\)/, $PRINT;
     }

     # Ende einer Kommentarzeile/blocks -> Zeilenumbruch 
     if ($_ =~/\*\//){
         print OUTFILE "\@*\n";
     }

 }

}
