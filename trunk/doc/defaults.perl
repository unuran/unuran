#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl spript extracting the default settings of the parameters
# of UNURAN methods
# 
# usage:  ./defaults.perl ../src/methods/*.c
#      (defaults.perl is in the directory /unuran/doc/
# Input:   c-files
# Output:  $OUTFILE contains docomentation in texinfo format
# 
# E. JANKA
# $Id$
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# This script searches for the following key word:
#
# =DEF
#
# usage: 
#    XXX.kkk = 44.444  /* =DEF METHOD wert some words of description
#                            maybe some lines                   */
#    this will produce following output:
#    "wert = 44.444   some words of description maybe some lines"
#    "METHOD" is the method using the parameter
#  
#
# =[A-Z,0-9]* Unknown `=...' key words will produce warnings  
#
# key words used by other scripts (no warnings):
# =METHOD
# =ROUTINES
# =END
# =REQUIRED
#
# (=>)     This string within a comment won't appear in
#          the texinfo output. 
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# comments possible 
$CENABLE = 0;
# $_ ist Kommentarzeile (bei $CENABLE = 1) 
$KOMMENT = 0;
# comment with more than just one line
$BLOCK   = 0;
# contains line to be printed
$PRINT = "";


# permitted key words (no warning)
@COMMAND =("=METHOD", "=ROUTINES", "=DEF", "=OPTIONAL", "=END");

# output file:
open(OUTFILE, ">qstart_paramdefaults.texi");


while($_ = <>)
{ 
    chomp; 

    # check for wrong key worde (e.g. typing error)
    if ($_ =~/^(\s*\/\*\s*|\s*)(=\w+?)\W(.*)/){
	$ERROR = 1;
	foreach $command (@COMMAND){
	    if ($2 eq $command){
		$ERROR = 0;
	    }
	}
	if ($ERROR == 1){
	    print "WARNING: unknown command: ", $2, "\n";
	}
    }

    # preset value of $PRINT
    $PRINT = "";

    # a new description block begins
    if($BLOCK == -1){
       if ($_ =~ /^\s*\/\*/ ){
	   $_ =~ s/\s*\/\*//;
	   $BLOCK = 1;
       }
       else {
	   $BLOCK = 0;
       }
   }

    # lines continuing a block of description
    if ($BLOCK == 1){
	$_ =~ /^\s*(.*?)\s*(\*\/)*\s*$/;
        $PRINT = join '',$1, "\n";
	if ($_ =~ /\*\//){
	    $BLOCK = -1; # description block ends
	}
    }

    # search lines with =DEF
    if (($value, $method, $wert, $description) = 
      $_ =~ /.*=(.*?);\s*\/\*\s*=DEF\s+(\w+)\s+(\w+)\s*(.*)/) {
        # are there more lines of description?
	$BLOCK = 0;
	if ($description !~ s/(.*?)\*\/\s*$/$1/){
	   $BLOCK = 1;
        }

    $PRINT = join  ' ', "\@*". $method, ": \@", $wert, "=", $value,
                        $description, "\n";
    }
 


    # output
    if ($PRINT !~ /^\s*$/){
	print OUTFILE join '', split /\(=>\)/, $PRINT;
    }

}











