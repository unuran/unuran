#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl spript extracts the description of what is required 
# by a method implemented in UNURAN to texinfo format.
# 
# usage:  ./required.perl ../src/methods/*.h
#      (required.perl is in the directory /unuran/doc/
# Input:   header-files
# Output:  $OUTFILE contains docomentation in texinfo format
# 
# E. JANKA
# $Id$
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# This script searches for the following key word:
#
# =REQUIRED
#    usage: =REQUIRED Method
#          (e.g =REQUIRED NINV)
#
#          IMMEDIATLY following comments (up to the first
#          non-comment-blank line) will be used as description
#          of the method and written to the output file. 
#          Comments bevore the function declaration or after
#          a blank line will be handled as internal information
#          
# =[A-Z,0-9]* Unknown `=...' key words will produce warnings  
#
# key words used by other scripts (no warnings):
# =DEF
# =METHOD
# =ROUTINES
# =END
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
@COMMAND =("=METHOD", "=ROUTINES", "=REQUIRED", "=DEF", "=END");

# output file:
open(OUTFILE, ">qstart_required.texi");


while($_ = <>)
{ 
    chomp; 

    # comments start with `/*' (in case of $CENABLE=1)
    if ($CENABLE == 1 && $_ =~/^\s*\/\*/){    
	$KOMMENT = 1;
    }

    # comments end with a blank line
    if ($BLOCK == 0 && $_ =~ /^\s*$/){
        $KOMMENT = 0;
        $CENABLE = 0;
    }
  
    # check for wrong key worde (e.g. typing error)
    if ($_ =~/^(\s*\/\*\s*|\s*)(=\w+?)\W/){
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
   
    #  requirements of a method are described (=REQUIRED)
    if ($_ =~/^(\s*\/\*\s*|\s*)=REQUIRED\s*(\w*)\s*/){
        print OUTFILE  "\@*", $2, ":\n";
        # comments possible 
        $CENABLE = 1;
        # already within a comment
	$KOMMENT = 1;
 	$BLOCK=1;
    }

 # output of comments
 if ($KOMMENT == 1){

     # comment -- single line  (begins with `/*', ends with `*/')
     if ( $_ =~ /^(\s*\/\*\s*)(.*?)(\s*\*\/)/ ){
	$PRINT = join  '',$2,"\n";
     }
     # comment -- more lines (but not the first line)
     elsif( $BLOCK == 1 && $_ =~ /^\s*(\*\/)?(.*?)\s*(\*\/)?$/ ){
	 $PRINT = join '', $2, " ";
     }
     # comment -- more lines (first line)
     elsif($_ =~ /^(\s*\/\*\s*)(.*)/ ){
       $PRINT = join '', $2, " ";
       $BLOCK = 1;
     }
     
     # print comments: 
     #  `(=>)' will be cut and lines with key words are ingored 
     if( $PRINT !~ /(\s*\/\*\s*|\s*)=\w+/ ) { 
        print OUTFILE join ' ', split /\(=>\)/, $PRINT;
     }

     # end of a comment (`*/') -> linebreak, $BLOCK=0 
     if ($_ =~/\*\//){
         print OUTFILE "\@*\n";
         $BLOCK = 0;
     }

 }  # --- if (KOMMENT == 1) ende ---

}




