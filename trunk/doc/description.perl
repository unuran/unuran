#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl spript extracts the description of a method
# fo the texinfo format.   
# 
# usage:  ./description.perl ../src/methods/*.h
#      (description.perl is in the directory /unuran/doc/
# Input:   header-files
# Output:  $OUTFILE contains docomentation in texinfo format
# 
# E. JANKA
# $Id$
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# This script searches for keywords =METHOD
# The keyword may only be preceded by white space or `/*'
#
# =METHOD
# usage:   =METHOD  name [longform]
#          name and description of a method within a header file
#
#          name       ... Name of method (ONE word!)
#          [longform] ... optional, long form of the name
#          Beispiel:
#             =METHODS NINV Numerical INVersion 
#                     
#          IMMEDIATLY following lines of comment (up to the first
#          non-comment-blank line) will be used as description
#          of the method and written to the output file. 
#
# -----------------------------------------------------------------------
#
# keywords used by other scripts:
#          
# =ROUTINES
#          C functions between =ROUTINES and =END beginning
#          with `unur_' and endingd with `(...)' will be 
#          documented in the output file.
#
#          IMMEDIATLY following lines of comment (up to the first
#          non-comment-blank line) will be used as description
#          of the method and written to the output file. 
#          Comments bevore the function declaration or after
#          a blank line will be handled as internal information
#          
# =END     Ends the block beginning with  =ROUTINES.
#          Function declarations not being within this block
#          will be ingored. =END must not be omitted.
#
# (=>)     This string within a comment won't appear in
#          the texinfo output. 
#
# =[A-Z,0-9]* Unknown `=...' key words will produce warnings  
#
# =ERRORCODE this key word is not in use
#
# =DEF
#    usage: 
#      XXX.kkk = 44.444  /* =DEF METHOD value some words of description
#                            maybe some lines                   */
#      this will produce following output:
#      "wert = 44.444   some words of description maybe some lines"
#      "METHOD" is the method using the parameter
#    The begin of the comment ("\*") has to be in the same line as
#    the 
#    The same holds for `METHOD' and `value'. The description may
#    begin in the following line. There also are consecutive blocks
#    of comment allowed (no blank line in between).
#  
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
# =INT
#   analog to =ROUTINES, also ends with =END
#   searches for internal routines
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

# region with declaration of routines found
$ROUTINES = 0;

# permitted key words (no warning)
@COMMAND =("=METHOD", "=ROUTINES", "=REQUIRED", "=DEF", "=END", "=INT");

# output file:
open(OUTFILE, ">qstart_method_description.texi");


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
   
    # name and description of a method (=METHOD)
    if ($_ =~/^(\s*\/\*\s*|\s*)=METHOD\s*(\w+)\s*(.*)/){
        # comments possible 
        $CENABLE = 1;
        # already within a comment
	$KOMMENT = 1;
        # formatting the output -- Header, node for new method
        print OUTFILE "\n\n\@node ", $2, "\n";
        print OUTFILE "\@subsection ", $2, " ", $3, "\n\n";
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
