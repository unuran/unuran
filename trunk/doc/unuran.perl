#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl spript generating the referernce of the methods/functions
# of UNURAN in the texinfo format.   
# 
# usage:  ./unuran.perl ../src/methods/*.h
#      (unuran.perl is in the directory /unuran/doc/
# Input:   header-files
# Output:  $OUTFILE contains docomentation in texinfo format
# 
# E. JANKA und G. TIRLER
# $Id$
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# This script searches for keywords beginning with `=...'
# The keywords may only be preceded by white space or `/*'
#
# =METHOD  name [longform]
#          name and description of a method within a header file
#
#          name       ... Name of method (ONE word!)
#          [longform] ... optional, long form of the name
#          Beispiel:
#             =METHODS NINV Numerical INVersion 
#                     
#          IMMEDIATLY following comments (up to the first
#          non-comment-blank line) will be used as description
#          of the method and written to the output file. 
#          
# =ROUTINES
#          C functions between =ROUTINES and =END beginning
#          with `unur_' and endingd with `(...)' will be 
#          documented in the output file.
#
#          IMMEDIATLY following comments (up to the first
#          non-comment-blank line) will be used as description
#          of the method and written to the output file. 
#          Comments bevore the function declaration or after
#          a blank line will be handled as internal information
#          
# =END     Ends the block beginning with  =ROUTINES.
#          Function declarations not being within this block
#          will be ingored. =END must not be omitted.
#
# =OPTIONAL, =REQUIRED
#          These key words are allowed BETWEEN =ROUTINES and
#          and =END and might be used for specification of
#          the subroutines.
#          Both will be ingored by this script;
#          There won't be a warning if used.
#
# (=>)     This string within a comment won't appear in
#          the texinfo output. 
#
# =[A-Z,0-9]* Unknown `=...' key words will produce warnings  
#
# =ERRORCODE not in use
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

# known data types 
@TYPES = ("UNUR_PAR", "UNUR_GEN", "struct", "void", "int", "double", "float", "long", "char", "short", "unsigned", "signed");

# permitted key words (no warning)
@COMMAND =("=METHOD", "=ROUTINES", "=REQUIRED", "=OPTIONAL", "=END");

# output file:
open(OUTFILE, ">qstart_function_reference.texi");


while($_ = <>)
{ 
    chomp; 

    # Region with description of the routines begins with =ROUTINES
    if ($_ =~/(^\s*\/\*\s*|^\s*)=ROUTINES/){
	$ROUTINES = 1;
        # begin of itemize environment -- a funtion - an item
	print OUTFILE "\@itemize \@minus\{\}\n";
    }

    # Region with description of the routines ends with =END
    if ($_ =~/(^\s*\/\*\s*|^\s*)=END/){
       # =END only possible if $ROUTINES = 1 (error)
       if ($ROUTINES = 0){
           print "\n\nERROR: =END before =ROUTINES\n\n";
       }
       $ROUTINES = 0;
       # last itemize environment must be finished
       print OUTFILE "\@end itemize\n";
    }

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
    if ($_ =~/^(\s*\/\*\s*|\s*)(=.*?)\W/){
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

   # Screening for function declarations
   foreach $type (@TYPES){

      # Screening for functions with allowed names
      # (beginning with `unur_'  (`_unur_...' is ignored (-> intern) )
      if ($_ =~ /_unur/ ){
          $CENABLE = 0;
          $KOMMENT = 0;
      }
      elsif ( $ROUTINES == 1 && $_ =~/^\s*($type.*)\s*\((.*\))\s*;/){
           $CENABLE = 1;
           $DECL = $1;   # string bevore the braces
           $FUNC = $2;   # string between the braces 
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
