#!/usr/bin/perl

use strict;

my $VERBOSE = 1;

# ----------------------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------------------
#
# File: make_stringparser.pl
#
##############################################################################

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname < template > Cfile

    Generates a stringparser. A template file is read from STDIN.
    It contains the auxilliary routines and the skeleton for the 
    string interpreter. The routine reads all relevant information
    from the corresponding header files and inserts the different
    cases. The output is written on STDOUT.
      
EOM

    exit;
}

##############################################################################
# Methods not supported by string input
#
our @No_String_Methods = ("UNIF");

##############################################################################
# top source directory
#

my $top_srcdir = $ENV{'srcdir'} ? $ENV{'srcdir'} : '.';
$top_srcdir .= "/../..";

##############################################################################
# Get all header files in methods directory
#
my $methods_dir = "$top_srcdir/src/methods";
opendir (METHDIR, "$methods_dir") or die "can't open directory $methods_dir";
my @methods_h_files = grep {/[.]h/ } readdir METHDIR;
closedir METHDIR;

# Get all header files in distributions directory
#
my $distr_dir = "$top_srcdir/src/distributions";
my $distr_h_file = "$distr_dir/unuran_distributions.h";

##############################################################################
# Read template C file from STDIN and insert C code for string interpreter 
#
my $code;
while ( <STDIN> ){

    unless (/^s*=INPUT\s*(\w*)/){
	print $_;
    }

    else {
	my $type = $1;    # INPUT type 

	if ( $type eq "distrinfo" ){
	    distr_info();
	}
	elsif ( $type eq "methodinfo" ){
	    $code .= xxx();
	    method_info();
	}
	else{
	    print "Error: unknown qualifier after =INPUT: -$type-\n";
	}
    }

}

print $code;

##############################################################################
#
# The end
#
exit (0);

##############################################################################

##############################################################################
#                                                                            #
# Subroutines                                                                #
#                                                                            #
##############################################################################

##############################################################################
#
# Simple test for unmatched parenthesis.
# It returns the number of opening parenthesis `(' minus the number
# of closing parenthesis `)' in argument $_.
#
sub unmatched_parenthesis {
    my $open = 0;   # number of unmatched `('

    ++$open while /\(/g;
    --$open while /\)/g;

    return $open;
} # end of unmachted_parenthesis()


##############################################################################
#
# Read header files for all methods and ...
# 
# 
#
sub xxx {
    my $code;

    $code = str2method();
    $code .= str2set();

    return $code;
}

##############################################################################
#
# Make subroutine for getting parameter object for method
#
sub str2method {

    # Header for subroutine
    my $code = <<EOS;

struct unur_par *
_unur_str2method( const char *method, UNUR_DISTR *distr )
\t /*----------------------------------------------------------------------*/
\t /* get default parameters for method                                    */
\t /*                                                                      */
\t /* parameters:                                                          */
\t /*   method ... string that contains method name                        */
\t /*   distr  ... pointer to distribution object                          */
\t /*                                                                      */
\t /* return:                                                              */
\t /*   default parameters (pointer to structure)                          */
\t /*                                                                      */
\t /* error:                                                               */
\t /*   return NULL                                                        */
\t /*----------------------------------------------------------------------*/
{
\t struct unur_par *par;

#ifdef UNUR_ENABLE_LOGGING
\t /* write info into log file */
\t _unur_str_debug_string(1,"nxi",method);
#endif

EOS

    $code .= "\t ";

    # Read all header files 
    foreach my $hfile (@methods_h_files) {

	# Read content of header file
	open H, "< $methods_dir/$hfile" or  die ("can't open file: $methods_dir/$hfile");
	my $content = '';
	while (<H>) { $content .= $_; } 
	close H;

	# We skip over all header files that do not correspond
	# to a method.
	next unless $content =~ /[^\n]\s*=METHOD\s+(\w+)/;

	# ID for method
	my $method = "\L$1";

	# print code
	$code .= "if ( !strcmp( method, \"$method\") ) {\n";
	$code .= "\t\t par = unur_$method\_new(distr);\n";
	$code .= "\t }\n\t else ";
    }

    # Otherwise ...
    $code .= "{\n\t\t par = NULL;\n";
    $code .= "\t\t fprintf(stderr, \"ERROR: Unknown method: %s\\n\",method);\n\t }\n\n";

    # End 
    $code .= "\t return par;\n";
    $code .= "} /* end of _unur_str2method() */\n\n";

    # Return result
    return $code;

} # end of str2method() 

##############################################################################
#
# Make subroutine for setting parameters
#
sub str2set {

    # Header for subroutine
    my $code = <<EOS;

int
_unur_str2set( UNUR_PAR *par, const char *key, const char *value )
\t /*----------------------------------------------------------------------*/
\t /* set parameters for method                                            */
\t /*                                                                      */
\t /* parameters:                                                          */
\t /*   par   ... pointer to parameter for building generator object       */
\t /*   key   ... string that contains key for parameter                   */
\t /*   value ... string that contains list of arguments for parameter     */
\t /*                                                                      */
\t /* return:                                                              */
\t /*   1 ... on success                                                   */
\t /*   0 ... on error                                                     */
\t /*                                                                      */
\t /* error:                                                               */
\t /*   return 0                                                           */
\t /*----------------------------------------------------------------------*/
{
\t int result;

#ifdef UNUR_ENABLE_LOGGING
\t /* write info into log file */
\t _unur_str_debug_string(1,key,value);
#endif

EOS

    $code .= "\t ";



    # Read all header files 
    foreach my $hfile (@methods_h_files) {

	# Read content of header file
	open H, "< $methods_dir/$hfile" or  die ("can't open file: $methods_dir/$hfile");
	my $content = '';
	while (<H>) { $content .= $_; } 
	close H;

	# We skip over all header files that do not correspond
	# to a method.
	next unless $content =~ /[^\n]\s*=METHOD\s+(\w+)/;

	# ID for method
	print STDERR "$1: " if $VERBOSE;
	my $method = "\L$1";

	# Remove all comments and empty lines ...
	$content =~ s {/\*.*?\*/} []gsx;
	$content =~ s /\n\s*\n/\n/gsx;

	# Split into lines ...
	my @lines = split /\n/, $content;

	# Get all set calls
	foreach my $l (@lines) {
	    next unless $l =~ /^\s*(\w*\s+)unur_$method\_set_(\w+)\s*\((.+)([^\s])\s*$/; 

	    # name of set command
	    my $command = $2;

	    # list of arguments
	    my $args = $3;

	    # Check syntax of set command
	    if ( $4 ne ';' ) {
		# last character must be a ';'
		die "Unknown syntax (terminating ';' missing) in $hfile:\n$l\n";
	    }
	    if ( $1 !~ /^int\s+$/ ) {
		# type must be 'int'
		die "Unknown syntax (function type) in $hfile:\n$l\n";
	    }
	    if ( unmatched_parenthesis($l) ) {
		# parenthesis must match
		die "Unknown syntax (umatched parenthesis) in $hfile:\n$l\n";
	    }

	    # print name of parameter
	    print STDERR "$command ";

	    # process list of args
	    $args =~ s/\)\s*$//;             # remove closing parenthesis
	    my @args_list = split /\,/, $args;

	    # first argement must be of type UNUR_PAR
	    $a = shift @args_list;
	    unless ($a =~ /UNUR_PAR/) {
		die "Unknown syntax (first argument not of type UNUR_PAR) in $hfile:\n$l\n";
	    }

	    # number of arguments
	    my $n_args = $#args_list+1;

	    $code .= "if ( !strcmp(key, \"$command\") ){\n";
###	    $code .= "\t\t result = unur_$method\_set_$command(par);\n";
	    $code .= "\t\t result = 0;\n";

	    if (0) {
	    # only parameter object passed
#	    if ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*\)/ ){
##		print "if ( !strcmp(key, \"$2\") ){\n";
##		# set parameter
#                print "\t\t\tcheck = unur_$method\_set_$2(par);\n";
##                # check for error
##		print "\t\t\tif ( ! check ){\n";
##		print "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
##		print "\t\t\t}\n";
##		print "\t\t}\n\t\telse ";
##	    }
	    # parameter object and a single value passed
#	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*(\w+)\s+(\w+)\s*\)/ ){
##		print "if ( !strcmp(key, \"$2\") ){\n";
##		# set parameter
#		print "\t\t\tcheck = unur_$method\_set_$2(par, ($3) dblvalue);\n";
#                # check for error
##		print "\t\t\tif ( ! check ){\n";
##		print "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
##		print "\t\t\t}\n";
##		print "\t\t}\n\t\telse ";
##	    }
	    # parameter object and two doubles passed
#	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*double\s+\w+\s*,\s*double\s+\w+\s*\)/ ){
##		print "if ( !strcmp(key, \"$2\") ){\n";
##		# set parameter
#		print "\t\t\tcheck = unur_$method\_set_$2(par, list[0], list[1]);\n";
##                # check for error
##		print "\t\t\tif ( ! check ){\n";
##		print "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
##		print "\t\t\t}\n";
##		print "\t\t}\n\t\telse ";
##	    }
	    # parameter object, size_of_list and double list passed
#	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*int\s+\w+\s*,\s*double\s+\*\w+\s*\)/ ){
##		print "if ( !strcmp(key, \"$2\") ){\n";
##		# set parameter
#                print "\t\t\tif (no_of_elem != 0 ){\n";
#		print "\t\t\t\tcheck = unur_$method\_set_$2(par, no_of_elem, list);\n";
#		print "\t\t\t}\n";
#		print "\t\t\telse{\n";
#		print "\t\t\t\tcheck = unur_$method\_set_$2(par, (int) dblvalue, NULL);\n";
#		print "\t\t\t}\n";
                # check for error
#		print "\t\t\tif ( no_of_elem != (int) dblvalue ){\n";
#		print "\t\t\t\tfprintf(stderr, \"WARNING: Check size of list when setting $2\\n\");\n";
#		print "\t\t\t}\n";
#		print "\t\t\tif ( ! check ){\n";
#		print "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
#		print "\t\t\t}\n";
#		print "\t\t}\n\t\telse ";
#	    }
	    }


	    $code .= "\t }\n\t else ";

	}

	# end of method
	print STDERR "\n";
    }

    # Otherwise ...
    $code .= "{\n\t\t result = 0;\n";
    $code .= "\t\t fprintf(stderr, \"ERROR: Unknown ....: %s\\n\",key);\n\t }\n\n";

    # End 
    $code .= "\t return 0;\n";
    $code .= "} /* end of _unur_str2set() */\n\n";

    # Return result
    return $code;

} # end of str2set() 




# #######################################################################################
# 
# UNURAN's input as string interface: 
# generates code for generating parameter object
#
# #######################################################################################
sub method_info{


    print "\t/* ************************ */\n";
    print "\t/* determine choosen method */\n";
    print "\t/* ************************ */\n";

    print "\tif ( !strcmp( key, \"method\") ){\n";
    print "\t\tpar = _unur_str2method(value,distr);\n";
    print "\t}\n";

##    print "\tif ( !strcmp( key, \"method\") ){\n";
##    print "\t\tif ( no_of_elem != 0 ){\n";
##    print "\t\t\tfprintf(stderr, \"SYNTAX ERROR: No list expected when setting method.\\n\");\n";
##    print "\t\t}\n\t\t";



    print "\n\t/* ****************************************** */\n";
    print "\t/*                                            */\n";
    print "\t/* set parameters depending on choosen method */\n";
    print "\t/*                                            */\n";
    print "\t/* ****************************************** */\n\n\t";


    foreach my $hfile (@methods_h_files){

	# search routines setting parameters -- method dependent
	my $method = "UNKNOWN"; # initialize method-name to empty string
        my $check  = 0;         # Variable for checking syntax

	open INFILE,  "< $methods_dir/$hfile" or  die ("can't open file: $methods_dir/$hfile");
	while ( $_ =  <INFILE> ){

	    if ( $_ =~ /^\s*=METHOD\s+(\w+)/ ){ # h-file with method description found
		$method = "\L$1";
		my $METHOD = "\U$method";

		# key is "method"
		print "if ( par->method == UNUR_METH_$METHOD && strcmp(key, \"method\") ){\n\t\t";
	    }

	    # only parameter object passed
	    if ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*\)/ ){
		$check = unmatched_parenthesis($_);
		print "if ( !strcmp(key, \"$2\") ){\n";
		# set parameter
                print "\t\t\tcheck = unur_$method\_set_$2(par);\n";
                # check for error
		print "\t\t\tif ( ! check ){\n";
		print "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
		print "\t\t\t}\n";
		print "\t\t}\n\t\telse ";
	    }
	    # parameter object and a single value passed
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*(\w+)\s+(\w+)\s*\)/ ){
		$check = unmatched_parenthesis($_);
		print "if ( !strcmp(key, \"$2\") ){\n";
		# set parameter
		print "\t\t\tcheck = unur_$method\_set_$2(par, ($3) dblvalue);\n";
                # check for error
		print "\t\t\tif ( ! check ){\n";
		print "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
		print "\t\t\t}\n";
		print "\t\t}\n\t\telse ";
	    }
	    # parameter object and two doubles passed
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*double\s+\w+\s*,\s*double\s+\w+\s*\)/ ){
		$check = unmatched_parenthesis($_);
		print "if ( !strcmp(key, \"$2\") ){\n";
		# set parameter
		print "\t\t\tcheck = unur_$method\_set_$2(par, list[0], list[1]);\n";
                # check for error
		print "\t\t\tif ( ! check ){\n";
		print "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
		print "\t\t\t}\n";
		print "\t\t}\n\t\telse ";
	    }
	    # parameter object, size_of_list and double list passed
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*int\s+\w+\s*,\s*double\s+\*\w+\s*\)/ ){
		$check = unmatched_parenthesis($_);
		print "if ( !strcmp(key, \"$2\") ){\n";
		# set parameter
                print "\t\t\tif (no_of_elem != 0 ){\n";
		print "\t\t\t\tcheck = unur_$method\_set_$2(par, no_of_elem, list);\n";
		print "\t\t\t}\n";
		print "\t\t\telse{\n";
		print "\t\t\t\tcheck = unur_$method\_set_$2(par, (int) dblvalue, NULL);\n";
		print "\t\t\t}\n";
                # check for error
		print "\t\t\tif ( no_of_elem != (int) dblvalue ){\n";
		print "\t\t\t\tfprintf(stderr, \"WARNING: Check size of list when setting $2\\n\");\n";
		print "\t\t\t}\n";
		print "\t\t\tif ( ! check ){\n";
		print "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
		print "\t\t\t}\n";
		print "\t\t}\n\t\telse ";
	    }
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters(.*)\)/ ){
		$check = unmatched_parenthesis($_);
		print "if ( !strcmp(key, \"$2\") ){\n";
		print "\t\t\tfprintf(stderr, \"Setting $2 for method $method not supported via this interface.\\n\");\n";
		print "\t\t}\n\t\telse ";
	    }
	    elsif ( $_ =~ /^\s*\w+\s*unur_($method)_set_(.*?)\s*\(\s*\w+.*\)/  ){
		$check = 1;
	    }

            if ( $check != 0 ){
                die "Unknown Syntax of set command in $hfile line $.: $_\n";
	    }

	} # end of while-loop
	close INFILE;
	
	# closing the method-block 
	if ( $method !~ /UNKNOWN/){
	    print "{\n\t\t\tfprintf(stderr, \"Unknown option for method $method: %s\\n\", key);\n\t\t}\n";
	    print "\t}\n\telse ";
	}

    } # end of foreach loop

    # case of unknown method
    print "if ( !strcmp(key, \"method\") && par->method == UNUR_METH_UNKNOWN){\n";
    print "\t\tfprintf(stderr, \"ERROR: No or unknown method defined.\\n\");\n\t}\n";

    close INFILE;
}







# ###########################################################################
# 
# input as string interface:
# generates code for generating distribution object
#
# ###########################################################################
sub distr_info{

    print "\t/* ------------------------------------------- */\n";
    print "\t/*                                             */\n";
    print "\t/* key = \"distribution\"                        */\n";
    print "\t/*                                             */\n";
    print "\t/* ------------------------------------------- */\n";
    print "\tif ( !strcmp(key, \"distr\") ){ \n\t\t";
    open INFILE, "< $distr_h_file" or  die ("can't open file: $distr_h_file");

    while ( <INFILE> ){

	# search for standard distribution
	if ( $_ =~ /^\s*=DISTR\s+(\w+)/ ){
	    print "if ( !strcmp(value, \"$1\") ){\n";
	    print "\t\t\tdistr = unur_distr_$1(list, no_of_elem);\n\t\t}\n\t\telse ";
	}
    }

    # Error -- in case of unknown distribution
    print "{\n\t\t\tfprintf(stderr,\"ERROR: Unknown distribution: %s\\n\",value);\n";
    print "\t\t\tbreak;\n\t\t}\n\t}\n";

    print "\t/* ------------------------------------------- */\n";
    print "\t/*                                             */\n";
    print "\t/* key = \"domain\"                              */\n";
    print "\t/*                                             */\n";
    print "\t/* ------------------------------------------- */\n";
    print "\telse if ( !strcmp( key , \"domain\") ){\n";
    print "\t  /* list must contain exactly two entries */\n";
    print "\t  if ( no_of_elem != 2 ){\n";
    print "\t    fprintf(stderr, \"SYNTAX ERROR: To set the domain use a list with exactly 2 entries -- Standard domain is used instead.\\n\");\n";
    print "\t  }\n";
    print "\t  else if ( unur_distr_is_cont(distr) ){\n";
    print "\t    unur_distr_cont_set_domain( distr, list[0], list[1]);\n";
    print "\t  }\n";
    print "\t  else if ( unur_distr_is_discr(distr) ){\n";
    print "\t    unur_distr_discr_set_domain( distr, list[0], list[1]);\n";
    print "\t  }\n";
    print "\t  else{\n";
    print "\t    fprintf(stderr, \"ERROR: Can't set domain for this type of distribution\\n\");\n";
    print "\t    break;\n";
    print "\t  }\n";
    print "\t}\n";
    print "\telse {\n";
    print "\t  fprintf(stderr, \"SYNTAX ERROR?: Unknown key will be ignored: %s\\n\", key);\n";
    print "\t}\n";


    close INFILE;

} # end of distr_info()

