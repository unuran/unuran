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
# Global variables

# Store unsupported set calls
my $unsupported;

##############################################################################
# Read template C file from STDIN and insert C code for string interpreter 
#
while ( <STDIN> ){

    unless (/^s*=INPUT\s*(\w*)/){
	print $_;
    }

    else {
	my $type = $1;    # INPUT type 

	if ( $type eq "distrinfo" ){
	    distr_info();
	}
	elsif ( $type eq "list_of_methods" ){
	    print make_list_of_methods();
	}
	elsif ( $type eq "list_of_par_sets" ){
	    print make_list_of_par_sets();
	}
	else{
	    print "Error: unknown qualifier after =INPUT: -$type-\n";
	}
    }
}

if ($unsupported) {
    print STDERR "\nUnsupported set commands:\n";
    print STDERR "$unsupported\n";
}

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
# Make subroutine for getting parameter object for method
#
sub make_list_of_methods {

    my $code .= "\t ";

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
    $code .= "\t\t fprintf(stderr, \"ERROR: Unknown method: %s\\n\",method);\n\t }\n";

    # Return result
    return $code;

} # end of make_list_of_methods()

##############################################################################
#
# Make subroutine for setting parameters in parameter objects
#
sub make_list_of_par_sets {

    my $code_unsupported;
    my $code = "\t ";

    # print info on screen
    print STDERR "Set commands for Methods:\n" if $VERBOSE;

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
	print STDERR "  $1: " if $VERBOSE;
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
	    print STDERR "$command " if $VERBOSE;

	    # process list of args
	    $args =~ s/\)\s*$//;             # remove closing parenthesis
	    my @args_list = split /\,/, $args;

	    # first argement must be of type UNUR_PAR
	    my $a = shift @args_list;
	    unless ($a =~ /UNUR_PAR/) {
		die "Unknown syntax (first argument not of type UNUR_PAR) in $hfile:\n$l\n";
	    }

	    # number of arguments
	    my $n_args = $#args_list+1;

	    # get type of arguments
	    my $type_args = "";
	    foreach my $a (@args_list) {
		my $t;
		# type of argument
		if    ($a =~ /double/)   { $t = 'd'; }
		elsif ($a =~ /int/)      { $t = 'i'; }
		elsif ($a =~ /unsigned/) { $t = 'u'; }
		elsif ($a =~ /char/)     { $t = 'c'; }
		else                     { $t = '?'; }

		# arrays are indicated by capital letters
		if ($a =~ /\*/ or $a =~ /\[/) {
		    # this is interprated as list
		    $t = "\U$t";
		}
		$type_args .= $t;
	    }

	    # make set calls
	    my $set;
	    my $args_comment;

	    # beginning of case
	    $set = "if ( !strcmp(key, \"$command\") ) {\n";
	    $args_comment = "\t /* n = $n_args; type = $type_args: $args*/\n";
	    $set .= "\t$args_comment";

	    # use keyword "void" when no argument is required
	    $type_args = "void" if $type_args eq "";

	    # we support the following cases:
	    #   void   ... no argument required
	    #   "i"    ... one argument of type int required
	    #   "u"    ... one argument of type unsigned required
	    #   "d"    ... one argument of type double required 
	    #   "dd"   ... two arguments of type double required 
	    #   "iD"   ... one argument of type int and a list of doubles required
	    #              (the first argument is considered as size of the double array
	    #   "Di"   ... a list of doubles and one argument of type int required
	    if ($type_args =~ /^(void|i|u|d|dd|iD|Di)$/) {
		my $type = $1;
		$set .= "\t\t result = _unur_str_par_set_$type(par,type_args,args,unur_$method\_set_$command);\n";
		$code .= $set;
		$code .= "\t }\n\t else ";
	    }

	    else {
		# cannot handle this set command
		$code_unsupported .= $args_comment;
		$unsupported .= "  unur_$method\_set_$command()\n"
	    }
	}

	# end of method
	print STDERR "\n" if $VERBOSE;
    }

    # Otherwise ...
    $code .= "{\n\t\t result = 0;\n";
    $code .= "\t\t fprintf(stderr, \"ERROR: Unknown ....: %s\\n\",key);\n\t }\n";

    # print info on screen
    print STDERR "\n" if $VERBOSE;

    # add comment on unsupported code into C file
    $code .= "\n\t /* Unsupported set commands: */\n $code_unsupported\n";

    # Return result
    return $code;

} # end of make_list_of_par_sets() 







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

##############################################################################
##############################################################################
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

