#!/usr/bin/perl

use strict;

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
usage: $progname srcdir < template > Cfile

    Generates a stringparser. A template file is read from STDIN.
    It contains the auxilliary routines and the skeleton for the 
    string interpreter. The routine reads all relevant information
    from the corresponding header files and inserts the different
    cases. The top source directory for searching the header files
    must be given as argument.
    The output is written on STDOUT.
      
EOM

    exit;
}

##############################################################################
# Methods not supported by string input
#
our @No_String_Methods = ("UNIF");

##############################################################################
# Read top source directory from argument list ...
#
my $top_srcdir = shift;
(usage and die) unless $top_srcdir;

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
	    method_info();
	}
	elsif ( $type eq "known_methods" ){
	    known_methods();
	}
	else{
	    print "Error: unknown qualifier after =INPUT: -$type-\n";
	}
    }
}

##############################################################################

##############################################################################
#                                                                            #
# Subrotuines                                                                #
#                                                                            #
##############################################################################

##############################################################################
#
# Report unknown in UNURAN set command
#
sub fatal_set {
    my $file = $_[0];
    my $line = $_[1];
    my $command = $_[2];

    die "Unknown Syntax of set command in $file:$line:\n$command\n";
} # end of fatal_set()

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
















# ##################################################################
#
# macro assignments: known methods -- integers
#
# ##################################################################
sub known_methods{

    #initialize counter
    my $counter = 0;

    foreach my $hfile (@methods_h_files){

	open INFILE, "< $methods_dir/$hfile" or  die ("can't open file: $methods_dir/$hfile");

	while ( <INFILE> ){
	    if ( $_ =~ /^\s*=METHOD\s+(\w+)/ ){
		print "#define \UUNUR_METH_$1\t($counter)\n";
		$counter++;
		break();
	    }
	}
    }
} # end of known_methods()


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
    print "\t\tif ( no_of_elem != 0 ){\n";
    print "\t\t\tfprintf(stderr, \"SYNTAX ERROR: No list expected when setting method.\\n\");\n";
    print "\t\t}\n\t\t";

    foreach my $hfile (@methods_h_files){
        # reset variable
	my $No_String_Method  = 0;
	open INFILE, "< $methods_dir/$hfile" or  die ("can't open file: $methods_dir/$hfile");


        # generate prameter object for specified method -- if possible
	while ( <INFILE> ){
	    if ( $_ =~ /^\s*=METHOD\s+(\w+)/ ){

                # check if method supports string context
		foreach my $method (@No_String_Methods){
		    if ( "\U$1" eq $method){
			$No_String_Method  = 1;
		    }
		}

		print "if ( !strcmp( value, \"\L$1\") ){\n";
		if ( $No_String_Method == 1 ){
		    # print "\t\t\tpar = unur_\L$1_new();\n";
		    print "\t\t\tfprintf(stderr, \"Method $1 not intended for usage within this string context.\\n\");\n";
		}
		else{
		    print "\t\t\tpar = unur_\L$1_new(distr);\n";
		    while ( $_ !~ /^\s*=UP\s+Methods_for_\w+/ ){
			$_ = <INFILE>;
		    }
		    $_ =~ /^\s*=UP\s+Methods_for_(\w+)/;
		}
		print "\t\t}\n\t\telse ";
	    }
	} # end of while-loop
	close INFILE;
    } # end of foreach-loop

    print "{\n\t\t\tpar->method = UNUR_METH_UNKNOWN;\n";
    print "\t\t\tfprintf(stderr, \"ERROR: Unknown method: %s\\n\",value);\n";
    print "\t\t\tbreak;\n\t\t}\n\t}\n";


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
		fatal_set($hfile,$.,$_);
#                die "Unknown Syntax of set command in $hfile line $.: $_\n";
#		exit (-1);
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




