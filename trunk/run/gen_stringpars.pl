#!/usr/bin/perl

use strict;

# ----------------------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------------------
#
# File: gen_stringpars.pl
#
# Description:
# Generates a stringparser from the file stringpars.c.base with
# information found within the source of UNURAN and writes it to
# the file stringpars.c
#

# Methods not supporting by string input
#
our @No_String_Methods = ("UNIF");

# where relevant information is stored:
#
our $Methinfopath = "../src/methods/";
our $Distinfofile = "../src/distributions/unuran_distributions.h";

# determine all h-files in $Methinfopath
#
opendir (METHDIR, $Methinfopath) or die "can't open directory $!";
our @all_h_files = grep {/[.]h/ } readdir METHDIR;
closedir METHDIR;


# from the Basefile and the info found in the source of UNURAN
# the file Outfile is generated
#
open Basefile, "< ./stringpars.c.base" or die "can't open file: $!";
open Outfile,  "> ./stringpars.c"      or die "can't open file: $!";


# go through the basefile and insert C code
while ( <Basefile> ){

    if ( $_ !~/^s*=INPUT/){
	print Outfile $_;
    }
    else{
	$_ =~ /^\s*=INPUT\s+(\w+)\s*/;
	if ( $1 eq "distrinfo" ){
	    distr_info();
	}
	elsif ( $1 eq "methodinfo" ){
	    method_info();
	}
	elsif ( $1 eq "known_methods" ){
	    known_methods();
	}
	else{
	    print "Error: unknown qualifier after =INPUT: -$1-\n";
	}
    }

}

close Basefile;
close Outfile;



# ##################################################################
#
# Exits with Error if number of opening and closing parenthesis
# in the argument doesn't coincede
#
# ##################################################################
sub parenthesis_check{
    my $open = 0;   # number of `('
    my $close = 0;  # number of `)'
    my $retval = 0; 

    while ( $_ =~ /\(/g){
	$open++;
    }
    while ( $_ =~ /\)/g){
	$close++;
    }

    if ( $open != $close){ # Syntax error
	$retval = 0;
    }

    return $retval;
}


# ##################################################################
#
# macro assignments: known methods -- integers
#
# ##################################################################
sub known_methods{

    #initialize counter
    my $counter = 0;

    foreach my $hfile (@all_h_files){

	open INFILE, "< $Methinfopath/$hfile" or  die ("can't open file: $!");

	while ( <INFILE> ){
	    if ( $_ =~ /^\s*=METHOD\s+(\w+)/ ){
		print Outfile "#define \UUNUR_METH_$1\t($counter)\n";
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

    print Outfile "\t/* ------------------------------------------- */\n";
    print Outfile "\t/*                                             */\n";
    print Outfile "\t/* key = \"distribution\"                        */\n";
    print Outfile "\t/*                                             */\n";
    print Outfile "\t/* ------------------------------------------- */\n";
    print Outfile "\tif ( !strcmp(key, \"distr\") ){ \n\t\t";
    open INFILE, "< $Distinfofile" or  die ("can't open file: $!");

    while ( <INFILE> ){

	# search for standard distribution
	if ( $_ =~ /^\s*=DISTR\s+(\w+)/ ){
	    print Outfile "if ( !strcmp(value, \"$1\") ){\n";
	    print Outfile "\t\t\tdistr = unur_distr_$1(list, no_of_elem);\n\t\t}\n\t\telse ";
	}
    }

    # Error -- in case of unknown distribution
    print Outfile "{\n\t\t\tfprintf(stderr,\"ERROR: Unknown distribution: %s\\n\",value);\n";
    print Outfile "\t\t\tbreak;\n\t\t}\n\t}\n";

    print Outfile "\t/* ------------------------------------------- */\n";
    print Outfile "\t/*                                             */\n";
    print Outfile "\t/* key = \"domain\"                              */\n";
    print Outfile "\t/*                                             */\n";
    print Outfile "\t/* ------------------------------------------- */\n";
    print Outfile "\telse if ( !strcmp( key , \"domain\") ){\n";
    print Outfile "\t  /* list must contain exactly two entries */\n";
    print Outfile "\t  if ( no_of_elem != 2 ){\n";
    print Outfile "\t    fprintf(stderr, \"SYNTAX ERROR: To set the domain use a list with exactly 2 entries -- Standard domain is used instead.\\n\");\n";
    print Outfile "\t  }\n";
    print Outfile "\t  else if ( unur_distr_is_cont(distr) ){\n";
    print Outfile "\t    unur_distr_cont_set_domain( distr, list[0], list[1]);\n";
    print Outfile "\t  }\n";
    print Outfile "\t  else if ( unur_distr_is_discr(distr) ){\n";
    print Outfile "\t    unur_distr_discr_set_domain( distr, list[0], list[1]);\n";
    print Outfile "\t  }\n";
    print Outfile "\t  else{\n";
    print Outfile "\t    fprintf(stderr, \"ERROR: Can't set domain for this type of distribution\\n\");\n";
    print Outfile "\t    break;\n";
    print Outfile "\t  }\n";
    print Outfile "\t}\n";
    print Outfile "\telse {\n";
    print Outfile "\t  fprintf(stderr, \"SYNTAX ERROR?: Unknown key will be ignored: %s\\n\", key);\n";
    print Outfile "\t}\n";


    close INFILE;

} # end of distr_info()



# #######################################################################################
# 
# UNURAN's input as string interface: 
# generates code for generating parameter object
#
# #######################################################################################
sub method_info{


    print Outfile "\t/* ************************ */\n";
    print Outfile "\t/* determine choosen method */\n";
    print Outfile "\t/* ************************ */\n";
    print Outfile "\tif ( !strcmp( key, \"method\") ){\n";
    print Outfile "\t\tif ( no_of_elem != 0 ){\n";
    print Outfile "\t\t\tfprintf(stderr, \"SYNTAX ERROR: No list expected when setting method.\\n\");\n";
    print Outfile "\t\t}\n\t\t";
    foreach my $hfile (@all_h_files){
        # reset variable
	my $No_String_Method  = 0;
	open INFILE, "< $Methinfopath/$hfile" or  die ("can't open file: $!");


        # generate prameter object for specified method -- if possible
	while ( <INFILE> ){
	    if ( $_ =~ /^\s*=METHOD\s+(\w+)/ ){

                # check if method supports string context
		foreach my $method (@No_String_Methods){
		    if ( "\U$1" eq $method){
			$No_String_Method  = 1;
		    }
		}

		print Outfile "if ( !strcmp( value, \"\L$1\") ){\n";
		if ( $No_String_Method == 1 ){
		    # print Outfile "\t\t\tpar = unur_\L$1_new();\n";
		    print Outfile "\t\t\tfprintf(stderr, \"Method $1 not intended for usage within this string context.\\n\");\n";
		}
		else{
		    print Outfile "\t\t\tpar = unur_\L$1_new(distr);\n";
		    while ( $_ !~ /^\s*=UP\s+Methods_for_\w+/ ){
			$_ = <INFILE>;
		    }
		    $_ =~ /^\s*=UP\s+Methods_for_(\w+)/;
		}
		print Outfile "\t\t}\n\t\telse ";
	    }
	} # end of while-loop
	close INFILE;
    } # end of foreach-loop

    print Outfile "{\n\t\t\tpar->method = UNUR_METH_UNKNOWN;\n";
    print Outfile "\t\t\tfprintf(stderr, \"ERROR: Unknown method: %s\\n\",value);\n";
    print Outfile "\t\t\tbreak;\n\t\t}\n\t}\n";


    print Outfile "\n\t/* ****************************************** */\n";
    print Outfile "\t/*                                            */\n";
    print Outfile "\t/* set parameters depending on choosen method */\n";
    print Outfile "\t/*                                            */\n";
    print Outfile "\t/* ****************************************** */\n\n\t";


    foreach my $hfile (@all_h_files){

	# search routines setting parameters -- method dependent
	my $method = "UNKNOWN"; # initialize method-name to empty string
        my $check  = 0;         # Variable for checking syntax

	open INFILE,  "< $Methinfopath/$hfile" or  die ("can't open file: $!");
	while ( $_ =  <INFILE> ){

	    if ( $_ =~ /^\s*=METHOD\s+(\w+)/ ){ # h-file with method description found
		$method = "\L$1";
		my $METHOD = "\U$method";

		# key is "method"
		print Outfile "if ( par->method == UNUR_METH_$METHOD && strcmp(key, \"method\") ){\n\t\t";
	    }

	    # only parameter object passed
	    if ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*\)/ ){
		$check = parenthesis_check($_);
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		# set parameter
                print Outfile "\t\t\tcheck = unur_$method\_set_$2(par);\n";
                # check for error
		print Outfile "\t\t\tif ( ! check ){\n";
		print Outfile "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
		print Outfile "\t\t\t}\n";
		print Outfile "\t\t}\n\t\telse ";
	    }
	    # parameter object and a single value passed
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*(\w+)\s+(\w+)\s*\)/ ){
		$check = parenthesis_check($_);
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		# set parameter
		print Outfile "\t\t\tcheck = unur_$method\_set_$2(par, ($3) dblvalue);\n";
                # check for error
		print Outfile "\t\t\tif ( ! check ){\n";
		print Outfile "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
		print Outfile "\t\t\t}\n";
		print Outfile "\t\t}\n\t\telse ";
	    }
	    # parameter object and two doubles passed
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*double\s+\w+\s*,\s*double\s+\w+\s*\)/ ){
		$check = parenthesis_check($_);
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		# set parameter
		print Outfile "\t\t\tcheck = unur_$method\_set_$2(par, list[0], list[1]);\n";
                # check for error
		print Outfile "\t\t\tif ( ! check ){\n";
		print Outfile "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
		print Outfile "\t\t\t}\n";
		print Outfile "\t\t}\n\t\telse ";
	    }
	    # parameter object, size_of_list and double list passed
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*int\s+\w+\s*,\s*double\s+\*\w+\s*\)/ ){
		$check = parenthesis_check($_);
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		# set parameter
                print Outfile "\t\t\tif (no_of_elem != 0 ){\n";
		print Outfile "\t\t\t\tcheck = unur_$method\_set_$2(par, no_of_elem, list);\n";
		print Outfile "\t\t\t}\n";
		print Outfile "\t\t\telse{\n";
		print Outfile "\t\t\t\tcheck = unur_$method\_set_$2(par, (int) dblvalue, NULL);\n";
		print Outfile "\t\t\t}\n";
                # check for error
		print Outfile "\t\t\tif ( ! check ){\n";
		print Outfile "\t\t\t\tfprintf(stderr, \"ERROR: while trying to set $2\\n\");\n";
		print Outfile "\t\t\t}\n";
		print Outfile "\t\t}\n\t\telse ";
	    }
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters(.*)\)/ ){
		$check = parenthesis_check($_);
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		print Outfile "\t\t\tfprintf(stderr, \"Setting $2 for method $method not supported via this interface.\\n\");\n";
		print Outfile "\t\t}\n\t\telse ";
	    }
	    elsif ( $_ =~ /^\s*\w+\s*unur_($method)_set_(.*?)\s*\(\s*\w+.*\)/  ){
		$check = 1;
	    }

            if ( $check == 1 ){
                print "Unknown Syntax of set command in $hfile line $.: $_\n";
		exit (-1);
	    }

	} # end of while-loop
	close INFILE;
	
	# closing the method-block 
	if ( $method !~ /UNKNOWN/){
	    print Outfile "{\n\t\t\tfprintf(stderr, \"Unknown option for method $method: %s\\n\", key);\n\t\t}\n";
	    print Outfile "\t}\n\telse ";
	}

    } # end of foreach loop

    # case of unknown method
    print Outfile "if ( !strcmp(key, \"method\") && par->method == UNUR_METH_UNKNOWN){\n";
    print Outfile "\t\tfprintf(stderr, \"ERROR: No or unknown method defined.\\n\");\n\t}\n";

    close INFILE;
}




