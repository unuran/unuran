#!/usr/bin/perl

# $ID$

# where relevant information is stored
$Methinfopath = "../src/methods/";
$Distinfofile = "../src/distributions/unuran_distributions.h";


open Outfile,  "> ./stringpars.c" or die "can't open file: $!";
open Basefile, "< ./stringpars.c.base" or die "can't open file: $!";

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

# #######################################################################################
#
# macro assignments: known methods -- integers
#
# #######################################################################################
sub known_methods{

    #initialize counter
    $counter = 0;

    # determine all h-files
    opendir (METHDIR, $Methinfopath) or die "can't open directory $!";
    @all_h_files = grep {/[.]h/ } readdir METHDIR;
    closedir METHDIR;

    foreach $hfile (@all_h_files){

	open INFILE, "< $Methinfopath/$hfile" or  die ("can't open file: $!");

	while ( <INFILE> ){
	    if ( $_ =~ /^\s*=METHOD\s+(\w+)/ ){
		print Outfile "#define \U$1\t($counter)\n";
		$counter++;
		break;
	    }
	}
    }
} # end of known_methods()


# #######################################################################################
# 
# for input as string: 
# generates code for generating distribution object
#
# #######################################################################################
sub distr_info{

    print Outfile "\tif ( !strcmp(key, \"distr\") ){ \n\t\t";
    open INFILE, "< $Distinfofile" or  die ("can't open file: $!");

    while ( <INFILE> ){

	# search for standard distribution
	if ( $_ =~ /^\s*=DISTR\s+(\w+)/ ){
	    print Outfile "if ( !strcmp(value, \"$1\") ){\n";
	    print Outfile "\t\t\tdistr = unur_distr_$1(list, no_of_elem);\n";

	    # obtain type of distribution
	    while ( $_ !~ /\s*=UP\s+/ ){
		$_ = <INFILE>;
	    }
	    $_ =~ /\s*=UP\s+Stddist_(\w+)/;
	    print Outfile "\t\t\ttype = UNUR_DISTR_$1;\n\t\t}\n\t\telse ";
	}
    }

    # Error -- in case of unknownd distribution
    print Outfile "{\n\t\t\tprintf(\"Unknown distribution!\\n\");\n";
    print Outfile "\t\t\tbreak;\n\t\t}\n\t}\n";



    close INFILE;

} # end of distr_info()



# #######################################################################################
# 
# for input as string: 
# generates code for generating parameter object
#
# #######################################################################################
sub method_info{


    # determine all h-files
    opendir (METHDIR, $Methinfopath) or die "can't open directory $!";
    @all_h_files = grep {/[.]h/ } readdir METHDIR;
    closedir METHDIR;


    print Outfile "\t/* ************************ */\n";
    print Outfile "\t/* determine choosen method */\n";
    print Outfile "\t/* ************************ */\n";
    print Outfile "\tif ( !strcmp( key, \"method\") ){\n";
    print Outfile "\t\tif ( no_of_elem != 0 ){\n";
    print Outfile "\t\t\tfprintf(stderr, \"Warning: List with method provided.\\n\");\n";
    print Outfile "\t\t}\n\t\t";
    foreach $hfile (@all_h_files){

	open INFILE, "< $Methinfopath/$hfile" or  die ("can't open file: $!");

	while ( <INFILE> ){
	    if ( $_ =~ /^\s*=METHOD\s+(\w+)/ ){

		print Outfile "if ( !strcmp( value, \"\L$1\") ){\n";
		print Outfile "\t\t\tmethod = $1;\n";
		if ( "\U$1" eq "UNIF"){
		    print Outfile "\t\t\tpar = unur_\L$1_new();\n";
		}
		else{
		    print Outfile "\t\t\tpar = unur_\L$1_new(distr);\n";
		}
		print Outfile "\t\t}\n\t\telse ";
	    }
	} # end of while-loop
	close INFILE;
    } # end of foreach-loop

    print Outfile "{\n\t\t\tmethod = UNKNOWN;\n";
    print Outfile "\t\t\tfprintf(stderr, \"Unknown method\\n\");\n";
    print Outfile "\t\t\tbreak;\n\t\t}\n\t}\n";


    print Outfile "\n\t/* ****************************************** */\n";
    print Outfile "\t/*                                            */\n";
    print Outfile "\t/* set parameters depending on choosen method */\n";
    print Outfile "\t/*                                            */\n";
    print Outfile "\t/* ****************************************** */\n\n\t";


    foreach $hfile (@all_h_files){

	# search routines setting parameters -- method dependent
	$method = "UNKNOWN"; # initialize method-name to empty string

	open INFILE,  "< $Methinfopath/$hfile" or  die ("can't open file: $!");
	while ( $_ =  <INFILE> ){

	    if ( $_ =~ /^\s*=METHOD\s+(\w+)/ ){ # h-file with method description found
		$method = "\L$1";
		$METHOD = "\U$method";

		# key is "method"
		print Outfile "if ( method == $METHOD && strcmp(key, \"method\") ){\n\t\t";
	    }

	    # only parameter object passed
	    if ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*\)/ ){
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		print Outfile "\t\t\tunur_$method\_set_$2(par);\n";
		print Outfile "\t\t}\n\t\telse ";
	    }
	    # parameter object and a single value passed
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*(\w+)\s+(\w+)\s*\)/ ){
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		print Outfile "\t\t\tunur_$method\_set_$2(par, ($3) dblvalue);\n";
		print Outfile "\t\t}\n\t\telse ";
	    }
	    # parameter object and two doubles passed
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*double\s+\w+\s*,\s*double\s+\w+\s*\)/ ){
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		print Outfile "\t\t\tunur_$method\_set_$2(par, list[0], list[1]);\n";
		print Outfile "\t\t}\n\t\telse ";
	    }
	    # parameter object, size_of_list and double list passed
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters\s*,\s*int\s+\w+\s*,\s*double\s+\*\w+\s*\)/ ){
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		print Outfile "\t\t\tunur_$method\_set_$2(par, no_of_elem, list);\n";    
		print Outfile "\t\t}\n\t\telse ";
	    }
	    elsif ( $_ =~ /unur_($method)_set_(.*?)\s*\(\s*UNUR_PAR\s+\*parameters(.*)\)/ ){
		print Outfile "if ( !strcmp(key, \"$2\") ){\n";
		print Outfile "\t\t\tfprintf(stderr, \"Setting $2 for method $method not supported via this interface.\\n\");\n";
		print Outfile "\t\t}\n\t\telse ";
	    }


	} # end of while-loop
	close INFILE;
	
	# closing the method-block 
	if ( $method !~ /UNKNOWN/){
	    print Outfile "{\n\t\t\tfprintf(stderr, \"Unknown option for method $METHOD: %s\\n\", key);\n\t\t}\n";
	    print Outfile "\t}\n\telse ";
	}

    } # end of foreach loop

    # case of unknown method
    print Outfile "if ( !strcmp(key, \"method\") && method == UNDEF){\n";
    print Outfile "\t\tfprintf(stderr, \"Error: No or unknown method defined.\\n\");\n\t}\n";

    close INFILE;
}
