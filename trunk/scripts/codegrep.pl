#!/usr/bin/perl

use File::Find;
use Cwd;

# Start the search for the files in these directories
@Startdirs = ('..');

# the keys are the beginning strings of the relevant files
%Disttypes = (c => 'cont', d => 'discr');

# Hash holding all informations keys are the distributions
%Distr = ();


# ---------------------------------------------------------
# Finds recursively all C-files containing the relevant
# informations. 

sub find_files{
    foreach my $type (keys %Disttypes){
	if ($File::Find::name =~ /\/($type\_([a-z]+)\.c)/ ){
	    extract_code($1, $2, $type);
       }
    }
} # end of find_files()


# ---------------------------------------------------------
# 

sub extract_code{
    my ($file, $dist, $type) = @_;

    $Distr{$dist}{"SRC"}  = $file;
    $Distr{$dist}{"TYPE"} = $Disttypes{$type};

    open FILE, $file or die "Can't open file $file\n";
    # search for beginning of parameter definitions
    for ( $_=<FILE>, $end = 0; $end == 0 ; $_=<FILE> ){
	if ( $_ =~ /^\s*\/\*\s+parameters\W*\*\// ){
          $end = 1;
	}
    }
    # read parameters
    for ( $no_of_params=0, $end=0; $end == 0; $_=<FILE> ){
	if ($_ =~ /#define\s+(\S+)\s+params\[(\d+)\]/){
	    $no_of_params++;
	    $Distr{$dist}{"PARAMS"}[$2] = $1;
	}
	elsif ($_ =~ /^\/\*/){
	    $end = 1;
	}
    }
    $Distr{$dist}{"NO_OF_PARAMS"} = $no_of_params;
    close FILE;

    open FILE, $file or die "Can't open file $file\n";  
    # search for beginning of code
    while (<FILE>){
	if ($_ =~ /^(\w+)\s*$/){
	    $functype = $1;
            $_ = <FILE>;
            if ($_ =~ /(_unur_(pmf|pdf)_(\w+))\((.*)\)/) {
		my $probfunctype = $2;
                my $probfuncname = $3;
		if ($probfuncname ne $dist){
		    print "Error -- should not happen\n";
		}
                # change to upper case letters
                $probfunctype =~ tr/a-z/A-Z/;

                $Distr{$dist}{"PROBFUNC_TYPE"} = $probfunctype;
      		$Distr{$dist}{$probfunctype}{"TYPE"} = $functype;
                $Distr{$dist}{$probfunctype}{"NAME"} = $1;
                $Distr{$dist}{$probfunctype}{"ARGS"} = $4; 
		$Distr{$dist}{$probfunctype}{"BODY"} = "\n";
                while ( $_ !~ /\/\*\s+end of _unur_\w+_\w+\(\)/ ){
		    $_ = <FILE>;
                    # ignore code lines setting double *params,
                    # blank lines and the opening and closing
                    # parentheses
                    if ( ($_ !~ /\sdouble\s+\*params\s=/) &&
                         ($_ !~ /^\s*$/) &&
                         ($_ !~ /^\{\s*/) &&
                         ($_ !~ /\/\*\s+end of _unur_\w+_\w+\(\)/) ){
			$Distr{$dist}{$probfunctype}{"BODY"} = join "",
			$Distr{$dist}{$probfunctype}{"BODY"}, $_;
		    }
                    # chech for kind of constant
                    if ($_ =~ /.*?(\w+CONSTANT)/){
			$Distr{$dist}{"CONSTANT"} = $1;
		    }

      		}
		$_ =~ /\/\*\s+end of _unur_\w+_(\w+)\(\)/;
		if ($probfuncname ne $1){
		    print "Error -- End of wrong function";
		}

	    }
            
	    
	}
 
    }
    close FILE;

} # end of extract_code()

# ---------------------------------------------------------


find (\&find_files, @Startdirs);

# write out to file
open (OUTFILE , ">distr_settings.data");

foreach $key (keys %Distr){
    foreach $key2 (keys %{ $Distr{$key}}){
        if ($key2 ne "PARAMS" && $key2 ne $Distr{$key}{"PROBFUNC_TYPE"}){
	    print OUTFILE "$key, $key2: $Distr{$key}{$key2} \n";
	    }
    }
    for (my $i=0; $i < $Distr{$key}{"NO_OF_PARAMS"}; $i++){
	print OUTFILE "$key, PARAM[$i]: $Distr{$key}{PARAMS}[$i]\n";
    }
    foreach $key3 (keys %{ $Distr{$key}{$Distr{$key}{"PROBFUNC_TYPE"}} } ){
	print OUTFILE "$key, ",
                  $Distr{$key}{"PROBFUNC_TYPE"}, ", $key3: ",
                  $Distr{$key}{$Distr{$key}{"PROBFUNC_TYPE"}}{$key3},"\n";
    }

}












