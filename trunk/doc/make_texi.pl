#!/usr/bin/perl

$VERBOSE = 0;

############################################################
# $Id$

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname ?????

EOM

    exit;
}

############################################################

use FileHandle;

############################################################

require "METHOD_to_texi.pl";
require "DISTRIBUTION_to_texi.pl";

############################################################

%section_TAGs = 
    (
     "=METHOD"       => { "scan"   => \&scan_METHOD,
			  "format" => \&format_METHOD },

     "=DISTRIBUTION" => { "scan"   => \&scan_DISTR,
			  "format" => \&format_DISTR }
     );

############################################################

# read directory name from argument list ...
$top_dir = shift;
(usage and die) unless $top_dir;

# header files in directory tree ...
# (files are stored in an associate array with key=filename and
# value=complete path of file.)
%header_files;

# search subtree for header files ...
open (FILES, "find $top_dir |") or die "cannot find directory \"$top_dir\".";
while (<FILES>) {
    chomp;
    next unless /^.*\/+(.*\.h)$/;
    next if /source/;        # we are not interested in these files
    next if /struct/;        # we are not interested in these files
    $header_files{$1} = $_;  # store file and path of file
}
close FILES;

# now scan all header files
foreach $file (sort keys %header_files) {
    scan_file($header_files{$file});
}

# format all sections
foreach my $section (keys %section_TAGs) {
    print STDERR "format $section\n" if $VERBOSE;
    &{$section_TAGs{$section}{"format"}}();
}

# print sections
print $texi_DISTRs;
print $texi_METHODs;

# end of job
exit 0;

############################################################

# scan given file ...
sub scan_file {
    my $file = $_[0];
    my $handle = new FileHandle;
    
    # open file ...
    open $handle, $file or die "cannot find file $file\n";

    # print info
    print STDERR "$file " if $VERBOSE;
    
    # search for section TAG
    while (<$handle>) {
	last if /^\s*=[A-Z]/;
    }

    # is there a section TAG ?
    unless (/^\s*(\/\*)?\s*(=[A-Z]+)/) {
	print STDERR "-- no doc --\n" if $VERBOSE; 
	return; 
    }
    my $line = $_;

    # have found tag
    my $tag = $2;
    print STDERR "  $tag: " if $VERBOSE;

    # scan if TAG is known
    if ($section_TAGs{$tag}) {
	&{$section_TAGs{$tag}{"scan"}}($file, $handle, $line);
    }
    else {
	print STDERR "invalid section TAG!!\n\n" if $VERBOSE;
    }

    # close file
    close $handle;

} # end of scan_file() 

############################################################

sub scan_END {
    # there is nothing to do here!
    return;
}

############################################################

sub scan_subsection {
    my $handle = $_[0];
    my $lines = $_[1];
    
    # clear lines buffer
    $$lines = '';
    
    # scan file
    while (<$handle>) {
	# next TAG ?
	if (/^\s*(\/\*)?\s*(=[A-Z]+)/) {
	    return $2;
	}

	# process input line 
	chomp;
	$_ =~ s/^\s+//;      # chop off heading and ...
	$_ =~ s/\s+$//;      # trailing blanks
	
	# add to lines 
	$$lines .= "$_\n";
    }

    # EOF
    return '';

} # end of scan_subsection() 

############################################################

sub scan_ROUTINES {
    my $entry = $_[0];    # pointer to TYPE entry

    # valid C data types 
    my @C_TYPES = ( "UNUR_PAR",
		    "UNUR_GEN",
		    "UNUR_DISTR",
		    "UNUR_FUNCT_CONT",
		    "struct",
		    "void",
		    "int",
		    "double",
		    "float",
		    "long",
		    "char",
		    "short",
		    "unsigned",
		    "signed" );

    # remove double blank lines
    $$entry =~ s/\n\s*\n\s*\n/\n\n/g;

    # simply comment markers: /** --> /*
    $$entry =~ s/\/\*+/\/\*/g;
    $$entry =~ s/\*+\//\*\//g;

    # remove blanks and new lines inside comments at markers
    $$entry =~ s/\/\*\s+/\/\*/g;
    $$entry =~ s/\s+\*\//\*\//g;

    # remove blank lines at begin and end of entry
    $$entry =~ s/^(\s*\n)+//g;
    $$entry =~ s/(\s*\n)+$//g;

    # if there are blocks of comments, separated by empty lines
    # then delete all but the first block of comments
#    $$entry =~ s/\*\/\s*\n\s*\n\s*\/\*[.\n]*\*\//\*\//g;

    # split into blocks
    my @blocks = split /\*\/\s*\n\s*\n/, $$entry ;

    # store processed text
    my $proc = '';

    # deftypefn block not closed
    my $defblock_open = 0;
    
    # process blocks
    my $fkt_block = '';
    foreach $block (@blocks) {
	# remove anyting that starts with an #
	$block =~ s/^\#.*$//mg;

	# remove all comment markers
	$block =~ s/\*\/\s*//g;
	$block =~ s/\s*\/\*//g;

	# skill over empty blocks
	next unless $block;

	# split into function declaration and its description
	(my $fkt, my $body) = split /\;/, $block, 2; 

	# add blanks around braces
	$fkt =~ s/\(/ \( /g;
	$fkt =~ s/\)/ \) /g;

	# get function type
	(my $fn_type, my $rest) = split /[\s]+/, $fkt, 2;

	# pointer ?
	my $pointer = (($fn_type =~ /(\*+)$/) ? $1 : '');

	# check function type
	my $type_ok = 0;
	foreach $type (@C_TYPES) {
	    if ($fn_type eq $type) {
		$type_ok = 1;
		last;
	    }
	}
	
	# if this is not a valid type, skip to next block.
	next unless $type_ok;

	# get function name
	(my $fn_name, $rest) = split /[\s\(]+/, $rest, 2;

	# pointer ?
	unless ($pointer) {
	    if ($fn_name =~ /^(\*+)/) {
		$pointer = $1;
		$fn_name =~ s/^(\*+)//;
	    }
	}

	# argument list of function
	$rest =~ s/\s*\)\s*$//;     # remove trailing brace ')'
	$rest =~ s/\s+(\*+)/$1 /g;  # move pointer '*' from argument name to type

	my @argslist = split /\s*\,\s*/, $rest;
	
	# write entry
	my $first = 1;
	$fkt_block .= (($defblock_open) ? "\@deftypefnx" : "\@deftypefn");
	$fkt_block .= " %%%Function%%% $fn_type$pointer $fn_name (";
	foreach $arg (@argslist) {
	    (my $type, my $arg_name) = split /\s+/, $arg, 2;
	    if ($first) { $first = 0; }
	    else { $fkt_block .= ", "; }
	    if ($arg_name) {
		$fkt_block .= "$type \@var\{$arg_name\}";
	    }
	    else {
		$fkt_block .= "$type";
	    }
	}
	$fkt_block .= ")\n";
	# description
	if ($body) {
	    $fkt_block .= "$body\n";
	    $fkt_block .= "\@end deftypefn\n";
	    $defblock_open = 0;
	    # for info file
	    my $fkt_string = $fkt_block;
	    $fkt_string =~ s/%%%Function%%%/Function/g;
	    $proc .= "\@ifinfo\n$fkt_string\@end ifinfo\n";
	    # for other output formats
	    $fkt_string = $fkt_block;
	    $fkt_string =~ s/%%%Function%%%/--/g;
	    $proc .= "\@ifnotinfo\n$fkt_string\@end ifnotinfo\n\n";
	    # clear block
	    $fkt_block = '';
	}
	else { 
	    $defblock_open = 1;
	}
    }

    die "last function without description" if $defblock_open;

    # store new lines
    $$entry = $proc;

    return;

} # end of scan_ROUTINES()

############################################################

