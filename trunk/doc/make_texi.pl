#!/usr/bin/perl
############################################################

$VERBOSE = 1;

############################################################
# $Id$
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname <top_src_dir>

Perl spript extracts documentation for unuran library from
header files and transform it into texinfo format.

The output is written on STDOUT.

For a detailed description see README.doc.

EOM
    exit;
}
############################################################

use FileHandle;

############################################################

require "DISTRIBUTION_to_texi.pl";
require "ERROR_to_texi.pl";
require "METHOD_to_texi.pl";
require "URNG_to_texi.pl";

############################################################

%section_TAGs = 
    (
     "=DISTR"   => { "scan"   => \&scan_DISTR,
		     "format" => \&format_DISTR },

     "=ERROR"   => { "scan"   => \&scan_section,
		     "format" => \&format_ERROR },

     "=METHOD"  => { "scan"   => \&scan_METHOD,
		     "format" => \&format_METHOD },

     "=URNG"    => { "scan"   => \&scan_section,
		     "format" => \&format_URNG },

     );

#...........................................................

%TAGs =
    (
     
     "=TYPE"        => { "scan" => \&scan_TYPE },

     "=REQUIRED"    => { "scan" => \&scan_do_nothing },
     
     "=OPTIONAL"    => { "scan" => \&scan_do_nothing },
     
     "=SEEALSO"     => { "scan" => \&scan_do_nothing },
     
     "=ABSTRACT"    => { "scan" => \&scan_do_nothing },
     
     "=DESCRIPTION" => { "scan" => \&scan_DESCRIPTION },
     
     "=ROUTINES"    => { "scan" => \&scan_ROUTINES },
     
     "=END"         => { "scan" => \&scan_END },
     );


############################################################

# list of all routines
$list_routines;

# scanned input
# it is stored in the form
# $in->{section TAG}->{section name}->{TAG}->{entry}
$in;

# texinfo output
$texi;

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
    print STDERR "format $section\n" if $VERBOSE > 1;
    &{$section_TAGs{$section}{"format"}}();
}

# make some modifications
foreach my $section (keys %section_TAGs) {
    print STDERR "transform $section\n" if $VERBOSE > 1;
    transform_special_strings(\($texi->{$section}));
}

# print sections
print $texi->{"=DISTR"};
print $texi->{"=METHOD"};
print $texi->{"=URNG"};
print $texi->{"=ERROR"};

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

sub scan_section {
    my $file = $_[0];
    my $handle = $_[1];
    my $line = $_[2];

    # get section TAG ($sect), short ($name) and long title ($title)
    chomp $line;          # remove trailing `\n'
    $line =~ s/^\s+//;    # trim blanks from beginning of line
    $line =~ s/\s+$//;    # trim blanks from end of line
    (my $sect, my $name, my $title) = split /\s+/, $line, 3;

    # print section on screen
    print STDERR "$sect: $name -- $title\n\t" if $VERBOSE;

    # store section
    $in->{$sect}->{$name}->{"=NAME"} = $title;

    # store file name
    $in->{$sect}->{$name}->{"=FILE"} = $file;

    # scan all subsections 
    my $this_TAG = "=END";  # add =END tag to section tag
    my $in_line;            # store scanned line

    while (my $next_TAG = scan_subsection($handle,\$lines)) {
	if ($this_TAG ne "=END") {
	    print STDERR "  $this_TAG" if $VERBOSE;
	    if ($TAGs{$this_TAG}) {
		# store subsection text
		$in->{$sect}->{$name}->{$this_TAG} .= $lines;
	    }
	    else {
		print STDERR "  invalid!!\n\n\t" if $VERBOSE;
	    }
	}
	# anything between =END and next tag is ignored
	$this_TAG = $next_TAG;
    }

    # close line on screen
    print STDERR "\n" if $VERBOSE;
  
    # scan and format all sections
    foreach $tag (keys %TAGs) {
	&{$TAGs{$tag}{"scan"}}( \($in->{$sect}->{$name}->{$tag}) );
    }

    return;
} # end of scan_section() 

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
		    "UNUR_URNG",
		    "UNUR_FUNCT_CONT",
		    "FILE",
		    "extern",
		    "struct",
		    "const",
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

	# "verbatim" block ?
	if ($block =~/^\s*==DOC/) {
	    # remove subTAG and copy text verbatim
	    $block =~ s/^\s*==DOC\s*//;
	    $proc .= "$block\n\n";
	    next;  # block
	}

	# split into function declaration and its description
	(my $fkt, my $body) = split /\;/, $block, 2; 

	# check function type
	my $type_ok = 0;
	foreach $type (@C_TYPES) {
	    if ("$fkt " =~ /^\s*$type /) {
		$type_ok = 1;
		last;
	    }
	}
	# if this is not a valid type, skip to next block.
	next unless $type_ok;

	# add blanks around braces
	$fkt =~ s/\(/ \( /g;
	$fkt =~ s/\)/ \) /g;

	# move pointer '*' from name to type
	$fkt =~ s/\s+(\*+)/$1 /g;

	# get argument list of function
	(my $fkt_decl, my $fn_args) = split /\s+[\(\)]\s+/, $fkt, 3; 
	my @argslist = split /\s*\,\s*/, $fn_args;
	
	# $fkt_decl should contain of at least two words
	unless ($fkt_decl =~ /^\s*(.*)\s+([\*\w+]+)\s*$/ ) {
	    die "type or name missing for function: '$fkt_decl'";
	}
	
	# get function type and name
	# the first part in $fkt_decl is the function type,
	# the last word is the function name
	my $fn_type = $1;
	my $fn_name = $2;

	# routine name must be unique
	die "Function defined twice: $fn_name" if $list_routines->{$fn_name};

	# store in list of all routines
	$list_routines->{$fn_name} = 1;

	# write entry
	my $first = 1;
	if (@argslist) {
	    # this is a function with arguments
	    # make anchor
	    $fkt_block .= "\@anchor{funct:$fn_name}\n";
	    # make texinfo tag
	    $fkt_block .= (($defblock_open) ? "\@deftypefnx" : "\@deftypefn");
	    $fkt_block .= " %%%Function%%% \{$fn_type\} $fn_name (";
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
	}
	else {
	    # this is a function does not have arguments
	    # maybe it is an variable
	    # make anchor
	    $fkt_block .= "\@anchor{var:$fn_name}\n";
	    # make texinfo tag
	    $fkt_block .= (($defblock_open) ? "\@deftypevarx" : "\@deftypevar");
	    $fkt_block .= " \{$fn_type\} $fn_name\n";
	    $fkt_block .= "\@findex $fn_name\n";
	}
	# description
	if ($body) {
	    $fkt_block .= "$body\n";
	    if (@argslist) {
		$fkt_block .= "\@end deftypefn\n"; }
	    else {
		$fkt_block .= "\@end deftypevar\n"; }
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
	
    die "last function without description: $fkt_block" if $defblock_open;

    # store new lines
    $$entry = $proc;

    return;

} # end of scan_ROUTINES()

############################################################

sub scan_TYPE {
    my $entry = $_[0];    # pointer to TYPE entry

    # we are interested in the very first word only ...
    if ( $$entry =~ /^\s*([A-Za-z]+)/ ) {
	$$entry = $1;
    }
    else {   # the string does not start with a word --> clear string
	$$entry = '';
    }
    return;
} # end if scan_TYPE()

############################################################

sub scan_DESCRIPTION {
    return;
} # end of scan_DESCRIPTION() 

############################################################

sub scan_END {
    # there is nothing to do here!
    return;
}

############################################################

sub scan_do_nothing {
    return;
} # end of scan_do_nothing()

############################################################

sub transform_special_strings {
    my $line = $_[0];

    # trim blanks
    $$line =~ s/[ \t\r\f]+\n/\n/g;

    # NULL --> @code{NULL}
    # TRUE --> @code{TRUE}
    # FALSE --> @code{FALSE}
    $$line =~ s/ (NULL|TRUE|FALSE)/ \@code\{$1\}/g;
    $$line =~ s/^(NULL|TRUE|FALSE)/\@code\{$1\}/g;

    # transform (\w+)\(\)   --> @command($1)
    my $first = "\n\@ifhtml\n\@ref\{funct:";
    my $middle = "\}\n\@end ifhtml\n\@ifnothtml\n";
    my $last = "\n\@end ifnothtml\n";

    $$line =~ s/\s+(\w+)\(\)([\.\,\;\:])\s*/$first$1,\@command\{$1\}$2$middle\@command\{$1\}$2$last/g;
    $$line =~ s/\s+(\w+)\(\)(\n|\s*)/$first$1,\@command\{$1\}$middle\@command\{$1\}$last/g;

    return;
} # end of transform_special_strings()

############################################################
