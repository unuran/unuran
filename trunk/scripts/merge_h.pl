#!/usr/bin/perl

############################################################
# $Id$

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname <file.h> 
      
Scans <file.h> and inserts all header file found in subtree 
rooted at current working directory.
Files with name containing "config" are ignored.
All comments and blank lines are removed.
The output is written on stdout.

EOM

    exit;
}

############################################################

use FileHandle;
use File::Find;

############################################################

# read master file name from argument list ...
$master_file = shift;
(usage and die) unless $master_file;

# header files in sub tree ...
%header_files;
# included header files ...
%header_included;

# search subtree for header files ...
open (FILES, "find . |");
while (<FILES>) {
    chomp;
    next unless /^.*\/+(.*\.h)$/;
    next if /config/;
    $header_files{$1} = $_;
}

# scan master file ...
scan_file ($master_file);

exit 0;

############################################################

# scan given file ...
sub scan_file {
    my $file = $_[0];
    my $handle = new FileHandle;
    
    print STDERR "scanning $file ...\n";

    # open file ...
    open  $handle, $file or die "cannot find file $file\n";
    
    # read file ...
    my $content = '';
    while (<$handle>) {
	$content .= $_;
    }
    close $handle;

    # remove all comments and empty lines ...
    $content =~ s {/\*.*?\*/} []gsx;
    $content =~ s /\n\s*\n/\n/gsx;

    # split into lines ...
    my @lines = split /\n/, $content;

    foreach $line (@lines) {
	unless ($line =~ /^\#include\s+<(.*)>/) {
	    print "$line\n";
	    next;
	}

	# have found a file to be included ...
	my $include_file = $1;
	print STDERR "$include_file  ";

	# we include header files only once ...
	if (defined( $header_included{$include_file} ) ) {
	    print STDERR "already included ... skip\n";
	    next;
	}
	
	# we do not include header files out of the subtree ...
	unless (defined( $header_files{$include_file} ) ) {
	    print STDERR "file not found in subtree ... ignore\n";
	    print "/*-----*/\n";
	    print "$line\n";
	    print "/*-----*/\n";
	    next;
	}
	    
	# have found header file ...
	print STDERR "to be inserted\n";
	$header_included{$include_file} .= 1;
	print "/*-----*/\n";
	print "/* `$include_file' */\n";

	# scan header file ...
	scan_file ($header_files{$include_file});

	print "/* end of `$include_file' */\n";
	print "/*-----*/\n";
    }

}

############################################################

############################################################
    
