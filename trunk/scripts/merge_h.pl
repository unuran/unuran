#!/usr/bin/perl

############################################################
# $Id$

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname <file.h> 
      
Scans <file.h> and inserts all header files found in subtree 
rooted at current working directory.
Other header files and those with names containing "config" 
are ignored and just #included. 
All comments and blank lines are removed.
The output is written on stdout.

EOM

    exit;
}

############################################################

sub h_file_header {
    print <<EOM;
/* file automatically generated by merge_h.pl                                */

/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unuran.h                                                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

EOM
}

############################################################

use FileHandle;

############################################################

# read master file name from argument list ...
$master_file = shift;
(usage and die) unless $master_file;

# header files in sub tree (except those containing "config") ...
# (files are stored in an associate array with key=filename and
# value=complete path of file.)
%header_files;
# included header files ...
%header_included;

# search subtree for header files ...
open (FILES, "find . |");
while (<FILES>) {
    chomp;
    next unless /^.*\/+(.*\.h)$/;
    next if /config/;
    $header_files{$1} = $_;  # store file and path of file
}
close FILES;

# insert file header ...
h_file_header;

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

	$header_included{$include_file} .= 1;
	
	# we do not include header files out of the subtree ...
	unless (defined( $header_files{$include_file} ) ) {
	    print STDERR "file not found in subtree ... #include\n";
	    print "\n$line\n\n";
	    next;
	}
	    
	# have found header file ...
	print STDERR "to be inserted\n";
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
    
