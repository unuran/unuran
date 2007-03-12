#!/usr/bin/perl
##############################################################################
#                                                                            #
#           UNURAN -- Universal Non-Uniform Random number generator          #
#                                                                            #
##############################################################################
#                                                                            #
#   FILE:    makedll_def.pl                                                  #
#                                                                            #
#   Read all UNURAN header files and create .def file for .dll               #
#                                                                            #
##############################################################################
#                                                                            #
#   Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold              #
#   Department of Statistics and Mathematics, WU Wien, Austria               #
#                                                                            #
#   This program is free software; you can redistribute it and/or modify     #
#   it under the terms of the GNU General Public License as published by     #
#   the Free Software Foundation; either version 2 of the License, or        #
#   (at your option) any later version.                                      #
#                                                                            #
#   This program is distributed in the hope that it will be useful,          #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#   GNU General Public License for more details.                             #
#                                                                            #
#   You should have received a copy of the GNU General Public License        #
#   along with this program; if not, write to the                            #
#   Free Software Foundation, Inc.,                                          #
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                   #
#                                                                            #
##############################################################################

use strict;

my $DEBUG = 0;

############################################################
# constants

my $DEP_file = ".dep-unuran_def";

############################################################

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname <file.h> 
      
Scans <file.h> and inserts all header files found in subtree 
rooted at current working directory.
Other header files and those with names containing "config" 
are ignored.
The function declaration and global variables are extracted
and .def file is created.
The output is written on stdout.

EOM

    exit -1;
}

############################################################

# dependencies
my $DEP = "unuran.def: ";

############################################################

use FileHandle;

############################################################
# global variables
	
# name of library
my $LIBRARY = "LIBRARY \"unuran\"\n";

# DATA section (to be set manually!!)
my $DATA = 
    "\tunur_errno\tDATA\n";


#    "\tunur_errno\tDATA\n".
#    "\tINFINITY\tDATA\n";

# EXPORTS section (functions; set automatically)
my $EXPORTS;
my @EXPORTS;

############################################################

# read master file name from argument list ...
my $master_file = shift;
(usage and die) unless $master_file;

# header files in sub tree (except those containing "config") ...
# (files are stored in an associate array with key=filename and
# value=complete path of file.)
my %header_files;
# included header files ...
my %header_included;

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
#h_file_header;

# scan master file ...
scan_file ($master_file,0);

# insert bottom of file
#h_file_bottom;

# write file
print $LIBRARY;
print "EXPORTS\n";
print $DATA;
foreach my $e (sort @EXPORTS) { $EXPORTS .= "\t$e\n"; }
print $EXPORTS;

# write dependencies
#open DEP, ">$DEP_file" or die "Cannot open file for writing: $DEP_file";
#print DEP "$DEP\n";
#close DEP;

exit 0;

############################################################

# scan given file ...
sub scan_file {
    my $file = $_[0];
    my $level = $_[1]+1;
    my $handle = new FileHandle;

    print STDERR "$file\n" if $DEBUG;
    
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

    foreach my $line (@lines) {
	unless ($line =~ /^\#include\s+[<\"](.*)[>\"]/) {
	    # remove all blank lines 
	    next if $line !~ /\w/;
	    # remove all (other) preprocessor directives
	    next if $line =~ /^\s*\#/;
	    # find all non-functions
	    next if $line !~ /\(.+/;
	    next if $line =~ /^\s*typedef\s+/;
	    # internal functions
	    next if $line =~ /[\s\*]_unur/;
	    # exported functions
	    $line =~ /.*?[\s\*]+(\w+)\s*\(/;
	    push @EXPORTS, $1;
	    next;
	}

	# have found a file to be included ...
	my @tmp = split /\//, $1;
	my $include_file = pop @tmp; 
	print STDERR "$include_file  " if $DEBUG;

	# we do not include header files out of the subtree ...
	unless (defined( $header_files{$include_file} ) ) {
	    print STDERR "file not found in subtree ... #include\n" if $DEBUG;
	    next;
	}
	    
	# we include header files only once ...
	if (defined( $header_included{$include_file} ) ) {
	    print STDERR "already included ... skip\n" if $DEBUG;
	    next;
	}

	$header_included{$include_file} .= 1;
	
	# have found header file ...
	print STDERR "to be inserted\n" if $DEBUG;

	# add dependency
        $DEP .= "$header_files{$include_file} ";

	# scan header file ...
	scan_file ($header_files{$include_file},$level);

    }

}

############################################################

############################################################
    
