#!/usr/bin/perl
############################################################
# $Id: merge_h.pl 2502 2005-07-26 14:16:49Z leydold $
############################################################

use strict;

my $DEBUG = 0;

############################################################

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname <file> 
      
Removes all comments and blank lines from C file <file>

EOM

    exit -1;
}

############################################################

# read file name from argument list ...
my $file = shift;
(usage and die) unless $file;
print "stripping $file ...\n";

# open file ...
open IN, $file or die "cannot open file $file\n";
    
# read file ...
my $content = '';
while (<IN>) { $content .= $_; }
close IN;

# remove all comments and empty lines ...
$content =~ s {/\*.*?\*/} []gsx;
$content =~ s /\n\s*\n/\n/gsx;

# print into file
open OUT, ">$file";
print OUT $content;
close OUT;

# end
exit 0;

############################################################
