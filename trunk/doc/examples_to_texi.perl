#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl spript writing the example testx.c to the texi-file   
# qstart_testx.texi
# 
# usage:  ./examples_to_texi.perl ../examples/testx.c
#      (examples_to_texi.perl is in the directory /unuran/doc/
# Input:   c-files (examples)
# Output:  $OUTFILE contains docomentation in texinfo format
# 
# E. JANKA
# $Id$
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

while($_ = <>)
{ 
    if ($ARGV ne $longname){
	$longname = $ARGV;
	$longname =~ /(\w*).c$/;
	$name = $1;
	close OUTFILE;
	open(OUTFILE, ">qstart_$name.texi");
    }

    $_ = join '@{', split /\{/, $_;
    $_ = join '@}', split /\}/, $_;
    print OUTFILE;
}











