#!/usr/bin/perl

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perl script generating a list of preset distributions
# available 
# 
# Call:  ./distr.perl ../src/distributions/*.c
# Input:  c-files with code of distributions 
# Output: description of distributions in texi-format
# 
# E.JANKA  September 2000
# $Id$
#
# This script scans c-files for description of distributions
# and produces texi-output
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# file for  output:
open(OUTFILE, ">qstart_distributions.texi");


while($_ = <>)
{ 
    chomp;


}











