#!/usr/bin/perl
############################################################
# $Id$
############################################################

use strict;

############################################################
# constants

############################################################

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname <file.tex> 
      
   Read TeX file that contains a single picture and
   converts it into EPS, PDF, PNG and TXT format.
   It is assumed that these pictures are produced by pstricks 
   commands.
   Each of these is written into the respective file <file.FORMAT>
   where FORMAT is "eps", "pdf", "png" or "txt".
   The path information is stripped and the files are written in
   working directory (i.e. the directory where the command 
   has been executed).

EOM
    exit -1;
}

############################################################

# read name of TeX file from argument list ...
my $tex_file = shift;
(usage and die) unless $tex_file;

# get other file names
my $base = $tex_file;
$base =~ s[^.*/]{};
$base =~ s/(.*)\..*$/$1/;
my $eps_file = "$base.eps";
my $pdf_file = "$base.pdf";
my $png_file = "$base.png";
my $txt_file = "$base.txt";

print STDERR "TeX file = $tex_file\n";
print STDERR "EPS file = $eps_file\n";
print STDERR "PDF file = $pdf_file\n";
print STDERR "PNG file = $png_file\n";
print STDERR "TXT file = $txt_file\n";
print STDERR "----------------------\n\n";

############################################################

# Read dimensions

my $xunit = 0;
my $yunit = 0;
my $x0 = 0;
my $y0 = 0;
my $x1 = 0;
my $y1 = 0;
open TEX, "$tex_file";
while (<TEX>) {
    if ( /\\psset\{xunit=(\d+)mm\,yunit=(\d+)mm\}/ ) {
	$xunit = $1;
	$yunit = $2;
    }
    if ( /\\begin\{pspicture\}\((-?[\d\.]+),\s*([-?\d\.]+)\)\((-?[\d\.]+),\s*([-?\d\.]+)\)/ ) {
	$x0 = $1; $y0 = $2;
	$x1 = $3; $y1 = $4;
    }
}
close TEX;

($xunit) or die "Cannot find units";
($x0) or die "Cannot find bounding box";

my $xlength = (5+($x1-$x0)*$xunit)."mm";
my $ylength = (5+($y1-$y0)*$yunit)."mm";

my $xsize = (15+($x1-$x0)*$xunit)."mm";
my $ysize = (15+($y1-$y0)*$yunit)."mm";

print STDERR "xunit = $xunit mm\n";
print STDERR "yunit = $yunit mm\n";

print STDERR "bounding box = ($x0,$y0) x ($x1,$y1)\n";
print STDERR "bounding box = ($x0,$y0) x ($x1,$y1)\n";
print STDERR "             = $xlength x $ylength\n";

print STDERR "----------------------\n\n";

############################################################

# Make temp texfile  
my $tmp_base = "tmp-$base-$$";
my $tmp_tex_file = "$tmp_base.tex";
my $tmp_dvi_file = "$tmp_base.dvi";
my $tmp_ps_file = "$tmp_base.ps";

print STDERR "make temp file $tmp_tex_file\n";

open TEX, ">$tmp_tex_file" or die "Cannot open file $tmp_tex_file";
print TEX
    "\\documentclass{article}\n",
    "\\usepackage{pstricks,pst-plot,pst-node}\n",
    "\\usepackage{amsfonts}\n",
    "\\usepackage{amssymb}\n",
    "%%\\usepackage[papersize={$xsize,$ysize},total={$xlength,$ylength},noheadfoot,verbose]{geometry}\n",
    "\\pagestyle{empty}\n",
    "\\parindent 0pt\n",
    "\\parskip 0pt\n",
    "\\newcommand{\\onehalf}{{\\textstyle\\frac{1}{2}}}\n",
    "\\begin{document}\n",
    "\\mbox{}\\vfill\n",
    "\\centering\n",
    "\\input{$tex_file}\n",
    "\\vfill\\mbox{}\n",
    "\\end{document}\n";
close TEX;

# Make dvi
system("latex $tmp_tex_file");

# Make eps
system("dvips -o $tmp_ps_file $tmp_dvi_file");
system("ps2epsi $tmp_ps_file $eps_file");

############################################################


############################################################
# Clear working space

#system("rm -v $tmp_base.*");

exit 0;

############################################################

