#!/usr/bin/perl
############################################################

use strict;
my $VERBOSE = 1;

############################################################
# $Id$
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname <top_src_dir>

Perl script extracts documentation for unuran library from
header files and transform it into texinfo format.

The output is written on STDOUT.

For a detailed description see README.

EOM
    exit;
}
############################################################

use FileHandle;

############################################################

my %node_TAGs = 
    (
     "=TOP"      => { "format" => \&texi_NODE },

     "=NODE"     => { "format" => \&texi_NODE },

     "=DISTR"    => { "format" => \&texi_NODE },

     "=METHOD"   => { "format" => \&texi_NODE },

     "=APPENDIX" => { "format" => \&texi_NODE },

     );

#...........................................................

my %TAGs =
    (
     
     "=UP"          => { "scan" => \&scan_UP },
     
     "=DESCRIPTION" => { "scan" => \&scan_do_nothing },
     
     "=ROUTINES"    => { "scan" => \&scan_ROUTINES },
     
     "=REQUIRED"    => { "scan" => \&scan_chop_blanks },

     "=OPTIONAL"    => { "scan" => \&scan_chop_blanks },

     "=SPEED"       => { "scan" => \&scan_chop_blanks },
     
     "=SEEALSO"     => { "scan" => \&scan_do_nothing },
     
     "=ABSTRACT"    => { "scan" => \&scan_do_nothing },
     
     "=REF"         => { "scan" => \&scan_REF },

     "=PDF"         => { "scan" => \&scan_PDF },
     "=PMF"         => { "scan" => \&scan_PDF },
     "=CONST"       => { "scan" => \&scan_PDF },
     "=CDF"         => { "scan" => \&scan_PDF },
     "=DOMAIN"      => { "scan" => \&scan_DOMAIN },
     "=FPARAM"      => { "scan" => \&scan_FPARAM },
     "=STDGEN"      => { "scan" => \&scan_STDGEN },

     "=END"         => { "scan" => \&scan_do_nothing },
     );

#############################################################

# list of nodes
my $LIST_nodes;
my @LIST_nodes_sorted;

# list of routines
my $LIST_routines;

# scanned input
# it is stored in the form
# $IN->{node name}->{TAG} = entry
my $IN;

# texinfo output
my $TEXI = "\@c automatically generated by `make_texi.pl'\n\n";

#############################################################

# read directory name from argument list ...
my $top_dir = shift;
(usage and die) unless $top_dir;

# header files in directory tree ...
# files are stored with key = filename and
# value = complete path of file.
my %header_files;
scan_dir($top_dir);

# now scan all header files
foreach my $file (sort keys %header_files) {
    scan_file($file);
}

# check node structure
check_node_structure();

# format texinfo output
texi_node("TOP",0);

# make some modifications
transform_special_strings(\($TEXI));

# write all output on STDOUT
print $TEXI;

# end of job
exit 0;

#############################################################
# get all header files in all sub directories ...
# 

sub scan_dir {
    my $dir = $_[0];

    # search subtree for header files ...
    open (FILES, "find $dir |") or die "cannot find directory \"$dir\".";

    while (<FILES>) {
	chomp;
	next unless /^.*\/+(.*\.d?h)$/;
	next if /source/;        # we are not interested in these files
	next if /struct/;        # we are not interested in these files
	$header_files{$1} = $_;  # store file and path of file
    }
    close FILES;
} # end of scan_dir()

#############################################################
# scan given file ...
#

sub scan_file {
    my $file = $_[0];
    my $file_handle = new FileHandle;
    
    # open file ...
    open $file_handle, $header_files{$file} or die "cannot find file $file\n";

    # scan file
    my $have_found_TAG = 0;
    while (1) {

	# search for node TAG
	while (<$file_handle>) {
	    last if /^\s*=[A-Z]/;
	}
	
	# is there a node TAG ?
	unless (/^\s*(\/\*)?\s*(=[A-Z]+)/) {
	    print STDERR "[$file] -- no doc --\n" if $VERBOSE and not $have_found_TAG;
	    close $file_handle or die "wrong file handle";
	    return; 
	}
	
	# prepare line with node TAG
	s/^\s*//;  s/\s*$//;    # trim blanks
	my $node_line = $_;     # store line
	
	# have found node TAG
	(my $node_type, my $node_name, my $node_title) = split /\s+/, $node_line, 3;

	# print info on screen
	print STDERR "[$file]  $node_type ($node_name):\t" if $VERBOSE;
	
	# is this a valid node TAG
	unless ($node_TAGs{$node_type}) {
	    print STDERR "\n\t[$file] invalid node TAG `$node_type' (skip over rest of file)\n";
	    close $file_handle or die "wrong file handle";
	    return; 
	}
	
	# node must be unique
	if ($LIST_nodes->{$node_name}) {
	    print STDERR "\n\t[$file] node `$node_name' used twice (skip over node)\n";
	    close $file_handle or die "wrong file handle";
	    return; 
	} 
	
	# store node name
	$LIST_nodes->{$node_name} = 1;
	$have_found_TAG = 1;

	# read till end of node
	my $node_text = $node_line."\n";
	while (<$file_handle>) {
	    if (/^\s*(\/\*)?\s*=EON/) {
		# End-Of-Node TAG
		last;
	    }
	    else {
		# add line to stored text
		s/^\s*//;  s/\s*$//;    # trim blanks
		$node_text .= $_."\n";  # store content of node
	    }
	}
	
	# now scan node text;
	scan_node($node_text, $file);
    }

} # end of scan_file() 

#############################################################

sub scan_node {
    my $node_text = $_[0];
    my $file = $_[1];

    # add =END TAG to text (for convinience)
    $node_text .= "\n=END";

    # split into lines
    my @lines = split /\n/, $node_text;

    # first line contains node definition
    my $node_line = shift @lines;

    # get node name
    (my $node_type, my $node_name, my $node_title) = split /\s+/, $node_line, 3;

    # store node 
    $IN->{$node_name}->{"=NODE_TYPE"} = $node_type;
    $IN->{$node_name}->{"=TITLE"} = $node_title;
    $IN->{$node_name}->{"=FILE"} = $file;

    # scan all TAGs (node sections)
    my $this_TAG = "=END";  # add =END tag to node TAG

    foreach my $l (@lines) {
	# next TAG ?
	if ($l =~ /^\s*(\/\*)?\s*(=[A-Z]+)\s*(.*)$/) {
	    # store next TAG
	    $this_TAG = $2;
	    # save rest of line
	    $l = $3;
	    print STDERR "  $this_TAG" if $VERBOSE and $this_TAG ne "=END";
	    unless ($TAGs{$this_TAG}) {
		print STDERR "\n\t[$file] $node_name: invalid TAG `$this_TAG'\n";
		last;
	    }
	}

	# append to stored lines
	# (except for =END TAG)
	unless ($this_TAG eq "=END") {
	    $IN->{$node_name}->{$this_TAG} .= $l."\n";
	}
    }

    # scan and format all node sections
    foreach my $tag (keys %TAGs) {
	&{$TAGs{$tag}{"scan"}} ($node_name,$tag);
    }

    # close line on screen
    print STDERR "\n" if $VERBOSE;

    # there must be an UP TAG
    die "UP missing for node $node_name" unless $IN->{$node_name}->{"=UP"};
    
} # end of scan_node()

#############################################################
# check node structure
#
sub check_node_structure {

    # we need a TOP node
    die "TOP node missing" unless $IN->{TOP};
    
    # add node `(dir)' to list of node
    $LIST_nodes->{"(dir)"} = 1;

    # check =UP nodes
    foreach my $n (keys %{$IN}) {
	unless ($LIST_nodes->{$IN->{$n}->{"=UP"}}) {
	    print STDERR "$n: invalid UP node `".$IN->{$n}->{"=UP"}."'\n";
	}
    }

    # sort nodes
    @LIST_nodes_sorted = sort nodes_by_order_tag keys %$LIST_nodes;

    # print node structure
    print_node("TOP","");

} # end of check_node_structure() 

sub print_node {
    my $node = $_[0];    # node for which all subnodes should be printed
    my $indent = $_[1];  # line indent

    # print node on screen
    print STDERR "$indent$node: (".$IN->{$node}->{"=TITLE"}.")\n" if $VERBOSE;

    # search for all subnodes
    foreach my $n (@LIST_nodes_sorted) {
	next unless $IN->{$n}->{"=UP"} eq $node;
	print_node($n,"$indent\t");
    }

} # end if print_node() 

#############################################################
# compare two nodes by their order TAG (lexicographically)
#

sub nodes_by_order_tag {
    $IN->{$a}->{"=ORDERING"} cmp $IN->{$b}->{"=ORDERING"};
}

#############################################################
# make texinfo output
#

sub texi_node {
    my $node = $_[0];    # node for which all subnodes should be printed
    my $level = $_[1];   # level for node.

    # make section type
    my $section = '';
    if ($IN->{$node}->{"=NODE_TYPE"} eq "=APPENDIX") {
      LEVEL: {
	  if ($level == 0) { # TOP node
	      $section = "\@top ";
	      last LEVEL; }
	  if ($level == 1) { 
	      $section = "\@appendix ";
	      last LEVEL; }
	  if ($level == 2) { 
	      $section = "\@appendixsec";
	      last LEVEL; }
	  if ($level == 3) { 
	      $section = "\@appendixsubsec ";
	      last LEVEL; }
	  if ($level >= 4) { 
	      $section = "\@appendixsubsubsec ";
	      last LEVEL; }
      }
    }
    else {
      LEVEL: {
	  if ($level == 0) { # TOP node
	      $section = "\@top ";
	      last LEVEL; }
	  if ($level == 1) { 
	      $section = "\@chapter ";
	      last LEVEL; }
	  if ($level == 2) { 
	      $section = "\@section ";
	      last LEVEL; }
	  if ($level == 3) { 
	      $section = "\@subsection ";
	      last LEVEL; }
	  if ($level >= 4) { 
	      $section = "\@subsubsection ";
	      last LEVEL; }
      }
    }

    # make title
    my $title = $IN->{$node}->{"=TITLE"};
    if ($IN->{$node}->{"=NODE_TYPE"} eq "=METHOD") {
	$title = "$node  --  ".$IN->{$node}->{"=TITLE"};
    }
    if ($IN->{$node}->{"=NODE_TYPE"} eq "=DISTR") {
	$title = "\@code{$node}  --  ".$IN->{$node}->{"=TITLE"}."\n";
	$title .= "\@findex unur_distr_$node\n";
    }

    # make menu
    my $menu;
    foreach my $n (@LIST_nodes_sorted) {
	next unless $IN->{$n}->{"=UP"} eq $node;
	$menu .= "* $n\:\: ".$IN->{$n}->{"=TITLE"}."\n";
    }
    if ($menu) {
	$menu = "\@menu\n".$menu."\@end menu\n\n";
    }

    # print header file name
    $TEXI .= "\@c -------------------------------------\n";
    $TEXI .= "\@c ".$IN->{$node}->{"=FILE"}."\n";
    $TEXI .= "\@c\n\n";

    # print node and section
    $TEXI .= "\@node $node\n";
    $TEXI .= "$section  $title\n\n";
    
    # print menu
    $TEXI .= $menu;

    # print REQUIRED, OPTIONAL, etc. of method
    if ($IN->{$node}->{"=NODE_TYPE"} eq "=METHOD") {
	if ($IN->{$node}->{"=REQUIRED"}) {
	    $TEXI .= "\@table \@i\n";
	    if ($IN->{$node}->{"=REQUIRED"}) {
		$TEXI .= "\@item Required:\n".$IN->{$node}->{"=REQUIRED"}."\n";
	    }
	    if ($IN->{$node}->{"=OPTIONAL"}) {
		$TEXI .= "\@item Optional:\n".$IN->{$node}->{"=OPTIONAL"}."\n";
	    }
	    if ($IN->{$node}->{"=SPEED"}) {
		$TEXI .= "\@item Speed:\n".$IN->{$node}->{"=SPEED"}."\n";
	    }
	    if ($IN->{$node}->{"=REF"}) {
		$TEXI .= "\@item reference:\n".$IN->{$node}->{"=REF"}."\n";
	    }
	    $TEXI .= "\@end table\n\@sp 1\n\n";
	}
    }

    # print PDF, domain, etc. distribution 
    if ($IN->{$node}->{"=NODE_TYPE"} eq "=DISTR") {
	if ($IN->{$node}->{"=PDF"} or $IN->{$node}->{"=PMF"}) {
	    $TEXI .= "\@table \@i\n";
	    if ($IN->{$node}->{"=PDF"}) {
		$TEXI .= "\@item PDF:\n".$IN->{$node}->{"=PDF"}."\n";
	    }
	    if ($IN->{$node}->{"=PMF"}) {
		$TEXI .= "\@item PMF:\n".$IN->{$node}->{"=PMF"}."\n";
	    }
	    if ($IN->{$node}->{"=CONST"}) {
		$TEXI .= "\@item constant:\n".$IN->{$node}->{"=CONST"}."\n";
	    }
	    if ($IN->{$node}->{"=CDF"}) {
		$TEXI .= "\@item CDF:\n".$IN->{$node}->{"=CDF"}."\n";
	    }
	    if ($IN->{$node}->{"=DOMAIN"}) {
		$TEXI .= "\@item domain:\n".$IN->{$node}->{"=DOMAIN"}."\n";
	    }
	    if ($IN->{$node}->{"=FPARAM"}) {
		$TEXI .= $IN->{$node}->{"=FPARAM"};
	    }
	    if ($IN->{$node}->{"=REF"}) {
		$TEXI .= "\@item reference:\n".$IN->{$node}->{"=REF"}."\n";
	    }
	    if ($IN->{$node}->{"=STDGEN"}) {
		$TEXI .= "\@item special generators:\n".$IN->{$node}->{"=STDGEN"}."\n";
	    }
	    $TEXI .= "\@end table\n\n";
	}
    }
	
    # print description
    $TEXI .= $IN->{$node}->{"=DESCRIPTION"};

    # print function reference
    if ($IN->{$node}->{"=ROUTINES"}) {
	$TEXI .= "\n\@subsubheading Function reference\n\n";
	$TEXI .= $IN->{$node}->{"=ROUTINES"}."\n\n";
    }

    # print `end of header file'
    $TEXI .= "\n\@c\n";
    $TEXI .= "\@c end of ".$IN->{$node}->{"=FILE"}."\n";
    $TEXI .= "\@c -------------------------------------\n";
	
    # search for all subnodes
    foreach my $n (@LIST_nodes_sorted) {
	next unless $IN->{$n}->{"=UP"} eq $node;
	
	# print all sub nodes
	texi_node($n,$level+1);
    }

} # end of texi_node() 

#############################################################
# scan reference
#

sub scan_REF {
    my $node_name = $_[0];   # name of node
    my $tag = $_[1];         # TAG (node section)

    # content of node
    my $entry = $IN->{$node_name}->{$tag};

    # empty ? 
    return unless $entry;

    # split into entries
    my @bibrefs = split /\]\s*\[/, $entry;

    # now print entries
    $IN->{$node_name}->{$tag} = "";
    foreach my $bre (@bibrefs) {
	$bre =~ s/[\[\]]//g;    # remove all square brackets
	$IN->{$node_name}->{$tag} .= transform_bibref($bre);
    }

} # end of scan_REF()

#############################################################
# texify string
#

sub texify_string {
    my $string = $_[0];

    # infinity --> oo 
    $string =~ s/(infinity)/\\infty/g;

    # replace '*' by space 
    $string =~ s/\*/\\, /g;

    # name of funktions
    $string =~ s/(exp|max|min|sqrt)/\\$1/g;

    # small greek letters
    $string =~ s/(alpha|beta|gamma|lambda|mu|nu|pi|phi|sigma|tau|theta|zeta)/\\$1/g;

    # capital greek letters
    $string =~ s/(Gamma)/\\$1/g;

    # <, <=, etc.
    $string =~ s/<=/\\leq/g;
    $string =~ s/>=/\\geq/g;

    # fractions
    $string =~ s/\\frac\{([^\}]+)\}\{([^\}]+)\}/\{$1\\over $2\}/g;

    return $string;

} # end of texify_string()



#############################################################
# scan domain for distribution
#

sub scan_DOMAIN {
    my $node_name = $_[0];   # name of node
    my $tag = $_[1];         # TAG (node section)

    # content of node
    my $entry = $IN->{$node_name}->{$tag};

    # empty ? 
    return unless $entry;

    # trim heading blanks
    $entry =~ s/^\s*//;

    # chop off trailing blanks
    $entry =~ s/\s+$//;

    # remove newlines
    $entry =~ s/\n+/ /g;

    # format tex output
    my $texentry = texify_string($entry);

    # return result
    $IN->{$node_name}->{$tag} = "\@iftex\n\@tex\n\$$texentry\$\n\@end tex\n\@end iftex\n";
    $IN->{$node_name}->{$tag} .= "\@ifnottex\n$entry\n\@end ifnottex\n";

} # end of scan_DOMAIN()

#############################################################
# scan PDF for distribution
#

sub scan_PDF {
    my $node_name = $_[0];   # name of node
    my $tag = $_[1];         # TAG (node section)
    my $texentry;            # content of node converted to TeX

    # content of node
    my $entry = $IN->{$node_name}->{$tag};

    # empty ? 
    return unless $entry;

    # trim heading blanks
    $entry =~ s/^\s*//;

    # chop off trailing blanks
    $entry =~ s/\s+$//;

    # remove newlines
    $entry =~ s/\n+/ /g;

    # text mode
    if ($entry =~ /\@text\{(.*)\}/) {
	$texentry = $1;
	$entry = $1;
    }

    else {
	# format tex output
	$texentry = "\$".texify_string($entry)."\$";

	# format other output
	$entry =~ s/\\over\s+/\//g;

	$entry =~ s/\\hbox{\s*(\w+)\s*}/ $1 /g;
	$entry =~ s/\\hfil+\\break/\n\n/g;

	$entry =~ s/\\frac\{([^\}]+[\s\+\-]+[^\}]+)\}\{([^\}]+[\s\+\-]+[^\}]+)\}/\($1\)\/\($2\)/g;
	$entry =~ s/\\frac\{([^\}]+)\}\{([^\}]+[\s\+\-]+[^\}]+)\}/$1\/\($2\)/g;
	$entry =~ s/\\frac\{([^\}]+[\s\+\-]+[^\}]+)\}\{([^\}]+)\}/\($1\)\/$2/g;
	$entry =~ s/\\frac\{([^\}]+)\}\{([^\}]+)\}/$1\/$2/g;

	$entry =~ s/\{/\(/g;
	$entry =~ s/\}/\)/g;
    }

    # return result
    $IN->{$node_name}->{$tag} = "\@iftex\n\@tex\n$texentry\n\@end tex\n\@end iftex\n";
    $IN->{$node_name}->{$tag} .= "\@ifnottex\n$entry\n\@end ifnottex\n";
	
} # end of scan_PDF()

#############################################################
# scan list of parameters for distribution
#

sub scan_FPARAM {
    my $node_name = $_[0];   # name of node
    my $tag = $_[1];         # TAG (node section)

    # content of node
    my $entry = $IN->{$node_name}->{$tag};

    # empty ? 
    return unless $entry;

    # trim heading blanks
    $entry =~ s/^\s*//;

    # chop off trailing blanks
    $entry =~ s/\s+$//;

    # remove blanks around `:'
    $entry =~ s/[ \t]*\:[ \t]*/\:/g;

    # split into lines
    my @lines = split /\n+/, $entry;

    # process lines
    my $out;
    my $texout; 
    my $flist;
    my $n_total = 0;
    my $n_optional = 0;
    my $opt_level;
    my $last_opt_level;

    foreach my $l (@lines) {
	# each line must start with a [\d+]
	next unless $l =~ /\[*\d+/;
	# split into columns
	my @cols = split /\:/, $l;
	die "\nwrong number of columns for =FPARAM: $#cols" if $#cols != 4;

	# get entries
	my $number  = $cols[0];
	$number     =~ s/(.*)(\d+).*/\[$2\]/;
	$opt_level  = $1;

	my $name    = $cols[1];
	my $cond    = $cols[2];
	my $default = $cols[3];
	my $type    = $cols[4];

	# append list of parameters
	if ($opt_level ne $last_opt_level) {
	    $last_opt_level = $opt_level;
	    $flist .= " [";
	}
	if ($n_total) {
	    $flist .= ", $name";
	}
	else {
	    $flist .= " $name";
	}

	# process
	$out    .= "\@item \@code{$number} \@tab $name \@tab $cond \@tab $default \@tab \@i{($type)}\n";

	$texout .= "\@item \@code{$number}\n";
	$texout .= "\@tab\@tex\$$name \$\@end tex\n";
	$texout .= "\@tab\@tex\$$cond \$\@end tex\n";
	$texout .= "\@tab $default\n";
	$texout .= "\@tab (\@i{$type})\n";
	
	++$n_total;
	++$n_optional if length $default;
    }
    my $n_required = $n_total - $n_optional;

    $last_opt_level =~ s/\[/ \]/g;
    $flist .= $last_opt_level;

    # make TeX output
    $texout =~ s/<=/\\leq/g;
    $texout =~ s/>=/\\geq/g;
    $texout =~ s/(alpha|beta|gamma|lambda|mu|nu|pi|phi|sigma|tau|theta|zeta)(\W)/\\$1$2/g;

    my $texout_header  = "\@iftex\n";
    $texout_header .= "\@item parameters $n_required ($n_total): \@r{$flist}\n\@sp 1\n";
    $texout_header .= "\@multitable {No.} {namex} {999999999999} {defaultx} {xxxxxxxxxxxxxxxxxxxxxxxx}\n";
    $texout_header .= "\@item No. \@tab name \@tab \@tab default\n";

    $IN->{$node_name}->{$tag} = $texout_header.$texout."\@end multitable\n\@end iftex\n";

    # make other output
    my $out_header  = "\@ifnottex\n";
    $out_header .= "\@item parameters $n_required ($n_total): $flist\n";
    $out_header .= "\@multitable {No.xx} {namexxx} {99999999999} {defaultx} {xxxxxxxxxxxxxxxxxxxxxxxx}\n";
    $out_header .= "\@item No. \@tab name \@tab \@tab default\n";

    $IN->{$node_name}->{$tag} .= $out_header.$out."\@end multitable\n\@end ifnottex\n";

} # end of scan_FPARAM()

#############################################################
# scan list of parameters for distribution
#

sub scan_STDGEN {
    my $node_name = $_[0];   # name of node
    my $tag = $_[1];         # TAG (node section)

    # content of node
    my $entry = $IN->{$node_name}->{$tag};

    # empty ? 
    return unless $entry;

    # trim heading blanks
    $entry =~ s/^\s*//;

    # chop off trailing blanks
    $entry =~ s/\s+$//;

    # split into lines
    my @lines = split /\n+/, $entry;

    my $out = "\@table \@code\n";

    # process lines
    foreach my $l (@lines) {

	# split into indicator and description
	(my $id, my $body) = split /\s+/, $l, 2;
	
	if ($body =~ /\[(\w+)\]/) {
	    # there are references
	    my @references;
	    while ($body =~ /\[(\w+)\]/) {
		$body =~ s /\[(\w+)\]\s*//;
		push @references, $1;
	    }
	    $out .= "\@item $id\n";

	    # tex output
	    my $texbody = texify_string($body);
	    $texbody =~ s/(\@tex)/\n$1/g;
	    $texbody =~ s/(\@end\s+tex)/\n$1\n/g;
	    $out .= "\@iftex\n";
	    $out .= $texbody;
	    foreach my $r (@references) {
		$out .= " [$r]";
	    }
	    $out .= "\n\@end iftex\n";

	    # remove tex marks
	    $body =~ s/\$//g;
	    $body =~ s/\@tex//g;
	    $body =~ s/\@end\s+tex//g;

	    $out .= "\@ifhtml\n";
	    $out .= "$body ";
	    foreach my $r (@references) {
		$out .= "\@ref{bib:$r, [$r]} ";
	    }
	    $out .= "\n\@end ifhtml\n";

	    $out .= "\@ifinfo\n";
	    $out .= "$body ";
	    foreach my $r (@references) {
		$out .= "[$r] ";
	    }
	    $out .= "\n\@end ifinfo\n";
	}
	
	else {
	    # there is no reference
	    $out .= "\@item $id\n";
	    $out .= "$body\n";
	}
    }

    $out .= "\@end table\n";
    
    # make other output
    $IN->{$node_name}->{$tag} = $out;

} # end of scan_STDGEN()

#############################################################
# chop off trailing blanks
#

sub scan_chop_blanks {
    my $node_name = $_[0];   # name of node
    my $tag = $_[1];         # TAG (node section)

    # content of node
    my $entry = $IN->{$node_name}->{$tag};

    # trim heading blanks
    $entry =~ s/^\s*//;

    # chop off trailing blanks
    $entry =~ s/\s+$//;

    # remove newlines
    $entry =~ s/\n+/ /g;

    # return result
    $IN->{$node_name}->{$tag} = $entry;

} # end of scan_chop_blanks()

#############################################################
# dummy routine
#

sub scan_do_nothing {
    my $node_name = $_[0];   # name of node
    my $tag = $_[1];         # TAG (node section)

    # nothing to do
    return;
} # end of scan_do_nothing()

#############################################################
# scan TAG (node section) =UP
#

sub scan_UP {
    my $node_name = $_[0];   # name of node
    my $tag = $_[1];         # TAG (node section);  not used

    # content of node
    my $entry = $IN->{$node_name}->{"=UP"};
    $entry =~ s/^\s*//;  # trim heading blanks

    # we have two entries which are separated by blanks
    # the first entry is the name of the UP node
    # the second entry is used to order within this upper node
    (my $upper_name, my $order) = split /\s+/, $entry, 2;

    # store upper node
    $IN->{$node_name}->{"=UP"} = $upper_name;

    # extract string for ordering (lexicographic ordering is used)
    if ($order =~ /\[(\w+)\]/) {
	$IN->{$node_name}->{"=ORDERING"} = $1;
    }
    else {  # use node name
	$IN->{$node_name}->{"=ORDERING"} = $node_name;
    }
	
} # end if scan_UP()

#############################################################
# scan TAG (node section) =ROUTINES
#

sub scan_ROUTINES {
    my $node_name = $_[0];   # name of node
    my $tag = $_[1];         # TAG (node section)

    # valid C data types 
    my @C_TYPES = ( "UNUR_PAR",
		    "UNUR_GEN",
		    "UNUR_DISTR",
		    "UNUR_URNG",
		    "UNUR_FUNCT_CONT",
		    "UNUR_FUNCT_CVEC",
		    "UNUR_VFUNCT_CVEC",
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

    # content of node
    my $entry = $IN->{$node_name}->{"=ROUTINES"};

    # trim heading blanks
    $entry =~ s/^\s*//;
    
    # remove double blank lines
    $entry =~ s/\n\s*\n\s*\n/\n\n/g;

    # simplify comment markers: /** --> /*
    $entry =~ s/\/\*+/\/\*/g;
    $entry =~ s/\*+\//\*\//g;

    # remove blanks and new lines inside comments at markers
    $entry =~ s/\/\*\s+/\/\*/g;
    $entry =~ s/\s+\*\//\*\//g;

    # remove blank lines at begin and end of entry
    $entry =~ s/^(\s*\n)+//g;
    $entry =~ s/(\s*\n)+$//g;

    # if there are blocks of comments, separated by empty lines
    # then delete all but the first block of comments
##  $entry =~ s/\*\/\s*\n\s*\n\s*\/\*[.\n]*\*\//\*\//g;

    # split into blocks
    my @blocks = split /\*\/\s*\n\s*\n/, $entry ;

    # store processed text
    my $proc = '';

    # deftypefn block closed
    my $defblock_open = 0;
    
    # process blocks
    my $fkt_block = '';
    foreach my $block (@blocks) {

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
	    next;  # next block
	}

	# split into function declaration and its description
	(my $fkt, my $body) = split /\;/, $block, 2; 

	# check function type
	my $type_ok = 0;
	foreach my $type (@C_TYPES) {
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
	die "Function defined twice: $fn_name" if $LIST_routines->{$fn_name};

	# store in list of all routines
	$LIST_routines->{$fn_name} = 1;

	# write entry
	my $first = 1;
	if (@argslist) {
	    # this is a function with arguments
	    # make anchor
	    $fkt_block .= "\@anchor{funct:$fn_name}\n";
	    # make texinfo tag
	    $fkt_block .= (($defblock_open) ? "\@deftypefnx" : "\@deftypefn");
	    $fkt_block .= " %%%Function%%% \{$fn_type\} $fn_name (";
	    foreach my $arg (@argslist) {
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
    $IN->{$node_name}->{"=ROUTINES"} = $proc;

    return;

} # end of scan_ROUTINES()

#############################################################
# format texinfo output for =NODE
#

sub texi_NODE {
    my $node = $_[0];

} # end of texi_TOP()

#############################################################
# transform special strings
#

sub transform_special_strings {
    my $line = $_[0];

    # trim blanks
    $$line =~ s/[ \t\r\f]+\n/\n/g;

    # @cc_start --> /*
    # @cc_stop  --> */
    $$line =~ s/\@cc_start/\/*/g;
    $$line =~ s/\@cc_stop/*\//g;

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

} # end of transform_special_strings()

#############################################################

sub transform_bibref {
    my $entry = $_[0];   # entry to be transformed

    # empty ? 
    return unless $entry;

    # remove newlines
    $entry =~ s/\n+/ /g;

    # trim heading blanks
    $entry =~ s/^\s*//;

    # chop off trailing blanks
    $entry =~ s/\s+$//;

    # split into ref and optional remark
    # which is separated by a colon
    (my $anchor, my $remark) = split /[\:\,]\s*/, $entry, 2;
    if ($remark) { 
	$remark =~ s/\,/\;/g;   # we cannot use a comma here
	$remark = ": $remark"; 
    }
    
    # output 
    my $entrywithlink = "\@ref{bib:$anchor,, [$anchor$remark]}";
    $entry = "[$anchor$remark]";

    # output
    my $output = 
	"\@ifhtml\n$entrywithlink\n\@end ifhtml\n"
	."\@ifnothtml\n$entry\n\@end ifnothtml\n";

    return $output;

} # end of transform_bibref()

#############################################################

