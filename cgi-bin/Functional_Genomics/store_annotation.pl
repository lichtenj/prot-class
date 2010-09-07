#!/usr/bin/perl -w

use strict;

use CGI;
my $cgi=new CGI;

print $cgi->header;
print $cgi->start_html('Annotation');

print $cgi->h1('Annotation Change');

# my $query = CGI->new;
my $annotation = $cgi->param('annotation');
my $source = $cgi->param('source');
my $file = $cgi->param('file');

$file =~ /(.+)\/.+\.annot$/;
my $dir = $1;

system("mkdir $dir");

open(OUT, ">$file") or print $cgi->p("Cannot open");
print OUT $annotation;
close OUT;

print $cgi->p("Annotation changed to:");
print $cgi->p($annotation);
print $cgi->p($file);

print '<a href="../../'.$source.'/">Return to Table</a>';

print $cgi->end_html;
exit (0);

