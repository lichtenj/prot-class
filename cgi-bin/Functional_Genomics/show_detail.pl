#!/usr/bin/perl -w
use strict;
use CGI;

my $cgi=new CGI;

my $seq = $cgi->param('seq');

print $cgi->header;
print $cgi->start_html('Upload');

print $cgi->h1("Details");

print $seq;

print $cgi->end_html;
exit (0);

