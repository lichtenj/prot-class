#!/usr/bin/perl -w

use strict;

use CGI;
my $cgi=new CGI;

# my $query = CGI->new;
my $userid = $cgi->param('userid');
my $source = $cgi->param('source');

print $cgi->header;
print $cgi->start_html('Annotation');
print $cgi->h1('Annotation Archival');

system('tar -cf /var/www/'.$source.'.tar /var/www/'.$source.'/*');

print '<a href="../../'.$source.'.tar">Export Annotations</a><BR><BR>';
print '<a href="../../'.$source.'/">Return to Table</a>';

print $cgi->end_html;
exit (0);

