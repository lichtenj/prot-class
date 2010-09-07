#!/usr/bin/perl -w

use strict;

use CGI;
my $cgi=new CGI;

# my $query = CGI->new;
my $userid = $cgi->param('userid');
my $search = $cgi->param('searchterm');

print $cgi->header;
print $cgi->start_html('Annotation');
print $cgi->h1('Relevant Annotations');

#print $userid.'<BR>';
#print $search.'<BR>';
my $command = 'ls /var/www/ProtClass_Annotations/'.$userid.'* | grep "'.$search.'" /var/www/ProtClass_Annotations/*.annot';

my $out = `$command`;

$out =~ /$userid\_([\d\w\_\.\|]+)\.annot\:(.+)/;

#print $out;
#print '<BR>';
print '<b>'.$1.'</b>';
print '<BR>';
print $2;
print '<BR><BR>';

print $cgi->end_html;
exit (0);

