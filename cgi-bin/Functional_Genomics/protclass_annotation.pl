#!/usr/bin/perl -w

use strict;

use CGI;
my $cgi=new CGI;

# my $query = CGI->new;
my $names = $cgi->param('file');
my $source = $cgi->param('source');

print $cgi->header;
print $cgi->start_html('Annotation');
print $cgi->h1('Annotation');
print $cgi->p('File: '.$names);
print $cgi->p('Time: '.$source);

open(IN, '/var/www/'.$source.'/ProtClass_Annotations/'.$names) or print $cgi->p("No Annotation File exists at the momemnt");
my @test = <IN>;
my $text = join("",@test);
print $cgi->p($text);

print '<form action="./store_annotation.pl" method="POST">
<input name=file type="hidden" value="/var/www/'.$source.'/ProtClass_Annotations/'.$names.'"/>
<input name=source type="hidden" value="'.$source.'"/>
<textarea name="annotation" cols=40 rows=4>'.$text.'</textarea><BR><input type="submit" value="Publish"></form>';

print $cgi->end_html;
exit (0);

