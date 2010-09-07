#!/usr/bin/perl -w

use strict;
use Archive::Tar;
use CGI;

my $cgi=new CGI;

# my $query = CGI->new;
my $file = $cgi->param('file');

print $cgi->header;
print $cgi->start_html('Upload');
#print $cgi->h1('Annotation Upload');

my $upload_dir = "/var/www/";
$file =~ /(.+)\.tar/;
my $subdir = $1;

open(UPLOADFILE, ">$upload_dir/$file" ) or die "$!";
binmode UPLOADFILE;
while (<$file>)
{
	print UPLOADFILE;
}
close UPLOADFILE;

#system('pwd');

my $tar = Archive::Tar->new;
$tar->read('/var/www/'.$file);
$tar->extract();

system('mv /var/www/cgi-bin/Functional_Genomics/var/www/'.$subdir.'/ /var/www/'.$subdir);

print $cgi->h1("Upload sucessful");
print 'Click here to see the <a href="/'.$subdir.'">results</a>';
#print $cgi->p($file);

print $cgi->end_html;
exit (0);

