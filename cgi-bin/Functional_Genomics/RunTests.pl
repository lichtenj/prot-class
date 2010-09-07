#!/usr/bin/perl 
#
# File	  : maketests.pl -- Make test files for Ohio University Bioinformatics Framework
# Author  : Ohio University
# Created : July 2006
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#
#     Run all the manifest files found in the specified directory
#
use strict;

my $spec = '..\\tests\*.man';
my $exef = '..\\bin\pmf_main.exe';

print "Testing exe: $exef\n";
print "Running tests: $spec\n";

foreach my $file (glob($spec))
{
	print "running $file\n";
	system("$exef $file");
}
