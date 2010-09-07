#!/usr/bin/perl -w
#
# File	  : gps_main.pl -- Gravity Persistant Signals main progam
# Author  : Ohio University 
# Created : July 2005
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#
use strict;
use gps;

my ($this,$class) = new gps();	# create GPS object 
$this->GetDefaults($ARGV[0]);	# load parameters from filename
my $result = $this->main();	    # run the experiment's logic
