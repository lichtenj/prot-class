#!/usr/bin/perl 
#
# File	  : make.pl -- Make utility for Ohio University Bioinformatics Framework
# Author  : Ohio University
# Created : July 2006
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#

use strict;
my $dir = '..\bin';
my $log = $dir . '\BioHio.Build.log';

# I've decided to force a "make all" just for simplicity
# you can always modify this script to just make the
# applications you're interested in.
my $make_all = 1;

open LOG, ">$log" or die $!;

print LOG 'Build Log for Bio Ohio'."\n";

# create the main GUI executables and main experiment exes
#
print "** building pmf ...\n";
print LOG `perlapp --force --exe $dir\\pmf.exe pmf_app.pl`, "\n";
print LOG `perlapp --force --exe $dir\\pmf_main.exe pmf_main.pl`, "\n";
print "**\n";

if ( $make_all )
{
    print "** building BioHio ...\n";
    print LOG `perlapp --force --exe $dir\\BioHio.exe main.pl`, "\n";
    print "**\n";

    print "** building hrgp...\n"; 
    print LOG `perlapp --force --exe $dir\\hrgp.exe hrgp_app.pl`, "\n";
    print LOG `perlapp --force --exe $dir\\hrgp_main.exe hrgp_main.pl`, "\n";
    print "**\n";
    
    print "** building gps ... \n";
    print LOG `perlapp --force --exe $dir\\gps.exe gps_app.pl`, "\n";
    print LOG `perlapp --force --exe $dir\\gps_main.exe gps_main.pl`, "\n";
    print "**\n";
    
    print "** building PromoterPredict ... \n";
    print LOG `perlapp --force --exe $dir\\PromoterPredict.exe PromoterPredict_app.pl`, "\n";
    print LOG `perlapp --force --exe $dir\\PromoterPredict_main.exe PromoterPredict_main.pl`, "\n";
    print "**\n";
    
    print "** building test executables ... \n";
    print LOG `perlapp --force --exe $dir\\MakeTests.exe maketests.pl`, "\n";
    print LOG `perlapp --force --exe $dir\\RunTests.exe runtests.pl`, "\n";
    print "**\n";
    
    # make a version in acceptable form from the date and time
    # make test version slightly different (one greater in seconds)
    # so that the test version can be installed overtop of a non-test version
    my ($version);
    {
        use POSIX;
        # need to keep seconds to one digit to fit into version string
        # version string is max:255.255.65535.65535, so the best we can
        # do is a new version every 10 seconds... should be enough
        my $S = int(POSIX::strftime( '%S' => localtime)/10);
        $version = POSIX::strftime( '%y.%m.%d.%H%M' => localtime). $S;
    }

    print "** building installer version $version ... \n";
    
    `"C:\\Program Files\\Caphyon\\Advanced Installer\\AdvancedInstaller.com" /edit hrgp.aip /SetVersion $version`;
    `"C:\\Program Files\\Caphyon\\Advanced Installer\\AdvancedInstaller.exe" /build hrgp.aip`;
    
    print "**\n";
    
    print "** building test cases ...\n";
    `$dir\\MakeTests.exe`;
}

print "**\n** Done\n";

