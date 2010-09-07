#!/usr/bin/perl -w

#

# File	  : hrgp_app.pl -- Main for HRGP application GUI

# Author  : Ohio University

# Created : June 2006

#

#     Copyright (c) 2006  Ohio University. All rights reserved.

#     This program is free software; you can redistribute it and/or

#     modify it under the same terms as Perl itself.

#

use strict;



#############################################################################

# subclass the application class to separate out specific functionality for

# this specific main app

package hrgp_app;

use application;



use vars qw(@ISA);

@ISA = qw(application);



# override this virtual function to include all the experiment type

# modules you want to "compile" into the exe program

sub Modules

{

  # use all the modules needed by this application and

  # return a list of modules used (so the caller can know

  # which default manifests to load, for instance)

  use hrgp;

  return (['hrgp']);

}



# override this virtual function to create the new experiment menu specific to

# the modules included in sub Modules()

sub ModulesMenu

{

  my $this = shift or die;

  return

  [

    [qw/command/, 'Hydroxyproline-rich Glycoprotein Analysis (HRGP)', -underline=>0, -command => sub{$this->FileNew('hrgp');}],

  ];

}

############################################################################

package main;



my ($app) = new hrgp_app;

$app->main();

