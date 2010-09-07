#!/usr/bin/perl -w
#
# File	  : gps_app.pl -- Main for GPS application GUI
# Author  : Ohio University
# Created : July 2005
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#
use strict;

#############################################################################
# subclass the application class to separate out specific functionality for
# this specific main app
package main_app;
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
  use gps;
  return (['gps']);
}

# override this virtual function to create the new experiment menu specific to
# the modules included in sub Modules()
sub ModulesMenu
{
  my $this = shift or die;
  return
  [
    [qw/command/, 'Gravity Persistent Signals Analysis (GPS)', -underline=>0, -command => sub{$this->FileNew('gps');}],
  ];
}
############################################################################
package main;

my ($app) = new main_app;
$app->main();
