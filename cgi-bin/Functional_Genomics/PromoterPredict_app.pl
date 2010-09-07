#!/usr/bin/perl -w
#
# File	  : main.pl -- Main for Promoter Prediction application GUI
# Author  : Ohio University
# Created : July 2006
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#

use strict;

package pp_app;
use application;

use vars qw(@ISA);
@ISA = qw(application);

sub Modules
{
  use PromoterPredict;
  return (['PromoterPredict']);
}

sub ModulesMenu
{
  my $this = shift or die;
  return
  [
    [qw/command/, 'Promoter Prediction Experiment Module (PromoterPredict)', -underline=>0, -command => sub{$this->FileNew('PromoterPredict');}],
  ];
}
############################################################################
package main;

my ($app) = new pp_app(@_);
$app->main();

sub ExePath
{
  return $app->ExePath();
}
