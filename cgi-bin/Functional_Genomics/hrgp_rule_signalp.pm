#!/usr/bin/perl -w
#
# File	  : hrgp_rule_signalp.pm
# Author  : Ohio University CNS
# Created : July 2005
# Version : <none>
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.

# First create a module derived from "manifest"
package hrgp;

use lib '/Perl/site/lib/';
use lib '../lib/';

use strict;
use Bio::Tools::Sigcleave;
#use AminoAcid 'GetHydropathy';


#######################################################################
# Function: rule_Signal_P
#
# Purpose: The rule for predicting signal peptide locations. The function
#          currently uses the BioPerl::Tools::sigcleave function, which
#          predicts protein signals. The default threshold value of 3.5
#          is estimated to correctly identify 95% of signal peptides and
#          correctly reject 95% of non-signal peptides. The cleavage site
#          is believed to be correctly predicted 75-80% of the time. The
#          number of candidate cleavage sites is saved in the $result hash.
#          The pretty_print string, which shows all predicted cleavage
#          sites, including score, length, and start, are also stored.
#
#          We need to investigate a way to output the highest scored
#          cleavage site, which is the first entry in the $result->
#          {SignalP}->{pretty_print} string.
#
# Param 1:
# Param 2:
#
# Return value: True if the sequence has a predicted signal peptide
#               location; false otherwise.
#######################################################################
sub rule_Signal_P
{
	my $this = shift or die;
	my $seq = shift or die;
	my $result = shift or die;
	my $threshold = $this->{config}->{SignalP}->{threshold};

	# get hash reference to the AgPeptide section of the config hash
	my $hr = $this->{config}->{SignalP};

	my $seq_str = substr($seq->seq, $hr->{'start_position'}, $hr->{'search_length'}); 

	# create the object that detects cleavages
	my $cleaver = new Bio::Tools::Sigcleave( -seq => $seq_str, -threshold => $threshold, -type => 'AMINO' );

	# save the cleaver in results if the count > 0
	$result->{SignalP}->{cleaver} = ($cleaver->result_count > 0) ? $cleaver : undef;
	$result->{SignalP}->{count}   = $cleaver->result_count;
	
	return defined $result->{SignalP}->{cleaver};
}

return 1;

