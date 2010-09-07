#!/usr/bin/perl -w
#
# File	  : hrgp_rule_classify.pm
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


#######################################################################
# Function: firstClassify
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub firstClassify {
 	my $this 	= shift or die;		# object
	my $seq 	= shift or die;		# current SeqIO
	my $result 	= shift or die;		# results duh
	my $hr		= $this->{config};	# easier typing

	# make sure that the result hash of class types
	# is initialized to '-' for the display
	$result->{has_biased_AA_chars}
		= $result->{has_ag_peptide_chars}
		= $result->{has_fasciclin_chars}
		= $result->{has_lrr_chars}
		= $result->{has_extensin_chars}
		= $result->{has_prp_chars}
		= $result->{passed_ad_hoc}
		= $result->{classification}
		= '-';
	
	# now proceed to check characterizations in order
	# to classify the sequence
	
	# Does this protein have any AGP characteristics?
	if ( ($result->{BiasedAA}->{percent} >= $hr->{BiasedAA}->{threshold}) )
	{
		$result->{has_biased_AA_chars} = 'yes';
	}
	
	if ( ($result->{AgPeptide}->{percent} >= $hr->{AgPeptide}->{threshold}) )
	{
		$result->{has_ag_peptide_chars} = 'yes';
	}
	
	if ( ($result->{Fasciclin}->{count}) && ($result->{Fasciclin}->{count} > 0) )
	{
		$result->{has_fasciclin_chars} = 'yes';
	}

	
	# Does this protein have LRRs? (yes if any rule matched, i.e. count > 0)
	{
		if ( $result->{LRR}->{count}  )
		{
			$result->{has_lrr_chars} = 'yes';
		}
	}

	# Does this protein have extensin characteristics? (yes if any rule matched, i.e. count > 0)
	{
		if ( $result->{Extensin}->{count}  )
		{
			$result->{has_extensin_chars} = 'yes';
		}
	}

	# Does this protein have PRP characteristics? (yes if any rule matched, i.e. count > 0)
	{	
		if ( $result->{PRP}->{count}  )
		{
			$result->{has_prp_chars} = 'yes';
		}
	}
	
	# we don't do any classifications for the ad hoc rules, but still
	# mark whether or not any of them passed
	{	
		#my $rule_count = keys %{$result->{PRP}->{pattern}};
		if ( $result->{AdHoc}->{count}  )
		{
			$result->{passed_ad_hoc} = 'yes';
		}
	}
		
}


my $classifications =
{
	"Classical AGP" => "(BAA composition and a signal peptide)",
	"AG-Peptide" => "(BAA composition and a signal peptide)",
	"Fasciclin" => "(fasciclin domain(s) and a signal peptide)",
	"LRR-extensin" => "(extensin motifs, LRRs, and a signal peptide)",
	"Extensin" => "(extensin motifs, no LRRs, and a signal peptide)",
	"chimeric" => "(AGP characteristics and extensin motifs)",
};

#######################################################################
# Function: finalShowalterClassify
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub finalShowalterClassify
{
 	my $this 	= shift or die;		# object
	my $seq 	= shift or die;		# current SeqIO
	my $result 	= shift or die;		# results duh
	my $hr		= $this->{config};	# easier typing
	my $counts  = shift or die;

	# these are boolean flags used to keep track of what we classify
	# a protein as (used for determining if a protein may be chimeric)
	my $looksLike =
	{
		classicalAGP => 0,
		agPeptideAGP => 0,
		fasciclinAGP => 0,
		extensin => 0,
	};

	# determine if it's likely a biased-amino acid protein (Classical AGP)
	if ($result->{has_biased_AA_chars} eq 'yes')
	{
		if ($result->{SignalP}->{count} > 0) {
			$result->{classification} = 'Classical AGP';
			$counts->{likelyClassicals}++;
			$looksLike->{classicalAGP} = 1;
		} else {
			#print "  UNLIKELY AGP (has BAA composition but no signal peptide)\n";
		}
	}
	
	# determine if it's likely an AG-peptide
	if ($result->{has_ag_peptide_chars} eq 'yes') {
		if ($result->{SignalP}->{count} > 0) {
			$result->{classification} = 'AG-Peptide';
			$counts->{likelyAgPeptides}++;
			$looksLike->{agPeptideAGP} = 1;
		} else {
			#print "  UNLIKELY AGP (has BAA composition but no signal peptide)\n";
		}
	}
	
	# determine if it's likely a fasciclin protein (Fasciclin AGP)
	if ($result->{has_fasciclin_chars} eq 'yes') {
		if ($result->{SignalP}->{count} > 0) {
			$result->{classification} = 'Fasciclin';
			$counts->{likelyFasciclins}++;
			$looksLike->{fasciclinAGP} = 1;
		} else {
			#print "  UNLIKELY FASCICLIN AGP (has fasciclin domain(s) but no signal peptide)\n";
		}
	}
	
	# determine if it's likely an extensin
	if ($result->{has_extensin_chars} eq 'yes') {
		if ($result->{SignalP}->{count} > 0) {
			if ($result->{has_lrr_chars} eq 'yes') {
				$result->{classification} = 'LRR-extensin';
				$counts->{likelyLrrExtensins}++;
				$looksLike->{extensin} = 1;
			} else {
				# regular extensin, or at least non-LRR extensin
				$result->{classification} = 'Extensin';
				$counts->{likelyOtherExtensins}++;
				$looksLike->{extensin} = 1;
			}
		} else {
			#print "  UNLIKELY EXTENSIN (has extensin motifs but no signal peptide)\n";
		}
	}

	if ( ($looksLike->{classicalAGP} && $looksLike->{extensin}) ||
	     ($looksLike->{agPeptideAGP} && $looksLike->{extensin})    ||
		 ($looksLike->{fasciclinAGP} && $looksLike->{extensin}) ) {
			$result->{classification} = 'chimeric';
			$counts->{possiblyChimeric}++;
	}
}

return 1;
