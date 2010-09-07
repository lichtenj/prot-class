#!/usr/bin/perl -w
#
# File	  : hrgp_rule_gpi.pm
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
use AminoAcid;

use Text::CSV;

# This-module specific requirements
use Bio::SeqIO;
use Bio::Tools::Sigcleave;
use Bio::Seq::RichSeq;


#######################################################################
# Function: rule_GPI_Anchor
#
# Purpose: The rule for predicting GPI anchor locations.
#
# Param 1:
# Param 2:
#
# Return value: True if the sequence has a predicted GPI anchor
#               location; false otherwise.
#######################################################################
sub rule_GPI_Anchor
{
	our $testno = 0 unless $testno; $testno++;
	my $this = shift or die;
	my $seq = shift or die;
	my $result = shift or die;
	my $len = length $seq->seq;
	my $id = $seq->display_id;

	# IDEALLY, this should be performed in the main while loop so each
	# function doesn't have to worry about it -- but how can we modify
	# a SeqIO object's sequence?
	# ALSO, how can we handle *'s in the sequence to mean "any amino acid" --
	# should we disallow databases of this type?
	
	# now see if the sequence ends with a '*'
	# by simply searching for the *
	# if it does, update $len by decerementing it by 1
	# otherwise, $len is the true length of the sequence
	my $numAsterisks = ($seq->seq =~ /\*/);
	$numAsterisks = 0 unless $numAsterisks;
	if ($numAsterisks > 0) {
		#print "this sequence has $numAsterisks *'s\n";
		$len--;
	}

	my $hr = $this->{config}->{GPI};
	my %bestMatch;
	$bestMatch{score} = 0;
	
	#our $gpicsv = new Text::CSV unless $gpicsv;
	#open CSVGPI, ">gpi_result.csv";
	
	#print "\nAnalyzing $id for GPI anchor...\n\n";
	my $minus = substr($seq->seq,length($seq->seq())-1,1) eq '*' ? 1 : 0;	
	#print "minus is $minus\n";
	
	# calculate the range for location of the omega sites
	# minus 3 for length of possible first omega site
	# minus 1 for trailing asterisk
	my $w_start = $len - ($hr->{max_phobic_region_length}+$hr->{max_philic_region_length}) - 3 - $minus;
	my $w_stop  = $len - ($hr->{min_phobic_region_length}+$hr->{min_philic_region_length}) - 3 - $minus;

	my $gpi = substr($seq->seq(), $w_start);
	my $sitelen = (length $gpi) - $minus;
	
	#print "Whole GPI anchor test site: $sitelen\n";
	#print $gpi,"\n";
	
	# find every possible omega site in the range
	my @w_sites;
	my $potentials = 0;
	#while ( $str =~ /([GASTVRNBQZLKMY][A-Z][GASCDNVLRH])/cgi )
	my $try = 0;

	my @gpi = split //, $gpi;	# array rep. of the whole GPI site
	my $siz = @gpi;				# size of whole GPI site
	for my $w ( 0..$siz)
	{
		my $tail_len = $siz - ($w) - 3 - $minus;
		my $tail_min = $hr->{min_phobic_region_length}+$hr->{min_philic_region_length};
		last if ( $tail_len < $tail_min ); # we've reached the smallest region allowed

		my $omega = "$gpi[$w]$gpi[$w+1]$gpi[$w+2]";

		#print "getting omega prob for $omega: ";
		my $prob = AminoAcid::GetOmegaProb($omega); # get the probability that this IS an omega site
		#print "prob = $prob\n";
		
		# if computed probability exceeds the threshold specified by the user,
		# proceed to compute hydrophocity regions
		if ( $prob && $prob >= $hr->{omega_site_threshold} )
		{
			#print "$w, $omega, $tail_len, $tail_min\n";
			for my $phil ( $hr->{min_philic_region_length} .. $hr->{max_philic_region_length} )
			{
				my $phob = $tail_len - $phil;
				if ( ($phob >= $hr->{min_phobic_region_length}) && ($phob <= $hr->{max_phobic_region_length}) )
				{
					my $cursor = $w+3;
					my ($phil_str,$phob_str) = ('','');
					my ($phil_hyd,$phob_hyd) = (0.0,0.0);;
					
					for (0..$phil-1)
					{
						$phil_str .= $gpi[$cursor];
						$phil_hyd += AminoAcid::GetHydropathy($gpi[$cursor]);
						$cursor++;
					}
					$phil_hyd /= ($phil);
					$phil_hyd = sprintf("%.3f",$phil_hyd);
					
					for (0..$phob-1)
					{
						$phob_str .= $gpi[$cursor];
						$phob_hyd += AminoAcid::GetHydropathy($gpi[$cursor]);
						$cursor++;
					}
					$phob_hyd /= ($phob);
					$phob_hyd = sprintf("%.3f",$phob_hyd);
					
					my $diff = sprintf("%.3f", $phob_hyd - $phil_hyd); # difference of average hydophobicity between the two regions
					# if the difference in average hydrophoobicity between the two regions exceeds
					# the specified threshold, we have a potential match
					if ( $diff >= $hr->{hydrophobicity_threshold} )
					{
						if ( ($diff +$prob) > $bestMatch{score} ) {
							$bestMatch{score} = $diff + $prob;
							$bestMatch{wStart} = $w_start + $w;
							$bestMatch{wRelative} = $w;
							$bestMatch{wSite} = $omega;
							$bestMatch{entireSite} = $gpi;
							$bestMatch{boundaryLeft} = $phob + $w +$w_start + 2;
							$bestMatch{boundaryRight} = $bestMatch{boundaryLeft} + 1;
							$bestMatch{philicSide} = $phil_str;
							$bestMatch{philicSideVal} = $phil_hyd;
							$bestMatch{phobicSide} = $phob_str;
							$bestMatch{phobicSideVal} = $phob_hyd;
							$bestMatch{prob} = $prob;
							$bestMatch{diff} = $diff;
							
							$potentials++;
							
							my @record = ($testno, $omega, $prob, $diff, "($phil_hyd -- $phil_str)", "($phob_hyd -- $phob_str)" );
							#$gpicsv->combine(@record);
						
							#print CSVGPI $gpicsv->string(),"\n";

							#printf "%4s, %6.3f %6.3f ($phil_hyd -- $phil_str), ($phob_hyd -- $phob_str)\n", $omega, $prob, $diff;
							#print "split: $phil, $phob, $diff, ($phil_hyd -- $phil_str), ($phob_hyd -- $phob_str)\n";
						}
					}
				}
			}
		}
		next;
	}
	
	if ($bestMatch{wStart}) {
		#print "\nThere were a total of $potentials potential sites\n";
		#print "   Best site had omega site probability $bestMatch{prob}, hydrophobic difference $bestMatch{diff}\n";
		#print "      entire GPI anchor site: $gpi\n";
		#print "      w site starts at position $bestMatch{wStart} in original sequence, site = $bestMatch{wSite}\n";
		#print "      hydrophobic barrier = $bestMatch{boundaryLeft} - $bestMatch{boundaryRight}\n\n";
		#print "      hydrophilic side = $bestMatch{philicSide} ($bestMatch{philicSideVal}), hydrophobicSide = $bestMatch{phobicSide} ($bestMatch{phobicSideVal})\n";
		
		$result->{GPI}->{count} = $potentials;
		$result->{GPI}->{matches} = \%bestMatch;
	}
}


return 1;

