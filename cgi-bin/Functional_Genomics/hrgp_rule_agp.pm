#!/usr/bin/perl -w
#
# File	  : hrgp_rule_agp.pm
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
# Function: rule_Biased_AA
#
# Purpose: The biased-amino-acid composition rule for AGPs. Determines 
#          if the passed sequence meets the percentage threshold for
#          biased amino acid composition and also the sequence length
#          restrictions.The percentage for each sequence (i.e., it meets
#          the length restrictions) is stored in the $result hash.
#
# Default values: threshold, 50%; no length restrictions
#
# Param 1:
# Param 2:
#
# Return value: True if the sequence meets the length restrictions and
#               biased AA composition requirements; false otherwise.
#######################################################################
sub rule_Biased_AA
{
	my $this = shift or die;
	my $seq = shift or die;
	my $result = shift or die;
	my $len = length($seq->seq);

	#
	# IDEALLY, this should be performed in the main while loop so each
	# function doesn't have to worry about it -- but how can we modify
	# a SeqIO object's sequence?
	# ALSO, how can we handle *'s in the sequence to mean "any amino acid" --
	# should we disallow databases of this type?
	#
	# now see if the sequence ends with a '*'
	# by simply searching for the *
	# if it does, update $len by decerementing it by 1
	# otherwise, $len is the true length of the sequence
	#
	#my $numAsterisks = ($seq->seq =~ /\*/);
	#$numAsterisks = 0 unless $numAsterisks;
	#if ($numAsterisks > 0) {
	#	#print "this sequence has $numAsterisks *'s\n";
	#	$len--;
	#}

	# get hash reference to the BiasedAA section of the config hash
	my $hr = $this->{config}->{BiasedAA};
	
	#if ($debug)
	#{
	#	my $id = $seq->display_id;
	#	print "\tParameters for rule_Biased_AA\n";
	#	print "\t\tid:  $id\n";
	#
	#	if ($hr->{short}) {
	#		print "\t\tmin length:  $hr->{short}\n";
	#	}
	#	if ($hr->{long}) {
	#		print "\t\tmax length:  $hr->{long}\n";
	#	}
	#	if ($hr->{amino_acids}) {
	#		print "\t\tAAs: $hr->{amino_acids}\n\n";
	#	}
	#}

	# convert the amino acids to uppercase
	$hr->{amino_acids} =~ tr/a-z/A-Z/;

	my %locations = $this->locate($seq->seq, 'BiasedAA');
	my $percent=0;
	if( %locations ){
		$percent = $locations{
			      (reverse sort {
				$locations{$a}<=>$locations{$b}
				} keys %locations) [0]};
	}
	
	$result->{'BiasedAA'}->{'percent'} = $percent;
	return ( ($percent >= $hr->{threshold}) ? 1:0);
}
#######################################################################

# Function: rule_Biased_AA_PRP

#

# Purpose: The biased-amino-acid composition rule for AGPs. Determines 

#          if the passed sequence meets the percentage threshold for

#          biased amino acid composition and also the sequence length

#          restrictions.The percentage for each sequence (i.e., it meets

#          the length restrictions) is stored in the $result hash.

#

# Default values: threshold, 50%; no length restrictions

#

# Param 1:

# Param 2:

#

# Return value: True if the sequence meets the length restrictions and

#               biased AA composition requirements; false otherwise.

#######################################################################

sub rule_BiasedPRP

{

	my $this = shift or die;

	my $seq = shift or die;

	my $result = shift or die;

	my $len = length($seq->seq);



	#

	# IDEALLY, this should be performed in the main while loop so each

	# function doesn't have to worry about it -- but how can we modify

	# a SeqIO object's sequence?

	# ALSO, how can we handle *'s in the sequence to mean "any amino acid" --

	# should we disallow databases of this type?

	#

	# now see if the sequence ends with a '*'

	# by simply searching for the *

	# if it does, update $len by decerementing it by 1

	# otherwise, $len is the true length of the sequence

	#

	#my $numAsterisks = ($seq->seq =~ /\*/);

	#$numAsterisks = 0 unless $numAsterisks;

	#if ($numAsterisks > 0) {

	#	#print "this sequence has $numAsterisks *'s\n";

	#	$len--;

	#}



	# get hash reference to the BiasedAA section of the config hash

	my $hr = $this->{config}->{BiasedPRP};

	

	#if ($debug)

	#{

	#	my $id = $seq->display_id;

	#	print "\tParameters for rule_Biased_AA\n";

	#	print "\t\tid:  $id\n";

	#

	#	if ($hr->{short}) {

	#		print "\t\tmin length:  $hr->{short}\n";

	#	}

	#	if ($hr->{long}) {

	#		print "\t\tmax length:  $hr->{long}\n";

	#	}

	#	if ($hr->{amino_acids}) {

	#		print "\t\tAAs: $hr->{amino_acids}\n\n";

	#	}

	#}



	# convert the amino acids to uppercase

	$hr->{amino_acids} =~ tr/a-z/A-Z/;



	my %locations = $this->locate($seq->seq, 'BiasedPRP');

	my $percent=0;

	if( %locations ){

		$percent = $locations{

			      (reverse sort {

				$locations{$a}<=>$locations{$b}

				} keys %locations) [0]};

	}

	

	$result->{'BiasedPRP'}->{'percent'} = $percent;

	return ( ($percent >= $hr->{threshold}) ? 1:0);

}



#######################################################################
# Function: rule_AG_Peptide
#
# Purpose: The AG-Peptide rule for AGPs, similar to rule_Biased_AA.
#          Determines if the passed sequence meets the percentage
#          threshold of biased amino acid composition and also the
#          sequence length restrictions. The percentage for each
#          candidate sequence (i.e., it meets the length restrictions)
#          is stored in the $result hash.
#
# Default values: threshold 35%; 50 <= sequence length <= 75
#
# Param 1:
# Param 2:
#
# Return value: True if the sequence meets the length restrictions and
#               biased AA composition requirements; false otherwise.
#######################################################################
sub rule_AG_Peptide
{
	my $this = shift or die;
	my $seq = shift or die;
	my $result = shift or die;
	my $len = length($seq->seq);

	#
	# IDEALLY, this should be performed in the main while loop so each
	# function doesn't have to worry about it -- but how can we modify
	# a SeqIO object's sequence?
	# ALSO, how can we handle *'s in the sequence to mean "any amino acid" --
	# should we disallow databases of this type?
	#
	# now see if the sequence ends with a '*'
	# by simply searching for the *
	# if it does, update $len by decerementing it by 1
	# otherwise, $len is the true length of the sequence
	#
	#my $numAsterisks = ($seq->seq =~ /\*/);
	#$numAsterisks = 0 unless $numAsterisks;
	#if ($numAsterisks > 0) {
	#	#print "this sequence has $numAsterisks *'s\n";
	#	$len--;
	#}

	# get hash reference to the AgPeptide section of the config hash
	my $hr = $this->{config}->{AgPeptide};

	#if ($debug) {
	#	my $id = $seq->display_id;
	#	print "\tParameters for rule_AG_Peptide\n";
	#	print "\t\tid:  $id\n";
	#
	#	if ($hr->{short}) {
	#		print "\t\tmin length:  $hr->{short}\n";
	#	}
	#	if ($hr->{long}) {
	#		print "\t\tmax length:  $hr->{long}\n";
	#	}
	#	if ($hr->{amino_acids}) {
	#		print "      AAs: $hr->{amino_acids}\n\n";
	#	}
	#}

	# return now if the sequence is too short
	# or too long
	if ( $hr->{short} && ($len < $hr->{short}) ) {
		$result->{AgPeptide}->{percent} = 0;
		return 0;
	}
	if ( $hr->{long} && ($len > $hr->{long}) ) {
		$result->{AgPeptide}->{percent} = 0;
		return 0;
	}

	# convert the amino acids to uppercase
	$hr->{amino_acids} =~ tr/a-z/A-Z/;

	my %locations = $this->locate($seq->seq, 'AgPeptide');
	my $percent=0;
	if( %locations ){
		$percent = $locations{
			(reverse sort {	$locations{$a}<=>$locations{$b}	} keys %locations) [0]};
	}
	
	$result->{AgPeptide}->{percent} = $percent;
	return ( ($percent >= $hr->{threshold}) ? 1:0);
}



#######################################################################
# Function: rule_Fasciclin
#
# Purpose: The rule for detecting fasciclins, a type of AGP.
#
# Param 1:
# Param 2:
#
# Return value: True if the sequence has characteristics of fasciclins;
#               false otherwise.
#######################################################################
sub rule_Fasciclin
{
	my $this = shift or die;
	my $seq = shift or die;
	my $result = shift or die;

	# get hash reference to the AgPeptide section of the config hash
	my $hr = $this->{config}->{Fasciclin};
	
	#if ( $debug ) {
	#	my $id = $seq->display_id;
	#	print "\tParameters for rule_HMM_Fasciclin\n";
	#	print "\t\tid:  $id\n";
	#}
	
	my $searchStr;
	if ($hr->{search_length}) {
		$searchStr = substr($seq->seq(), 0, $hr->{search_length});
	} else {
		$searchStr = $seq->seq();
	}
	
	my $h1Count = 0;
	if ($hr->{h1_motif}) {
		my @h1MotifMatches = ($searchStr =~ /($hr->{h1_motif})/g);
		$h1Count = @h1MotifMatches;
		if ($h1Count > 0) {
			$result->{Fasciclin}->{'pattern'}->{$h1MotifMatches[0]} = 1;
			#print "     match: $h1MotifMatches[0]\n";
			
			if ($hr->{agp_region_motif}) {
				@h1MotifMatches = ($seq->seq() =~ /($hr->{agp_region_motif})/g);

				my $count1 = @h1MotifMatches;
				if ($count1 > 0) {
					#print "I found $count1 AGP region: $1\n";
				} else {
					$h1Count = 0;
				}
			}
		}
	}
	
	my $h2Count = 0;
	#if ($hr->{h2_motif}) {
	#	my @h2MotifMatches = ($searchStr =~ /($hr->{h2_motif})/g);
	#	$h2Count = @h2MotifMatches;
	#	if ($h2Count > 0) {
	#		$result->{Fasciclin}->{'pattern'}->{$h2MotifMatches[0]} = 1;
	#		#print "     match: $h2MotifMatches[0]\n";
	#		
	#		if ($hr->{agp_region_motif}) {
	#			@h2MotifMatches = ($seq->seq() =~ /($hr->{agp_region_motif})/g);
	#
	#			my $count1 = @h2MotifMatches;
	#			if ($count1 > 0) {
	#				#print "I found $count1 AGP region: $1\n";
	#			} else {
	#				$h2Count = 0;
	#			}
	#		}
	#	}
	#}

	# return the number of fasciclin domain motifs found
	$result->{Fasciclin}->{count} = $h1Count + $h2Count;
	return ( (($h1Count + $h2Count) > 0 ? 1:0) );
}

sub locate
{
	my $this = shift or die;
	my $seq = shift or die;
	my $rule_type = shift or die;

	my $h;
	my $len = length $seq;

        if($rule_type eq "Window"){
		$h = $this->{config}->{general}->{window_parameters};
	} else {
		$h = $this->{config}->{$rule_type};
		return if ( $h->{short} && ($len < $h->{short}) );
		return if ( $h->{long} && ($len > $h->{long}) );
	}

	my $window = $h->{window} || $len;
	my %locations;
	for(my $i = 0; $i< $len - $window +1; $i++)
	{
		my $winstr = substr($seq, $i, $window);
		my $count = eval '( $winstr =~ tr/' . $h->{amino_acids} . '// )';		
		$count = 0 unless $count;
		my $percent = $count / $window * 100;
		next if ( $percent < $h->{threshold} );
		$locations{$i} = $percent;
		# print $i, ' ', $percent, "\n";
	}
	return %locations;
}

return 1;

