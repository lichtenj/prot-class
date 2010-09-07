#!/usr/bin/perl -w

#

# File	  : hrgp_rule_lrr_extensin_prp_ad_hoc.pm

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

# Function: rule_LRR

#

# Purpose: The rule for detecting LRRs. Scans the sequence for

#          repetitive motifs in the config hash. If a motif occurs in

#          the sequence at least a certain number of times (or constitutes

#          a given percentage of the sequence), this sequence considered

#          to have extensin characteristics. The results of each motif and

#          threshold pair are stored in the $result hash.

#

# Param 1: reference to calling object

# Param 2: the current sequence (SeqIO object)

# Param 3: reference to the current result record

#

# Return value: True if the sequence passed at least one motif-threshold

#               test; false otherwise.

#######################################################################

sub rule_LRR

{

	my $this = shift or die;		# object reference

	my $seq = shift or die;			# current sequence

	my $result = shift or die;		# current result



	# get hash reference to the LRR section of the config hash

	my $hr = $this->{config}->{LRR};



	return $this->_rule_pattern_test($seq,$result,'LRR');

}





#######################################################################

# Function: rule_Extensin

#

# Purpose: The rule for detecting extensins. Scans the sequence for

#          repetitive motifs in the config hash. If a motif occurs in

#          the sequence at least a certain number of times (or constitutes

#          a given percentage of the sequence), this sequence considered

#          to have extensin characteristics. The results of each motif and

#          threshold pair are stored in the $result hash.

#

# Param 1: reference to calling object

# Param 2: the current sequence (SeqIO object)

# Param 3: reference to the current result record

#

# Return value: True if the sequence passed at least one motif-threshold

#               test; false otherwise.

#######################################################################

sub rule_Extensin

{

	my $this = shift or die;		# object reference

	my $seq = shift or die;			# current sequence

	my $result = shift or die;		# current result



	# get hash reference to the Extensin section of the config hash

	my $hr = $this->{config}->{Extensin};



	return $this->_rule_pattern_test($seq,$result,'Extensin');

}





#######################################################################

# Function: rule_PRP

#

# Purpose: The rule for detecting PRPs. Scans the sequence for

#          repetitive motifs in the config hash. If a motif occurs in

#          the sequence at least a certain number of times (or constitutes

#          a given percentage of the sequence), this sequence considered

#          to have extensin characteristics. The results of each motif and

#          threshold pair are stored in the $result hash.

#

# Param 1: reference to calling object

# Param 2: the current sequence (SeqIO object)

# Param 3: reference to the current result record

#

# Return value: True if the sequence passed at least one motif-threshold

#               test; false otherwise.

#######################################################################

sub rule_PRP

{

	my $this = shift or die;		# object reference

	my $seq = shift or die;			# current sequence

	my $result = shift or die;		# current result



	# get hash reference to the PRP section of the config hash

	my $hr = $this->{config}->{PRP};

	

	return $this->_rule_pattern_test($seq,$result,'PRP');

}





#######################################################################

# Function: rule_AdHoc

#

# Purpose: Ad hoc regular expressions entered by user. Scans the sequence

#          for repetitive motifs in the config hash. If a motif occurs in

#          the sequence at least a certain number of times (or constitutes

#          a given percentage of the sequence), this sequence considered

#          to have extensin characteristics. The results of each motif and

#          threshold pair are stored in the $result hash.

#

# Param 1: reference to calling object

# Param 2: the current sequence (SeqIO object)

# Param 3: reference to the current result record

#

# Return value: True if the sequence passed at least one motif-threshold

#               test; false otherwise.

#######################################################################

sub rule_AdHoc

{

	my $this = shift or die;		# object reference

	my $seq = shift or die;			# current sequence

	my $result = shift or die;		# current result



	# get hash reference to the ad hoc section of the config hash

	my $hr = $this->{config}->{AdHoc};

	

	return $this->_rule_pattern_test($seq,$result,'AdHoc');

}





#######################################################################

# Function: _rule_pattern_test

#

# Purpose: Generic function which implements the regular expression

#          search for LRRs, extensins, PRPs, and ad hoc regular

#          expressions.

#

# Param 1: reference to calling object

# Param 2: the current sequence (SeqIO object)

# Param 3: reference to the current result record

# Param 4: the type of regular expression test (LRR, Extensin, etc.)

#

# Return value: True if the sequence passed at least one motif-threshold

#               test; false otherwise.

#######################################################################

sub _rule_pattern_test

{

	my $this = shift or die;		# object reference

	my $seq = shift or die;			# current sequence

	my $result = shift or die;		# current result

	my $type = shift or die;		# 'Extensin' or 'PRP'

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

	my $numAsterisks = ($seq->seq =~ /\*/);

	$numAsterisks = 0 unless $numAsterisks;

	if ($numAsterisks > 0) {

		#print "this sequence has $numAsterisks *'s\n";

		$len--;

	}

	

	

	# get hash reference to the $type section of config hash

	my $hr = $this->{config}->{$type};

	my $mt = $this->{config}->{general}->{motif_threshold};

	

	$result->{$type}->{pattern} = {};	# result rules



	# Look for all patterns by counting the number of times

	# each pattern is found in the sequence

	my $count_all = 0;

	foreach my $pattern_and_thresh (@{$hr->{pattern}})

	{

		my ($pattern,$thresh);

		if ($pattern_and_thresh)

		{

			# retrieve pattern and number from array element

			$pattern_and_thresh =~ /([^ \t]+)(?:[ \t]+(\d+))*/i;



			# if number is missing then $num == general motif threshold

			($pattern, $thresh) = ($1, $2?$2:$mt );	



			#if ($pattern_and_thresh =~ /\%/gi) {

			#	#print "pattern: $pattern\nthresh: $thresh PERCENT\n";

			#	$result->{$type}->{pattern}->{percent}->{$pattern} = 1;

			#} else {

			#	#print "pattern: $pattern\nthresh: $thresh COUNT\n";

			#	$result->{$type}->{pattern}->{percent}->{$pattern} = undef;		

			#}



			#print "    $result->{$type}->{pattern}->{percent}->{$pattern}\n";

		}

		

		# if no pattern or threshold == 0, then don't bother

		if ( $pattern && $thresh )

		{

			my @matches = ( $seq->seq() =~ /$pattern/gi );

			my $count = @matches;





			if ($pattern_and_thresh =~ /\%/gi) {

				#print "testing for $pattern percentage of $thresh\n";



				my $percentage = $count*length($pattern)*100.0/$len;

				#print "percentage = $percentage\n";

				if ($percentage >= $thresh) {

					# store the pattern and increment count_all

					$count_all += $count;

					$result->{$type}->{pattern}->{$pattern} = $percentage;

				} else {

					$result->{$type}->{pattern}->{$pattern} = 0;

				}

			} else {

				#print "testing for $pattern count of $thresh\n";

				if ($count >= $thresh) {

					# store the pattern and increment count_all

					$count_all += $count;

					$result->{$type}->{pattern}->{$pattern} = $count;

				} else {

					$result->{$type}->{pattern}->{$pattern} = 0;

				}

			}

		}

	}

	$result->{$type}->{count} = $count_all;

	return $count_all;

}





return 1;



