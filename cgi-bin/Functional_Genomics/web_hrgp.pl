#!/usr/bin/perl -w
#
# Author  	: Jens Lichtenberg
# Created 	: December 2008
# Version 	: 1.0
#
#     Copyright (c) 2008  Ohio University Genomics Facility. All rights reserved.
#     This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.
#
use lib '/Perl/site/lib/';
use lib '../lib/';
use lib '/var/www/cgi-bin/Functional_Genomics';

use hrgp_rule_agp;
use hrgp_rule_lrr_extensin_prp_ad_hoc;
use hrgp_rule_signalp;
use hrgp_rule_gpi;
use hrgp_rule_classify;
use hrgp_display_misc;

use strict;
#use base ("hrgp_manifest");
use manifest;
use AminoAcid;

use Text::CSV;

# This-module specific requirements
use Bio::SeqIO;
use Bio::SeqIO::fasta;  #gud
use Bio::SeqIO::swiss;  #gud
use Bio::Tools::Sigcleave;
use Bio::Seq::RichSeq;

my $timestamp = shift or die;

#
# Default Parameters
#
my $FILE="8.39.2-12.15.2008.log";
my $BAA_OPT="No";
my $BAA_SHORT="50";
my $BAA_LONG="50";
my $BAA_THRESHOLD="50";
my $BAA_WINDOW="50";
my $BAA_AA="50";
my $AgPEP_OPT="No";
my $AgPEP_SHORT="50";
my $AgPEP_LONG="50";
my $AgPEP_THRESHOLD="50";
my $AgPEP_WINDOW="50";
my $AgPEP_AA="50";
my $FASC_OPT="No";
my $FASC_H1M="50";
my $FASC_LENGTH="50";
my $FASC_MOTIF="50";
my $EXT_OPT="No";
my $EXT_0="50";
my $EXT_1="50";
my $PRP_OPT="No";
my $PRP_0="50";
my $PRP_1="50";
my $SigP_OPT="No";
my $SigP_START="50";
my $SigP_LENGTH="50";
my $SigP_THRESHOLD="50";
my $GPI_OPT="No";
my $GPI_MIN_PHIL="50";
my $GPI_MAX_PHIL="50";
my $GPI_MIN_PHOB="50";
my $GPI_MAX_PHOB="50";
my $GPI_OMEGA="50";
my $GPI_HYDRO="50";
my $BPRP_OPT="No";
my $BPRP_SHORT="50";
my $BPRP_LONG="50";
my $BPRP_THRESHOLD="50";
my $BPRP_WINDOW="50";
my $BPRP_AA="50";
my $submit="Run";

#
# Main Function
#
#open(LOG ">$timestamp.log");

#close LOG;

	print "Starting HRGP Analysis".'<BR>';	
	print "Timestamp: $timestamp".'<BR>';

	system("ls ../".$timestamp);
	
	print '<BR><BR>';

	open(PARAMS, "../".$timestamp."/parameters.csv");
	while(my $records = <PARAMS>)
	{
		chomp($records);
		$records =~ s/\,/\=\"/;
		$records = '$'.$records.'";';
		#print $records.'<BR>';
		eval($records);
	}

	my $input = "./dat/ath1.pep";
	# try to open input file of sequences to analyze
	if ( ! -f $input ) 
	{
		print "$! : $input\n";
		return;
	}
	else
	{
		print "File opened".'<BR>';
		print LOG "File opened\n";
	}

	# create an input file object to read in the sequences to analyze
	my $in = Bio::SeqIO->new( -file => $input, -format => 'fasta' );
	if ( ! $in ) {
		print "Could not open input sequence file: $input\n";
		return;
	}

	# create an output file handle in swiss format
	# without a file handle, uses stdout
	my $verboseOut = Bio::SeqIO->newFh( -format => 'swiss' );
	if (! $verboseOut ) {
		print "Could not open swiss SeqIO for verbose output\n";
		return;
	}
	
	# a useful hash for keeping counts of various items
	my $count =
	{
		seq_count => 0,	# sequence loop counter
		hrgps => 0,		# proteins that passed at least one of our tests
		baa => 0,		# proteins that passed BiasedAA
		pep => 0,		# proteins that passed AgPeptide
		fas => 0,		# proteins that passed Fasciclin
		#lrr => 0,		# proteins that passed LRR
		ext => 0,		# proteins that passed Extensin
		prp => 0,		# proteins that passed Proline Rich
                
                baa_prp => 0,		# proteins that passed Proline Rich
		ad_hoc => 0,	# # proteins that passed an ad hoc rule
		sig => 0,		# proteins that passed Signalp
		gpi => 0, 		# proteins that passed GPI Anchor
		baa_pep => 0,	# proteins that passed both BiasedAA and AgPeptide
		baa_fas => 0,		# proteins that passed both BiasedAA and Fasciclin
		pep_fas => 0,		# proteins that passed both AgPeptide and Fasciclin
		baa_pep_fas => 0,	# proteins that passed all three of BiasedAA, AgPeptide, and Fasciclin
	};

	#
	# this is the main while loop -- it analyzes every sequence in the input file
	# 
	while ( my $seq = $in->next_seq() )
	{
		if( $seq->seq =~ /\*$/ ){
		  $seq = $seq->trunc(1, length($seq->seq)-1);
		}
		$count->{seq_count}++;
	
		my $current = default_result(); 	# record the results of this loop
		my $keeper = 0;				# number of tests this sequence passed

		# Exclude/Include any sequences based on keywords
#		if ( $this->{config}->{Annotations}->{'include_test'} )
#		{
#		  #print "testing\n";
#			if ( $this->rule_ExcludeKeywords($seq, $current) )
#			{
#			  next;
#			}
#			elsif(   $this->rule_IncludeAnnotation($seq,$current)
#				  || $this->rule_IncludeKeywords($seq, $current) )
#			{
#			 ++$keeper;
#			}
#		}

		# Characterize Biased AA Proteins
		if ( $BAA_OPT eq "Yes" && rule_Biased_AA($seq, $current) ) 
		{
			++$keeper;
			$count->{baa}++;
#			$passedBiasedAA{$seq->display_id} = 1; # for Schultz comparison
		}
		
#		# Characterize AG Peptides
#		if ( $this->{config}->{AgPeptide}->{include_test} && $this->rule_AG_Peptide($seq, $current) ) {
#			++$keeper;
#			$count->{pep}++;
#			$passedAgPeptide{$seq->display_id} = 1; # for Schultz comparison
#			
#			# if this protein also passed biasedAA, add it to %passedBiasedAA_and_AgPeptide
#			# and remove it from %passedBiasedAA and %passedAgPeptide
#			if ($passedBiasedAA{$seq->display_id}) {
#				$passedBiasedAA_and_AgPeptide{$seq->display_id} = 1; # for Schultz comparison
#				delete($passedBiasedAA{$seq->display_id});
#				delete($passedAgPeptide{$seq->display_id});
#			}
#		}
#		
#		# Characterize Fasciclins
#		if ( $this->{config}->{Fasciclin}->{include_test} && $this->rule_Fasciclin($seq, $current) ) {
#			++$keeper;
#			$count->{fas}++;
#			$passedFasciclin{$seq->display_id} = 1; # for Schultz comparison
#		}		
#
#		# Characterize LRRs
#		#if ( $this->{config}->{LRR}->{include_test} && $this->rule_LRR($seq, $current) ) {	
#		#	++$keeper;
#		#	$count->{lrr}++;
#		#}
#
#
#		# Characterize Extensins
#		if ( $this->{config}->{Extensin}->{include_test} && $this->rule_Extensin($seq, $current) ) {	
#			++$keeper;
#			$count->{ext}++;
#		}
#
#		# Characterize PRPs
#		if ( $this->{config}->{PRP}->{include_test} && $this->rule_PRP($seq, $current) ) {
#			++$keeper;
#			$count->{prp}++;
#		}
#		# Characterize Biased AA PRPs
#
#		if ( $this->{config}->{BiasedPRP}->{include_test} && $this->rule_BiasedPRP($seq, $current) ) {
#
#			++$keeper;
#
#			$count->{baa_prp}++;
#
#			$passedBiasedAA_PRP{$seq->display_id} = 1; # for Schultz comparison
#
#		}
#		# ad hoc rules (regular expressions)
#		if ( $this->{config}->{AdHoc}->{include_test} && $this->rule_AdHoc($seq, $current) ) {
#			++$keeper;
#			$count->{ad_hoc}++;
#		}
#		
#		
#		# only continue if seq passed at least one characterization
#		# test or we're testing all regardless
#		next unless ($keeper || $this->{config}->{general}->{test_all});
#		
#		# it's a keeper so remember the id and
#		# continue with some more tests
#		$current->{id} = $seq->display_id;
#		$current->{length} = length($seq->seq);
#
#		# Predict Signal Peptide Site
#		if ( $this->{config}->{SignalP}->{include_test} && $this->rule_Signal_P($seq, $current) ) {
#			++$keeper;
#			$count->{sig}++;
#			
#			# these counts are just for temporary comparisons with Schultz's results
#			if ($passedBiasedAA{$seq->display_id}) {
#				$passedBiasedAA_and_sigP{$seq->display_id} = 1;
#			} elsif ($passedAgPeptide{$seq->display_id}) {
#				$passedAgPeptide_and_sigP{$seq->display_id} = 1;
#			} elsif ($passedBiasedAA_and_AgPeptide{$seq->display_id}) {
#				$passedBiasedAA_and_AgPeptide_and_sigP{$seq->display_id} = 1;			
#			}
#			
#			if ($passedFasciclin{$seq->display_id}) {
#				$passedFasciclin_and_sigP{$seq->display_id} = 1;
#			}
#		}
#
#		# Predict GPI Anchor Site
#		if ( $this->{config}->{GPI}->{include_test} && $this->rule_GPI_Anchor($seq, $current) ) {
#			++$keeper;
#			$count->{gpi}++;
#		}
#
#		if ( $keeper ) {
#			$count->{hrgps}++;
#			
#			# Make first high-level classification
#			$this->firstClassify($seq, $current);
#	
#			# make final Showalter-style classification
#			$this->finalShowalterClassify($seq, $current, $finalCounts);
#			
#			$count->{baa_pep}++ if ( $current->{BiasedAA}->{percent} && $current->{AgPeptide}->{percent} );
#			$count->{baa_fas}++ if ( $current->{BiasedAA}->{percent} && $current->{Fasciclin}->{percent} );
#			$count->{pep_fas}++ if ( $current->{AgPeptide}->{percent} && $current->{Fasciclin}->{percent} );
#			$count->{baa_pep_fac}++ if ( $current->{BiasedAA}->{percent} && $current->{AgPeptide}->{percent} && $current->{Fasciclin}->{percent} );
#
#			# display test results for this sequence
#			$this->displayResult($verboseOut, $seq, $current);
#			
#			# write sequence to fasta file if user selected this option
#			if ($this->{config}->{general}->{write_fasta_file}) {
#				$fastaOut->write_seq($seq);
#			}
#			
#			#print STDERR $seq->display_id . "\n";
#		}
#
	} # end of main while loop -- every sequence has been analyzed



#######################################################################
# Function: default_result
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub default_result
{
	return
	{
		'id' => '',
		'length' => '',

		'BiasedAA' =>
		{
			'percent' => 0,
		},
		'AgPeptide' =>
		{
			'percent' => 0,
		},
		#'LRR' =>
		#{
		#	'flag' => '-',
		#	'count' => 0,
		#	'pattern' => {},
		#},
		'Extensin' =>
		{
			'flag' => '-',
			'count' => 0,
			'pattern' => {},
		},
		'PRP' =>
		{
			'flag' => '-',
			'count' => 0,
			'pattern' => {},
		},
		'Fasciclin' =>
		{
			'count' => 0,
			'pattern' => {},
		},
		'SignalP' =>
		{
			'cleaver' => undef,
			'count' => 0,
		},
		'GPI' =>
		{
			'count' => 0,
			'matches' => {},
		}
	}
}

