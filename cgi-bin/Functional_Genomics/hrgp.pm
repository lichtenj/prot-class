#!/usr/bin/perl -w
#
# File	  : hrgp.pm -- Derived class for Bioinformatics Experiments
# Author  : Ohio University CNS
# Created : July 2005
# Version : <none>
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.

#############################################################################
# subclass the manifest class to separate out specific functionality
package hrgp_manifest;
use manifest;
use vars qw(@ISA);
@ISA = qw(manifest);
sub new
{
  my( $class ) = shift;
  my( $this ) = $class->SUPER::new( @_ );
  return $this;
}

#############################################################################
package hrgp;

use hrgp_rule_agp;
use hrgp_rule_lrr_extensin_prp_ad_hoc;
use hrgp_rule_signalp;
use hrgp_rule_gpi;
use hrgp_rule_classify;
use hrgp_display_misc;

use lib '/Perl/site/lib/';
use lib '../lib/';

use strict;
use base ("hrgp_manifest");
use manifest;
use AminoAcid;

use Text::CSV;

# This-module specific requirements
use Bio::SeqIO;
use Bio::SeqIO::fasta;  #gud
use Bio::SeqIO::swiss;  #gud
use Bio::Tools::Sigcleave;
use Bio::Seq::RichSeq;

 $SIG{CHLD} = 'IGNORE';
#######################################################################
# These temporary global arrays contain all of Schultz's
# AGPs, by classification, for comparing with her results
#######################################################################
my (%passedBiasedAA, %passedBiasedAA_and_sigP,
	%passedAgPeptide, %passedAgPeptide_and_sigP,
	%passedBiasedAA_and_AgPeptide, %passedBiasedAA_and_AgPeptide_and_sigP,
	%passedFasciclin, %passedFasciclin_and_sigP, %passedBiasedAA_PRP);

my %schultzClassicals = ('At5g64310', '1', 'At2g22470', '1', 'AL161596',  '1', 'At5g10430', '1', 'At1g35230', '1',
						 'At5g14380', '1', 'At5g65390', '1', 'At2g14890', '1', 'At4g09030', '1', 'At3g01700', '1',
						 'At5g18690', '1', 'At2g47930', '1', 'At3g06360', '1');
				   
my %schultzAgPeptides = ('At3g13520', '1', 'At4g26320', '1', 'At5g56540', '1', 'At5g11740', '1', 'At2g46330', '1',
						 'At3g61640', '1', 'At1g55330', '1', 'At5g53250', '1', 'At3g57690', '1', 'At5g40730', '1');
				   
my %schultzLysRichs  =  ('At2g23130', '1', 'At4g37450', '1', 'At1g68725', '1');
				   
my %schultzFasciclins = ('At5g55730', '1', 'At4g12730', '1', 'At2g24450', '1', 'At3g46550', '1', 'At4g31370', '1',
						 'At2g20520', '1', 'At2g04780', '1', 'At2g45470', '1', 'At1g03870', '1', 'At3g60900', '1',
						 'At5g03170', '1', 'At5g60490', '1', 'At5g44130', '1', 'At3g12660', '1', 'At3g52370', '1',
						 'At2g35860', '1', 'At5g06390', '1', 'At3g11700', '1', 'At1g15190', '1', 'At5g40940', '1',
						 'At5g06920', '1');

#######################################################################
# These are Schultz's GPI anchors
#######################################################################

# final counts for Showalter-classified proteins
my $finalCounts =
{
	likelyClassicals => 0,
	likelyAgPeptides => 0,
	likelyFasciclins => 0,
	likelyLrrExtensins => 0,
	likelyOtherExtensins => 0,
	possiblyChimeric => 0
};


#######################################################################
# Function: dprint
#
# Purpose: Debug routine that prints its argument only if
#          global $debug == 1.
#
# Param 1: String or list of items to print to the screen
#
# Return value:
#######################################################################
my $debug = 0;	# debug flag
sub dprint {
	print @_ if $debug;
}

#######################################################################
# HRGP constructor
sub new
{
  my($class) = shift;
  my($self) = $class->SUPER::new(@_);

  $self->{csv} = new Text::CSV;

  bless $self, $class;  # gud
  return ( $self, $class );
}

########################################################################
# override this function for HRGP specific functionality
#
sub setup
{
  my $this = shift or die;
  $this->SUPER::setup();
  # configure the text widget to override default
  $this->{txt}->configure(-wrap => 'none');
}

########################################################################
# override this function to save notes first
#
sub Close
{
  my $this = shift or die;
  SaveNotes('notes.backup');
  $this->SUPER::Close();
}

#######################################################################
# Function: hrgp_main
#
# Purpose: The main function for the HRGP module. Processes
#          each sequence in the input classifies it, and
#          displays the result in a table.
#
# Param 1: config file name (used to create the config hash)
# Param 2: log file for status messages and sequence results
sub main
{
	dprint "Starting hrgp_main\n";
	
	my $this = shift or die "No this in hrgp main";
	my $input = $this->{config}->{general}->{input};
	
	if ( not $this->{config}->{general}->{input} ) {
		print "Experiment has no input file specification.";
		return;
	}
	#
	# insert code to save the ad hoc rules to a custom rules file
	# 
	
	#
	# we need to decide how we're going to facilitate the opening of
	# the custom rules file
	#
	
	#
	# next several lines of code check the experiment settings and open the
	# input file; it also creates the csv output file and fasta output
	# file if one or both of those options are selected
	#
	
	#
	# validate the user's input to the settings pane
	#
	validate_input($this);
	
	# try to open input file of sequences to analyze
	if ( ! -f $input ) {
		print "$! : $input\n";
		return;
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
	
	# if user wants to save characterized sequences to a fasta file, create the file
	my $fastaOut;
	if ($this->{config}->{general}->{write_fasta_file}) {
		my $fileName = $this->{filename} . '.fasta';
		$fastaOut = Bio::SeqIO->new(-file => ">$fileName", -format => 'fasta' );
		if (! $fastaOut ) {
			print "Could not open output fasta file: $fileName\n";
			return;
		}
	}
	
	# if user wants to save a csv file, create it and write the header
	if ($this->{config}->{general}->{write_csv_file}) {
		my $fileName = $this->{filename} . '.csv';
		if ( ! open( CSVFILE, ">$fileName" ) ) {
			print "Could not open CSV file: $fileName\n";
			return;
		}
		
		my $local = localtime(time);
		print CSVFILE "HRGP Experiment Results for $local\n\n";
		print CSVFILE join ",", @{$this->GetFieldNames()}, "\n";
	}
	
	# print the experiment settings
	$this->displaySettings();

	# print header record for screen display
	if (!$this->{config}->{general}->{verbose}) {
		print join "\t", @{$this->GetFieldNames()}, "\n";
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
		my $keeper = 0;						# number of tests this sequence passed

		# Exclude/Include any sequences based on keywords
		if ( $this->{config}->{Annotations}->{'include_test'} )
		{
		  #print "testing\n";
			if ( $this->rule_ExcludeKeywords($seq, $current) )
			{
			  next;
			}
			elsif(   $this->rule_IncludeAnnotation($seq,$current)
				  || $this->rule_IncludeKeywords($seq, $current) )
			{
			 ++$keeper;
			}
		}

		# Characterize Biased AA Proteins
		if ( $this->{config}->{BiasedAA}->{include_test} && $this->rule_Biased_AA($seq, $current) ) {
			++$keeper;
			$count->{baa}++;
			$passedBiasedAA{$seq->display_id} = 1; # for Schultz comparison
		}
		
		# Characterize AG Peptides
		if ( $this->{config}->{AgPeptide}->{include_test} && $this->rule_AG_Peptide($seq, $current) ) {
			++$keeper;
			$count->{pep}++;
			$passedAgPeptide{$seq->display_id} = 1; # for Schultz comparison
			
			# if this protein also passed biasedAA, add it to %passedBiasedAA_and_AgPeptide
			# and remove it from %passedBiasedAA and %passedAgPeptide
			if ($passedBiasedAA{$seq->display_id}) {
				$passedBiasedAA_and_AgPeptide{$seq->display_id} = 1; # for Schultz comparison
				delete($passedBiasedAA{$seq->display_id});
				delete($passedAgPeptide{$seq->display_id});
			}
		}
		
		# Characterize Fasciclins
		if ( $this->{config}->{Fasciclin}->{include_test} && $this->rule_Fasciclin($seq, $current) ) {
			++$keeper;
			$count->{fas}++;
			$passedFasciclin{$seq->display_id} = 1; # for Schultz comparison
		}		

		# Characterize LRRs
		#if ( $this->{config}->{LRR}->{include_test} && $this->rule_LRR($seq, $current) ) {	
		#	++$keeper;
		#	$count->{lrr}++;
		#}


		# Characterize Extensins
		if ( $this->{config}->{Extensin}->{include_test} && $this->rule_Extensin($seq, $current) ) {	
			++$keeper;
			$count->{ext}++;
		}

		# Characterize PRPs
		if ( $this->{config}->{PRP}->{include_test} && $this->rule_PRP($seq, $current) ) {
			++$keeper;
			$count->{prp}++;
		}
		# Characterize Biased AA PRPs

		if ( $this->{config}->{BiasedPRP}->{include_test} && $this->rule_BiasedPRP($seq, $current) ) {

			++$keeper;

			$count->{baa_prp}++;

			$passedBiasedAA_PRP{$seq->display_id} = 1; # for Schultz comparison

		}
		# ad hoc rules (regular expressions)
		if ( $this->{config}->{AdHoc}->{include_test} && $this->rule_AdHoc($seq, $current) ) {
			++$keeper;
			$count->{ad_hoc}++;
		}

		# only continue if seq passed at least one characterization
		# test or we're testing all regardless
		next unless ($keeper || $this->{config}->{general}->{test_all});
		
		# it's a keeper so remember the id and
		# continue with some more tests
		$current->{id} = $seq->display_id;
		$current->{length} = length($seq->seq);

		# Predict Signal Peptide Site
		if ( $this->{config}->{SignalP}->{include_test} && $this->rule_Signal_P($seq, $current) ) {
			++$keeper;
			$count->{sig}++;
			
			# these counts are just for temporary comparisons with Schultz's results
			if ($passedBiasedAA{$seq->display_id}) {
				$passedBiasedAA_and_sigP{$seq->display_id} = 1;
			} elsif ($passedAgPeptide{$seq->display_id}) {
				$passedAgPeptide_and_sigP{$seq->display_id} = 1;
			} elsif ($passedBiasedAA_and_AgPeptide{$seq->display_id}) {
				$passedBiasedAA_and_AgPeptide_and_sigP{$seq->display_id} = 1;			
			}
			
			if ($passedFasciclin{$seq->display_id}) {
				$passedFasciclin_and_sigP{$seq->display_id} = 1;
			}
		}

		# Predict GPI Anchor Site
		if ( $this->{config}->{GPI}->{include_test} && $this->rule_GPI_Anchor($seq, $current) ) {
			++$keeper;
			$count->{gpi}++;
		}

		if ( $keeper ) {
			$count->{hrgps}++;
			
			# Make first high-level classification
			$this->firstClassify($seq, $current);
	
			# make final Showalter-style classification
			$this->finalShowalterClassify($seq, $current, $finalCounts);
			
			$count->{baa_pep}++ if ( $current->{BiasedAA}->{percent} && $current->{AgPeptide}->{percent} );
			$count->{baa_fas}++ if ( $current->{BiasedAA}->{percent} && $current->{Fasciclin}->{percent} );
			$count->{pep_fas}++ if ( $current->{AgPeptide}->{percent} && $current->{Fasciclin}->{percent} );
			$count->{baa_pep_fac}++ if ( $current->{BiasedAA}->{percent} && $current->{AgPeptide}->{percent} && $current->{Fasciclin}->{percent} );

			# display test results for this sequence
			$this->displayResult($verboseOut, $seq, $current);
			
			# write sequence to fasta file if user selected this option
			if ($this->{config}->{general}->{write_fasta_file}) {
				$fastaOut->write_seq($seq);
			}
			
			#print STDERR $seq->display_id . "\n";

			#
			#Find blast similarities
			#

			my $filename = $seq->display_id;
			$filename =~ s/\|/\_/g;
 			open(TMP, ">/mnt/ramdisk/".$filename.".LOG");
 			print TMP '>'.$filename."\n";
 			print TMP $seq->seq();
 			close TMP;
 

 			my @childs;
 			my $pid = fork();
			if($pid)
			{
				#this is the parent
				push(@childs,$pid);
			}
			elsif($pid == 0)
			{
				#this is the child
				system('signalp -t euk -short /mnt/ramdisk/'.$filename.'.LOG > /var/www/SigP_Results/'.$filename.'.sigp');
				exit(0);
			}
			else
			{
				die "Couldn't fork\n";
			}
			
			$pid = fork();
			if($pid)
			{
				#this is the parent
				push(@childs,$pid);
			}
			elsif($pid == 0)
			{
				#this is the child
				my $blast_threshold = 0.00000001;
				if($this->{config}->{BLAST}->{THRESHOLD})
				{
					$blast_threshold = $this->{config}->{BLAST}->{THRESHOLD};
				}

 				system("/opt/ncbi-2.2.23/bin/blastp -query /mnt/ramdisk/".$filename.".LOG -db /opt/ncbi-2.2.23/bin/TAIR9_pep -task blastp -outfmt 7 -evalue ".$blast_threshold." -out /var/www/Blast_Results/".$filename.".blast");
				exit(0);
			}
			else
			{
				die "Couldn't fork\n";
			}

			$pid = fork();
			if($pid)
			{
				#this is the parent
				push(@childs,$pid);
			}
			elsif($pid == 0)
			{
				#this is the child
				system('perl gpi_test.pl /mnt/ramdisk/'.$filename.'.LOG '.$filename);
				system('mv '.$filename.'.gpi /var/www/GPI_Results/'.$filename.'_gpi.html');
				exit(0);
			}
			else
			{
				die "Couldn't fork\n";
			}

			foreach(@childs)
			{
				waitpid($_,0);
			}

			$pid = fork();
			if($pid == 0)
			{
				$SIG{CHLD} = 'IGNORE';
				exec('perl pattern.pl /mnt/ramdisk/'.$filename.'.LOG /var/www/Patterns/'.$filename.'.pattern; rm /mnt/ramdisk/'.$filename.'.LOG');
			}
		
 			
		}
	} # end of main while loop -- every sequence has been analyzed
	
	# calculate total agps
	my $total_agps = $count->{baa} + $count->{pep} + $count->{fas};
	
	# since double hits are counted twice, subtract these variables
	# from the total... also triple hits
	$total_agps -= $count->{baa_pep};
	$total_agps -= $count->{baa_fas};
	$total_agps -= $count->{pep_fas};
	$total_agps -= $count->{baa_pep_fas} * 3;
	
	my $total_hrgps = $count->{hrgps};
	
	{
		my $baa = $this->{config}->{BiasedAA}->{include_test};
		my $pep = $this->{config}->{AgPeptide}->{include_test};
		my $fas = $this->{config}->{Fasciclin}->{include_test};
		#my $lrr = $this->{config}->{LRR}->{include_test};
		my $ext = $this->{config}->{Extensin}->{include_test};
		my $prp = $this->{config}->{PRP}->{include_test};
                my $baa_prp = $this->{config}->{BiasedPRP}->{include_test};

		my $ad_hoc = $this->{config}->{AdHoc}->{include_test}; 
		my $sig = $this->{config}->{SignalP}->{include_test};
		my $gpi = $this->{config}->{GPI}->{include_test};
		
		print "\n\n Experiment Summary \n";
		print "--------------------\n";
		print "    $count->{seq_count} proteins analyzed\n";
		print "\n    Total unique predicted HRGPs: $total_hrgps (each protein counted only once)";
		print "\n\n    $total_agps proteins had characteristics of at least one type of AGP" if $baa || $pep || $fas; 
		print "\n        $count->{baa} proteins passed the biased-amino-acid (BAA) test" if $baa;
		print "\n        $count->{pep} proteins passed the AG-peptide (PEP) test" if $pep;
		print "\n        $count->{fas} proteins passed the fasciclin (FAS) test" if $fas;
		print "\n            $count->{baa_pep} proteins passed both the BAA and PEP tests" if $baa && $pep;
		print "\n            $count->{baa_fas} proteins passed both the BAA and FAS tests" if $baa && $fas;
		print "\n            $count->{pep_fas} proteins passed both the PEP and FAS tests" if $pep && $fas;
		print "\n            $count->{baa_pep_fas} proteins passed all AGP tests (BAA, PEP, and FAS)" if $baa && $pep && $fas;
		#print "\n\n    $count->{lrr} proteins had LRR characteristics" if $lrr;					
		print "\n\n    $count->{ext} proteins had extensin characteristics" if $ext;
		print "\n\n    $count->{prp} proteins had PRP characteristics" if $prp;
                
                print "\n\n    $count->{baa_prp} proteins had Biased PRP characteristics" if $baa_prp;
		print "\n\n    $count->{ad_hoc} proteins passed an ad hoc rule" if $ad_hoc;
		print "\n\n    $count->{sig} proteins had a predicted signal peptide" if $sig;
		print "\n\n    $count->{gpi} proteins had a predicted GPI anchor" if $gpi;
	}

	print "\n\n-----------------------------------------------------------------------------\n";
	
	# print the final stats based on Showalter-style classification
	print "$finalCounts->{likelyClassicals} likely classical AGPs\n";
	print "$finalCounts->{likelyAgPeptides} likely AG-peptides\n";
	print "$finalCounts->{likelyFasciclins} likely fasciclin AGPs\n";
	#print "$finalCounts->{likelyLrrExtensins} likely LRR-extensins\n";
	print "$finalCounts->{likelyOtherExtensins} likely other extensins\n";
	print "$finalCounts->{possiblyChimeric} possibly chimeric AGP-extensins\n";
	
	print "-----------------------------------------------------------------------------\n";	
	
	# compare the results with Schultz if the user selected this option
	if ($this->{config}->{general}->{compare_with_schultz}) {
		compareWithSchultz();
	}
	
	close CSVFILE;
} # end of main routine



#######################################################################
# Function: compareWithSchultz
#
# Purpose: This function compares the results of our program
#          with those of Schultz and presents the user with
#          a table of the comparisons.
#
# Param 1: 
#
# Return value: none
#######################################################################
sub compareWithSchultz {
	print "\nDisplaying results in a format helpful for comparing with Schultz\n";
	 
	my $size = keys(%passedBiasedAA_and_sigP);
	my $rawId;

	print "\n$size protein(s) passed only the BiasedAA test AND had a signal peptide\n";
	print "    Biased-AA\t\tSchultz\n";
	print "    ---------\t\t-------\n";
	foreach my $key (keys(%passedBiasedAA_and_sigP)) {
		print "    $key\t\t";
		
		if ($key =~ /(A.*)\.[123456789]*/) {
			$rawId = $1;
			if ($schultzClassicals{$rawId}) {
				print "Classical\n";
				$schultzClassicals{$rawId} = 'found';
			} elsif ($schultzAgPeptides{$rawId}) {
				print "AG-Peptide\n";
				$schultzAgPeptides{$rawId} = 'found';
			} elsif ($schultzLysRichs{$rawId}) {
				print "Lys-Rich\n";
				$schultzLysRichs{$rawId} = 'found';
			} elsif ($schultzFasciclins{$rawId}) {
				print "Fasciclin\n";
				$schultzFasciclins{$rawId} = 'found';
			} elsif ($rawId eq 'At4g40090') {
				print "- * Schultz used AGP3 fragment, GenBank id AL161596; ours is the complete AGP3, Genbank id AF082300\n";
			} else {
				print "-\n";
			}
		} else {
			#print "couldn't match regexpr\n";
		}
	}
	
	$size = keys(%passedAgPeptide_and_sigP);
	print "\n\n$size protein(s) passed only the AgPeptide test AND had a signal peptide\n";
	print "    AG-Peptide\t\tSchultz\n";
	print "    ----------\t\t-------\n";
	foreach my $key (keys(%passedAgPeptide_and_sigP)) {
		print "    $key\t\t";
		
		if ($key =~ /(A.*)\.[123456789]*/) {
			$rawId = $1;
			if ($schultzClassicals{$rawId}) {
				print "Classical\n";
				$schultzClassicals{$rawId} = 'found';
			} elsif ($schultzAgPeptides{$rawId}) {
				print "AG-Peptide\n";
				$schultzAgPeptides{$rawId} = 'found';
			} elsif ($schultzLysRichs{$rawId}) {
				print "Lys-Rich\n";
				$schultzLysRichs{$rawId} = 'found';
			} elsif ($schultzFasciclins{$rawId}) {
				print "Fasciclin\n";
				$schultzFasciclins{$rawId} = 'found';
			} else {
				print "-\n";
			}
		} else {
			#print "couldn't match regexpr\n";
		}
	}

	$size = keys(%passedBiasedAA_and_AgPeptide_and_sigP);
	print "\n\n$size protein(s) passed both the BiasedAA test AND AgPeptide test AND had a signal peptide\n";
	print "    Both\t\tSchultz\n";
	print "    ----\t\t-------\n";
	foreach my $key (keys(%passedBiasedAA_and_AgPeptide_and_sigP)) {
		print "    $key\t\t";
		
		if ($key =~ /(A.*)\.[123456789]*/) {
			$rawId = $1;
			if ($schultzClassicals{$rawId}) {
				print "Classical\n";
				$schultzClassicals{$rawId} = 'found';
			} elsif ($schultzAgPeptides{$rawId}) {
				print "AG-Peptide\n";
				$schultzAgPeptides{$rawId} = 'found';
			} elsif ($schultzLysRichs{$rawId}) {
				print "Lys-Rich\n";
				$schultzLysRichs{$rawId} = 'found';
			} elsif ($schultzFasciclins{$rawId}) {
				print "Fasciclin\n";
				$schultzFasciclins{$rawId} = 'found';
			} else {
				print "-\n";
			}
		} else {
			#print "couldn't match regexpr\n";
		}
	}
	
	$size = keys(%passedFasciclin_and_sigP);
	print "\n\n$size protein(s) passed the Fasciclin test AND had a signal peptide\n";
	print "    Both\t\tSchultz\n";
	print "    ----\t\t-------\n";
	foreach my $key (keys(%passedFasciclin_and_sigP)) {
		print "    $key\t\t";
		
		if ($key =~ /(A.*)\.[123456789]*/) {
			$rawId = $1;
			if ($schultzClassicals{$rawId}) {
				print "Classical\n";
				$schultzClassicals{$rawId} = 'found';
			} elsif ($schultzAgPeptides{$rawId}) {
				print "AG-Peptide\n";
				$schultzAgPeptides{$rawId} = 'found';
			} elsif ($schultzLysRichs{$rawId}) {
				print "Lys-Rich\n";
				$schultzLysRichs{$rawId} = 'found';
			} elsif ($schultzFasciclins{$rawId}) {
				print "Fasciclin\n";
				$schultzFasciclins{$rawId} = 'found';
			} else {
				print "-\n";
			}
		} else {
			#print "couldn't match regexpr\n";
		}
	}
	
	
	print "\n\nSchultz classical AGPs not found by the program (signal peptide required)\n";
	foreach my $key (keys(%schultzClassicals)) {
		if ($schultzClassicals{$key} ne 'found') {
			print "    $key";
			if ($key eq 'AL161596') {
				print " <-- Schultz used AGP3 fragment; our database has the complete AGP3, At4g40090.1, from GenBank id AF082300";
			}
			print "\n";
		}
	}
	
	print "\nSchultz AG-peptides not found by the program (signal peptide required)\n";
	foreach my $key (keys(%schultzAgPeptides)) {
		if ($schultzAgPeptides{$key} ne 'found') {
			print "    $key\n";
		}
	}
	
	print "\nSchultz lys-rich proteins not found by the program (signal peptide required)\n";
	foreach my $key (keys(%schultzLysRichs)) {
		if ($schultzLysRichs{$key} ne 'found') {
			print "    $key\n";
		}
	}
	
	print "\nSchultz fasciclins not found by the program (signal peptide required)\n";
	foreach my $key (keys(%schultzFasciclins)) {
		if ($schultzFasciclins{$key} ne 'found') {
			print "    $key\n";
		}
	}
}



sub rule_IncludeAnnotation
{
  my $this = shift or die;
  my $seq = shift or die;
  my $t = $seq->desc();

  # get hash reference to the Annotations section of the config hash
  my $hr = $this->{config}->{Annotations};
  #print "rule $hr->{annotation_expression} => $t\n";
  # Check whether to exclude based on keywords
  if ( $hr->{annotation_expression} )
  {
	#print "$t => $hr->{annotation_expression}\n";
	if ( $t =~ /$hr->{annotation_expression}/i )
	{
	  return 1;
	}
  }
  return 0;   # no inclusionary keywords found
}


#{
  # Regular expression to recognize a gene locus id
  # our $generx = '\w\w[0-9]?[0-9]g[0-9]{5}(\.\d)?';
  #my $generx_ath_os = '\w\w[0-9]?[0-9]g[0-9]{5}(\.\d)?';
  #my $generx_yeast =  'Y[A-Z][LR]\d\d\d(-\w)?';
  
  #my $generx_poplar =  '[0-9]{6}';
  #my $generx = "($generx_ath_os|$generx_yeast|$generx_poplar)";
	my $generx = '.+';	
  # this function will mak a key out of just
  # the locus id part of what is passed to it
  sub MakeKey
  {
	my $k = shift;
	$k =~ s/\|/\_/g;
	if ( $k =~ /($generx)/ )
	{
	  $k = $1;
	}
	return $k;
  }

  sub tagit
  {
	my $this = shift or die;
	my $text = shift or die;
	my $position = shift;
	
    $text->tagConfigure('protein_link',
                     -foreground=>'red',
                     -underline=>1,
                     );
    $text->tagBind('protein_link', "<Button-1>", sub{$this->do_link();} );

	$text->tagBind('protein_link', "<Any-Enter>", sub
                {
                  my $this = shift;
                  $this->configure(-cursor => 'top_left_arrow')
                  });
	$text->tagBind('protein_link', "<Any-Leave>", sub
                {
                  shift->configure(-cursor => 'xterm') }
                );
	
	$this->add_tags( $text, $generx, 'protein_link', $position );
  }

  # scope control for dialog widget
  {
	my $dialog = undef;

	sub do_link
	{
	  my $this = shift or die;
	
	  # retrieve the text that is under the cursor, to the line end
	  my $start = 'current wordstart';
	  my $stop  = 'current wordend';
	  
	   # get the first word that matches the gene regular expression and fix the
	  # word to match the pattern of locus IDs and 
	  my $gene = $this->{txt}->get( 'current wordstart', 'current lineend');
	  

	  # fix the gene id format
	  $gene = MakeKey($gene);
          
	  # get the sequence based on the $gene id, also this function will return an
	  # modified gene ID which includes the version number if the original wasn't found
	  ($gene, my $seq) = $this->GetSeq($gene);
	
	  # fix the format again because the getseq may have changed it
	  $gene = MakeKey($gene);
	
	  my $button;
	  {
		if (!$dialog) # only initialize once
		{
		  $dialog = application::MW()->DialogBox(
                              -buttons => [ "Details",
                                            "Web SignalP",
                                            "Web Big-PI",
                                            "Web Search",
                                            "Cancel" ] );
		}
		$dialog->configure(-title=>"Get information for $gene");
	  
		## add widgets to the dialog box with $dialog->Add()
		#$dialog->add("Label", -text => "Get information for $gene")->pack();
		
		$dialog->bind('<Escape>',            sub { $dialog->{'selected_button'} = 'Cancel' });
		$dialog->protocol('WM_DELETE_WINDOW',sub { $dialog->{'selected_button'} = 'Cancel' });
	
		while(1)
		{
		  # later, when you need to display the dialog box
		  $button = $dialog->Show();
	  
		  if ($button eq "Details")
		  {
			$this->Repeats($gene, $seq);
		  }
                  elsif ($button eq "Web SignalP")
		  {
			$this->WebCheckSignalP($gene, $seq);
		  }
		  elsif ($button eq "Web Big-PI")
		  {
			$this->WebCheckGPI($gene, $seq);
		  }
                  #elsif ($button eq "Web Search")
		  #{
		  #	$this->CallLink($gene);
		  #}
		  else
		  {
			last;
		  }
		}
	  }  
	}
  }
  
  #use Win32::OLE;
  #my $IE  = Win32::OLE->new("InternetExplorer.Application") || die "Could not start Internet Explorer.Application\n";

#  sub CallLink
#  {
#	use URI;
#	my $this = shift or die;
#	my $gene = shift or die;
# 
#	my $url = URI->new('http://www.arabidopsis.org/servlets/Search');
# 
#	$url->query_form(
#				 'sub_type' => 'gene',
#				 'name' => $gene,
#				 'method' => '1',
#				 'search_action' => 'detail',
#				 'type' => 'general',
#				 ); 
# 
#	my $cmd = "rundll32 url.dll,FileProtocolHandler $url";
#	
#	# if we are forcing the user to use IE then...
#	if ( 0 ) # $this->{UseIE} )
#	{
#	  $IE->Navigate(ou::UnFixFileName($url));
#	  $IE->{visible} = 1;
#	}
#	# else just pretend the file was clicked on
#	else
#	{
#	  system("$cmd");
#	}
#	#
#	#
#	#die $cmd;
#	#my $p = Proc::Background->new($cmd);
#  }
#}

{
  use LWP::UserAgent;
  use HTTP::Request::Common;
  use HTML::TokeParser;
  use HTML::TreeBuilder 3;  # make sure our version isn't ancient
  my $browser; 
  
  sub WebCheckGPI
  {
	my $this = shift or die;
	my $gene = shift or die;
	my $seq = shift or die;

	my $URL = 'http://mendel.imp.ac.at/gpi/cgi-bin/gpi_pred.cgi';
	application::MW()->Busy(-recurse=>1);
	$browser = LWP::UserAgent->new( ) unless $browser;
        my $to_tr = $seq->seq();
	my %parms;
	$parms{Sequence} = $to_tr;
	$parms{"LSet"} = 'metazoa';
	# send the variables and get the results
	my $doc = do_POST( $URL, \%parms );
        $doc =~ s/<IMG[^>]*>//g;
        $doc =~ s/<A [^>]*>//g;
        $doc =~ s/<\/[Aa].*>//g;
	open TEMPHTML, ">temp.html";
        print TEMPHTML $doc;
        close TEMPHTML;
        application::MW()->Unbusy;
        system("start temp.html");

	#use HTML::FormatText;
	#my $html = parse_htmlfile($doc);
	#my $formatter = HTML::FormatText->new(leftmargin => 0, rightmargin => 50);
	#print $formatter->format($html);
	
	## build a tree from the results
	#my $root = HTML::TreeBuilder->new;
	#$root->parse($doc);
	#$root->eof( );  # done parsing for this tree
	# 
	##get the link out of $doc
	#my $page;
  	#my ($link) = extractLink($doc);
	#$page = $browser->request(GET $link);
	#$page = $browser->request(GET $link); #do 2 times because of calculation and redirect
	# 
	#$root->parse($page->content);
	#$root->eof( );  # done parsing for this tree
	# 
	## find all the text with a particular tag
	#my @found = $root->find_by_tag_name('PRE');
	# 
	#my $seq_nn;
	#my $seq_hmm; # 2 sequeneces for nn, hmm in the html
	#{
	#  $seq_hmm = (pop @found)->as_text;
	#  $seq_nn = (pop @found)->as_text;
	#  my @splinks = extractLink($page->content); #signalp website html links
	#  my $siteAddr = 'http://www.cbs.dtu.dk';
	#  GetImg($siteAddr.$splinks[1],'1.gif');
	#  GetImg($siteAddr.$splinks[3],'2.gif');  #2nd and 4th are 2 gif links;
	# 
	#  if ( -f '1.gif' && -f '2.gif' )
  	#  {
	#	# at this point we have the text from the webpage and two images
	#	# so build a dialog to diplay the results
	#	
	#	my $dialog = application::MW()->DialogBox( -title => "Sigal P and GPI web check for $gene",
	#								  -buttons => [ "Ok", "Cancel" ] );
	#	#2 image file to display in the dialog
	#	my $img1 = $dialog->Photo(-file => '1.gif');
	#	my $img2 = $dialog->Photo(-file => '2.gif');
	#	
	#	  #$dialog->Label(-text => "SignalP-nn result:")->grid; 
	#	  $dialog->Label(-text => "SignalP-nn result:\n\n".$seq_nn)->grid
	#		($dialog->Label(-image => $img1)); 
	#	  #$dialog->Label(-text => "SignalP-hmm result:")->grid; 
	#	  $dialog->Label(-text => "SignalP-hmm result:\n\n".$seq_hmm)->grid
	#		($dialog->Label(-image => $img2)); 
	#	
	#	$dialog->bind('<Escape>',            sub { $dialog->{'selected_button'} = 'Cancel' });
	#	$dialog->protocol('WM_DELETE_WINDOW',sub { $dialog->{'selected_button'} = 'Cancel' });
	#  
	#	# later, when you need to display the dialog box
	#	my $button = $dialog->Show();
	#	if ($button eq "Ok")
	#	{
	#	  
	#	}
	#	elsif ($button eq "Cancel")
	#	{
	#	}
	#	else
	#	{
	#	  # this shouldn't happen
	#	}
	#  }
	#}  
  }
  
  sub WebCheckSignalP
  {
	my $this = shift or die;
	my $gene = shift or die;
	my $seq = shift or die;
  
	my $URL = 'http://www.cbs.dtu.dk/cgi-bin/nph-webface';
	
	application::MW()->Busy(-recurse => 1); 

	$browser = LWP::UserAgent->new( ) unless $browser;
  
	my %parms;
	
	my $to_tr = $seq->seq();
	
	$parms{SEQPASTE} = $to_tr;
	$parms{configfile} = '/usr/opt/www/pub/CBS/services/SignalP-3.0/SignalP.cf';
	$parms{graphmode} = 'gif';
	$parms{orgtype} = 'euk';
	$parms{format} = 'summary';
	#$parms{SEQSUB} = '';
	$parms{method} = 'nn+hmm';
	$parms{trunc} = '70';
	
	# send the variables and get the results
	my $doc = do_POST( $URL, \%parms );
	
	# build a tree from the results
	my $root = HTML::TreeBuilder->new;
	$root->parse($doc);
	$root->eof( );  # done parsing for this tree
  
	#get the link out of $doc
	my $page;
	my ($link) = extractLink($doc);
	$page = $browser->request(GET $link);
	$page = $browser->request(GET $link); #do 2 times because of calculation and redirect
  
	$root->parse($page->content);
	$root->eof( );  # done parsing for this tree
  
	# find all the text with a particular tag
	my @found = $root->find_by_tag_name('PRE');
  
	my $seq_nn;
	my $seq_hmm; # 2 sequeneces for nn, hmm in the html
	{
	  $seq_hmm = (pop @found)->as_text;
	  $seq_nn = (pop @found)->as_text;
	  my @splinks = extractLink($page->content); #signalp website html links
	  my $siteAddr = 'http://www.cbs.dtu.dk';
	  GetImg($siteAddr.$splinks[1],'1.gif');
	  GetImg($siteAddr.$splinks[3],'2.gif');  #2nd and 4th are 2 gif links;
  
	  if ( -f '1.gif' && -f '2.gif' )
	  {
		# at this point we have the text from the webpage and two images
		# so build a dialog to diplay the results

		my $dialog = application::MW()->DialogBox( -title => "Sigal P web check for $gene",
							   -buttons => [ "Ok", "Cancel" ] );
                
		#2 image file to display in the dialog
		my $img1 = $dialog->Photo(-file => '1.gif');
		my $img2 = $dialog->Photo(-file => '2.gif');
		
		  #$dialog->Label(-text => "SignalP-nn result:")->grid; 
		  $dialog->Label(-text => "SignalP-nn result:\n\n".$seq_nn)->grid
			($dialog->Label(-image => $img1)); 
		  #$dialog->Label(-text => "SignalP-hmm result:")->grid; 
		  $dialog->Label(-text => "SignalP-hmm result:\n\n".$seq_hmm)->grid
			($dialog->Label(-image => $img2)); 
		
		$dialog->bind('<Escape>',            sub { $dialog->{'selected_button'} = 'Cancel' });
		$dialog->protocol('WM_DELETE_WINDOW',sub { $dialog->{'selected_button'} = 'Cancel' });
	  
		# later, when you need to display the dialog box
		my $button = $dialog->Show();
		if ($button eq "Ok")
		{
		  
		}
		elsif ($button eq "Cancel")
		{
		}
		else
		{
		  # this shouldn't happen
		}
	  }
	}
        application::MW()->Unbusy;
  }

  sub do_POST
  {
	#application::MW()->Busy(-recurse => 1); 
	# Parameters:
	#  the URL,
	#  an arrayref or hashref for the key/value pairs,
	#  and then, optionally, any header lines: (key,value, key,value)
	$browser = LWP::UserAgent->new( ) unless $browser;
	
	#print 
	my $resp = $browser->post(@_);

	#application::MW()->Unbusy(-recurse => 1); 

	return ($resp->content, $resp->status_line, $resp->is_success, $resp)
	  if wantarray;
	return unless $resp->is_success;
	return $resp->content;
  }
  
  # extract link from html code
  sub extractLink
  {
  # for links in the html
	use HTML::Parse;
	use HTML::Element;
	my ($html) = @_;
	my $parsed_html = HTML::Parse::parse_html($html);
	my @ret;
	my $link;
	#return @{ $parsed_html->extract_links()->[0] };
	for (@{ $parsed_html->extract_links( ) }) {
	  $link = $_->[0];
	  push @ret, $link;
	#$link =~ s/=wait$/=none/i;
	  #print "$link\n";
	}
	return @ret;
  }
  
  # download an image from internet address
  sub GetImg
  {
	  use HTTP::GetImages;
	  use File::Copy;
	  
	  my ($rtimg, $lxfile) = @_;
	  $_ = new HTTP::GetImages (
		  dir  => '.',
		  todo => [$rtimg,],
	 );
  
	  my $hash = $_->imgs_as_hash;
	  my ($tmp) = keys %{$hash}; #only one in the array
	  #print $tmp."\n";
	  $tmp =~ s/(\.\/)//g; #tmp is "./.gifgif"
	  copy($tmp,$lxfile) or die "File cannot be copied.";
	  return 1;
  }
} # browser scope control

########################################################################
# open the page files for the user annotations but be sure to save any
# notes that may have still exist from a previous version (CheckOldNotes)
# note this functionality is currently only used by HRGP so it should
# be "virtualized"
#
{
  my $notes = {};
  #dbmopen(%$notes, 'notes', 0666) or die "$!";
  #if ( CheckOldNotes() )
  #{
  #    print "Saving notes\n";
  #    SaveNotes('notes.backup');
  #}
  
  ########################################################################
  # this function gives access to the static local hash
  sub Keyword
  {
	my $k = shift;
        return $notes->{"k:$k"};
  }
  
  ########################################################################
  # this function gives access to the static local hash
  sub Comment
  {
	my $k = shift;
	return $notes->{"c:$k"};
  }

  ########################################################################
  sub CheckOldNotes
  {
    # if the notes file hasn't already been initialized
    # by reading in the old annotations and keywords
    # then do it now
    if ( (-f 'annotations.dir' && -f 'annotations.pag') || (-f 'keywords.dir' && -f 'keywords.pag') )
    {
      if ( -f 'annotations.dir' && -f 'annotations.pag' )
      {
        {
          my $annotations = {};
          dbmopen(%$annotations,'annotations', 0666) or die "$!";
          foreach my $a (keys %{$annotations})
          {
              if ( $a )
              {
                  my $val = $annotations->{$a};
                  $val =~ s/\n//g;
                  $val =~ s/[\s\t\n]+$//;
                  $notes->{"c:$a"} = $val if $val;
              }
          }
          dbmclose (%$annotations);
          $annotations = {};
        }
        foreach my $f(('annotations.dir','annotations.pag'))
        {
          rename($f, "$f.bak") or warn "Couldn't rename $f: $!\n";
        }
      }
      if ( -f 'keywords.dir' && -f 'keywords.pag' )
      {
        {
          my $keywords = {};
          dbmopen(%$keywords,'keywords', 0666) or die "$!";
          
          foreach my $k (keys %$keywords)
          {
              my $val = $keywords->{$k};
              $val =~ s/\n//g;
              $val =~ s/[\s\t\n]+$//;
              $notes->{"k:$k"} = $val if $val;
          }
          dbmclose %$keywords;
          $keywords = {};
        }
        foreach my $f(('keywords.dir','keywords.pag'))
        {
          rename($f, "$f.bak") or warn "Couldn't rename $f: $!\n";
        }
      }
      return 1; # we did find old files so return 1
    }
    return 0; # no old files found
  }
  
  sub SaveNotes
  {
    my $file = shift or die;
    
    if ( open( OUT, ">$file") )
    {
      use Data::Dumper;
      my $dumper = Data::Dumper->new([$notes],['$notes']);
      $dumper->Indent(1);
      $dumper->Purity(1);
      $dumper->Deepcopy(1);
      my $sav = $dumper->Dump();
      print OUT "#!Biohio_notes\n";
      print OUT 'hrgp'->serial_comment();
      print OUT $sav;
      print OUT "__DATA__\n";
      close(OUT);
    }
  }
  
  sub RestoreNotes
  {
    my $file = shift or die;
    if ( open (IN , "<$file") )
    {
      my $bangline = <IN>;
      # if this is not a valid file
      if ( not ( $bangline =~ /\#\!Biohio_notes/ ) )
      {
        print "Not a valid Bio Ohio notes File: $file\n";
        return ;
      }
      
      my $read;
      {
        local $/ = undef;
        $read = <IN>;
      }
      eval $read;
      close(IN);
    }
  }

  #######################################################################
  # Function: rule_ExcludeKeywords
  #
  # Purpose: 
  #
  # Param 1: 
  # Param 2: 
  # Param 3:
  #
  # Return value: 
  #######################################################################
  sub rule_ExcludeKeywords
  {
	  my $this = shift or die;
	  my $seq = shift or die;
	  my $result = shift or die;
	  my $len = length($seq->seq);
  
	  # get hash reference to the Annotations section of the config hash
	  my $hr = $this->{config}->{Annotations};
	  
	  # Check whether to exclude based on keywords
	  if ( $hr->{exclude_keywords} )
	  {
		  my $data_keywords = $notes->{"k:".MakeKey($seq->display_id)};  # gud
		  
		  if ( $hr->{exclude_keywords} && $data_keywords )
		  {
			  my @exclude_keys = split(/[ \t;,]+/, $hr->{exclude_keywords} );
			  foreach my $k ( @exclude_keys )
			  {
				  if ( $data_keywords =~ /$k/i )
				  {
					  return 1; # exclusionary keyword found
				  }
			  }
		  }
	  }
	  return 0;   # no exclusionary keywords found
  }
  
  sub rule_IncludeKeywords
  {
	  my $this = shift or die;
	  my $seq = shift or die;
	  my $result = shift or die;
	  my $len = length($seq->seq);
  
	  # get hash reference to the Annotations section of the config hash
	  my $hr = $this->{config}->{Annotations};
	  
	  # Check whether to exclude based on keywords
	  if ( $hr->{include_keywords} )
	  {
		  #my $data_keywords = $manifest::notes->{"k:$seq->display_id"};
		  my $data_keywords = $notes->{"k:".MakeKey($seq->display_id)};  # gud
		  
		  if ( $hr->{include_keywords} && $data_keywords )
		  {
			  my @include_keys = split(/[ \t;,]+/, $hr->{include_keywords} );
			  foreach my $k ( @include_keys )
			  {
				  if ( $data_keywords =~ /$k/i )
				  {
					  return 1; # inclusionary keyword found
				  }
			  }
		  }
	  }
	  return 0;   # no inclusionary keywords found
  }

  ######################################################################
  # scope for the repeats dialog
  {
	my $dialog = undef; 
	
	# display the sequence and repeats in dialog
	sub Repeats
	{
	  my $this = shift or die;
	  my $gene = shift or die;
	  my $seq = shift or die;
	  
	  my $keys = $notes->{"k:$gene"};
	  my $comm = $notes->{"c:$gene"};
	  
	  my $frq = Freqs($seq, 2, 30 );
	
	  if ( !$dialog ) # we only want to initialize once
	  {
		$dialog = application::MW()->DialogBox(-buttons=>["Ok","Cancel"])
	  }
	  $dialog->configure( -title => "Details for $gene" );
	
	  # add widgets to the dialog box with $dialog->Add()
	  $dialog->add('Label',
					-text => 'Keywords:',
					-takefocus => 0,
					-anchor => 'e',
					-justify => 'left',
				   )
	  ->grid(-row=>1,-column=>0,-sticky=>'ew');
	  $dialog->add("Entry",
							  -width => 100,
							  -textvariable=> \$keys,
							  #-height => 1,
				  )
	  ->grid(-row=>1,-column=>1,-sticky=>'ew');
	
	  $dialog->add('Label',
					-text => 'Comments:',
					-takefocus => 0,
					-anchor => 'e',
					-justify => 'left',
				   )
	  ->grid(-row=>2,-column=>0,-sticky=>'ew');
	  $dialog->add("Entry",
								  -width => 100,
								  -textvariable=>\$comm,
								  #-height => 5,
								  #-scrollbars => 'osoe',
								  )
	  ->grid(-row=>2,-column=>1,-sticky=>'ew');
	
	  my $top = $dialog->Scrolled('Text',
					  -width => 100,
					  -height => 10,
					  -scrollbars => 'osoe',
					  )
	  ->grid(-row=>3,-column=>0, -columnspan => 2, -sticky=>'ew');

	  my $mid = $dialog->Scrolled('Text',
					  -width => 100,
					  -height => 10,
					  -scrollbars => 'osoe',
					  )
	  ->grid(-row=>4,-column=>0, -columnspan => 2, -sticky=>'ew');
	
	  my $bot = $dialog->Scrolled('Text',
					  -width => 100,
					  -height => 10,
					  -scrollbars => 'osoe',
					  )
	  ->grid(-row=>5,-column=>0, -columnspan => 2, -sticky=>'ew');
	  
	  $this->configure_links($top);
	  $this->configure_links($mid);          
	  $this->configure_links($bot);
	  
	  $top->insert('end', ">".$seq->display_id." ", 'comment');
	  #print Data::Dumper->Dump([$seq]);
	  $top->insert('end', "$seq->{primary_seq}->{desc}\n", 'comment');
	  $top->insert('end', $seq->seq()."\n\n");


          # middle panel for AGP detection rules or a general window
# -------------------------------------------------------------------------------
          if( !$this->{config}->{BiasedAA}->{include_test} and
              !$this->{config}->{AgPeptide}->{include_test} and
              !$this->{config}->{general}->{window_visibility} )
          {
            $mid->insert('end', "     To use the window feature, please turn on \"Window Visibility\" on the left panel.");
          }
          else {          
            my $str = sprintf("%10s, %7s, %s\n", 'Rule', 'Window%', 'Sequence');
            $mid->insert('end', $str);
          }
  
            foreach my $rule_type (qw/ BiasedAA AgPeptide Window /) {
              next if($rule_type eq 'Window' and !$this->{config}->{general}->{window_visibility});
              next if($rule_type eq 'BiasedAA' and !$this->{config}->{BiasedAA}->{include_test});
              next if($rule_type eq 'AgPeptide' and !$this->{config}->{AgPeptide}->{include_test});
              
              my $h;
              if( $rule_type eq 'Window' ){
                $h = $this->{config}->{general}->{window_parameters};
              } else {
                $h = $this->{config}->{$rule_type};
              }

              if($rule_type eq 'BiasedAA'){
                $top->tagConfigure($rule_type, -background => 'lightgreen');
              }
              elsif($rule_type eq 'AgPeptide'){
                $top->tagConfigure($rule_type, -underline => 1);
              }
              else {
                $top->tagConfigure($rule_type, -background => 'lightblue');
              }
              
              my %locations = $this->locate($seq->seq, $rule_type);
              unless(%locations){
                my $str = sprintf("%10s  (Undetected)\n", $rule_type);
                $mid->insert('end', $str);
              }
              foreach my $key ( reverse sort
                               {$locations{$a}<=>$locations{$b}} keys %locations ){
                my $win_seq = substr($seq->seq, $key, $h->{window} || length($seq->seq) );
                my $str = sprintf("%10s  %7.2f  ", $rule_type, $locations{$key});
                $mid->insert('end', $str);
                $mid->insert('end', $win_seq, $win_seq);
                $mid->insert('end', "\n");
                
                $this->add_tags($top, $win_seq, $rule_type) if $h->{window};
                
                #########################################
                $mid->tagBind($win_seq, "<Any-Enter>",
                              sub         {
                                                $top->tagConfigure($win_seq,-background=>'red');
                                                $mid->tagConfigure($win_seq,-background=>'red');
                                                shift->configure(-cursor => 'top_left_arrow');
                                                $this->add_tags($top, $win_seq, $win_seq);
                                          }
                                          );
                ##########################################
                $mid->tagBind($win_seq, "<Any-Leave>",
                              sub         {   
                                                $top->tagConfigure($win_seq, -background=>'white');
                                                $mid->tagConfigure($win_seq, -background=>'white');
                                                shift->configure(-cursor => 'xterm');
                                                $this->rem_tags($top, $win_seq, $win_seq);
                                          }
                                          ); 
              }
            }
          

# -------------------------------------------------------------------------------

	  my $str = sprintf("%7s,%7s,%8s, sequence\n",
						'length',
						'count',
						'percent',
						);
	  $bot->insert('end',$str);
	  foreach my $seq ( sort {$frq->{$b}->{pct} <=> $frq->{$a}->{pct}} keys %$frq )
	  {
		my $str = sprintf("%7d,%7d,%8.2f,",
						  length($seq),
						  $frq->{$seq}->{frq},
						  $frq->{$seq}->{pct},
						  );
		$bot->insert('end',"$str ");
		$bot->insert('end',$seq, $seq);
                $bot->insert('end',"\n");
                
                #repeats to be highlighted by default
                if ($seq eq 'AP' || $seq eq 'PA' || $seq eq 'SP' || $seq eq 'TP' ||
                    $seq eq 'SPPP' || $seq eq 'SPPPP' || $seq eq 'SPPPPP' ) {
                  $bot->tagConfigure($seq, -background=>'yellow');
                  $top->tagConfigure($seq, -background=>'yellow');
                  $this->add_tags( $top, $seq, $seq );
                } else {
                  $bot->tagConfigure($seq, -background=>'white');
                  $top->tagConfigure($seq, -background=>'white');
                }
                
		#########################################
		$bot->tagBind($seq, "<Button-1>",
                              sub {
                                my $color = $bot->tagCget($seq, -background);
                                if($color eq 'yellow'){
                                  $top->tagConfigure($seq, -background=>'white');
                                  $bot->tagConfigure($seq, -background=>'white');
                                  $this->rem_tags($top, $seq, $seq);
                                } else {                                
                                  $bot->tagConfigure($seq, -background=>'yellow');
                                  $top->tagConfigure($seq, -background=>'yellow');
                                  $this->add_tags($top, $seq, $seq);
                                }
			      } );
		##########################################
		#$bot->tagBind($seq, "<Any-Leave>", sub
		#			  {   
		#				$top->tagConfigure($seq, -background=>'white');
		#				$bot->tagConfigure($seq, -background=>'white');
		#				shift->configure(-cursor => 'xterm');
		#				$this->rem_tags( $top, $seq, $seq );
		#			  }
		#			  );
	  }
	
	  $dialog->bind('<Escape>',            sub { $dialog->{'selected_button'} = 'Cancel' });
	  $dialog->protocol('WM_DELETE_WINDOW',sub { $dialog->{'selected_button'} = 'Cancel' });
	
	  # later, when you need to display the dialog box
	  my $button = $dialog->Show();
	  if ($button eq "Ok")
	  {
		# don't count empty strings
		$comm =~ s/[\s\t\n]+$// if $comm;
		$keys =~ s/[\s\t\n]+$// if $keys;
		
		$notes->{"c:$gene"} = ($comm?$comm:undef);
		$notes->{"k:$gene"} = ($keys?$keys:undef);
	  }
	  elsif ($button eq "Cancel")
	  {
	  }
	  else
	  {
		# this shouldn't happen
	  }
	}
  }
  
  sub Freqs
  {
    my $freqs = undef;
    my ($s,$min,$max) = @_;
    my $str = $s->seq();
    
    {
      local $/ = '*';
      chomp($str);
    }
  
    my $id =  $s->display_id;
    my $slen = length($str);
    for my $len ($min..$max)
    {
      for my $pos ( 0..$slen )
      {
        if ( $len <= $slen-$pos )
        {
          my $word = substr($str,$pos,$len);
          $freqs->{$word}++;
          #print "word:$word == $freqs->{$word}\n";
        }
      }
    }
    my @len_sorted = reverse sort {length($a) cmp length($b)} keys %$freqs;
    # discard keys that were only seen once
    #my $ret;  
    my %ret;
  
    #$ret .= sprintf("%7s,%7s,%7s,%s\n", 'length', 'count', 'percent', 'sequence' );
    foreach my $k ( @len_sorted )
    {
      if ( $freqs->{$k} > 1 )
      {
        my $l = length($k);
        my $f = $freqs->{$k};
        my $pcnt = ($l*$f)*100/$slen;
        
        if ( $pcnt > 1 )
        {
          $ret{$k}->{pct} = $pcnt;
          $ret{$k}->{frq} = $freqs->{$k};
          #$ret{seq} = $k;
          #$ret{pct} = $pcnt;
          #$ret{frq} = $freqs->{$k};
          
          #$ret .= sprintf("%7d,%7d,%6.2f %s\n", $l, $freqs->{$k}, $pcnt, $k );
        }
      }
    }
    return \%ret; # ($ret ? $ret : 'freqs here');
  }
}

sub GetSeq
{
  my $this = shift or die;
  my $gene = shift or die;
  
  my $input = $this->{'config'}->{'general'}->{'input'};

  # create an input file object to read in the sequences to analyze
  my $in = Bio::SeqIO->new( -file => $input, -format => 'fasta' );
  if ( ! $in )
  {
    print "Could not open input sequence file: $input\n";
    return;
  }
  
  # create an output file handle in swiss format
  while ( my $seq = $in->next_seq() )
  {
    if ( $seq->display_id =~ /$gene/ )
    {
      if( $seq->seq =~ /\*$/ ){
	$seq = $seq->trunc(1, length($seq->seq)-1);
      }      
      return ($seq->display_id,$seq) if wantarray;
      return $seq;
    }
  }
  return undef;
}

return 1;
