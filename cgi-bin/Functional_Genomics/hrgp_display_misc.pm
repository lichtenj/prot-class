#!/usr/bin/perl -w
#
# File	  : hrgp_display_misc.pm
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

use Text::CSV;
use manifest;

# This-module specific requirements
use Bio::SeqIO;
use Bio::Tools::Sigcleave;
use Bio::Seq::RichSeq;


#######################################################################
# Function: displaySettings
#
# Purpose: Displays all experiment settings (parameters, motif rules,
#          etc.)
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub displaySettings
{
	my $this = shift or die;
	my $hr = $this->{config};
	
	print "Experiment Settings", ($this->{filename} ? " for $this->{filename}\n" : "\n" );

	print "    General Settings\n";
	print "        Input File: $hr->{general}->{input}\n";
	print "        Input Format: $hr->{general}->{format}\n";
	print "        Write CSV File: $hr->{general}->{write_csv_file}\n";
	print "        Write Fasta File: $hr->{general}->{write_fasta_file}\n";	
	print "        Verbose: $hr->{general}->{verbose}\n";
	print "        Default Motif Threshold: $hr->{general}->{motif_threshold}%\n";
	print "        Test All: $hr->{general}->{test_all}\n";
	print "        Compare with Schultz: $hr->{general}->{compare_with_schultz}\n\n";	

	print "    Biased Amino Acid Protein Settings\n";
	if ( $hr->{BiasedAA}->{include_test}) {
		print "        min length: $hr->{BiasedAA}->{short}\n";
		print "        max length: $hr->{BiasedAA}->{long}\n";
		print "        threshold: $hr->{BiasedAA}->{threshold}\%\n";
		print "        window: $hr->{BiasedAA}->{window}\n";		
		print "        amino acids: $hr->{BiasedAA}->{amino_acids}\n\n";		
	}
	else {
		print "        Test not included in experiment\n\n";
	}

	print "    AG-Peptide Settings\n";
	if ( $hr->{AgPeptide}->{include_test}) {
		print "        min length: $hr->{AgPeptide}->{short}\n";
		print "        max length: $hr->{AgPeptide}->{long}\n";
		print "        threshold: $hr->{AgPeptide}->{threshold}\%\n";
		print "        window: $hr->{AgPeptide}->{window}\n";		
		print "        amino acids: $hr->{AgPeptide}->{amino_acids}\n\n";
	}
	else {
		print "        Test not included in experiment\n\n";
	}

	print "    Fasciclin Settings\n";
	if ( $hr->{Fasciclin}->{include_test}) {
		print "        consensus h1 motif: $hr->{Fasciclin}->{h1_motif}\n";
		#print "        consensus h2 motif: $hr->{Fasciclin}->{h2_motif}\n";		
		print "        search length (from beginning of sequence): $hr->{Fasciclin}->{search_length}\n";
		print "        AGP glycomodule region motif: $hr->{Fasciclin}->{agp_region_motif}\n\n";

	}
	else {
		print "        Test not included in experiment\n\n";
	}

	#print "    LRR Motifs and Threshold Counts/Percentages\n";
	#if ( $hr->{LRR}->{include_test}) {
	#	my $i = 0;
	#	foreach my $pattern_and_thresh (@{$hr->{LRR}->{pattern}}) {
	#		if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/)
	#		{
	#			print "        LRR[$i]: $pattern_and_thresh";
	#			
	#			if ($pattern_and_thresh =~ /\%/gi)
	#			{
	#				print " (percent)\n";
	#			} else {
	#				print " (count)\n";
	#			}
	#		}	
	#		++$i;
	#	}
	#	print "\n";
	#} else {
	#	print "        Test not included in experiment\n\n";
	#}



	print "    Extensin Motifs and Threshold Counts/Percentages\n";
	if ( $hr->{Extensin}->{include_test}) {
		my $i = 0;
		foreach my $pattern_and_thresh (@{$hr->{Extensin}->{pattern}}) {
			if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/) {
				print "        Ext[$i]: $pattern_and_thresh";
				if ($pattern_and_thresh =~ /\%/gi) {
					print " (percent)\n";
				} else {
					print " (count)\n";
				}
			}

			++$i;
		}
		print "\n";
	} else {
		print "        Test not included in experiment\n\n";
	}
	
	print "    PRP Motifs and Threshold Counts/Percentages\n";
	if ( $hr->{PRP}->{include_test}) {
		my $i = 0;
		foreach  my $pattern_and_thresh (@{$hr->{PRP}->{pattern}}) {
			if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/ ) {
				print "        Prp[$i]: $pattern_and_thresh";
				if ($pattern_and_thresh =~ /\%/gi) {
					print " (percent)\n";
				} else {
					print " (count)\n";
				}
			}
			++$i;
		}
		print "\n";
	} else {
		print "        Test not included in experiment\n\n";
	}
	
	print "    Ad Hoc Regular Expressions and Threshold Counts/Percentages\n";
	if ( $hr->{AdHoc}->{include_test}) {
		my $i = 0;
		foreach  my $pattern_and_thresh (@{$hr->{AdHoc}->{pattern}}) {
			if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/ ) {
				print "        AdHoc[$i]: $pattern_and_thresh";
				if ($pattern_and_thresh =~ /\%/gi) {
					print " (percent)\n";
				} else {
					print " (count)\n";
				}
			}
			++$i;
		}
		print "\n";
	} else {
		print "        No ad hoc tests included in experiment\n\n";
	}
	
		print "    Signal Peptide Prediction Settings\n";
	if ( $hr->{SignalP}->{include_test}) {
		print "        start position: residue $hr->{SignalP}->{start_position}\n";
		print "        search length (from start position): $hr->{SignalP}->{search_length} residues\n";
		print "        threshold: $hr->{SignalP}->{threshold}\n\n";
	}
	else {
		print "        Test not included in experiment\n\n";
	}

	print "    GPI Anchor Prediction Settings\n";
	if ( $hr->{GPI}->{include_test}) {
		print "        min philic region length: $hr->{GPI}->{min_philic_region_length}\n";
		print "        max philic region length: $hr->{GPI}->{max_philic_region_length}\n";
		print "        min phobic region length: $hr->{GPI}->{min_phobic_region_length}\n";
		print "        max phobic region length: $hr->{GPI}->{max_phobic_region_length}\n";
		print "        omega site probability threshold: $hr->{GPI}->{omega_site_threshold}\n";
		print "        hydrophobicity drop-off threshold: $hr->{GPI}->{hydrophobicity_threshold}\n\n";
	}
	else {
		print "        Test not included in experiment\n\n";
	}
}


#######################################################################
# Function: displayResult
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub displayResult
{
	my $this = shift or die;
	my $out = shift or die;		# seqIO object
	my $seq = shift or die;		# seq  object
	my $current = shift or die; 

	# retrieve array reference of field values
	my $values = $this->GetFieldValues($current);
	
	# write result to CSV file if this option was selected
	if ($this->{config}->{general}->{write_csv_file})
	{
		$this->{csv}->combine(@$values);	
		print CSVFILE $this->{csv}->string(),"\n";
	}

	print "VERBOSE	\n" if $this->{config}->{general}->{verbose};
	
	# if verbose output wasn't selected, then display concise output
	if ( ! $this->{config}->{general}->{verbose} )
	{
		# replace zeros with dashes
		my @fields = map {$_=!$_?'-':$_} @{$values};
		
		# little fix for databases like rice that have more in the primary_id
		# field than just locus id... just use the locus id part
		$fields[0] = MakeKey($fields[0]);

		# print an asterisk if the sequence has an annotation
		my $anno = (Keyword($fields[0]) || Comment($fields[0]));
		my $annotation_indicator = ( $anno ? '* ' : '  ' );
		
		print $annotation_indicator, join("\t",@fields), "\n";
	}
	else
	{
		# else print out the long display
		my $len = length($seq->seq);
		print "\n", '='x70, "\n";
		
		#my $key = MakeKey($seq->display_id);
		my $key = $seq->display_id;
		
		# print keywords and comment
		my $comm = Comment($key);
		my $keys = Keyword($key);
		my $call = $current->{classification};

		#my $annotation_indicator = ( $anno ? '* ' : '  ' );
		print "Keywords  : $keys\n" if $keys;
		print "Comments  : $comm\n" if $comm;
		print "Prediction: $call\n" if $call;
		
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
		
		# print the sequence to "out" which is a SeqIO object
		# this will display the seq in swiss format
		print $out $seq, "\n";
		
		if ( $current->{has_biased_AA_chars} eq 'yes' )
		{
			PrintComposition($this,$seq,'BiasedAA',$current);
		}
                
                if ( $current->{has_biased_PRP_chars} eq 'yes' )

		{

			PrintComposition($this,$seq,'BiasedPRP',$current);

		}

		if ( $current->{has_ag_peptide_chars} eq 'yes' )
		{
			PrintComposition($this,$seq,'AgPeptide',$current);
		}
		if ( $current->{has_fasciclin_chars} eq 'yes' )
		{
			PrintMotifs($this, $seq, 'Fasciclin', $current);
		}
		#if ( $current->{has_lrr_chars} eq 'yes' )
		#{
		#	PrintMotifs($this,$seq,'LRR',$current);
		#}
		if ( $current->{has_extensin_chars} eq 'yes' )
		{
			PrintMotifs($this,$seq,'Extensin',$current);
		}
		if ( $current->{has_prp_chars} eq 'yes' )
		{
			PrintMotifs($this,$seq,'PRP',$current);
		}
		if ( $current->{passed_ad_hoc} && ($current->{passed_ad_hoc} eq 'yes') )
		{
			PrintMotifs($this,$seq,'AdHoc',$current);
		}
		if ( defined $current->{SignalP}->{cleaver} )
		{
			print ' ' x 5, '-' x 65, "\n";
			print "     $this->{class} Signal Peptide Prediction\n";
			my $str = $current->{SignalP}->{cleaver}->pretty_print;
			$str =~ s/\n/\n     /g;
			print "     " . $str;
		}
		if ($current->{GPI}->{count})
		{
			PrintGpiAnalysis($this, $seq, $current);
		}
	}
}

#######################################################################
# Function: MakePattern
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub MakePattern
{
	my $str = shift or die;
	my $rex = shift or die;
	my $rx = '('.join('|',@$rex).')';

	# create a pattern of '^' coresponding to
	# the matched segments of the sequence
	$str =~ s/($rx)/('^' x length $1)/egi;
	$str =~ s/[^\^]/ /g;

	return $str;
}

#######################################################################
# Function: PrintMotifs
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub PrintMotifs
{
	my $this 	= shift or die;		# class
	my $seq 	= shift or die;		# Bio::Seq object
	my $type 	= shift or die;
	my $result = shift or die;

	# display the motifs that got hits
	my @motifs = keys %{$result->{$type}->{'pattern'}};
	my @hits = grep { $result->{$type}->{'pattern'}->{$_} } @motifs;
	
	print ' ' x 5, '-' x 65, "\n";
	{
		my $pad = sprintf("%-4s", uc($type) );
		print "     $type: Motif Information\n";
		print "     Count: $result->{$type}->{count} total occurrences across all motifs\n";
		#print "     Percentage: $result->{$type}->{pattern}->{$pattern} %\n";
		print "     Matches: " . '( '.join(", ", @hits  ).' )', "\n";
	}

	my $rex = \@hits;

	my $str = $seq->seq;
	my $pat = MakePattern( $str, $rex );
	
	{
		my $len = length($str);
		
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
		
		my $str_linepos;
		my $res_linepos;
		my $i=0;
		my $j=0;
	
		while($j<$len)
		{
			print ' 'x5;
			for (; $i <= $len; $i += 10)
			{
				my $block = substr($str,$i,10);
				my $bl = length($block);
				
				print( $block, " ");
				$str_linepos += $bl+1;
				if( ($i+$bl)%60 == 0 || $bl<10 ) 
				{
					print("\n");
					$i += 10;
					last;
				}
			}
			print ' 'x5;
			for (; $j <= $len; $j += 10)
			{
				my $block = substr($pat,$j,10);
				my $bl = length($block);
				
				print( $block, " ");
				$res_linepos += $bl+1;
				if( ($j+10)%60 == 0 || $bl<10 ) 
				{
					print("\n");
					$j+=10;
					last;
				}
			}
		}
	}
}

#######################################################################
# Function: MakePatternComp
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub MakePatternComp
{
	my $str = shift or die;
	my $aas = shift or die;  # amino acids

	# transliterate the string
	eval '( $str =~ tr/' . $aas . '/\^/ )';
	$str =~ s/[^\^]/ /g;
	return $str;
}

#######################################################################
# Function: PrintComposition
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub PrintComposition
{
	my $this 	= shift or die;		# class
	my $seq 	= shift or die;		# Bio::Seq object
	my $type 	= shift or die;
	
	my $result = shift or die;
	
	my $comp = $this->{config}->{$type}->{amino_acids};
	
	print ' ' x 5, '-' x 65, "\n";
	print "     $type Amino Acid Composition\n";
	print "     Counted Amino Acids: $comp\n";
	print "     Match percent: $result->{$type}->{percent}\n\n";

	my $str = $seq->seq;
	my $pat = MakePatternComp( $str, $comp );
	
	{
		my $len = length($str);
		
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
		
		my $str_linepos;
		my $res_linepos;
		my $i=0;
		my $j=0;
	
		while($i<$len)
		{
			print ' 'x5;
			for (; $i <= $len; $i += 10)
			{
				my $block = substr($str,$i,10);
				my $bl = length($block);
				
				print( $block, " ");
				$str_linepos += $bl+1;
				if( ($i+$bl)%60 == 0 || $bl<10 ) 
				{
					print("\n");
					$i += 10;
					last;
				}
			}
			print ' 'x5;
			for (; $j <= $len; $j += 10)
			{
				my $block = substr($pat,$j,10);
				my $bl = length($block);
				
				print( $block, " ");
				$res_linepos += $bl+1;
				if( ($j+10)%60 == 0 || $bl<10 ) 
				{
					print("\n");
					$j+=10;
					last;
				}
			}
		}
	}
}

#######################################################################
# Function: PrintGpiAnalysis
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub PrintGpiAnalysis
{
	my $this 	= shift or die;		# class
	my $seq 	= shift or die;		# Bio::Seq object	
	my $result = shift or die;

	print ' ' x 5, '-' x 65, "\n";
	print "     GPI Anchor Site Predictions\n\n";

	print "     Best site had omega site probability $result->{GPI}->{matches}->{prob}, hydrophobic difference $result->{GPI}->{matches}->{diff}\n";
	print "        entire GPI anchor site: $result->{GPI}->{matches}->{entireSite}\n";
	print "        w site starts at position $result->{GPI}->{matches}->{wStart} in original sequence, site = $result->{GPI}->{matches}->{wSite}\n";
	print "        hydrophobic barrier = $result->{GPI}->{matches}->{boundaryLeft} - $result->{GPI}->{matches}->{boundaryRight}\n";
	print "        hydrophilic side = $result->{GPI}->{matches}->{philicSide} ($result->{GPI}->{matches}->{philicSideVal}), hydrophobicSide = $result->{GPI}->{matches}->{phobicSide} ($result->{GPI}->{matches}->{phobicSideVal})\n";
	
	# print the entire anchor site with omega site and barrier labelled
	print "\n      Components of the GPI Anchor Site [pre-omega, omega, hydrophilic, hydrophobic]\n";
	my $preOmegaSite = substr($result->{GPI}->{matches}->{entireSite}, 0, $result->{GPI}->{matches}->{wRelative});
	print "           $preOmegaSite $result->{GPI}->{matches}->{wSite} $result->{GPI}->{matches}->{philicSide} $result->{GPI}->{matches}->{phobicSide}\n";
}

#######################################################################
# Function: GetFieldValues
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
# Retrieve an array of field values for the current record, In display format
# i.e. dashes or 'yes'es or values
sub GetFieldValues
{
	my $this = shift or die;
	my $result = shift or die;
	my $hr = $this->{config};

	my @rec;
	
	push @rec, $result->{id};

	# print BiasedAA results
	# Note: percent will only be recorded if the protein passed min and max length requirements
	# We still have to check here if the threshold was met
	if ($this->{config}->{BiasedAA}->{include_test})
	{
		my $val = sprintf("%d", $result->{BiasedAA}->{percent});
		if ($val != 0)
		{
			my $str = $val . '%';
			push @rec, $str;
		}
		else
		{
			push @rec, $val;
		}
	}
	if ($this->{config}->{BiasedPRP}->{include_test})

	{

		my $val = sprintf("%d", $result->{BiasedPRP}->{percent});

		if ($val != 0)

		{

			my $str = $val . '%';

			push @rec, $str;

		}

		else

		{

			push @rec, $val;

		}

	}
	# print AG-Peptide results
	if ($this->{config}->{AgPeptide}->{include_test}) {
		my $val = sprintf("%d", $result->{AgPeptide}->{percent});
		if ($val != 0) {
			my $str = $val . '%';
			push @rec, $str;
		} else {
			push @rec, $val;
		}
	}
	
	# print fasciclin results
	if ($this->{config}->{Fasciclin}->{include_test}) {
		push @rec, $result->{Fasciclin}->{count};
	}
	
	# print LRR rule results
	#if ($this->{config}->{LRR}->{include_test}) {
	#	foreach my $pattern_and_thresh (@{$hr->{LRR}->{pattern}}) {
	#		my ($pattern, $thresh); 
	#		if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/) {
	#			# retrieve pattern and number from array element
	#			$pattern_and_thresh =~ /([^ \t]+)(?:[ \t]+(\d+))*/;
	#
	#			# if number is missing then $num == 1
	#			($pattern, $thresh) = ($1,$2?$2:1);				
	#			
	#			my $val = sprintf("%d", $result->{LRR}->{pattern}->{$pattern});
	#			#if ($val != 0) {
	#			#	if ($pattern_and_thresh =~ /\%/gi) {
	#			#		# result is percentage
	#			#		my $str = $val . '%';
	#			#		push @rec, $str;
	#			#	} else {
	#			#		# result is number of occurrences
	#			#		push @rec, $val;
	#			#	}
	#			#} else {
	#				push @rec, $val;
	#			#}
	#		}
	#	}
	#}
	
	# print extensin rule results
	if ($this->{config}->{Extensin}->{include_test}) {
		foreach my $pattern_and_thresh (@{$hr->{Extensin}->{pattern}}) {
			my ($pattern, $thresh); 
			if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/) {
				# retrieve pattern and number from array element
				$pattern_and_thresh =~ /([^ \t]+)(?:[ \t]+(\d+))*/;

				# if number is missing then $num == 1
				($pattern, $thresh) = ($1,$2?$2:1);				
				
				my $val = sprintf("%d", $result->{Extensin}->{pattern}->{$pattern});
				#if ($val != 0) {
				#	#if ($pattern_and_thresh =~ /\%/gi) {
				#	#	# result is percentage
				#	#	my $str = $val . '%';
				#	#	push @rec, $str;
				#	#} else {
				#	#	# result is number of occurrences
				#		push @rec, $val;
				#	#}
				#} else {
					push @rec, $val;
				#}
			}
		}
	}

	# print PRP rule results
	if ($this->{config}->{PRP}->{include_test}) {
		foreach my $pattern_and_thresh (@{$hr->{PRP}->{pattern}}) {
			my ($pattern, $thresh); 
			if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/) {
				# retrieve pattern and number from array element
				$pattern_and_thresh =~ /([^ \t]+)(?:[ \t]+(\d+))*/;

				# if number is missing then $num == 1
				($pattern, $thresh) = ($1,  $2?$2:1);				
				
				my $val = sprintf("%d", $result->{PRP}->{pattern}->{$pattern});
				#if ($val != 0) {
				#	if ($pattern_and_thresh =~ /\%/gi) {
				#		# result is percentage
				#		my $str = $val . '%';
				#		push @rec, $str;
				#	} else {
				#		# result is number of occurrences
				#		push @rec, $val;
				#	}
				#} else {
					push @rec, $val;
				#}
			}
		}
	}
	
	# print ad hoc rule results
	if ($this->{config}->{AdHoc}->{include_test}) {
		foreach my $pattern_and_thresh (@{$hr->{AdHoc}->{pattern}}) {
			my ($pattern, $thresh); 
			if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/) {
				# retrieve pattern and number from array element
				$pattern_and_thresh =~ /([^ \t]+)(?:[ \t]+(\d+))*/;

				# if number is missing then $num == 1
				($pattern, $thresh) = ($1,  $2?$2:1);				
				
				my $val = sprintf("%d", $result->{AdHoc}->{pattern}->{$pattern});
				#if ($val != 0) {
				#	if ($pattern_and_thresh =~ /\%/gi) {
				#		# result is percentage
				#		my $str = $val . '%';
				#		push @rec, $str;
				#	} else {
				#		# result is number of occurrences
				#		push @rec, $val;
				#	}
				#} else {
					push @rec, $val;
				#}
			}
		}
	}
	
	# print Signal Peptide prediction results
	if ($this->{config}->{SignalP}->{include_test}) {
		push @rec, $result->{SignalP}->{count} ? $result->{SignalP}->{count} : 0; 
	}
	
	# print GPI Anchor predicition results
	if ($this->{config}->{GPI}->{include_test}) {
		push @rec, $result->{GPI}->{count}; 
	}
	
	# print final classification, keywords, and comments
	{
		my ($call,$key,$comm);
		
		if ( $result->{classification} && ($result->{classification} ne '-') )
		{
			$call = $result->{classification};
		}
		
		# gud: for rice the id needs MakeKey before the comments can be retrieved; cf. later:
		# MakeKey: "little fix for databases like rice that have more in the primary_id
		# field than just locus id... just use the locus id part"
		$key = Keyword(MakeKey($result->{id}));  # gud
		$comm= Comment(MakeKey($result->{id}));  # gud
		# print "Key: ".$key."Comment: ".$comm."\n\n";

		push @rec, $call || '';
		push @rec, $key || '';
		push @rec, $comm || '';
	}
	
	return \@rec;
}

#######################################################################
# Function: GetFieldNames
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
# Retrieve an array of field names for the current record
sub GetFieldNames
{
	my $this = shift or die;
	my $hr = $this->{config};

	my @fields;
	
	# the first field was changed from ID to 'Protein ID'
	# so that it is wider and closer to the actual size of
	# the ID so the fields line up better without an extra tab
	# which mis-aligns the csv file
	push @fields, 'Protein ID';
	
	if ($this->{config}->{BiasedAA}->{include_test}) {
		push @fields, 'BAA';
	}
        if ($this->{config}->{BiasedPRP}->{include_test}) {

		push @fields, 'BPRP';

	}
	if ($this->{config}->{AgPeptide}->{include_test}) {
		push @fields, 'AgPep';
	}

	if ($this->{config}->{Fasciclin}->{include_test}) {
		push @fields, 'Fas';
	}

	# print LRR rule names
	#if ($this->{config}->{LRR}->{include_test}) {
	#	my $lrrcount = 0;
	#	foreach my $pattern_and_thresh (@{$hr->{LRR}->{pattern}}) {
	#		my ($pattern, $thresh); 
	#		if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/) {
	#			push @fields, "LRR$lrrcount";
	#		}
	#		$lrrcount++;
	#	}
	#}

	# print extensin rule names
	if ($this->{config}->{Extensin}->{include_test}) {
		my $extcount = 0;
		foreach my $pattern_and_thresh (@{$hr->{Extensin}->{pattern}}) {
			my ($pattern, $thresh); 
			if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/) {
				push @fields, "Ext$extcount";
			}
			$extcount++;
		}
	}

	# print PRP rule names
	if ($this->{config}->{PRP}->{include_test}) {
		my $prpcount = 0;
		foreach my $pattern_and_thresh (@{$hr->{PRP}->{pattern}}) {
			my ($pattern, $thresh); 
			if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/) {
				push @fields, "Prp$prpcount";
			}
			$prpcount++;
		}
	}
	
	# print ad hoc rule names
	if ($this->{config}->{AdHoc}->{include_test}) {
		my $adhoc_count = 0;
		foreach my $pattern_and_thresh (@{$hr->{AdHoc}->{pattern}}) {
			my ($pattern, $thresh); 
			if ($pattern_and_thresh && $pattern_and_thresh !~ /^[ \t]*$/) {
				push @fields, "AdHoc$adhoc_count";
			}
			$adhoc_count++;
		}
	}
	
	if ($this->{config}->{SignalP}->{include_test}) {
		push @fields, 'SigP';
	}
	
	if ($this->{config}->{GPI}->{include_test}) {
		push @fields, 'GPI'; 
	}
	
	push @fields, 'Class';
	push @fields, 'Keywords';
	push @fields, 'Comments';
	push @fields, 'Blast Results';
	return \@fields;
}

#######################################################################
# Function: config_hash
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub config_hash
{
	return
	{
	  'general' =>
		{
			'verbose' => 0,
			'format' => 'FASTA',
			'input' => '../dat/HRGP/ATH1.pep',
			'write_csv_file' => 1,
			'write_fasta_file' => 1,
			'motif_threshold' => 1,
			'test_all' => 0,
			'compare_with_schultz' => 0,
			'window_visibility' => 0,
			'window_parameters' => {
						'window' => 20,
						'threshold' => 50,
						'amino_acids' => 'PAST'
					       }
	  },
	  'BiasedAA' =>
		{
			'include_test' => 1,
			'short' => '',
			'long' => '',
			'amino_acids' => 'PAST',
			'window' => '',
			'threshold' => 50,
	  },
	  'AgPeptide' =>
		{
			'include_test' => 1,
			'short' => 50,
			'long' => 75,
			'amino_acids' => 'PAST',
			'window' => '',
			'threshold' => 35,
	  },
	  'Fasciclin'  =>
		{
			'include_test' => 1,
			'h1_motif' => '[MALIT]T[VILS][FLCM][CAVT][PVLIS][GSTKRNDPEIV]+[DNS][DSENAGE]+[ASQM]',
			#'h2_motif' => '[LVGI][AIVLSTP][IVT][YFHQILV][AKTVQGHER][TIVFS][DNSETAI]',
			'search_length' => '',
			'agp_region_motif' => '[^P]+P[^P]+P'
	  },
	#  'LRR'	=>
	#  {
	#		'include_test' => 1,
	#		'pattern' =>
	#		[
	#			'L..L.L..[NC].L',,
	#			' ',
	#			' ',
	#			' ',
	#			' ',
	#		]
	#  },
	  'Extensin' =>
		{
			'include_test' => 1,
			'pattern' =>
			[
				'SPPPP',
				'SPPPPKH',
				'SPPPPKK',
				'SPPPPKKPYYPP',
				'SPPPPSP',
				'SPPPPSPKYYYK',
				'SPPPPSPSPPPP',
				'SPPPPYYYH',
				'SPPPPYYYK',
				'SPPPPTPVYK',
				'SPPPPVYK',
				'SPPPPVYSPPPP',
				'SPPPPVHSPPPPVA',
				'SPPPPVK',
				'SPPPPVKSPPPP',
				'SPPPPVKP'
			]
	  },
	  'PRP' =>
		{
			'include_test' => 1,
			'pattern' =>
			[
				'PP[VHTA][YTEP]K',
				'PPYVPP',
				'PPTPRP',
				'KPP[QVTA]KPP',
				'AKPSYPPTYK',
				'KPTLSPPVYT',
				'PPV[YHE][KT]',
				'PPPKIEHPPPVPVYK',
				'KKPCPP',
				'P[VI]YK',
				'(PPV){2}',
				'PPRHKPP',
				'PPVHK',
				'PPVHKPP',
				'TPVHKPP',
				'PPIHKPP',
				'PPVVKPP',
				'TPVYK',
				'PPVYKPP',
				'PPVYK',
				'KPPVYK',
				'PPVYKPP',
				'PPVEK',
				'KPPVEK',
				'KPPIYKPP',
				'PPIYK',
				'PPTYK',
				'KPPTPKP',
				'PPTEKPP',
				'PPHEKPP',
				'PPHVKPP',
				'KPPGAK',
				'KPPATKPP',
				'PPYVPP',
				'PPPATPPP',
				'PPTPRP',
				'KPPQK',
				'KPPVKPP',
				'PPVKPP',
			]
	  },
	  'SignalP' =>
		{
			'include_test' => 1,
			'threshold' => '2.9',
			'start_position' => 0,
			'search_length' => 70, 
	  },
	  'GPI' =>
		{
			'include_test' => 1,
			'max_phobic_region_length' => 24,
			'min_phobic_region_length' => 8,
			'max_philic_region_length' => 12,
			'min_philic_region_length' => 8,
			'omega_site_threshold' => 0.1, #2.5,
			'hydrophobicity_threshold' => 0.2, #2.5,
	  },
	  'Annotations' =>
		{
			'include_test' => 1,
			'exclude_keywords' => 'Junk; Unlikely',
			'include_keywords' => '',
			'annotation_expression' => ''
	  },
	  'AdHoc' =>
		{
			'include_test' => 1,
			'pattern' =>
			[
				' ',
				' ',
				' ',
				' ',
				' ',
			]
	  },
	}
}


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


#######################################################################
# Function: min
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub min {
	my $one = shift;
	my $two = shift;
	return $one < $two ? $one : $two;
}


#######################################################################
# Function: max
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub max {
	my $one = shift;
	my $two = shift;
	return $one > $two ? $one : $two;
}


#######################################################################
# Function: sortkeys
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub sortkeys 
{
  {
	'include_test'			=> 0,
	
	# for general
	'input'					=> 0,
	'write_csv_file'		=> 1,
	'write_fasta_file'		=> 2,
	'verbose'				=> 3,
	'motif_threshold'		=> 4,
	'test_all'				=> 5,
	'compare_with_schultz'	=> 6,

	# for biased AA and AG-peptide
	'short'			=> 1,
	'long'			=> 2,
	'threshold'		=> 3,
	'window'		=> 4,	
	'amino_acids'	=> 5,

	# for fasciclin
	'h1_motif'			=> 1,
	#'h2_motif'			=> 2,
	#'search_length'		=> 4,
	#'agp_region_motif'	=> 5,

	# for signal peptide
	'start_position'	=> 1,
	'search_length'		=> 2,
	'threshold'			=> 3,

	# for GPI anchor
	'min_philic_region_length'			=> 1,
	'max_philic_region_length'			=> 2,
	'min_phobic_region_length'			=> 3,
	'max_phobic_region_length'			=> 4,
	'omega_site_threshold'	=> 5,
	'hydrophobicity_threshold'	=> 6,

	# groupings
	'general'       => 0,
	'BiasedAA'      => 1,
	'AgPeptide'     => 2,
	'Fasciclin'     => 3,
	#'LRR'			=> 4,
	'Extensin'      => 5,
	'PRP'           => 6,
	'SignalP'       => 7,
	'GPI'           => 8,
	'AdHoc'			=> 9,
	'Annotations'	=> 10,
  }
};


#######################################################################
# Function: sort_func
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub sort_func
{
	my ($this,$a,$b) = @_;
	
  my $k1 =  defined sortkeys()->{$a} ? sortkeys()->{$a} : 99999;
  my $k2 =  defined sortkeys()->{$b} ? sortkeys()->{$b} : 99999;
  return $k1 <=> $k2;
}


#######################################################################
# Function: type_func
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
sub type_func
{
  my $this = shift or die;
  my $key  = shift;
  
  my $types =
  {
    'input' => 'browsefile',
    'verbose' => 'checkbox',
    'include_test' => 'checkbox',
	'write_csv_file' => 'checkbox',
	'write_fasta_file' => 'checkbox',	
    'overwrite' => 'checkbox',
	'test_all' => 'checkbox',
	'compare_with_schultz' => 'checkbox',
    'format' => 'hidden',
    'window_parameters' => 'settings',
    'window_visibility' => 'checkbox'
  };
  return $types->{$key} ? $types->{$key} : '';
}


#######################################################################
# Function: validate_input
#
# Purpose: 
#
# Param 1: 
# Param 2: 
# Param 3: 
#
# Return value: 
#######################################################################
# this function examines the input config file and determines if it contains
# any errors, such as invalid amino acids or spaces in regular expressions
sub validate_input {
	my $this = shift or die "No this in hrgp validateInput";

##### this is debugged and works	
	#
	# General Settings
	#
	# check the file name
	#my $fn = $this->{config}->{general}->{input};
	
	# check that motif threshold is a number and is not blank
	my $val = ($this->{config}->{general}->{motif_threshold} =~ /[0-9]+/);
	if ($val > 0) {
		#
	} else {
		die "Default motif threshold must be an integer and not left blank.\n\nPlease correct this setting and try again.\n";
	}
	
	
	
	#
	# BiasedAA
	#
	# check that short and long are numbers or are blank
	if ( !($this->{config}->{BiasedAA}->{short}) || ($this->{config}->{BiasedAA}->{short} =~ /[0-9]+/)) {
	} else {
		die "BiasedAA minimum length must be an integer or left blank.\n\nPlease correct this setting and try again.\n";
	}
	if ( !($this->{config}->{BiasedAA}->{long}) || ($this->{config}->{BiasedAA}->{long} =~ /[0-9]+/)) {
	} else {
		die "BiasedAA maximum length must be an integer or left blank.\n\nPlease correct this setting and try again.\n";
	}
	
	# check that threshold is a number and is not blank
	if ( !($this->{config}->{BiasedAA}->{threshold} =~ /([0-9]+)|([0-9]+\.[0-9]+)/) ) {
		die "BiasedAA threshold must be an integer or floating point number.\n\nPlease correct this setting and try again.\n";
	}
	# The above reg-expr is matching things like 29.9.6, which it shouldn't!

	# check that amino acids are legal amino acids, contain no spaces or punctuation, and is not blank
	if ( !($this->{config}->{BiasedAA}->{amino_acids} =~ /[ABCDEFGHIKLMNOPQRSTVWXYZ]+/i) ) {
		die "BiasedAA amino acids must be a sequence of one or more amino acids (e.g. PAST).\n\nPlease correct this setting and try again.\n";
	}
	# I can't get this regular expression right...
	
	
	
	
	#
	# AgPeptide
	#
	# check that short and long are numbers or are blank
	
	# check that threshold is a number and is not blank
	
	# check that amino acids are legal amino acids, contain no spaces or punctuation, and is not blank	
	
	
	
	#
	# Fasciclin
	#
	# check that motif is a valid regular expression with no spaces and only legal amino acids
	
	# check that search length is a number and is not blank
	
	# check that agp region motif is a valid regular expression with no spaces and only legal amino acids
	
	
	#
	# Extensin
	#
	# check that each pattern is a valid regular expression with no spaces and only legal amino acids
	
	
	#
	# PRP
	#
	# check that each pattern is a valid regular expression with no spaces and only legal amino acids
	
	
	#
	# SignalP
	#
	# check that start position is a number and is not blank
	
	# check that search length is a number and is not blank
	
	# check that threshold is a number and is not blank
	
	
	#
	# GPI
	#
	# check that min and max philic and min and max phobic are numbers and are not blank
	
	# check that omega site threshold is a number and is not blank
	
	# check that hydrophobicity threshold is a number and is not blank

	return 1;
}



return 1;

