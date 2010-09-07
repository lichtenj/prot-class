#!/usr/bin/perl -T
#
# File	  : Dump.pl -- append form data to file o.dat.
# Author  : Jens Lichtenberg, Ohio University
# Created : August 2005
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#
#		This script does two things:
#			1) Let's the user know the request was submitted and the result page
#				 will be emailed to them
#			2) Runs the tests and formats the results to HTML
use HTML::Entities ();
use strict;
use CGI;
use CGI::Carp qw ( fatalsToBrowser ); 
use Data::Dumper;
use File::Basename;

# We fork off the execution of the tests and don't care
# what happens to the child from this point.
$SIG{CHLD} = 'IGNORE';

# Make our path more secure
$ENV{'PATH'} = '/usr/bin:/bin';

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday)=localtime(time());
my $timestamp = $hour.".".$min.".".$sec."-".($mon+1).".".$mday.".".($year+1900);

my $cgi_query = CGI->new();

# Where the results will be saved
my $fsdirectory = $ENV{'DOCUMENT_ROOT'} . '/' . $timestamp;
if ($fsdirectory =~ m/^([a-zA-Z0-9\.\_\-\/\s]+)$/) 
{ 
	$fsdirectory = $1; 
}
else 
{
	die "Invalid File System Directory: " . $fsdirectory;
}

my $directory = $timestamp;

# Upload the user's input file
# Only allow certain characters in the filename for security
my $safe_filename_characters = "a-zA-Z0-9_.-";
my $users_fasta = $cgi_query->param('FILE');
if($users_fasta eq 'other')
{
	$users_fasta = $cgi_query->param('UPFILE');
	# Make sure the upload was ok
	if ( !$users_fasta )  
	{  
 		print $cgi_query->header ( );  
 		print "There was a problem uploading the FASTA file.";  
 		exit(1);
 	}
 	
 	# Seperate the components of the filename
	my ( $name, $path, $extension ) = fileparse ( $users_fasta, '\..*' );  
	$users_fasta = $name . $extension;

	# Untaint the input field by matching against safe characters
	if ( $users_fasta =~ /^([$safe_filename_characters]+)$/ )  
	{  
 	$users_fasta = $1;  
	}  
	else  
	{  
 		die "Filename contains invalid characters";  
	}  
	
	my $upload_filehandle = $cgi_query->upload("UPFILE");

	my $upload_dir = "/var/www/UserFASTAUploads";
	open (UPLOADFILE, ">$upload_dir/$users_fasta" ) or die "$!";  
	binmode UPLOADFILE;  
 
	while ( <$upload_filehandle> )  
	{  
 		print UPLOADFILE;  
	}  
 
	close UPLOADFILE;
	
	$users_fasta = $upload_dir . '/' . $users_fasta;
}





my $pid = fork();
if($pid)
{
	# The parent
	print $cgi_query->header();
	print $cgi_query->start_html('Functional Genomics Toolkit');
	print "Conducting the analysis. Your results will be available at the following location: <br />\n";
	print "http://" . $ENV{'SERVER_NAME'} . '/' . $directory . "<br />\n";
	print "You will receive an email notification upon completion of the analysis.";
	print $cgi_query->end_html();
	exit(0);
}
elsif($pid == 0)
{
	# The child
	close STDOUT;
	close STDIN;
	mkdir($fsdirectory, 0774);

	my $this->{config} = {
  		'general' => {
  			'test_all' => 0,
			'verbose' => 0,
			'input' => $users_fasta,
			'write_csv_file' => 1,
			'window_visibility' => 0,
			'format' => 'FASTA',
			'motif_threshold' => 1,
			'write_fasta_file' => 1,
			'window_parameters' => {
				'threshold' => 50,
  				'window' => 20,
  				'amino_acids' => 'PAST'
			},
			'compare_with_schultz' => 0
		},
		'AgPeptide' => {
			'short' => $cgi_query->param('AgPEP_SHORT'),
			'threshold' => $cgi_query->param('AgPEP_THRESHOLD'),
			'window' => $cgi_query->param('AgPEP_WINDOW'),
			'amino_acids' => $cgi_query->param('AgPEP_AA'),
			'long' => $cgi_query->param('AgPEP_LONG'),
			'include_test' => $cgi_query->param('AgPEP_OPT')
	  	},
	  	'AdHoc' => {
	   		'pattern' => [
	     		' ',
	      		' ',
	      		' ',
	      		' ',
	      		' '
	    	],
	    	'include_test' => 1
	  	},
	  	'BLAST' => {
	   		'THRESHOLD' => $cgi_query->param('BLAST'),
	  	},
	  	'Annotations' => {
	    	'include_keywords' => '',
	    	'annotation_expression' => '',
	    	'exclude_keywords' => 'Junk; Unlikely',
	    	'include_test' => 1
	  	},
	  	'Extensin' => {
	    	'pattern' => [
	      		$cgi_query->param('EXT_0'),
	      		$cgi_query->param('EXT_1'),
	    	],
	    	'include_test' => $cgi_query->param('EXT_OPT')
	  	},
	  	'GPI' => {
	    	'omega_site_threshold' => $cgi_query->param('GPI_OMEGA'),
	    	'min_philic_region_length' => $cgi_query->param('GPI_MIN_PHIL'),
	    	'min_phobic_region_length' => $cgi_query->param('GPI_MIN_PHOB'),
	    	'hydrophobicity_threshold' => $cgi_query->param('GPI_HYDRO'),
	    	'max_philic_region_length' => $cgi_query->param('GPI_MAX_PHIL'),
	    	'max_phobic_region_length' => $cgi_query->param('GPI_MAX_PHOB'),
	    	'include_test' => $cgi_query->param('GPI_OPT')
	  	},
	  	'BiasedAA' => {
	    	'short' => $cgi_query->param('BAA_SHORT'),
	    	'threshold' => $cgi_query->param('BAA_THRESHOLD'),
	    	'window' => $cgi_query->param('BAA_WINDOW'),
	    	'amino_acids' => $cgi_query->param('BAA_AA'),
	    	'long' => $cgi_query->param('BAA_LONG'),
	    	'include_test' => $cgi_query->param('BAA_OPT')
	  	},
	  	'Fasciclin' => {
	    	'search_length' => $cgi_query->param('FASC_LENGTH'),
	    	'agp_region_motif' => $cgi_query->param('FASC_MOTIF'),
	    	'h1_motif' => $cgi_query->param('FASC_H1M'),
	    	'include_test' => $cgi_query->param('FASC_OPT')
	  	},
	  	'SignalP' => {
	    	'search_length' => $cgi_query->param('SigP_LENGTH'),
	    	'threshold' => $cgi_query->param('SigP_THRESHOLD'),
	    	'start_position' => $cgi_query->param('SigP_START'),
	    	'include_test' => $cgi_query->param('SigP_OPT')
	  	},
	  	'PRP' => {
	    	'pattern' => [
	      		$cgi_query->param('PRP0'),
	      		$cgi_query->param('PRP1'),
	    	],
	   		'include_test' => $cgi_query->param('PRP_OPT')
		},
		'BiasedPRP' => {
			'include_test' => $cgi_query->param('BPRP_OPT'),
			'short' => $cgi_query->param('BPRP_SHORT'),
			'long' => $cgi_query->param('BPRP_LONG'),
			'threshold' => $cgi_query->param('BPRP_THRESHOLD'),
			'window' => $cgi_query->param('BPRP_WINDOW'),
			'amino_acids' => $cgi_query->param('BPRP_AA')
		}
	};
	my $email = HTML::Entities::decode($cgi_query->param('USEREMAIL'));
	if ($email =~ /^((\w|\-|\.)+\@((\w|\-|\_)+\.)+[a-zA-Z]{2,})$/) 
	{ 
		$email = $1;
	} 
	else 
	{ 
		$email = 'genomics@shiningnight.net'; 
	} 
	
	$0 = 'Functional genomic tests for ' . $email;
	
	# Write out the user specified parameters
	open(PARAM, ">$fsdirectory/parameters.csv");
	foreach my $field ( $cgi_query->param() )
	{
		my $value = "";
		if($cgi_query->param($field) =~ /\w/)
		{
			$value = $cgi_query->param($field);		
		}
		print PARAM $field.",".$value."\n";
	}
	close PARAM;
	
	# Writing Parameter Hash
	open(OUT, ">./Functional_Genomics/hrgp.cfg");
   	print OUT "#!manifest\n";
	# write the actual configuration to the file
	my $dumper = Data::Dumper->new([$this->{config}],['this->{config}']);
	$dumper->Indent(1);
	$dumper->Purity(1);
	$dumper->Deepcopy(1);
	my $sav = $dumper->Dump();
	print OUT $sav;
	# add a perl imbedded data tag in case anyone ever
	# trys to concatenate onto this file... it'll still work
	print OUT "__DATA__\n";
	close(OUT);
	
	chdir("./Functional_Genomics");
	my $output = `perl ./hrgp_main.pl`;
	my @sets = split(/\n\n/, $output);
	$sets[9] =~ s/Blast Results\t/Blast Results/;
	my $table = '<TABLE border="1">';
	my @tmp = split(/\n/,$sets[9]);

	my $userid = $email;
	$userid =~ s/\@/\_/g;

	foreach my $line (@tmp)
	{
		if($line =~ /Keywords/)
		{
			$line =~ s/\t/\<\/TD\>\<TD\>/g;
			$table .= '<TR><TD>'.$line.'</TD>';
			$table .= '<TD>SignalP 3.0</TD>';
			$table .= '<TD>Big-PI Prediction</TD>';
			$table .= '<TD>Patterns</TD>';
			$table .= '</TR>';
		}
		else
		{
			$table .= '<TR valign="top">';
			my @fields = split(/\t/,$line);
			my $gene = shift(@fields);
			$gene =~ s/\ +//;
			my $filename = $gene;
			$filename =~ s/\|/\_/g;
			if ($filename =~ m/^([a-zA-Z0-9\.\_\-\/\s]+)$/) 
			{ 
				$filename = $1; 
			}
			else 
			{
				die "Invalid Filename: " . $filename;
			}
			$table .= '<TD id="'.uc($gene).'">'.$gene.'</TD>';
			my $comments = pop(@fields);
			foreach my $field (@fields)
			{
				$table .= '<TD>'.$field.'</TD>';
			}

			if( ! -e $fsdirectory.'/ProtClass_Annotations/'.$userid.'_'.$filename.'.annot' )
			{
				open(ANNOT, ">$fsdirectory/ProtClass_Annotations/".$userid.'_'.$filename.".annot");
				print ANNOT $gene."\n";
				close ANNOT;
			}
			
			my $command = '../cgi-bin/Functional_Genomics/protclass_annotation.pl?file='.$userid.'_'.$filename.'.annot&source='.$timestamp;
			$table .= '<TD>[<a href="'.$command.'">LINK</a>]</TD>';
			$table .= '<TD>[<a href="../Blast_Results/'.$filename.'.blast">LINK</a>]<BR>';
			
			open(BIN, '/var/www/Blast_Results/'.$filename.'.blast') or print ERROR_LOG "Cannot open " . $filename . " of blast results\n";
			my $blast_hash;
			while(my $line = <BIN>)
			{
				if($line =~ /\#/){next;}
				my @tmp = split(/\s/,$line);
				if(! $blast_hash->{$tmp[1]})
				{
					$table .= '<a href="#'.$tmp[1].'">'.$tmp[1].'</a><BR>';
					$blast_hash->{$tmp[1]} = 1;
				}
			}
			close BIN;

			$table .= '</TD>';
			open(SIG, '/var/www/SigP_Results/'.$filename.'.sigp') or print ERROR_LOG "Cannot open ". $filename ." of sigp results\n";
			my @sigtmp = <SIG>;
			my $sigline = $sigtmp[-1];
			my @sigs = split(/\s+/, $sigline);
			$table .= '<TD>[<a href="../SigP_Results/'.$filename.'.sigp">'.$sigs[-2].'</a>]</TD>';

			#GPI
			open(GPI, '/var/www/GPI_Results/'.$filename.'_gpi.html') or print ERROR_LOG "Cannot open " . $filename ." of gpi results\n";
			my @tmp = <GPI>;
			my $doc = join("",@tmp);
			my $score = $doc;
			$score =~ /PValue\s+\=\s+([\.\w\-]+)/;
			my $pval = $1;
			$score =~ /Total\ Score\<\/B\>\.+\:([\.\:\ \<\>\-\(\)\w\=\/]+)/;
			$score = $1;
			$score =~ /\<B\>\s*([\-\d\.]+)[\<\/FONT\>]?\<\/B\>/;
			my $gpi_score = $1;			
			$table .= '<TD>[<a href="../GPI_Results/'.$filename.'_gpi.html">'.$gpi_score.' ('.$pval.')</a>]</TD>';

			$table .= '<TD>[<a href="../Patterns/'.$filename.'.pattern">LINK</a>]</TD>';

			$table .= '</TR>';
		}
	}
	
	$table .= '</TABLE>';
	$sets[9] = $table;
	open(LOG, ">$fsdirectory/index.html");
	my $output = join('<BR><BR>',@sets);
	$output =~ s/\n/\<BR\>/g;
	print LOG $output;
	my $ann_download = '../cgi-bin/Functional_Genomics/annotation_download.pl?userid='.$userid.'&source='.$timestamp;
	print LOG 'Export Analysis <a href="'.$ann_download.'">here</a>.';
	print LOG '<form name="annot_search" action="../cgi-bin/Functional_Genomics/search_annot.pl" method="POST">';
	print LOG '<br/>Search Annotations: <input type="text" name="searchterm" />';
	print LOG '<input type="hidden" name="userid" value="'.$userid.'" />';
	print LOG '<input type="submit" value="GO" />';
	print LOG '</form>';
	close LOG;

	open(MESSAGE, ">email.log") or die "Cannot open";
        print MESSAGE 'Your complete results are now available at: http://'.$ENV{'SERVER_NAME'} . '/' . $directory;
        print MESSAGE "\n\n";
        print MESSAGE "\nPlease DO NOT respond to this email. In case of questions please contact ".'lichtenj@ohio.edu'.".\n\n";
	close MESSAGE;

	system("mv hrgp.cfg.csv $fsdirectory/results.csv");
	system("mv hrgp.cfg.fasta $fsdirectory/results.fasta");

	system('mail '.$email.' -s "Functional Genomics Analysis Results" < email.log');	
	
	exit(0);
}
else
{
	print "Could not start analysis.  System resources may be low at the moment, try again later or notify an administrator.";
	exit(1);
}