#!/usr/bin/perl -w

#
# File	  : PMF.pm -- Protein Mass Spectrometry Bioinformatics Experiment
# Author  : Ohio University 
# Created : July 2005
# Version : <none>
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#

use strict;
package pmf;

use base ("manifest");
use manifest;

# This-module specific requirements
use LWP::UserAgent;
use File::Glob;
use ou;  # some utility functions

########################################################################
# PMF constructor, See perltoot, Tom's Object Oriented Tutorial for Perl
#                  type perldoc perltoot, or google perltoot
#
sub new
{
  my($class) = shift;
  my($self) = $class->SUPER::new(@_);
  $self->{UseIE} = 1;		# force internet explorer (or not)
  return ( $self, $class ); 
}

########################################################################
# override this virtual function for PMF specific functionality
#
sub setup
{
  my $this = shift or die;
  $this->SUPER::setup();
}

########################################################################
# Main experiment logic for this class, this function is called by
# the experimen main: pmf_main.pl
#
sub main
{
	my $this = shift or die;
	
	# store settings from the base class into variables for [easier] local use
	my $ms_spec 	= $this->{config}->{General}->{MS_Input};
	my $msms_spec 	= $this->{config}->{General}->{Tandem_MS_Input};
	my $species 	= $this->{config}->{General}->{Taxonomy};
	my $weight 		= $this->{config}->{General}->{Mass};
	
	# retrieve the output specification from the user config
	my $output_spec = $this->{config}->{General}->{Output_Directory};
	$output_spec = ou::UnFixFileName($output_spec);
	
	# check for valid output directory
	if ( $output_spec && ! -d $output_spec )
	{
		if ( ! mkdir($output_spec) )
		{
			$this->ErrorMessage("$output_spec : $!");
			return;
		}
	}
	
	# Some similar parameters have different formats for different
    # database searches so retrieve a hash specifying which format
    # for which type such as PP for protein prospector (see below)
	my $taxonomy 	= Parameters()->{SPECIES};
	my $mass 		= Parameters()->{Mass};
	
	# No reason to continue without filespecs
	if ( ! $ms_spec )
	{
		$this->ErrorMessage('Missing a file spec for MS.');
		return;
	}

	# load the filespecs into file arrays for processeing
	# and check that the specs look good
	my @ms_files = File::Glob::bsd_glob($ms_spec) if $ms_spec;
	if (@ms_files == 0)
	{
		$this->ErrorMessage('MS Filespec did not find any files.');
        return;
	}
	
	# the ms files are the basis for the wells that are searched
	# so the ms/ms files must correspond eg well count == number of files
	my $number_wells = (@ms_files);
	printf "Processing $number_wells wells.\n\n";
	
	my $well_count = 0;

	# main loop processes each file from the ms spec
	# and finds each cooresponding well file from the
	# msms spec. All the processing is done for the current well
	foreach my $ms_file ( sort @ms_files )
	{
		# make sure the file has a long filename so we
		# can parse out the well information (from the name)
		$ms_file = ou::UnFixFileName($ms_file);
		
		######################################################
		# process ms_file for current well here and retrieve results
		# filename in the variable
		my $htm;
		{
			# First read the data file into a scalar or
			# continue on to next file
			my $data = ReadPeaks($ms_file);
			next unless $data;

			# submit to the search websites, passing GUI parameters
			# but only if the option [in config] is set
			if($this->{config}->{Mass_Spectrometry}->{MS_Fit})
			{
				our $!=undef;
				my $tax = $taxonomy->{$species}->{PP};
				my $m = $mass->{$weight}->{PP};
				
				$htm = $this->ms($m, $tax, $ms_file, $data, $output_spec);
			}
		}

		$well_count++;
		PrintLinkLine( $htm, $well_count, $number_wells, 'MS Fit' ) if $htm;			
	}
}

########################################################################
# Print a link 'nuff said... but not
# this function adds delimiters to the string so that the
# TK text widget can recognize it easily, someday this may not
# be neccesary, with different GUI library, also show percent done
sub PrintLinkLine
{
	my $file = shift or die;
	my $well_count = shift || 0;
	my $number_wells = shift || 1;
	my $type = shift || '';

	# calculate percent with default blank
	my $percent = '';	
	$percent = int($well_count/$number_wells*100) if $number_wells;
	
	# clean up the link
	my $link = File::Spec->rel2abs(ou::UnFixFileName($file));

	# show the link line
	printf "%3.0f%% %s Results: <%s>\n", $percent, $type, $link;
}

########################################################################
# Submits configuration data to the MS Fit website
sub ms
{
	my $this 	= shift or die;	
	my $mass 	= shift or die;	# parent_mass_convert parameter from gui
	my $tax 	= shift or die;	# taxonomy (pedantic species)
	my $file 	= shift or die;	# peaklist filename
	my $data 	= shift or die;	# contents of peaklist file
	my $output	= shift or die; # output directory

	# submit to ms site
	{
		my $url    = $this->parameter()->{MS}->{PP}->{url};
		my $submit = $this->parameter()->{MS}->{PP}->{submit};

		# first make sure some vital parameters are set
		{
			$submit->{output_type} = 'HTML';
			
			$submit->{species} = $tax;
			$submit->{parent_mass_convert} = $mass;
	
			$submit->{data} = $data;
			$submit->{max_reported_hits} = $this->{config}->{General}->{Number_Hits};
		}
		
		# post the parameters to the website and get the results, get the
        # results and store them in the variable $html
		my $html;
		{
			# build the output filespec based on the passed output directory
			# rely on a volume being present or not in the output spec
			my (undef,undef,$file) = File::Spec->splitpath( "$file.ms.htm" );
			my $htm = File::Spec->catpath( '', $output, $file );
			
			# don't repost to website if not overwrite and exists
			if ( -f $htm && !($this->{config}->{General}->{Overwrite}) )
			{
				# include this code if you want verbose
				#my $fixed = ou::UnFixFileName($htm);
				#print "File exists... not overwriting: $fixed\n";
			}
			else
			{
				$this->do_POST( $htm, $url, $submit );
			}
			ms_fix_links($htm);
			return $htm;
		}
	}
}

########################################################################
# convert all action= links in the html text
# from relative to absolute
sub ms_fix_links
{
	# file what needs fixing
	my $htm = shift or die;
	
	# open the parameter file and load into variable
	local $/ = undef;
	open IN, "<$htm" or die $!;
	my $html_text = <IN>;
	close IN;

	# base to which all relatives are converted to absolute
	my $base = 'http://prospector.ucsf.edu/ucsfbin4.0/msfit.cgi';
	
	# all the action urls from the html text
	my @rels = ( $html_text =~ m/action="(.*)"/gi );

	# convert each to relative
	foreach my $rel (@rels)
	{
		my $abs = URI->new_abs( $rel, $base );
		$html_text =~ s/$rel/$abs/gi
	}
	
	# write it back
	open OUT, ">$htm" or die $!;
	print OUT $html_text;
	close OUT;
}

########################################################################
# the basic LWP wrapper function for form actions
{
    # create this object ex-functio (static) to save time
    # because it it used repeatedly for web form action
    my $browser = LWP::UserAgent->new( );

	# Parameters to do_POST are:
	#  the URL, an arrayref or hashref for the key/value pairs,
	#  and then, optionally, any header lines: (key,value, key,value)
	sub do_POST
	{
		my $this = shift or die;
		my $file = shift or die; # for recievin' results

		my $resp = $browser->post(@_ ,'Content_Type' => 'form-data');
		
		if ( $resp->is_success )
		{
			if ( open(RESP, ">$file") )
			{
				print RESP $resp->content;
				close RESP;
				return 1;
			}
		}	
		return 0;
	}
}

########################################################################
# The main virtual function for overridding default parameters
#
sub config_hash
{
	return
	{
		'General' =>
		{
			'Number_Hits' => 50,
			'Taxonomy' => 'M_MUSCULUS',
			'Mass' => 'Average',
			'MS_Input' => '../dat/ms/test*.pkl',
			'Output_Directory' => '../results',
			'Overwrite' => 1,
		},
		'Mass_Spectrometry' =>
		{
			'MS_Fit' => 1,
		},
	};
}

########################################################################
# the Following functions: sortkeys() and sort_func() are used to
# determine the order that the config settings appear in the
# Manifest configuration screen.  If the function sort_func() does
# not exist then the default sort from Manifest.pm is used.
#
# The value (number) returned from this hash will be used to determine
# how the key (parameter name) is sorted in the user interface display
sub sortkeys 
{
  return
  {
	# General before MS, PTM_CFG,MSMS,
	'General'  				          => 0,
	'Mass_Spectrometry'		     	  => 1,
	'Posttranslational_Modifications' => 2,
	'Tandem_Mass_Spectrometry'	      => 3,
	
	'MS_Input'  				      => 5,
	'Output_Directory'			      => 6,
	'Tandem_MS_Input' 			      => 7,
	'Gel_Layout'				      => 8,
	
	'Taxonomy'				          => 9,
	'Mass'					          => 10,
	'Number_Hits'				      => 11,
	'Overwrite'				          => 12,
	'Detailed_Output'			      => 13,
	
	#MS and PTM
	'MS_Fit'  				          => 1,
	'Mascot'    				      => 2,
	'Aldente'				          => 3,
	'Profound'				          => 4,
	
	#Tandem_MS
	'MS_Tag'				          => 1,
	'Mascot_Ion_Search' 			  => 2,	
  }
};

########################################################################
# determine the sort of keys $a and $b to an arbitrary max
sub sort_func
{
	my ($this,$a,$b) = @_;	
	my $k1 =  defined sortkeys()->{$a} ? sortkeys()->{$a} : 99999;
	my $k2 =  defined sortkeys()->{$b} ? sortkeys()->{$b} : 99999;
	return $k1 <=> $k2;
}

########################################################################
# return the type of the entry control to be used for a particular
# variable, this function doesn't indicate that the variable is used
# or even defined, but may be in the future
sub type_func
{
  my $this = shift or die;
  my $key  = shift;
  
  my $types =
  {
	'Overwrite' => 'checkbox',
	'Taxonomy' => 'dropdown',
	'Mass' => 'dropdown',
	'Detailed_Output' => 'checkbox',
	'launch_spreadsheet' => 'checkbox',
	'MS_Fit' => 'checkbox',
	'pmfLOCAL' => 'checkbox',
	'Mascot' => 'checkbox',
	'Profound' => 'checkbox',
	'Aldente' => 'checkbox',
	'pepMAPPER' => 'checkbox',
	'MS_Tag' => 'checkbox',
	'Mascot_Ion_Search' => 'checkbox',
	'Gel_Layout' => 'browsefile',
	'MS_Input' => 'browsefile',
	'Output_Directory' => 'browsedirectory',
	'Tandem_MS_Input' => 'browsefile',
  };
  return $types->{$key} ? $types->{$key} : '';
}

# virtual function which returns a list of values
# for the given variable
sub list_func
{
	my $this = shift or die;
	my $var  = shift or die;

	my $data = Parameters(); # $xml->XMLin($file);
  
	if ( $var eq 'Mass' )
	{    
	  $! = undef; # ignore errors from previous function
	  my $m_temp = $data->{Mass};
	  my @masstype = keys %$m_temp;
	  return \@masstype;
	}
	if ( $var eq 'Taxonomy' )
	{
	  $! = undef; # ignore errors from previous function
	  my $tax = $data->{SPECIES};
	  my @taxtype = keys %$tax;
	  return \@taxtype;
	}
	return undef;	
}

########################################################################
# retrieve the peaklist from a file, and return the memory
sub ReadPeaks
{
	my $file = shift or die;
	my $data;
	if ( ou::Open( *DAT, $file, '<' ) )
	{
		while ( <DAT> )
		{
			if( /([\.\+\-\d]+).*/ )
			{
				$data .= "$1\n";
			}
		}
		close(DAT);
	}
	return $data;
}

########################################################################
# This function parameter() stores the default settings for submission
# to the web site for PMF searching
#
sub parameter
{
	return
	{
		'MS' =>
		{
			'PP' =>
			{
				'execute' => 1,
				'reference' => '[Clauser1995]',
				'url' => 'http://prospector.ucsf.edu/ucsfbin4.0/msfit.cgi',
				'submit' => {
					'detailed_report' => 1,
					'data_source' => 'Data Paste Area',
					'mod_AA' => 'Peptide N-terminal Gln to pyroGlu',
					'accession_nums' => '',
					'output_filename' => 'lastres1',
					'search_type' => 'identity',
					'high_pi' => '10',
					'prot_low_mass' => 1000,
					'min_parent_ion_matches' => 1,
					'report_title' => 'MS-Fit',
					'cterm' => 'Free Acid',
					'tolerance_units' => 'ppm',
					'sort_type' => 'Score Sort',
					'instrument_name' => 'MALDI-TOF',
					'stop_button' => 1,
					'cys' => 'acrylamide',
					'user1_name' => 'Phosphorylation of S, T and Y',
					'dna_frame_translation' => 3,
					'parent_mass_systematic_error' => 0,
					'enzyme' => 'Trypsin',
					'parent_contaminant_masses' => '',
					'nterm' => 'Hydrogen',
					'low_pi' => '3',
					'input_program_name' => 'msfit',
					'data_format' => 'PP M/Z Charge',
					'full_pi_range' => 1,
					'output_type' => 'XML',
					'mowse_pfactor' => '0.4',
					'names' => '',
					'search_name' => 'msfit',
					'mowse_on' => 1,
					'met_ox_factor' => '1',
					'seqdb_dir' => '/seqdb/',
					'species_names' => '',
					'database' => 'NCBInr.2006.02.16',
					'min_matches' => 4,
					'comment' => 'Magic Bullet digest',
					'missed_cleavages' => 1,
					'parent_mass_tolerance' => 50,
					'prot_high_mass' => 100000,
					'input_filename' => 'lastres1',
				},
			},
		},			
	};
}

########################################################################
# This function is used as a crosswalk table so that certain parameters
# (the ones a user can specify) have a common format even when different
# websites require different formats. The first level is the field name
# the second level is the common format, and the third level is the
# site specific mapping.  
sub Parameters
{
	return
	{
	  'SPECIES' => {
		'SALMONELLA' => {
		  'PROFOUND' => 'All-taxa',
		  'MASCOT' => '. . . . . . Salmonella',
		  'ALDENTE' => '1|All',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'SALMONELLA'
		},
		'H_SAPIENS' => {
		  'PROFOUND' => 'Homo-sapiens',
		  'MASCOT' => '. . . . . . . . . . . . . . . . Homo sapiens (human)',
		  'ALDENTE' => '9606|Homo sapiens',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'HOMO SAPIENS'
		},
		'M_MUSCULUS' => {
		  'PROFOUND' => 'Mus-musculus',
		  'MASCOT' => '. . . . . . . . . . . . . . . . . . Mus musculus (house mouse)',
		  'ALDENTE' => '10090|Mus musculus',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'MUS MUSCULUS'
		},
		'C_ELEGANS' => {
		  'PROFOUND' => 'Caenorhabditis-elegans',
		  'MASCOT' => '. . . . . . Caenorhabditis elegans',
		  'ALDENTE' => '1|All',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'CAENORHABDITIS ELEGANS'
		},
		'ALL' => {
		  'PROFOUND' => 'All-taxa',
		  'MASCOT' => 'All entries',
		  'ALDENTE' => '1|All',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'All'
		},
		'S_CEREVISIAE' => {
		  'PROFOUND' => 'Saccharomyces-cerevisiae',
		  'MASCOT' => '. . . . . . Saccharomyces Cerevisiae (baker\'s yeast)',
		  'ALDENTE' => '4932|Saccharomyces cerevisiae',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'SACCHAROMYCES CEREVISIAE'
		},
		'B_SUBTILIS' => {
		  'PROFOUND' => 'Bacillus-subtilis',
		  'MASCOT' => '. . . . . . Bacillus subtilis',
		  'ALDENTE' => '1423|Bacillus subtilis',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'BACILLUS SUBTILIS'
		},
		'D_MELANOGASTER' => {
		  'PROFOUND' => 'Drosophila',
		  'MASCOT' => '. . . . . . Drosophila (fruit flies)',
		  'ALDENTE' => '7227|Drosophila melanogaster',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'DROSOPHILA MELANOGASTER'
		},
		'R_NORVEGICUS' => {
		  'PROFOUND' => 'Rattus',
		  'MASCOT' => '. . . . . . . . . . . . . . . . Rattus',
		  'ALDENTE' => '10116|Rattus norvegicus',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'RATTUS NORVEGICUS'
		},
		'UNCLASSIFIED' => {
		  'PROFOUND' => 'Unclassified',
		  'MASCOT' => '. . unclassified',
		  'ALDENTE' => '-1|Unclassified',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'UNREADABLE'
		},
		'A_THALIANA' => {
		  'PROFOUND' => 'Arabidopsis-thaliana',
		  'MASCOT' => '. . . . . . Arabidopsis thaliana (thale cress)',
		  'ALDENTE' => '3702|Arabidopsis thaliana',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'ARABIDOPSIS THALIANA'
		},
		'O_SATIVA' => {
		  'PROFOUND' => 'Oryza-sativa',
		  'MASCOT' => '. . . . . . Oryza sativa (rice)',
		  'ALDENTE' => '1|All',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'ORYZA SATIVA'
		},
		'MAMMALS' => {
		  'PROFOUND' => 'Mammalia',
		  'MASCOT' => '. . . . . . . . . . . . Mammalia (mammals)',
		  'ALDENTE' => '40674|Mammalia',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'MAMMALS'
		},
		'RODENT' => {
		  'PROFOUND' => 'Rodentia',
		  'MASCOT' => '. . . . . . . . . . . . . . Rodentia (Rodents)',
		  'ALDENTE' => '9989|Rodentia',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'RODENT'
		},
		'E_COLI' => {
		  'PROFOUND' => 'Escherichia-coli',
		  'MASCOT' => '. . . . . . Escherichia coli',
		  'ALDENTE' => '562|Escherichia coli',
		  'pepMAPPER' => 'swiss',
		  'PP' => 'ESCHERICHIA COLI'
		}
	  },
	  'Mass' => {
		'Average' => {
		  'MASCOT' => 'Average',
		  'ALDENTE' => '2|Average',
		  'pepMAPPER' => 'av',
		  'PP' => 'average'
		},
		'Monoisotopic' => {
		  'MASCOT' => 'Monoisotopic',
		  'ALDENTE' => '1|Monoisotopic',
		  'pepMAPPER' => 'mo',
		  'PP' => 'monoisotopic'
		}
	  }
	};	
}

########################################################################
# this scope holds functions dealing with inks in the text control
{
	my $tag_rx = '(\{[^\}]+\}|<[^>]+>)';
	sub tagit
	{
		my $this = shift or die;
		my $text = shift or die;
		
		$text->tagConfigure('file_link',
					-foreground=>'red',
					-underline=>1,
					);
		$text->tagBind('file_link', "<Button-1>", sub { $this->do_link() } );
		
		$text->tagBind('file_link', "<Any-Enter>", sub
			   {
				 my $this = shift;
				 $this->configure(-cursor => 'top_left_arrow')
				 });
		$text->tagBind('file_link', "<Any-Leave>", sub
			   {
				 shift->configure(-cursor => 'xterm') }
			   );
		$this->add_tags( $text, $tag_rx, 'file_link' );
	}

	# respond to a click on the link in the GUI...
	# this function could be added to the base class eventually
	sub do_link
	{
		my $this = shift or die;
		
		# indicate that the app is busy by presenting
		# an hour glass
		application::MW()->Busy(-recurse => 1); 

		# retrieve the line under the click to certain delimiters
		my $start = $this->{txt}->search(-backwards, -regex, '[<\{]', 'current' );
		my $end   = $this->{txt}->search(-forwards, -regex, '[>\}]', 'current' );
		my $cmd = $this->{txt}->get( "$start + 1 chars", $end);
		my $filename = ou::FixFileName($cmd);

		# if running windows do it this way 'cause it's [much] quicker
		# we only support windows in PMF
		if ( $^O eq 'MSWin32' )
		{
			# if we are forcing the user to use IE then...
			# we may want to make this a user option
			if ( $this->{UseIE} )
			{
				use Win32::OLE;
				my $IE  = Win32::OLE->new("InternetExplorer.Application")
					|| die "Could not start Internet Explorer.Application\n";
				$IE->Navigate(ou::UnFixFileName($filename));
				$IE->{visible} = 1;
			}
			# else just pretend the file was clicked on
			else
			{
				system("start $filename");
			}
		}
		
		application::MW()->Unbusy(-recurse => 1);
	}
}

return 1;