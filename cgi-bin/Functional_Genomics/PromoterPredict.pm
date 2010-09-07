#!/usr/bin/perl -w
#
# File	  : PromoterPredict.pm -- Derived class for Bioinformatics Experiments
# Author  : Ohio University CNS
# Created : July 2006
# Version : <none>
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.

# First create a module derived from "manifest"
package PromoterPredict; # 

use lib "/Perl/site/lib";
use lib '../lib/';

use strict;
use base ("manifest");
use manifest;

# This-module specific requirements
use Bio::SeqIO;
use Bio::SeqIO::fasta;

use PromoterToolWrapper;

########################################################
# constructor for this class
#
sub new
{
  my($class) = shift;
  my(%params) = @_;
  my($self) = manifest->new(@_);
  bless($self,$class);
  $self->{class} = $class;
  return ( $self, $class );
}

########################################################################
# override this function for GPS specific functionality
#
sub setup
{
  my $this = shift or die;
  $this->SUPER::setup();
}

sub main
{
	# first parameter is a hashref of the config
	my $this = shift or die;

	my $input = $this->{config}->{input};
	
	# make sure the file has a long filename so we
	# can parse out the well information (from the name)
	$input = ou::UnFixFileName($input);
	
	if ( not -f $input )
	{
	  print "$!\n";
	  return;
	}
	
	my $data  = ReadSequences($input);

	if ($this->{config}->{TSSG}->{include_test})
    {
      print "\t" . '*'x60 . "\n";
      print "Running TSSG promoter prediction tool on $input\n\n";
      
      my $promoter = PromoterToolWrapper->new("TSSG",
               $this->{config}->{TSSG}->{url},
                          $data, "tssg.html");
      $promoter->predict();
      PrintLinkLine( "tssg.html", 0, 0, 'TSSG' );
      #$promoter->show_result();
	}
	
	if ($this->{config}->{PROMOTER20}->{include_test})
    {
      # tom trying to implement in more generic way
      # (without using the wrapper)... not working yet 
      if (0)
      {
        my $url = $this->{config}->{PROMOTER20}->{url};
        my $submit = $this->{config}->{PROMOTER20}->{submit};
        #$submit->{'SEQSUB'} = $input;
        $submit->{seqpaste} = $data or die;

        my $htm = 'promoter20.html';
		$this->do_POST( $htm, $url, $submit );
        PrintLinkLine( $htm, 0, 0, 'PROMOTER20' );
      }
      else
      {
        print "\t" . '*'x60 . "\n";
        print "Running Promoter 2.0 promoter prediction Tool  on $input\n\n";
  
        my $promoter = PromoterToolWrapper->new("PROMOTER20",
                 $this->{config}->{PROMOTER20}->{url},
                        $data, "promoter20.html");
        $promoter->predict();
        PrintLinkLine( "promoter20.html", 0, 0, 'PROMOTER20' );
        #$promoter->show_result();
      }
	}
	return "\n";
}

sub ReadSequences
{
	my $file = shift or die;
	my $data;
	if ( ou::Open( *DAT, $file, '<' ) )
	{
      local $/ = undef;
      $data = <DAT>;
      close(DAT);
	}
	else
	{
		printf "Can't open input file $file : $!";
	}
	return $data;
}

sub config_hash
{
    return
    {
      'input' => '',
      
      'TSSG' =>
      {
        'include_test' => 1,
        'url' => 'http://www.softberry.com/cgi-bin/programs/promoter/tssg.pl',   
        'submit' =>
        {
          'FILE' => '',
          'SEQ_DATA' => '',
          'PL4_subgroup' => 'promoter',
          'PL4_topic' => 'tssg',
          'PL4_group' => 'programs'
        },
      },

      'PROMOTER20' =>
      {
          'include_test' => 1,
          'url' => 'http://www.cbs.dtu.dk/cgi-bin/nph-webface',
          #services/promoter/',
          'submit' =>
          {
            'SEQPASTE' => '',
            'configfile' => '/usr/opt/www/pub/CBS/services/Promoter-2.0/Promoter.cf',
            'SEQSUB' => '', # local file name
          },
          #{
          ##'SEQSUB' => '',# local file name
          #},
      },
    };
}

sub type_func
{
  my $this = shift or die;
  my $key  = shift;
  
  my $types =
  {
      'include_test' => 'checkbox',
      'input' => 'browsefile',
  };
  return $types->{$key} ? $types->{$key} : '';
}

sub sortkeys 
{
  {
      'input' => 0,
      'TSSG' => 1,
      'PROMOTER20' => 2,
      'include_test'  => 0,
  }
};

sub sort_func
{
    my ($this,$a,$b) = @_;
    
    my $k1 =  defined sortkeys()->{$a} ? sortkeys()->{$a} : 99999;
    my $k2 =  defined sortkeys()->{$b} ? sortkeys()->{$b} : 99999;
    return $k1 <=> $k2;
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


return 1;
