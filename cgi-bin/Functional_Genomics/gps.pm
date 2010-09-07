#!/usr/bin/perl -w
#
# File    : gps.pm -- Derived class for Bioinformatics Experiments
# Author  : Ohio University CNS
# Created : July 2005
# Version : <none>
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.

# First create a module derived from "manifest"
package gps; # 

use lib "/Perl/site/lib";
use lib '../lib/';

use strict;
use base ("manifest");
use manifest;

use LWP::UserAgent;
use HTTP::Request::Common qw(POST);

# This-module specific requirements
use Bio::SeqIO;
use Bio::SeqIO::fasta;

use MotifFinder;
use MotifToolWrapper;

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

#sub gps_main
#{
#   # first parameter is a hashref of the config
#   my $this = shift or die;
#
#   my $input = $this->{config}->{input};
#   my $length = $this->{config}->{'Branch_Bound'}->{length};
#   
#   -f $input or return "$!\n";
#   
#   print "Running Motif Finder on $input\n";
#   
#   my ($motif) = MotifFinder->new(fasta_file=>$input, len=>$length);
#   if (!defined($motif))
#   {
#       return "Cannot make the motif object, try to check the input parameters\n";
#   }
#   $motif->BBMS_search();
#   return $motif->show_motif();
#}

sub main
{
    # first parameter is a hashref of the config
    my $this = shift or die;

    my $input = $this->{config}->{input};

    #return Data::Dumper->Dump([$this->{config}]);
    
    # whatever this fuction returns will be displayed
    # in the ui window
    #return "Called gps Main\n";
    
    -f $input or return "$!\n";
    if($this->{config}->{Branch_Bound}->{include_test})
    {
        print "Running Branch And Bound Algorithm on $input\n";
        
        my $length = $this->{config}->{'Branch_Bound'}->{length};
        my $motif = MotifFinder->new(fasta_file=>$input, len=>$length);
        # my $motif = MotifFinder->new(fasta_file=>$input, len=>$length);
        if (!defined($motif))
        {
          print "\tCannot make the motif object, try to check the input parameters\n";
        }
        else
        {
          print "Starting search...\n";
          $motif->BBMS_search();

          print "Showing results...\n";
          $motif->show_motif();
        }
    }
    if ($this->{config}->{AlignACE}->{include_test}) {

        print "\t" . '*'x60 . "\n";
        print "\t" . '*'x60 . "\n";
        print "Running AlignACE Motif Tool on $input\n";
        
        my $data = ReadSequences($input);
        my $alignACE = MotifToolWrapper->new("AlignACE",
                 $this->{config}->{AlignACE}->{url},
                 $this->{config}->{AlignACE}->{submit},
                         $data, "AlignACEmotif.html");
        $alignACE->motif_search();
        $alignACE->show_motif();
    }
    #if ($this->{config}->{MEME}->{include_test}) {
    #    print "\t" . '*'x60 . "\n";
    #    print "\t" . '*'x60 . "\n";
    #    print "Running MEME Motif Tool  on $input\n";
    #
    #    my $data = ReadSequences($input);
    #    my $meme = MotifToolWrapper->new("MEME",
    #             $this->{config}->{MEME}->{url},
    #             $this->{config}->{MEME}->{submit},
    #                     $input, "MEMEmotif.html");
    #    $meme->motif_search();
    #    $meme->show_motif();
    #    
    #}
    if ($this->{config}->{MEME}->{include_test}) {
      print "\t" . '*'x60 . "\n";
      print "\t" . '*'x60 . "\n";
      print "Running MEME Motif Tool  on $input\n";
    
      # my $data = ReadSequences($input);
      open MEMEPROC, "meme.bin.exe $input |" or die "cannot pipe from meme: $!";
      while(<MEMEPROC>) {
        chomp;
        print $_."\n";
      }
    }
}

sub ReadSequences
{
    my $file = shift or die;
    my $data;
    if ( ou::Open( *DAT, $file, '<' ) )
    {
        while ( <DAT> )
        {
            $data .= $_;
        }
        close(DAT);
    }
    else
    {
        printf "Can't open input file $file : $!";
    }
    return $data;
}

#sub ShowAlignAceResults
#{
#   my $file = shift or die;
#   my $data;
#   my $showlines = 0;
#   if ( ou::Open( *DAT, $file, '<' ) )
#   {
#       while ( <DAT> )
#       {
#           # start showing lines but skip the first
#           # text line
#           if ( s/<\/TT><BR><TT>//gi )
#           {
#               s/\&nbsp;/ /gi;
#               print;
#           }
#       }
#       close(DAT);
#   }
#   else
#   {
#       printf "Can't open AlignAce result file $file : $!";
#   }
#   return $data;
#}


#Submits configuration data to mass spec websites
#sub AlignAce
#{
#   my $this = shift or die;
#   my $data = shift or die;
#   
#   # submit to AlignAce site
#   {
#       my $url    = $this->{config}->{AlignAce}->{url};
#       my $submit = $this->{config}->{AlignAce}->{submit};
#
#       # store the data for submission to the site
#       $submit->{seq} = $data;
#
#       print "\n" . '='x60 . "\n";
#       print "Motif searching...\n";
#       print "URL : $url\n";
#   
#       my $htm = 'motif.htm';
#       $this->do_POST( $htm, $url, $submit );
#       ShowAlignAceResults($htm);
#       return $htm;  
#   }
#}


sub config_hash
{
    return
    {
        'input' => '../dat/GPS/gps5.fa',

        'Branch_Bound' =>
        {
        'include_test' => 1,
        'length' => 15,
        },

        'AlignACE' =>
        {
        'include_test' => 1,
        'url' => 'http://atlas.med.harvard.edu/cgi-bin/alignace.pl',
        'submit' =>
        {
            'xp' => '15',                                       # number of sites to expect
            'fname' => '',                                  # file name when submtting whole files
            'name' => 'description',                # sequence description
            'gc' => '0.38',                                 # fractional background GC content
            #'seq' => 'Fasta paste here',
            'nc' => '10',                                       # number of columns to align
        },
        },

        'MEME' =>
        {
        'include_test' => 1,
        #'reference' => '[the MEME paper]',
        'url' => 'http://srs.nchc.org.tw/emboss-bin/emboss.pl',
        'submit' =>
        {
            'minw'      => '8',
            'page'      => '80',
            'distance'  => '1e-3',
            'model'     => 'oops',
            'timer'     => '0',
            'maxsize'   => '100000',
            'shorten'   => 'yes',
            'trace'     => 'no',
            'zprint'    => 'no',
            'outfile'   => 'outfile',
            'fastaprint'    => 'no',
            'palindromes'   => 'no',
            'wthree'    => 'yes',
            'maxw'      => '57',
            '_email'    => '',
            'maxsites'  => '0.',
            'llprint'   => 'no',
            'cfive'     => 'yes',
            'align'     => 'yes',
            'v'         => 'no',
            'cons'      => '',
            'ponly'     => 'no',
            '_app'      => 'meme',
            'b'         => '-1.0',
            'prior'     => 'dirichlet',
            'minsites'  => '0.',
            'allprint'  => 'no',
            'w'         => '0',
            'cthree'    => 'yes',
            'status'    => 'no',
            'brief'     => 'yes',
            'seed'      => '0',
            'chi'       => '1.0',
            'maxiter'   => '50',
            '_action'   => 'run',
            'prob'      => '1.0',
            'wprint'    => 'no',
            'spfuzz'    => '-1.0',
            'nsites'    => '0.',
            'nucleic'   => 'yes',
            'seqfrac'   => '1.0',
            'protein'   => 'no',
            'adj'       => 'root',
            'nmotifs'   => '1',
            'spmap'     => 'uni',
            'startsprint'   => 'no',
            'ntype'     => 'pair',
        },
        },
        
    };
}

#sub do_POST
#{
#   my $this = shift or die;
#   my $file = shift or die; # for results
#
#   my $browser = LWP::UserAgent->new( );
#
#   # Parameters:
#   #  the URL,
#   #  an arrayref or hashref for the key/value pairs,
#   #  and then, optionally, any header lines: (key,value, key,value)
#   $this->{browser} = LWP::UserAgent->new( ) unless $this->{browser};
#   my $resp = $browser->post(@_ ,'Content_Type' => 'form-data');
#   
#   if ( $resp->is_success )
#   {
#       if ( open(RESP, ">$file") )
#       { 
#           print RESP $resp->content;
#           close RESP;
#           return 1;
#       }
#       else
#       {
#           print "do_POST $file: $!\n";
#       }
#   }
#   else
#   {
#       print "Website submission failed\n";
#   }
#   return 0;
#}

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
      'Branch_Bound' => 1,
      'AlignACE' => 2,
      'MEME' => 3,
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


return 1;
