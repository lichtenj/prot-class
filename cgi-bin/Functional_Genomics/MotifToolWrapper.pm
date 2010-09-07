# File	  : MotifToolWrapper.pm -- TFBS prediction class
# Author  : Jinfei Zhang, Ohio University
# Created : November 2005
# 
# Version : <none>
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#Description:
#To wrap up motif discovery web tools: AlignACE and MEME


package MotifToolWrapper;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
use HTML::Parse;
use LWP::Simple;

use strict;
my %fields = (
    class       => undef,
    tool_name   => undef,
    url         => undef,
    submit      => undef,
    input       => undef,
    output      => undef,
	     
    #parsed motif information
);

sub new
{
    $fields{'class'} = shift or die;
    $fields{'tool_name'} = shift or die;
    $fields{'url'} = shift or die;
    $fields{'submit'} = shift or die;
    $fields{'input'}=shift or die;
    $fields{'output'}=shift or die;
    return bless \%fields;
}

sub motif_search
{
    my $self = shift;
    if ($self->{'tool_name'} =~ m/AlignACE/) {
	$self->{'submit'}->{'seq'}=$self->{'input'};
    }
    elsif ($self->{'tool_name'} =~ m/MEME/) {
	$self->{'submit'}->{'datafile'}=[$self->{'input'}];
	    #[ou::FixFileName($self->{'input'}];
    }
    else {
	print "unknow tool\n";
	return 0;
    }

    $self->submit_search();
}

sub submit_search
{
    my $self = shift;
    my $agent = LWP::UserAgent->new( );

    my $resp = $agent->post($self->{'url'},$self->{'submit'} ,
			    'Content_Type' => 'form-data');
    
    if ( $resp->is_success )
    {
	if ( open(RESP, ">". $self->{'output'}) )
	{ 
	    print RESP $resp->content;
	    close RESP;
	    return 1;
	}
	else
	{
	    print $self->{'tool_name'} . " search " .
		$self->{'output'} . ": $!\n";
	}
    }
    else
    {
	print $self->{'tool_name'} . "Website submission failed\n";
    }
    return 0;
}

sub show_motif
{
    my $self = shift;
    if( $self->{'tool_name'} =~ m/AlignACE/) {
	if ( ou::Open( *DAT, $self->{'output'}, '<' ) )
	{
	    while ( <DAT> )
	    {
		# start showing lines but skip the first
		# text line
		if ( s/<\/TT><BR><TT>//gi )
		{
		    s/\&nbsp;/ /gi;
		    print;
		}
	    }
	    close(DAT);
	}
	else
	{
	    printf "Can't open " . $self->{'tool_name'} . " result file " . $self->{'output'} . ": $!\n";
	}
    }
    elsif ($self->{'tool_name'} =~ m/MEME/) {
	my $data_link = ou::FixFileName(extract_filename($self->{'output'}));
	
	if ( getstore($data_link, $self->{'output'}))
	{
	    if ( open(RESP, "<".$self->{'output'}) )
	    {
		my $text ="";
		while(<RESP>) {
		    $text .=$_;
		}
		if ($text eq "") {
		    print "No results are returned\n";
		}
		else {
		    (my $tmp,my $rtn) = split(/\<xmp\>/,$text);
		    if(!$rtn) {
			print "No results are returned\n";
		    }
		    else {
			($rtn,$tmp) = split(/\<\/xmp\>/,$rtn);
			print $rtn;
		    }
	        }
		close RESP;
	    }
	    else
	    {
		print $self->{'tool_name'} . " output file "
		    . $self->{'output'} . ": $!\n";
	    }
	}
	else
	{
	    print $self->{'tool_name'} . "Website submission failed\n";
	}
    }
}
sub extract_filename
{
    my $data_file = shift or die;
    
    if ( open(DAT, $data_file) )
    {
	local $/=undef;
	my $html = <DAT>;
	my @link;
	close(DAT);
	
	my $parsed_html = HTML::Parse::parse_html($html);
	
	for (@{ $parsed_html->extract_links( ) })
	{
	    push(@link, $_->[0]);
	}
	
	return $link[1];
    }
    else
    {
	print "First parsing pass failed\n";
    }
}
