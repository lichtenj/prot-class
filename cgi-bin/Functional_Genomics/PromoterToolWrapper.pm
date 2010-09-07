# File	  : PromoterToolWrapper.pm -- promoter prediction class
# Author  : Ohio University
# Created : July 2006
# 
# Version : <none>
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#Description:
#To wrap up promoter prediction web tools: TSSG and Promoter 2.0


package PromoterToolWrapper;
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
    #$fields{'submit'} = {};#shift or die;
    $fields{'input'}=shift or die;
    $fields{'output'}=shift or die;
    return bless \%fields;
}

sub predict
{
    my $self = shift;
    $self->{'submit'} = {}; 
    if ($self->{'tool_name'} =~ m/TSSG/)
    {
        #$self->{'submit'}->{'FILE'}=
        $self->{'submit'}->{'PL4_group'}   ='programs'; #hidden input
        $self->{'submit'}->{'PL4_subgroup'}='promoter'; #hidden input
        $self->{'submit'}->{'PL4_topic'}   ='tssg'; #hidden input
        $self->{'submit'}->{'SEQ_DATA'}   =$self->{'input'};
    }
    elsif ($self->{'tool_name'} =~ m/PROMOTER20/)
    {
        $self->{'submit'}->{'configfile'}=
            "/usr/opt/www/pub/CBS/services/Promoter-2.0/Promoter.cf";
        #$self->{'submit'}->{'SEQPASTE'}=$self->{'input'};
        $self->{'submit'}->{'SEQSUB'}=$self->{'input'};#"promoter.fa";
    }
    else
    {
        print "unknown tool\n";
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
	    if($self->{'tool_name'} =~ m/PROMOTER20/) {
		my $tmp=$resp->content;
		#print $tmp . "\n";
		$tmp =~ s/(\n|\r)//g;
		#print $&;
		#$tmp=$&;
		#print $tmp . "\n";
		$tmp =~ /http.(\s|\S)+wait">/;
                if($&=='') {
                  print RESP $resp->content;
                  close RESP;
                  return 1;
                }
		my $newurl = substr($&,0,length($&)-2);
                #print $newurl . "\n";
                $resp=$agent->post($newurl);
                if($resp->is_success) {
                        print RESP $resp->content;
                        close RESP;
#print $resp->content;
                        return 1;
                }
                else {
                  	print $self->{'tool_name'} . 
                             " Website submission failed\n";
                        close RESP;
                        return 0;
                }
            }
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
	print $self->{'tool_name'} . " Website submission failed\n";
    }
    return 0;
}

sub show_result
{
    my $self = shift;
 	if ( ou::Open( *DAT, $self->{'output'}, '<' ) )
	{
	    while ( <DAT> )
	    {
		print;
		
	    }
	    close(DAT);
	}
	else
	{
	    printf "Can't open " . $self->{'tool_name'} . " result file " . $self->{'output'} . ": $!\n";
	}
    return ;
    if( $self->{'tool_name'} =~ m/PROMOTER/) {
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
