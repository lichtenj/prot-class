# File	  : MotifFinder.pm -- TFBS prediction class
# Author  : Jinfei Zhang, Ohio University
# Created : Aug 2005
# Version : <none>
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.

package MotifFinder;
use Bio::SeqIO;
use Carp;
use strict;

my $debug_flag=0;

my %fields = (
    name            => undef,
    len             => undef,
    seq_list        => undef,
    seq_numb        => undef,
    motif           => undef,
    consensus_index => undef,
    type            => undef,
    alphabet_size   => undef,
    base_index      => undef,
    motif_positions => undef,
    opt_score       => undef,
);

my $base= {
    "RNA" => [ qw ( A C U G) ],
    "DNA" => [ qw ( A C T G) ],
};

# $seqobj->alphabet(); ==>dna, rna, protein

sub new
{
    print "constructing motif finder\n";
    
    # this stores the "class" name in the fields
    $fields{'name'} = shift or die;
        
    # retrieve named parameters... this checks for even sized parameter
    # lists, extra parameters are allowed but ignored for now
    my %param = @_;
    
    # now... check that we have all parameters which we need
    
    if(!defined($param{'len'}))
    {
        print STDERR "need parameter for the motif length: 'len'\n";
        return undef;
    }
    $fields{'len'}=$param{'len'};

    if (defined($param{'debug'}))
    {
        $debug_flag=$param{'debug'};
    }

    print "fasta file: $param{'fasta_file'}\n";
    
    if (defined($param{'fasta_file'}))
    {
        my $in  = Bio::SeqIO->new(-file => $param{'fasta_file'},
				  -format => 'Fasta');
        my $i=0;
        while ( my $seq = $in->next_seq() )
        {
            if($seq->alphabet eq 'dna')
            {
                $fields{'seq_list'}[$i]=$seq->seq(); 
                $fields{'seq_list'}[$i]=~tr/ACTG/actg/;
            		$fields{'seq_list'}[$i]=~tr/actg//cd;
                $i++;
            }
            else
            {
                print STDERR "all of the seqences should be DNA\n";
                return undef;
            }
        }   
        $fields{'seq_numb'}=$i;
    }
    $fields{'type'}='dna';
    $fields{'base_index'} = {'a'=>1, 'c'=>2, 't'=>3, 'g'=>4};
    $fields{'alphabet_size'} = 4;

    # return the class fields and the class name, which is also stored
    # in the fields.
    #return bless (\%fields, $fields{'name'});
    return bless \%fields;
}

#$seqobj->seq();              # string of sequence
#$seqobj->subseq(5,10);       # part of the sequence as a string


#########################################################################
# Branch and Bound Median String Problem
sub BBMS_search
{
    my $self=shift;
    my @motif_index;
    my $motif_index_ref;

#    print $fields->{"seq_list"}[2],"\n";

    for (my $i=0; $i < $self->{'len'}; $i++)
    {
        $motif_index[$i]=1;
    }
    $self->{'opt_score'}=$self->{'len'}*$self->{'seq_numb'};
    $motif_index_ref=\@motif_index;
    my $level=1;
    while ($level>0)
    {
        #print "level is $level\n";
        if ($level<$self->{'len'})
        {
            my ($dist) = $self->total_distance(\@motif_index,$level);
            #my ($dist) = $self->total_distance($motif_index_ref,$level);
            if ($dist > $self->{'opt_score'})
            {
                #($motif_index_ref, $level) = $self->bypass(\@motif_index,$level);
                ($motif_index_ref, $level) = $self->bypass($motif_index_ref,$level);
                @motif_index = @$motif_index_ref;
            }
            else
            {
                #($motif_index_ref, $level) = $self->next_node(\@motif_index,$level);
		($motif_index_ref, $level) = $self->next_node($motif_index_ref,$level);
                @motif_index = \@$motif_index_ref;
            }
        }
        else
        {
            my ($dist, $ref) = $self->total_distance(\@motif_index,$level);
            if ($dist<$self->{'opt_score'})
            {
                $self->{'opt_score'}=$dist;
                @{$self->{'consensus_index'}}=@motif_index;
                @{$self->{'motif_positions'}}=@$ref;
            }
            #($motif_index_ref, $level)=$self->next_node(\@motif_index,$level);
            ($motif_index_ref, $level)=$self->next_node($motif_index_ref,$level);
            @motif_index=\@$motif_index_ref;
        }
    }
    #print "@{$self->{'motif_positions'}} \n";
    $self->build_motif();
     
    return $self->consensus_string();
}

sub build_motif
{
    my $self=shift;
    for (my $i =0; $i<$self->{'seq_numb'}; $i++)
    {
        $self->{'motif'}[$i]=substr($self->{'seq_list'}[$i], 
		    $self->{'motif_positions'}[$i],
		    $self->{'len'});
    }
    return @{$self->{'motif'}};
}

sub show_motif
{
    my $self = shift or die;

    my $str = $self->consensus_string();
    if( !defined($self->{'motif'}))
    {
        print STDERR "the motif has not yet created\n";
        return undef;
    }
    
    my $rtn;
    $rtn =  "=========================================================\n";
    $rtn.=  "            The consensus string: $str\n";
    $rtn.=  "                          Length: $self->{'len'}\n";
    $rtn.=  "Score (minimum Hamming distance): $self->{'opt_score'}\n";
    $rtn.=  "=========================================================\n";
    $rtn.=  "The motif string matrix is \n";
    for (my $i=0; $i<$self->{'seq_numb'}; $i++)
    {
        $rtn.= "\t(start location: " . ($self->{'motif_positions'}[$i]+1).")\t";
        $rtn.= $self->{'motif'}[$i] . "\n";
    }
    $rtn .= "\n";
    print $rtn;
    return $rtn;
}

sub draw_logo {
    my $self = shift;
    if( !defined($self->{'motif'})) {
	print STDERR "the motif has not yet created\n";
	return undef;
    }
    my $picto = Bio::Graphics::Pictogram->new(-width=>"800",
					      -height=>"500",
					      -fontsize=>"60",
					      -plot_bits=>1,
					      -background=>{
						  'A'=>0.25,
						  'C'=>0.18,
						  'T'=>0.32,
						  'G'=>0.25},
					      -color=>{'A'=>'red',
						       'G'=>'blue',
						       'C'=>'green',
						       'T'=>'magenta'});
    my @seq;
    my $i=1;
    foreach my $seq (@{$self->{'motif'}}) {
	my $seqobj = Bio::Seq->new( -display_id => "id$i",
				    -seq => $seq);
	$i++;
	push @seq, $seqobj;
    }
    my $svg = $picto->make_svg(\@seq);
    
    return $svg->xmlify."\n";
}    
sub bypass
{
    my $self=shift;
    my ($indexRef, $level)= @_;
    for (my $j=$level; $j>=1;$j--)
    {
        if (${$indexRef}[$j-1]<$self->{'alphabet_size'})
        {
            ${$indexRef}[$j-1]+=1;
	    #return ($indexRef, $level);
	    return ($indexRef, $j);
        }
	#----
	#else {
	    

        #	else {
        #	    return ($indexRef, 0);
        #	}
    }
    return ($indexRef, 0);
}

sub next_node
{
    my $self = shift or die;
    my ($indexRef, $level)= @_;
    #print @$indexRef,"---",$level,"--next_node\n";
    if ($level<$self->{'len'})
    {
        ${$indexRef}[$level]=1;
        return ($indexRef,$level+1);
    }
    else
    {
        for (my $j=$self->{'len'}; $j>=1;$j--)
        {
            if (defined(${$indexRef}[$j-1])
				&& ${$indexRef}[$j-1]<$self->{'alphabet_size'})
            {
                ${$indexRef}[$j-1]+=1;
                return ($indexRef,$j);
            }
        }
    }
    return ($indexRef,0);
}

sub total_distance
{
    my $self = shift or die;
    my ($indexRef, $level)= @_;
    my $opt_dist=0;#$self->seq_numb;
    my $dist=0;
    my $seq_dist=0;
    my @start_positions;
    for (my $row=0; $row<$self->{'seq_numb'};$row++)
    {
        my $opt_seq_dist=$self->{'len'};
        #	print "length of Row $row is ",  length($self->{'seq_list'}[$row]), "\n";
        #	print "content of Row $row is ",  $self->{'seq_list'}[$row], "\n";
        my $col_bound=length($self->{'seq_list'}[$row])-$level+1;
        for (my $col=0; $col <$col_bound;$col++) {
        my $seq_dist=0;
        for(my $k=0; $k<$level;$k++)
        {
            # this "if defined" added by Tom... to prevent a warning
            # question: what does the "if" statement mean? (extra {}?)
            my $val = $self->{'base_index'} {substr($self->{'seq_list'}[$row],$col+$k,1)};
            if ( defined ( $val ) && defined ($indexRef->[$k]) )
            {
                if( $indexRef->[$k] != $val )
                {
                    $seq_dist+=1;
                }
            }
        }
        if($opt_seq_dist>$seq_dist)
        {
            $opt_seq_dist=$seq_dist;
            #if($level==$self->{'len'})
            {
                $start_positions[$row]=$col;
            }
        }
    } 
    $opt_dist+=$opt_seq_dist;
    }
    return ($opt_dist,\@start_positions);
}

sub consensus_string
{
    my $self=shift;
    my $rtn="";
    for (my $i=0; $i<$self->{'len'}; $i++)
    {
        my $index=$self->{'consensus_index'}[$i];
        while (my ($key, $val)= each(%{$self->{'base_index'}}))
        {
            #	print "$key => $val\n" ;
			if ( defined($val) && defined($index) )
			{
				if($val==$index)
				{
					$rtn.=$key;
				}
			}
        }
    }
    return $rtn;
}

1;
