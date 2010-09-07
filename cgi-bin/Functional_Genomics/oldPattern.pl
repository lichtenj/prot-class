use strict;
use Bio::SeqIO;

my $file = shift or die;
my $output = shift or die;

my $inseq = Bio::SeqIO->new( -file => $file, -format => 'fasta' );
my $in = $inseq->next_seq;

my $seq = $in->seq();
my $min_word_length = 2;
my $threshold = 2;
my %hash;

for(my $length = $min_word_length; $length < length($seq); $length++)
{
        for(my $i = 0; $i + $length < length($seq); $i++)
        {
                if($hash{substr($seq,$i,$length)})
                {
                        $hash{substr($seq,$i,$length)}++;
                }
                else
                {
                        $hash{substr($seq,$i,$length)} = 1;
                }
       	}
       	foreach my $key (keys %hash)
			{
				if($hash{$key} == 1)
				{
					delete($hash{$key})
				}
			}
}

open(OUT, ">$output");

foreach my $aa (sort {$hash{$b} <=> $hash{$a}} keys %hash)
{
   if($hash{$aa} < $threshold)
   {
   	last;
	}
	else
	{
		print OUT $aa.','.$hash{$aa}."\n";
	}
}

close OUT;
