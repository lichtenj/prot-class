use strict;
use Bio::SeqIO;
use threads;
use threads::shared;
use Thread::Semaphore;

my $file = shift or die;
my $output = shift or die;

my $inseq = Bio::SeqIO->new( -file => $file, -format => 'fasta' );
my $in = $inseq->next_seq;

my $seq = $in->seq();
my $threshold = 2;
my $min_word_length = 2;
my %hash :shared;
my $counter :shared;
my $semaphore = Thread::Semaphore->new(10);

sub thread_search
{	

	my $length = $_[0];
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
	{
		lock($counter);
		$counter++;
		cond_signal($counter);
	}
	$semaphore->up();
}

$counter = $threshold;
for(my $length = $min_word_length; $length < length($seq); $length++)
{
	$semaphore->down();
	my $thr = threads->create('thread_search', $length);
	$thr->detach();
}

{
 	lock($counter);
	cond_wait($counter) until $counter == length($seq);
}
open(OUT, ">$output");
{

	foreach my $aa (sort {$hash{$b} <=> $hash{$a}} keys %hash)
	{
		print OUT $aa.','.$hash{$aa}."\n";
	}

}
close OUT;
