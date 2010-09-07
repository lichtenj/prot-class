use strict;
use Bio::SeqIO;
use LWP::UserAgent;

my $input = shift or die;
my $gene = shift or die;

my $inseq = Bio::SeqIO->new( -file => $input, -format => 'fasta' );
my $in = $inseq->next_seq;

WebCheckGPI($in, $gene);

sub WebCheckGPI
{
	my $seq = shift or die;
	my $gene = shift or die;

	my $URL = 'http://mendel.imp.ac.at/gpi/cgi-bin/gpi_pred.cgi';
	my $browser = LWP::UserAgent->new( );
	my $to_tr = $seq->seq();

	my %parms;
	$parms{Sequence} = $to_tr;
	$parms{"LSet"} = 'metazoa';

	my $doc = do_POST( $URL, \%parms );
	$doc =~ s/<IMG[^>]*>//g;
	$doc =~ s/<A [^>]*>//g;
	$doc =~ s/<\/[Aa].*>//g;

	open(OUT, ">".$gene.".gpi") or die "Cannot open";
	print OUT $doc;
	close OUT;
}

sub do_POST
{
	my $browser = LWP::UserAgent->new( );

	my $resp = $browser->post(@_);

	return ($resp->content, $resp->status_line, $resp->is_success, $resp) if wantarray;
	return unless $resp->is_success;
	return $resp->content;
}
