=head1 NAME

AminoAcid - Amino Acid Information Repository

=head1 SYNOPSIS

    use AminoAcid;

    my $sequence = "MAGTAPASTTSAP";
    my $amino = "M";
    
    print AminoAcid::CalcpI($sequence);
    print AminoAcid::GetHydropathy($amino);
    print AminoAcid::monMW($sequence);
    print AminoAcid::avgMW($sequence);
    print AminoAcid::digest($sequence, $enz);

=head1 DESCRIPTION

The properties module supplies functions to determine the isoelectric point (pI Value)
of a given amino acid sequence, the hydropathy of a specific amino acid, the monoisotopic molecular weight and the average molecular
weight of a given sequence. It is also able to enzymatically digest a protein sequence into peptide fragments

=head1 FUNCTIONS

The following functions are exported by <AminoAcid>:

=over

=cut

package AminoAcid;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(GetHydropathy);

use Bio::SeqIO;

my $aa = AminoAcids();
my $h_mono = 1.008725;
my $oh_mono = 17.00274;

my $h_avg = 1.00794;
my $oh_avg = 17.00734;

=item monMW()

Calculates the monoisotopic molecular weight of an amino acid sequence, by adding the specific monoisotopic weight
of an amino acid to the base mass of 18.0105 da.

Arguments:

  $seq - Amino Acid Sequence

=cut

sub monMW
{
 my $seq = shift or die;
 my $mass = 18.0105;
 my $char;
 my @chars = split (//, $seq);
 foreach $char (@chars)
 {
  $mass = $mass + AminoAcids()->{$char}->{mon};
 }
 return ($mass);
}

=item avgMW()

Calculates the average molecular weight of an amino acid sequence, by adding the specific average weight
of an amino acid to the base mass of 18.0152 da.

Arguments:

  $seq - Amino Acid Sequence

=cut

sub avgMW {
 my $seq = shift or die;
 my $mass = 18.0152;
 my $char;
 my @chars = split (//, $seq);
 foreach $char (@chars)
 {
  $mass = $mass + AminoAcids()->{$char}->{avg};
 }
 return ($mass);
}

=item GetMass()

Return the molecular weight of an amino acid sequence - Old version

Arguments:

  $seq - Amino Acid sequence
  $type - Monoisotopic or Average molecular weight, average by default

=cut

sub GetMass
{
	my $seq = shift or die;
	my $type = shift || 'average';
	
	$type = $type=~/mono/i ? 'mon' : 'avg';
	
	my ($mass);
	while ($seq =~ m/([A-Z])/gi)
	{
		$mass += $aa->{$1}->{$type};
	}
	$mass += $type eq 'avg' ? $h_avg : $h_mono;
	$mass += $type eq 'avg' ? $oh_avg: $oh_mono;
	
	return $mass;
}


=item digest()

Retrieve the hydrophobicity of the amino acid

Arguments:

  $amino - Amino Acid

=cut

sub GetHydropathy
{
	my $amino  = shift;
	return 0 unless $aa->{$amino}->{hyd};
}

#
# According to I, Robot this will turn into artificial intelligence since it
# seems to be disconnected from the rest of the module
#

{
	my $probs;
	sub GetOmegaProb
	{
		my $tri_string = uc(shift) or die;
		$probs = OmegaProbTable() unless $probs;
		$tri_string =~ /(.)(.)(.)/;
		my $one = $probs->[0]->{$1};
		my $two = $probs->[1]->{$2};
		my $tre = $probs->[2]->{$3};
		
		my $result = 100 * $one * $two * $tre;
		#print "$result = 100 x $one x $two x $tre\n";
		return $result;	
	}
}

=item digest()

Digests a peptide sequence. It returns the monoisotopic mass of each fragment.

Arguments:

  $seq - Amino Acid Sequence
  $enz - Cleavage Enzyme, if not specified Trypsin is chosen as default

=cut

sub digest
{
	my $seq = shift or die;
	my $enz = shift || 'Trypsin';
	my @fragments;

	my $hr = enzymes();
	my $expr = $hr->{$enz}->{expr};
	
	my $start=0;
	while ($seq =~ m/$expr/gi)
	{
		my $end = pos($seq)-1;
		my $frag = substr($seq,$start,($end-$start));
		my $mass = GetMass($frag);
		
		if ( length $frag > 2 )
		{		
			print "$mass, $frag\n";
			#print "$mass\n";
		}
		#print "$mass\n";
		
		#push @fragments, substr($seq,$start,($end-$start));
		
		#my $frag = substr($seq,$start,($p-$start));
		$start = pos($seq) = $end;
	}
	# last fragment
	#push @fragments, substr($seq,$start);
	
	#print join "\n",@fragments;
}

=item AminoAcids()

Provides the Properties of each Amino Acids. It gives information about the molecular formula of the amino acid,
the monoisotopic and average molecular weight, the name and three letter code, the hydrophobicity and isoelectric
charge on c and n terminal and its side chain (if existant).

The hydrophobicities given are the "Scaled" values from computational log(P)
determinations by the "Small Fragment Approach" (see, "Development of
Hydrophobicity Parameters to Analyze Proteins Which Bear Post- or
Cotranslational Modifications" Black, S.D. and Mould, D.R. (1991)
Anal. Biochem. 193, 72-82).
The equation used to scale raw log(P) values to the scaled values given is
as follows: Scaled Parameters = (Raw Parameters + 2.061)/4.484 .

Arguments:


=cut

sub AminoAcids
{
	return
	{
		'A' => {
			'frm' => 'C3H5ON',
			'avg' => '71.0788',
			'tla' => 'Ala',
			'mon' => '71.03711',
			'nam' => 'Alanine',
			'hyd' => 0.616,
			'pKC' => 2.35,
			'pKN' => 9.87,
		},
		'C' => {
			'frm' => 'C3H5ONS',
			'avg' => '103.1388',
			'tla' => 'Cys',
			'mon' => '103.00919',
			'nam' => 'Cysteine',
			'hyd' => 0.680,
			'pKC' => 1.71,
			'pKN' => 8.33,
			'pKside' => 10.78,
		},
		'D' => {
			'frm' => 'C4H5O3N',
			'avg' => '115.0886',
			'tla' => 'Asp',
			'mon' => '115.02694',
			'nam' => 'AsparticAcid',
			'hyd' => 0.028,
			'pKC' => 1.88,
			'pKN' => 9.6,
			'pKside' => 3.65
		},
		'E' => {
			'frm' => 'C5H7O3N',
			'avg' => '129.1155',
			'tla' => 'Glu',
			'mon' => '129.04259',
			'nam' => 'GlutamicAcid',
			'hyd' => 0.043,
			'pKC' => 2.19,
			'pKN' => 9.67,
			'pKside' => 4.25
		},
		'F' => {
			'frm' => 'C9H9ON',
			'avg' => '147.1766',
			'tla' => 'Phe',
			'mon' => '147.06841',
			'nam' => 'Phenylalanine',
			'hyd' => 1.00,
			'pKC' => 2.58,
			'pKN' => 9.24,
		},
		'G' => {
			'frm' => 'C2H3ON',
			'avg' => '57.0519',
			'tla' => 'Gly',
			'mon' => '57.02146',
			'nam' => 'Glycine',
			'hyd' => 0.501,
			'pKC' => 2.34,
			'pKN' => 9.60,
		},
		'H' => {
			'frm' => 'C6H7ON3',
			'avg' => '137.1411',
			'tla' => 'His',
			'mon' => '137.05891',
			'nam' => 'Histidine',
			'hyd' => 0.165,
			'pKC' => 1.78,
			'pKN' => 8.97,
			'pKside' => 5.97,
		},
		'I' => {
			'frm' => 'C6H11ON',
			'avg' => '113.1594',
			'tla' => 'Ile',
			'mon' => '113.08406',
			'nam' => 'Isoleucine',
			'hyd' => 0.943,
			'pKC' => 2.27,
			'pKN' => 9.62,
		},
		'K' => {
			'frm' => 'C6H12ON2',
			'avg' => '128.1741',
			'tla' => 'Lys',
			'mon' => '128.09496',
			'nam' => 'Lysine',
			'hyd' => 0.283,
			'pKC' => 2.20,
			'pKN' => 8.90,
			'pKside' => 10.28,
		},
		'L' => {
			'frm' => 'C6H11ON',
			'avg' => '113.1594',
			'tla' => 'Leu',
			'mon' => '113.08406',
			'nam' => 'Leucine',
			'hyd' => 0.943,
			'pKC' => 2.36,
			'pKN' => 9.60,
		},
		'M' => {
			'frm' => 'C5H9ONS',
			'avg' => '131.1926',
			'tla' => 'Met',
			'mon' => '131.04049',
			'nam' => 'Methionine',
			'hyd' => 0.738,
			'pKC' => 2.28,
			'pKN' => 9.21,
		},
		'N' => {
			'frm' => 'C4H6O2N2',
			'avg' => '114.1038',
			'tla' => 'Asn',
			'mon' => '114.04293',
			'nam' => 'Asparagine',
			'hyd' => 0.236,
			'pKC' => 2.02,
			'pKN' => 8.80,
		},
		'P' => {
			'frm' => 'C5H7ON',
			'avg' => '97.1167',
			'tla' => 'Pro',
			'mon' => '97.05276',
			'nam' => 'Proline',
			'hyd' => 0.711,
			'pKC' => 1.99,
			'pKN' => 10.60,
		},
		'Q' => {
			'frm' => 'C5H8O2N2',
			'avg' => '128.1307',
			'tla' => 'Gln',
			'mon' => '128.05858',
			'nam' => 'Glutamine',
			'hyd' => 0.251,
			'pKC' => 2.17,
			'pKN' => 9.13,
		},
		'R' => {
			'frm' => 'C6H12ON4',
			'avg' => '156.1875',
			'tla' => 'Arg',
			'mon' => '156.10111',
			'nam' => 'Arginine',
			'hyd' => 0.000,
			'pKC' => 2.18,
			'pKN' => 9.09,
			'pKside' => 13.2,
		},
		'S' => {
			'frm' => 'C3H5O2N',
			'avg' => '87.0782',
			'tla' => 'Ser',
			'mon' => '87.03203',
			'nam' => 'Serine',
			'hyd' => 0.359,
			'pKC' => 2.21,
			'pKN' => 9.15,
			'pKside' => 12.87,
		},
		'T' => {
			'frm' => 'C4H7O2N',
			'avg' => '101.1051',
			'tla' => 'Thr',
			'mon' => '101.04768',
			'nam' => 'Threonine',
			'hyd' => 0.450,
			'pKC' => 2.15,
			'pKN' => 9.12,
			'pKside' => 12.98,
		},
		'V' => {
			'frm' => 'C5H9ON',
			'avg' => '99.1326',
			'tla' => 'Val',
			'mon' => '99.06841',
			'nam' => 'Valine',
			'hyd' => 0.825,
			'pKC' => 2.32,
			'pKN' => 9.62,
		},
		'W' => {
			'frm' => 'C11H10ON2',
			'avg' => '186.2132',
			'tla' => 'Trp',
			'mon' => '186.07931',
			'nam' => 'Tryptophan',
			'hyd' => 0.878,
			'pKC' => 2.38,
			'pKN' => 9.39,
		},
		'Y' => {
			'frm' => 'C9H9O2N',
			'avg' => '163.1760',
			'tla' => 'Tyr',
			'mon' => '163.06333',
			'nam' => 'Tyrosine',
			'hyd' => 0.880,
			'pKC' => 2.20,
			'pKN' => 9.11,
			'pKside' => 10.07,
		},
	}
}

=item enzymes()

Describes the properties of cleavage enzymes. It provides cleavage sites and c and n terminal additions.

=cut

sub enzymes
{
	return
	{
		'Arg-C' => {
			'expr' => '[R][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Subtilisin' => {
			'expr' => '[^RHK][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Hydroxylamine' => {
			'expr' => '[N][G]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Cathepsin B' => {
			'expr' => '[R][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Cathepsin D' => {
			'expr' => '[LF][^VAG]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Chymotrypsin' => {
			'expr' => '[YWFL][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'CNBr_HSer' => {
			'expr' => '[M][A-Z]',
			'cterm' => 'H2O2-SCH3',
			'nterm' => 'H'
		},
		'Bromelain' => {
			'expr' => '[KAY][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Glu-C_Bic' => {
			'expr' => '[E][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Lys-N' => {
			'expr' => '[A-Z][K]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Thermolysin' => {
			'expr' => '[LFIVMA][^P]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Asp-N' => {
			'expr' => '[A-Z][D]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Glu-C_Phos' => {
			'expr' => '[ED][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Lys-C' => {
			'expr' => '[K][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Pepsin' => {
			'expr' => '[LF][^VAG]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Trypsin' => {
			'expr' => '[KR][^P]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Clostripain' => {
			'expr' => '[R][P]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Proteinase K' => {
			'expr' => '[YWF][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Papain' => {
			'expr' => '[RK][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Elastase' => {
			'expr' => '[AVLIGS][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'Cathepsin G' => {
			'expr' => '[YWF][A-Z]',
			'cterm' => 'OH',
			'nterm' => 'H'
		},
		'CNBr_HSerLac' => {
			'expr' => '[M][A-Z]',
			'cterm' => 'O-SCH3',
			'nterm' => 'H'
		},
	}
}

=item OmegaProbTable()



Arguments:



=cut

sub OmegaProbTable
{
		my @aa_permis;

	  $aa_permis[0]->{"A"} = "0.005";
	  $aa_permis[0]->{"R"} = "0.01";#ND
	  $aa_permis[0]->{"N"} = "0.19";
	  $aa_permis[0]->{"D"} = "0.13";
	  $aa_permis[0]->{"C"} = "0.04";
	  $aa_permis[0]->{"Q"} = "0.01";#0 selon Kodukula
	  $aa_permis[0]->{"E"} = "0.0";
	  $aa_permis[0]->{"G"} = "0.12";
	  $aa_permis[0]->{"H"} = "0.0";#ND
	  $aa_permis[0]->{"I"} = "0.005";#absent
	  $aa_permis[0]->{"L"} = "0.0";
	  $aa_permis[0]->{"K"} = "0.01";
	  $aa_permis[0]->{"M"} = "0.00";
	  $aa_permis[0]->{"F"} = "0.000";#absent
	  $aa_permis[0]->{"P"} = "0.000";
	  $aa_permis[0]->{"S"} = "0.430";
	  $aa_permis[0]->{"T"} = "0.010";
	  $aa_permis[0]->{"W"} = "0.000";
	  $aa_permis[0]->{"Y"} = "0.000";
	  $aa_permis[0]->{"V"} = "0.005";
	  $aa_permis[0]->{"X"} = "0.000"; # aa non défini
	  $aa_permis[0]->{"B"} = "0.320";  # B = N or D    B=Asx =(0.8+0.4)/2
	  $aa_permis[0]->{"Z"} = "0.010"; # Z = Q or E    Z=Glx =(0.1+0.0)/2
	  $aa_permis[0]->{"U"} = "0.04";  # U = sélénocystéine (très rare : 74 dans tout swissprot. Associé à la Cystéine =C)

	  $aa_permis[1]->{"A"} = "0.180";
	  $aa_permis[1]->{"R"} = "0.02";#ND
	  $aa_permis[1]->{"N"} = "0.04";#+
	  $aa_permis[1]->{"D"} = "0.02";
	  $aa_permis[1]->{"C"} = "0.01";
	  $aa_permis[1]->{"Q"} = "0.01";#ND
	  $aa_permis[1]->{"E"} = "0.01";#+
	  $aa_permis[1]->{"G"} = "0.20";#+
	  $aa_permis[1]->{"H"} = "0.02";#ND
	  $aa_permis[1]->{"I"} = "0.01";#absent
	  $aa_permis[1]->{"L"} = "0.01";#ND
	  $aa_permis[1]->{"K"} = "0.00";#ND
	  $aa_permis[1]->{"M"} = "0.0";
	  $aa_permis[1]->{"F"} = "0.01";#absent
	  $aa_permis[1]->{"P"} = "0.01";
	  $aa_permis[1]->{"S"} = "0.29";
	  $aa_permis[1]->{"T"} = "0.02";
	  $aa_permis[1]->{"W"} = "0.005";
	  $aa_permis[1]->{"Y"} = "0.0";#ND
	  $aa_permis[1]->{"V"} = "0.04";#ND
	  $aa_permis[1]->{"X"} = "0.0"; # aa non défini
	  $aa_permis[1]->{"B"} = "0.06"; # B = N or D    B=Asx =(0.1+0.4)/2
	  $aa_permis[1]->{"Z"} = "0.02";  # Z = Q or E    Z=Glx =(0.0+0.0)/2
	  $aa_permis[1]->{"U"} = "0.01";  # U = sélénocystéine (très rare : 74 dans tout swissprot. Associé à la Cystéine =C)

	  $aa_permis[2]->{"A"} = "0.3";
	  $aa_permis[2]->{"R"} = "0.02";#ND
	  $aa_permis[2]->{"N"} = "0.0";#ND
	  $aa_permis[2]->{"D"} = "0.01";
	  $aa_permis[2]->{"C"} = "0.0";
	  $aa_permis[2]->{"Q"} = "0.005";#ND
	  $aa_permis[2]->{"E"} = "0.0";
	  $aa_permis[2]->{"G"} = "0.17";
	  $aa_permis[2]->{"H"} = "0.005";
	  $aa_permis[2]->{"I"} = "0.02";#absent
	  $aa_permis[2]->{"L"} = "0.06";#ND selon Kodukula
	  $aa_permis[2]->{"K"} = "0.005";#ND
	  $aa_permis[2]->{"M"} = "0.10";#ND
	  $aa_permis[2]->{"F"} = "0.01";#absent
	  $aa_permis[2]->{"P"} = "0.0";
	  $aa_permis[2]->{"S"} = "0.20";
	  $aa_permis[2]->{"T"} = "0.04";
	  $aa_permis[2]->{"W"} = "0.0";
	  $aa_permis[2]->{"Y"} = "0.0";#ND
	  $aa_permis[2]->{"V"} = "0.01";
	  $aa_permis[2]->{"X"} = "0.0"; # aa non défini
	  $aa_permis[2]->{"B"} = "0.01"; # B = N or D    B=Asx =(0.1+0.0)/2
	  $aa_permis[2]->{"Z"} = "0.005"; # Z = Q or E    Z=Glx =(0.0+0.0)/2
	  $aa_permis[2]->{"U"} = "0.0";  # U = sélénocystéine (très rare : 74 dans tout swissprot. Associé à la Cystéine =C)
		
		return \@aa_permis;
}

=item CalcpI()

Calculates the isoelectric point of an amino acid sequence. It adds the pk c and n terminal
values of an amino acids together and in certain cases adds the pk of a side chain.

Arguments:

  $seq - Amino Acid Sequence

=cut

sub CalcpI
{
 my $seq = shift or die;
 my $char;
 my @chars = split (//, $seq);
 my %count;
 my $pka_side;
 my $cr;
 my $c;
 my $cn;
 my $charge;
 my $first = $chars[0];
 my $last = $chars[$#chars];
 my $ph;
 my $dph;
 my $reached = 0;
 my $pI = 0;
 
 #print AminoAcid::AminoAcids->{$first}->{pKC};
 
 my $pk_c = 0 - AminoAcids()->{$first}->{pKC};
 my $pk_n = AminoAcids()->{$first}->{pKN};
 
 foreach $char (@chars)
 {
  $count{$char}++;
 }
 
# print "$first, $pk_n, $last, $pk_c \n";

# begin pH loop from 0 to 14

 for ($ph = 0; $ph<14.5; $ph = $ph+0.1)
 {
  $charge = 0;
  
  $c = 10**($pk_n - $ph);
  $cr = $c/(1+$c);
  $cn = (+1) * $cr;
  # printf ("N -> %.2f ",$cn);
  $charge = $charge + $cn;
   
  $c = 10**($pk_c+$ph);
  $cr = $c/(1+$c);
  $cn = (-1) * $cr;
  # printf ("C -> %.2f \n",$cn);
  $charge = $charge + $cn;
  
  foreach $char (keys (%count)) {
     
   if ($char =~ m/K|R|H/)
   {
    $pka_side = AminoAcids()->{$char}->{pKside};
   }
   else {
    $pka_side = 0 - (AminoAcids()->{$char}->{pKside} ? AminoAcids()->{$char}->{pKside} : 0);
   }
  
   if ($pka_side > 0)
   {
    $c = 10**($pka_side-$ph);
    $cr = $c/(1+$c);
    $cn = (+1) * $count{$char} * $cr;
   }
   elsif ($pka_side < 0)
   {
    $c = 10**($pka_side+$ph);
    $cr = $c/(1+$c);
    $cn = (-1) * $count{$char} * $cr;
   }
   else {
    $cn = 0;
   }
   $charge = $charge + $cn;
  }
  
  if (($charge <= 0) && ($reached == 0))
  {
   $pI = $ph - 0.05;
   $reached = 1;
  }     
 } 

 return ($pI);
}

return 1;

=back

=head1 BUGS

<AminoAcid> currently only supports a subset of all possible amino acid properties. Also <Properties>
should support nucleotide input sequences and translate them or reject them.

=head1 AUTHOR

Jens Lichtenberg <lichtenj@ohio.edu>
Tom Conley <conleyt@ohio.edu>

=head1 COPYRIGHT

Copyright (c) 2006. CIDDS. All rights reserved.

=cut