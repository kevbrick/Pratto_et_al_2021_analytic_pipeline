use strict;
use Bio::SeqIO;
use Math::Round;
use List::Util qw(sum shuffle);
use POSIX qw(ceil);
use Statistics::Descriptive;

my $fa  = $ARGV[0];
my $win = $ARGV[1]?$ARGV[1]:25;

my %bendability = ('AAA'=>0.1,'CAG'=>9.6,'AAC'=>1.6,'CCA'=>0.7,'AAG'=>4.2,
'CCC'=>5.7,'AAT'=>0.0,'CCG'=>3.0,'ACA'=>5.8,'CGA'=>5.8,'ACC'=>5.2,
'CGC'=>4.3,'ACG'=>5.2,'CTA'=>7.8,'ACT'=>2.0,'CTC'=>6.6,'AGA'=>6.5,
'GAA'=>5.1,'AGC'=>6.3,'GAC'=>5.6,'AGG'=>4.7,'GCA'=>7.5,'ATA'=>9.7,
'GCC'=>8.2,'ATC'=>3.6,'GGA'=>6.2,'ATG'=>8.7,'GTA'=>6.4,'CAA'=>6.2,
'TAA'=>7.3,'CAC'=>6.8,'TCA'=>10.0,'TTT'=>0.1,'CTG'=>9.6,'GTT'=>1.6,
'TGG'=>0.7,'CTT'=>4.2,'GGG'=>5.7,'ATT'=>0.0,'CGG'=>3.0,'TGT'=>5.8,
'TCG'=>5.8,'GGT'=>5.2,'GCG'=>4.3,'CGT'=>5.2,'TAG'=>7.8,'AGT'=>2.0,
'GAG'=>6.6,'TCT'=>6.5,'TTC'=>5.1,'GCT'=>6.3,'GTC'=>5.6,'CCT'=>4.7,
'TGC'=>7.5,'TAT'=>9.7,'GGC'=>8.2,'GAT'=>3.6,'TCC'=>6.2,'CAT'=>8.7,
'TAC'=>6.4,'TTG'=>6.2,'TTA'=>7.3,'GTG'=>6.8,'TGA'=>10.0);

## Also make a randomized set of bvals
my %bRand = %bendability;
@bRand{ keys %bRand } = shuffle values %bRand;
#
# my @k1 = keys(%bendability);
#
# for my $n(0..$#k1){
# 	print join("\t",$n,$k1[$n],$bendability{$k1[$n]},$rand{$k1[$n]})."\n";
# }
#
# die();

my $faIn = Bio::SeqIO->new(-file => $fa , '-format' => 'Fasta');

my $seqNum = 0;

while ( my $oSeq = $faIn->next_seq() ) {

	$seqNum++;
	my (@arr,@arrR);
	my $pos = round($win/2);

	my $s = uc($oSeq->seq);

	my $out;
	while ($s =~ s/^([GATCNX]{3})//){
		my $val = $bendability{$1};
		$val = $val?$val:0;
		push @arr, $val;

    ##Randoms
		my $valR = $bRand{$1};
		$valR = $valR?$valR:0;
		push @arrR, $valR;

	}

	my $nW = round($win/2);
	$pos = 0;
	for (my $i = $nW; $i < $#arr-$nW; $i+=$nW){
		my $d  = Statistics::Descriptive::Full->new();
		$d->add_data(@arr[$i-$nW..$i+$nW]);

		my $r  = Statistics::Descriptive::Full->new();
		$r->add_data(@arrR[$i-$nW..$i+$nW]);

		$pos += $nW;
	  print join("\t",$seqNum,$pos,$d->mean,$r->mean)."\n"
	}

}
