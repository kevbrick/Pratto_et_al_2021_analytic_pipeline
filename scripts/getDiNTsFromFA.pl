use strict;
use Bio::SeqIO;

my $fa   = new Bio::SeqIO(-file => $ARGV[0], -format => 'fasta');
my $diNT = $ARGV[1]?$ARGV[1]:"CG";

print STDERR "Finding all $diNT\'s in $ARGV[0] ...\n";

while (my $s = $fa->next_seq){
	my $seq = $s->seq;

	my $len = $s->length;
	my $chrom = $s->display_id;

	my $nextpos = 0;

	for (my $n = 0; $n < $len; $n += 1000000){

		my $nTo     = ($n+1000000);
		my $subseq  = $s->subseq($n+1,$nTo>$len?$len:$nTo);

		while ($subseq =~ s/^(.*?)($diNT)//i){
			my $pos     = $nextpos + length($1);
			print join("\t",$chrom,$pos,$pos+2)."\n";
			$nextpos = $pos+2;
		}

		if ($nTo<$len){
			my $edgeseq = $s->subseq($nTo,$nTo+1);
			if ($edgeseq eq $diNT){
				print join("\t",$chrom,$nTo-1,$nTo+1)."\n";
			}
		}

		$nextpos += length($subseq);
	}
}
