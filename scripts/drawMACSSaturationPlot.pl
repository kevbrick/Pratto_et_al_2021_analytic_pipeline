use strict; 

my ($inTab, $outTab, $sample)  = @ARGV;

my $tf = 'satCurve_'.rand().rand().'.tab';

open TMP, '>', $outTab;
print TMP join("\t","reads","pc","hs")."\n";

open my $IN, '-|', 'sort -k1n,1n '.$inTab;

while (<$IN>){
	chomp;
	next if ($_ =~ /\stotal\s*$/);
	$_ =~ /^\s*(\d+).+\.N(\d+)_([\d\.]+)pc.+$/;
	my ($HS,$N,$pc) = ($1,$2,$3*100);
	print TMP join("\t",$N,$pc,$HS)."\n";
}

close TMP;
close $IN;

my $R = $tf.'.R';

makeRScript($R,$outTab,$sample);

system('R --no-save <'.$R);

sub makeRScript{
	my ($sName,$data,$sampleName) = @_;
	open RS, '>', $sName;
	print RS 'source("accessoryFiles/scripts/R/satCurveHS.R")'."\n";
	print RS 'satCurveHS(fIN = "'.$data.'", sampleName = "'.$sampleName.'")'."\n";
	close RS;
}
