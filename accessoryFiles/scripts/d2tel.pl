use strict;
use Math::Round;
use Getopt::Long;
use POSIX;

GetOptions ('gName=s' => \(my $gName),
            'g=s'     => \(my $fai),
            'i=s'     => \(my $bed));

unless ($fai && -e $fai){
    if ($gName){
        $fai = $ENV{GENOMES}."/".$gName."/genome.fa.fai";
    }
}

die ("Specify IDX file") unless ($fai && (-e $fai));

my %csLens;
open FAI, $fai;
while (<FAI>){
    chomp;
    my @F = split(/\t/,$_);
    $csLens{$F[0]} = $F[1];
}
close FAI;

open IN, $bed;
while (<IN>){
    chomp;
    my @F = split(/\t/,$_);
    my $pDist = $F[1];
    my $qDist = $csLens{$F[0]} - $F[2];
    print join("\t",@F,$pDist,$qDist)."\n";
}
close IN;
