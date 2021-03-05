#!/usr/bin/perl

use strict;

my $inFA = $ARGV[0];
my $out = $ARGV[1]?$ARGV[1]:'out.qp.bed';

my $cmd = "/home/kevbrick/quadparser2/quadparser2 -n -d ".$inFA." stdout";

open my $PIPE, '-|', $cmd; 

system('mkdir -p ./.tmpQP');
my $tf = './.tmpQP/qptmp'.rand().rand().'.bed';
open TMP, '>', $tf;

my $chrom;
while (<$PIPE>){
	chomp; 
	if ($_ =~ /^\s*(chr\S+)\s*.*$/){
		$chrom = $1; 
	}else{
		my ($coord, $type, $seq) = split("\t",$_);
		my ($start,$end) = split("-",$coord);
		
		my $G += ($seq =~ tr/[Gg]//);
		my $C += ($seq =~ tr/[Cc]//);
		
		my $strand = ($G>$C)?"+":"-";
		
		next unless ($start =~ /^\d+$/);
		next unless ($end =~ /^\d+$/);
		
		print TMP join("\t",$chrom, $start, $end, $seq, $type, $strand)."\n" if ($chrom);
	}
}

system('sort -k1,1 -k2n,2n -k6,6 '.$tf.' >'.$out);
system('rm '.$tf);
