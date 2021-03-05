#!/usr/bin/perl
use strict; 

system('awk \'{close(f);f=$1}{print > f".bed"}\' '.$ARGV[0]);

my $catCmd = 'cat ';
open CS, $ARGV[1];

while (<CS>){
	chomp; 
	my @F = split(/\t/,$_);
	if (-e "$F[0]\.bed"){
		system("sort -k1,1 -k2n,2n -k3n,3n $F[0]\.bed >$F[0].sorted.bed");
		$catCmd .= " $F[0].sorted.bed";
	}else{
		print STDERR "No data for $F[0]\.bed ...\n";
	}
}

close CS; 

system($catCmd);
