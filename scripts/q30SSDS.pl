#!/usr/bin/perl

use strict; 

while (<>){
	my @F = split(/\t/,$_);
	my @Q = split(/_/,$F[3]);
	print join("\t",@F) if ($Q[0] >= 30 && $Q[1] >= 30);
}
