#!/usr/bin/perl
use strict;
use Getopt::Long;
use Math::Round;
use List::Util qw(sum);

GetOptions ('s=s'          => \(my $statFile),
	    'g=s'          => \(my $genome),
            'run+'         => \(my $runMe),
            'n=i'          => \(my $ncells = 300),
            't=i'          => \(my $totTime = 7500),
            'justgetname+' => \(my $justGetName),
            'x=s'          => \(my $excludeCS),
            'bad+'         => \(my $getBADmodel),
            'bgpath=s'     => \(my $bgPath  = 'bgs'),
            'oripath=s'    => \(my $oriPath = 'ori'));

die('Please provide a genome ... (--g)') unless ($genome);
####
my (@RT,@ORI,%RTnames);
doLS($bgPath,$oriPath);
####

my ($time,$pcRep,$R2,$score,$oriPerMb,$recycle,$repTime,$bg,$ori,$useStr,$recycling);

if ($getBADmodel){
	($time,$pcRep,$R2,undef,undef,undef,$score,$oriPerMb,$recycle,$repTime,$bg,$ori,$useStr,$recycling,undef,undef,undef) = split(/\t/,`grep -v RMSE $statFile |tail -n5000 |head -n1`);
}else{
	($time,$pcRep,$R2,undef,undef,undef,$score,$oriPerMb,$recycle,$repTime,$bg,$ori,$useStr,$recycling,undef,undef,undef) = split(/\t/,`grep -v RMSE $statFile |head -n1`);
}

$bg  =~ s/^\S+\///;
$ori =~ s/^\S+\///;

my $rtBG  = $RTnames{$bg};
die("BG : $rtBG") unless (-e $rtBG);

my $oriBG  = $RTnames{$ori};
die("ORI [$ori]: $oriBG") unless (-e $oriBG);

my $name = $bg;
$name =~ s/^(.+?)\..+bedgraph/$1/;
$name =~ s/^RT_//;
$name =~ s/_(hg38|mm10).forModel.bedgraph//;

if ($justGetName){
	print $name;
	exit();
}

my $excl = " ";
if ($excludeCS){
	$excl = " -v $excludeCS";
}

my $cmd;
if ($getBADmodel){
	$cmd = 'Rscript '.$ENV{"RTSCRIPTS"}.'/run_allCSmodel.R -u '.$name.' -i BADMODEL_'.$name.' -t '.$totTime.' -c '.$ncells.' -g '.$genome.' -n '.$oriPerMb.' -b '.$rtBG.' -o '.$oriBG.' -e '.($time*10).' -f '.$pcRep.' -m TRUE -s '.$useStr.' '.$excl;
}else{
	$cmd = 'Rscript '.$ENV{"RTSCRIPTS"}.'/run_allCSmodel.R -u '.$name.' -i '.$name.' -t '.$totTime.' -c '.$ncells.' -g '.$genome.' -n '.$oriPerMb.' -b '.$rtBG.' -o '.$oriBG.' -e '.($time*10).' -f '.$pcRep.' -m TRUE -s '.$useStr.' '.$excl;
}

unless ($runMe) {
	print $cmd."\n";
}else{
	print STDERR $cmd."\n";
	system($cmd)
}

######################################################################################
sub getFile{
	my ($f,$d,$type) = @_;

	my $ret = '';
	$f =~ s/\.(bedgraph|bed)//;
	if ($type eq 'ori'){
 		for my $n(@ORI){
			if ($n =~ /$f.+(bed|bedgraph)$/){
				return $n;
			}
		}
	}

	print STDERR "[$type] ; $f ... NOT FOUND!!!\n";

	return;
}

######################################################################################
sub doLS{
	my ($pathBG,$pathOri) = @_;
	open my $PIPE, '-|' , 'ls '.$pathBG.'/*bedgraph '.$pathOri.'/*bedgraph';
	while (<$PIPE>){
		chomp;
		#print STDERR $_."\n";
		my $realfile = $_;
		my $fName    = $realfile;
		$fName =~ s/^.+\///;
		push @RT, $realfile;

		$RTnames{$fName} = $realfile;
	}
	close $PIPE;

}

# sub getBG{
# 		my ($name,$bgPath) = @_;
# 		opendir(DIR,$bgPath);
# 		while (my $fBG = readdir(DIR)){
# 			if ($fBG =~ /$name/){
# 				return $bgPath."/".$fBG;
# 			}
# 		}
#
# 		return 0;
# }
