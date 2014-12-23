#!/usr/bin/perl -w
# usage : perl generate_tracks_for_pacbio_reads.pl  blast_output.out  genome_size  >  output.bed

use strict;
use warnings;
use POSIX;

my ($megab, $genomesz)=@ARGV;
open MEGAB,"$megab";

while (<MEGAB>)
{
	my @array=split /\s+/,$_;
	my $readlength=0;
	if($array[0]=~/([0-9]+)_([0-9]+)$/)
	{
		$readlength=$2-$1;
	}
	my $strand='+';
	my $readstart=$array[6];
	my $readend=$array[7];
	my $refstart=$array[8];
	my $refend=$array[9];
	
	my $realrefst=$refstart-$readstart;
	my $realrefe=$refend+($readlength-$readend);

	my $blockst2=$readstart+1;
	my $blockst3=$readlength-1;
	
	if($array[8]>$array[9])
	{
		$strand='-';
		$readstart=$readlength-$array[7];
		$readend=$readlength-$array[6];
		$refstart=$array[9];
		$refend=$array[8];
		$realrefst=$refstart-$readstart;
		$realrefe=$refend+$array[6];
		$blockst2=$readstart+1;
		$blockst3=$readlength-1;
	}

	if($realrefst<0)
	{$realrefst=0;$blockst2=$refstart;}
	if($realrefe>$genomesz)
	{$realrefe=$genomesz;$blockst3=$realrefe-1;}
	
	my $blocksz=$refend-$refstart-1;
	#print STDOUT "$array[1]\t$realrefst\t$realrefe\t$array[0]\t1000\t$strand\t$refstart\t$refend\t0,0,0\n";
	print STDOUT "$array[1]\t$realrefst\t$realrefe\t$array[0]\t1000\t$strand\t$realrefst\t$realrefe\t0,0,0\t3\t1,$blocksz,1\t0,$blockst2,".($realrefe-$realrefst-1)."\n";
}
close MEGAB;
