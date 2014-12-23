#!/usr/bin/perl -w

use strict;
use warnings;
use POSIX;

my ($ctfile)=@ARGV;
open CTFILE, $ctfile;
my $prevpos=0;

my $zeroCovCount=1;
my $chr;
my $startPos;

$_ = <CTFILE>;
my @array=split /\s+/,$_;
$chr=$array[0];
$startPos=$array[1];
$prevpos = $startPos;

while(<CTFILE>)
{
	my @array=split /\s+/,$_;
	$chr=$array[0];
	my $pos=$array[1];
	if($pos > ($prevpos+1)){
		print STDOUT "$chr\t".($startPos-1)."\t".$prevpos."\t".$zeroCovCount."\n";
    $zeroCovCount++;
    $startPos = $pos;
	}
	$prevpos=$pos;
}
print STDOUT "$chr\t".($startPos-1)."\t".$prevpos."\t".$zeroCovCount."\n";
close CTFILE;
