#!/usr/bin/perl -w

use strict;
use locale;
use warnings;
use POSIX;

my ($fasta)=@ARGV;
my %hash=();

open FASTA, $fasta;
while(<FASTA>)
{
	if ($_ =~ /(^>(\S+))/)	
	{
		chomp $_;
		#my $header=$1;
		my $readname = $1;
		#my $pos = $2;
		#my $length = $3;
		
		my $seq = <FASTA>;
		chomp $seq;
		my $length=length $seq;
		
		$hash{$length}{$readname}=$seq;
	}
}

close FASTA;

foreach my $length (sort {$b <=> $a} keys %hash)
{
	foreach my $readname (keys %{$hash{$length}})
	{	
		#print "\t$readname\n";
		print STDOUT "$readname\t$length\n$hash{$length}{$readname}\n";
		#print STDOUT "$readname\t$length\n$hash{$length}{$readname}\n$seq\n";
		#print STDOUT "$readname\t$length\n";
	}
}
