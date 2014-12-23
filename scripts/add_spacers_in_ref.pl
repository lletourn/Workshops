#!/usr/bin/perl -w

use strict;
use warnings;
use POSIX;


my ($fasta, $gapsz)=@ARGV;
open FA, $fasta;
print STDOUT ">chr1\n";

my $first=1;
while(<FA>)
{
	if($_=~/^>/)
	{
    if(!$first) {
  		for(my $c=1;$c<=$gapsz;$c++)
	  	{
		  	print STDOUT "N";
  		}
		  print STDOUT "\n";
    }
    $first = 0;
	}
	else
	{	print STDOUT $_;}
}
# Add last spacer
for(my $c=1;$c<=$gapsz;$c++)
{
  print STDOUT "N";
}
print STDOUT "\n";

close FA;
