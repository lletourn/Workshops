#!/usr/bin/perl -w

use strict;
use warnings;
use POSIX;

# Jessica Wasserscheid, Genome Quebec Innovation Center, 2006
# Takes a fasta file containing sequences on multi-lines
# and remove any new line to end up with each sequence on one line

my ($multifasta) = @ARGV;
open F, "<$multifasta";

my @array = <F>;
my $i = 0;
close F;

for (my $j=0; $j < scalar @array; $j++)
{
	if (defined $array[$j])
	{
		if($array[$j] =~ /^>/)
		{
      if($j > 0) {
        print "\n";
      }
			print STDOUT "$array[$j]";
		}
		else
		{
			chomp $array[$j];
			print STDOUT $array[$j];
		}
		#$string .= "\n";
	}
}

close F;
