#!usr/bin/env perl

# Last modified Oct 28, 16:40 CST

use strict;
use warnings;
use diagnostics;
use Getopt::Long;

my $sequence;
my $infile_fh;
my $outfile_fh;
my $help;
my $three;

GetOptions (	"infile=s" => \$infile_fh,
				"outfile=s" => \$outfile_fh,
				"sequence=s" => \$sequence,
				"three" => \$three,
				"help" => \$help);

if ($help){
print <<"HELP";

--------------------------------------------------------------------------------
                              fastaptamer_distance
--------------------------------------------------------------------------------

Usage: fastaptamer_distance [-h] [-i INFILE] [-o OUTFILE] [-s] [-t]

  [-h]         = Help screen.
  [-i INFILE]  = Input file from fastaptamer_enrich; required.
  [-o OUTFILE] = Output file, plain-text with tab-separated values. 
  [-s]         = Sequence to compare edit distance to.
  [-t]         = fastaptamer_enrich file with three populations

fastaptamer_distance will calculate the edit distance from a defined sequence for 
each sequence in a fastaptamer_enrich file. Default behavior is designed for 
enrich file with 2 populations, use -t for three population enrich file.

HELP
exit;
}

open (INPUT, '<', $infile_fh) or die "No input file specified.\nSee help [-h] for correct usage.\n";
open (OUTPUT, '>', $outfile_fh) or die "No output file specified.\nSee help [-h] for correct usage.\n";
unless ($sequence) { die "No sequence specified.\nSee help [-h] for correct usage.\n" }


while (<INPUT>){

		if ($_ =~ /^Sequence/){
			if (defined ($three)){
				print OUTPUT "Edit Distance\tSequence\tLength\tRank (x)\tReads (x)\tRPM (x)\tRank (y)\tReads (y)\tRPM (y)\tRank (z)\tReads (z)\tRPM (z)\tEnrichment (y/x)\tEnrichment (z/y)\tEnrichment (z/x)\n";
			}
			elsif (undefined ($three)){
				print OUTPUT "Edit Distance\tSequence\tLength\tRank (x)\tReads (x)\tRPM (x)\tRank (y)\tReads (y)\tRPM (y)\tEnrichment (y/x)\n";
			}
		}
		elsif ($_ =~ /(\S+)(.*)/){
			my $distance = levenshtein ( $sequence, $1 );
			print OUTPUT "$distance\t$1$2\n";
		}

}
		
sub levenshtein {
	my ($s1, $s2) = @_;
	my ($len1, $len2) = (length $s1, length $s2);
	
	return $len2 if ($len1 == 0);
	return $len1 if ($len2 == 0);
	
	my %mat;
	
	for (my $i = 0; $i <= $len1; ++$i){
		for (my $j = 0; $j <= $len2; ++$j){
			$mat{$i}{$j} = 0;
			$mat{0}{$j} = $j;
		}
		$mat{$i}{0} = $i;
	}
	
	my @ar1 = split(//, $s1);
	my @ar2 = split(//, $s2);

	for (my $i = 1; $i <= $len1; ++$i){
		for (my $j = 1; $j <= $len2; ++$j){
			my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;
			$mat{$i}{$j} = min([$mat{$i-1}{$j} + 1,
			$mat{$i}{$j-1} + 1,
			$mat{$i-1}{$j-1} + $cost]);
		}
	}
    return $mat{$len1}{$len2};
}

sub min
{
    my @list = @{$_[0]};
    my $min = $list[0];

    foreach my $i (@list)
    {
        $min = $i if ($i < $min);
    }

    return $min;
}