#!/usr/bin/perl

use diagnostics;
use warnings;
use strict;
use Statistics::Distributions;

#############################################################################
my $mydir = "codeml_output";
my @species = ('AF', 'AP', 'AF-AP');

my $pamlindir  = "24-pamlout/$mydir";
my $outdir  = "27-paml_results/$mydir";
mkdir ($outdir, 0777);
#

my %tests;
# 	test name       	       (alt model,   null model)
#$tests{"SITE-TEST 12"}		= [("Sa_M2a", 	"Sn_M1a"		)];
#$tests{"SITE-TEST 78"}		= [("Sa_M8", 	"Sn_M7"			)];
$tests{"BRANCH-TEST"}		= [("Ba", 		"n_M0"			)];
#$tests{"BRANCH-TEST 2"}		= [("Ba", 		"Bn"			)];
#$tests{"CLADE-TEST"}		= [("CmC", 		"Sn_M1a"		)];
#$tests{"CLADE-TEST 2"}		= [("CmC", 		"CmC_M2a_rel"	)];
#$tests{"BRANCHSITE-TEST"}	= [("BSa", 		"BSn"			)];
###

my %res;
foreach my $test (keys %tests) {
	$res{$test} = "gene";
	foreach my $spec (sort @species) {
		$res{$test} .= "\t${$tests{$test}}[0]\_$spec\t${$tests{$test}}[1]\_$spec\tdf_$spec\tP_$spec";
	}
}

print "\n";
opendir (DIR, $pamlindir) or die "can't opendir $pamlindir: $!";
while (defined(my $genedir = readdir(DIR))) {
	if ($genedir =~ /^AF/) {
		print "Harvesting PAML results from ... $genedir     \r";
		my $thedir = "$pamlindir/$genedir";
		foreach my $test (keys %tests) {
			$res{$test} .= "\n$genedir";
		}
		foreach my $spec (sort @species) {
			my $file = "$thedir/$genedir.$spec.results.txt";
			open(IN, $file);
			while(<IN>) {
				chomp($_);
				my @data = split /\t/, $_;
				my $altLik = $data[3];
				my $nulLik = $data[4];
				my $df = $data[5];
				my $chisqprob;
				if ($altLik ne 'NA' && $nulLik ne 'NA' && $data[0] ne 'test') {
					my $delta = 2 * ($altLik - $nulLik);
					if ($delta < 0) {
						$delta = 0;
					}
					$chisqprob = Statistics::Distributions::chisqrprob ($df, $delta);
				} else {
					$chisqprob = "NA";
				}	
				foreach my $test (keys %tests) {
					if ($data[0] eq $test) {
						$res{$test} .= "\t$altLik\t$nulLik\t$df\t$chisqprob";
					}
				}
			}
			close IN;
		}
	}
}
closedir DIR;

foreach my $test (keys %tests) {
	my $test2 = $test;
	$test2 =~ s/ //g;
	open (OUT, ">$outdir/$test2.txt");
	print OUT "$res{$test}";
	close OUT;
}

exit;
