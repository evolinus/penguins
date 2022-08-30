#!/usr/bin/perl

use diagnostics;
use warnings;
use strict;

#############################################################################
my $mydir = "codeml_output";

my $indir  = "24-pamlout/$mydir";
my $outdir  = "27-paml_results/$mydir";
mkdir ($outdir, 0777);

my @treebranches = ('AF', 'AP', 'AF-AP');
#

my $outfile = "$outdir/dNdSw.txt";
open (OUT, ">$outfile");


my %values;
print OUT "gene\tomegafix";
foreach my $treebranch (sort @treebranches) {
	print OUT "\ttot_$treebranch\tN_$treebranch\tS_$treebranch\tw_$treebranch\tdS_$treebranch\tdN_$treebranch";
}

print "\n";
opendir (DIR, $indir) or die "can't opendir $indir: $!";
while (defined(my $genedir = readdir(DIR))) {
	if ($genedir =~ /^AF/) {
		print "Harvesting dN, dS and dN/dS from ... $genedir     \r";
		my $file = "$indir/$genedir/$genedir.results2.txt";
		print OUT "\n$genedir";
		my $omega = 0;
		my $branch = 0;
		open (IN, $file);
		while(<IN>) {
			chomp($_);
			my @data = split /\t/, $_;
			my $siz = scalar @data;
			if ($siz > 0) {
				if ($omega == 1 && $branch == 0 && $data[0] ne 'branch') {
					print OUT "\t$data[0]";
				}
				if ($omega == 1 && $branch == 1) {
					#print OUT "\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\t$data[8]";
					$values{$data[0]} = "\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]";
					#print "\n$data[0]\t$values{$data[0]}";
				}
				if ($data[0] eq 'omegafix') {
					$omega = 1;
				}
				if ($data[0] eq 'branch') {
					$branch = 1;
				}
			}
		}
		
		foreach my $treebranch (sort @treebranches) {
			my $thesevalues = $values{$treebranch};
			print OUT "$thesevalues";
		}	
		
		close IN;
	}
}
print "\n\n";
closedir DIR;
close OUT;

exit;
