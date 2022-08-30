#! /perl/bin/perl

use diagnostics;
use warnings;
use strict;
#############################################################################

my $mydir  = "star20-AF-EH-GS-FG-GC-NN-EG-PC-HL-TA-EM-CC-SM-CB-PL-PA-PP-AP-EC-OH";

#my $stops = "no";
my $wantstops = "no";

my @myspecies = split (/-/, $mydir);
shift(@myspecies);

my $indir = "16-translatorx_out/$mydir";
my $outdir = "19-alignments_nostops/$mydir";
mkdir ($outdir, 0777);

my $stoppedgenes = 0;
my $allgenes = 0;
opendir (DIR, $indir) or die "can't opendir $indir: $!";
while (defined(my $file = readdir(DIR))) {
	if ($file =~ /fasta/) {
		my $infile = "$indir/$file";
		my $outfile = "$outdir/$file";
		$outfile =~ s/fasta.//;
		my @allgeneseqs = slurp ($infile);
		my $nseqs = scalar @allgeneseqs;
		my $stops = 0;
		for (my $i=0; $i<$nseqs-1; $i+=2) {
			my $geneseqname = $allgeneseqs[$i];
			my $geneseq = $allgeneseqs[$i+1];
			my $seqlength = length $geneseq;
			for (my $i=0; $i<$seqlength-2; $i+=3) {
				my $codon = substr($geneseq,$i,3);
				$codon = uc $codon;
				if ($codon eq 'TAG' || $codon eq 'TAA' || $codon eq 'TGA') {
					$stops++;
				}
			}
		}
		if ($stops > 0) {
			$stoppedgenes++;
			print "\n$file contains internal stop codon(s)";
		}
		if ($stops == 0 || ($stops > 0 && $wantstops eq 'no')) {
			open (OUT, ">$outfile");
			for (my $i=0; $i<$nseqs-1; $i+=2) {
				my $geneseqname = $allgeneseqs[$i];
				my $geneseq = $allgeneseqs[$i+1];
				my $seqlength = length $geneseq;
				for (my $i=0; $i<$seqlength-2; $i+=3) {
					my $codon = substr($geneseq,$i,3);
					$codon = uc $codon;
					if ($codon eq 'TAG' || $codon eq 'TAA' || $codon eq 'TGA') {
						substr($geneseq,$i,3) = "NNN";
					}
				}
				print OUT "$geneseqname\n$geneseq\n";
			}
		}
		$allgenes++;
		close OUT;
	}	
}
print "\n\nAlignments with internal stop codons: $stoppedgenes of $allgenes\n\nFinito, buon appetito                                       \n\n";
exit;

###############################################################################
#	Slurps in all data
#	Receives FILENAME as argument; outputs an array of data
#
sub slurp {
	my ($fileName) = @_;
	my @fileData = ();
	open(GET_DATA, $fileName);
	while(<GET_DATA>) {
		chomp($_);
		push (@fileData, $_) ;
	}
	close GET_DATA;
	return @fileData;
}	




