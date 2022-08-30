#! /perl/bin/perl

use diagnostics;
use warnings;
use strict;
#############################################################################

my $mydir  = "star20-AF-EH-GS-FG-GC-NN-EG-PC-HL-TA-EM-CC-SM-CB-PL-PA-PP-AP-EC-OH";

my @myspecies = split (/-/, $mydir);
shift(@myspecies);

my $indir  = "10-alignments/$mydir";
my $outdir = "12-alignments_parsed/$mydir";
mkdir ($outdir, 0777);
open (LOG, ">12-alignments_parsed/$mydir.ortho_names.txt");

my %orthonames;
foreach my $species (@myspecies) {
	print LOG "$species\_newname\t$species\_originalname\t";
	my $infile = "0-$species/$species.cds_renaming_log.txt";
	open(IN, $infile);
	while(<IN>) {
		chomp($_);
		my @genenames = split /\t/, $_;
		$orthonames{$genenames[1]} = $genenames[0];
	}
	close IN;	
}	

opendir (DIR, $indir) or die "can't opendir $indir: $!";
while (defined (my $file = readdir(DIR))) {
	#if ($file =~ /.fasta$/) {
	if ($file =~ /.aln$/ || $file =~ /.fasta$/) {
		print "Processing $file       \r";
		my $infile  = "$indir/$file";
		$file =~ s/aln/fasta/;
		my $outfile = "$outdir/$file";
		open (OUT, ">$outfile");
		open (IN, $infile) or die print "Cannot open file \"$infile\"\n\n";
		my $c=0;
		print LOG "\n";
		while (my $line=<IN>) {
			chomp($line);
			if ($line =~ /^>/) {
				my $name = $line;
				$name =~ s/>//;
				if (exists $orthonames{$name}) {
					print LOG "$name\t$orthonames{$name}\t";
				} else {
					print LOG "NA\tNA\t";
				}
				$name =~ /(.*)_[0-9]+/;
				$name = $1;
				$name =~ s/_R_//;
				$c++;
				if ($c == 1) {
					print OUT ">$name\n";
				}
				if ($c > 1) {
					print OUT "\n>$name\n";
				}
			} else {
				print OUT "$line";
			}
		}
	close OUT;
	}
}
close LOG;
print "\n";
exit;
