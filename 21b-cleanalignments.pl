#!/usr/bin/perl 

use diagnostics;
use strict;
use warnings;
######################################################################

my $mydir = "prankete";
my $mycleandir = $mydir."_clean";

my $maxgapped = 0.20;  # 0-1; maximum fraction of sequences that can have a gap at a certain position
my @musthave = ('AP', 'AF');  # name of sequences that must always have no gap
my $minlength = 150;  # minimum length of alignment after degapping (should be multiple of 3)

my $dir = "21-alignments_clean_mastrolino";

my $indir = "$dir/$mydir";
my $outdir = "$dir/$mycleandir";
mkdir ($outdir, 0777);

open (LOG, ">$dir/$mydir.cleanerlog.txt");
print LOG "maximum fraction of sequences that can have a gap at a certain position = ".$maxgapped;
print LOG "\nname of sequences that must always have no gap:";
print LOG "\nminimum length of alignment after degapping (bp) = ".$minlength;
foreach my $must (@musthave) {
	print LOG " $must";
}
print LOG "\n\ngene\talignmentlength\tnewalignmentlength\ttaken";

opendir (DIR, $indir) or die "can't opendir $indir: $!";
while (defined (my $file = readdir(DIR))) {
	if ($file =~ /.fasta$/) {
		my $infile = $indir."/".$file;	 
		my $outfile = $outdir."/".$file;
		open (IN, $infile) or die print "Cannot open file \"$infile\"\n\n";
		$/ = "\n>";
		my $gene = $file;
		$gene =~ s/\.fasta//;
		#print "\n$gene\n";
	
		my %sequences;
		my %newsequences;
		my $originallength = 0;
		while (<IN>) {
			(my $species, my $seq) = split /\n/, $_;
			chomp $species;
			$species =~ s/>//;
			my @seqq = split "", $seq;
			$sequences{$species} = [@seqq];
			$newsequences{$species} = "";
			$originallength = length $seq;
		}		
		close IN;
		
		my $nseqs = keys %sequences;
		my $newlength = 0;
		for (my $i=0; $i<$originallength; $i++) {
			my $countgaps = 0;
			my $musts = 0;
			foreach my $key (sort(keys(%sequences))) {
				#print " $key\t$i\t${$sequences{$key}}[$i]\r";
				if (${$sequences{$key}}[$i] eq '-') {
					$countgaps++;
				}
			}	
			foreach my $must (@musthave) {
				if (${$sequences{$must}}[$i] eq '-') {
					$musts++;
				}
			}
			my $fractiongapped = $countgaps/$nseqs;
			if ($fractiongapped <= $maxgapped && $musts == 0) {
				$newlength++;
				foreach my $key (sort(keys(%newsequences))) {
					$newsequences{$key} .= ${$sequences{$key}}[$i];
				}	
			}
		}
		
		print LOG "\n$gene\t$originallength\t$newlength";
		if ($newlength < $minlength) {
			print LOG "\tno";
		} else {
			print LOG "\tyes";
			open (OUT, ">$outfile");
			foreach my $key (sort(keys(%newsequences))) {
				print OUT ">$key\n$newsequences{$key}\n"
			}
			close OUT;
		}
		
	}
}	
close DIR;
close LOG;

exit;
