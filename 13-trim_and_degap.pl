#!/usr/local/bin/perl

use diagnostics;
use warnings;
use strict;
#############################################################################

my $mydir  = "star20-AF-EH-GS-FG-GC-NN-EG-PC-HL-TA-EM-CC-SM-CB-PL-PA-PP-AP-EC-OH";

my @myspecies = split (/-/, $mydir);
shift(@myspecies);
my $focal = $myspecies[0]; 

my $howmany = scalar (@myspecies);

my $tooldir = "1000-tools";
my $indir = "12-alignments_parsed/$mydir";
my $degappedaligndir = "14-translatorx_in/$mydir";
my $prankoutdir = "16-translatorx_out/$mydir";
mkdir ($degappedaligndir, 0777);
mkdir ($prankoutdir, 0777);
open (PRANK, ">15-prankcommands_$mydir.sh");

my $file;
my %myseqs;
opendir (DIR, $indir) or die "can't opendir $indir: $!";
while (defined($file = readdir(DIR))) {
	if ($file =~ /fasta/) {
		my $name = $file;
		print "Processing ... $name         \r";
		$name =~ s/mafft.//;
		open (OUT, ">$degappedaligndir/$name");

		my @data = slurp ("$indir/$file");
		my $nseqs = scalar (@data);
		for (my $i=0; $i<$nseqs; $i+=2) {
			my $species = $data[$i];
			my $seq = $data[$i+1];
			$species =~ s/>//;
			$myseqs{$species} = $seq;
		}
		
		my $focalseq = $myseqs{$focal};
		my $length = length $focalseq;
		my @splitseq = split (//, $focalseq);
		my $start = 0;
		for (my $i=0; $i<$length; $i++) {
			if ($splitseq[$i] eq "-") {
				$start++;
			} elsif ($splitseq[$i] ne "-") {
				last;
			}
		} 
		my $end = $length - 1;
		for (my $i=$length-1; $i>0; $i--) {
			if ($splitseq[$i] eq "-") {
				$end--;
			} elsif ($splitseq[$i] ne "-") {
				last;
			}
		} 
		my $subseqlength = $end - $start + 1 - 3;  # -3 is to remove eventual stop codons
		
		my $species;
		# trim sequences so that they and end where the focal sequence does
		my %mytrimmedseqs;
		foreach $species (sort(keys(%myseqs))) {
			$mytrimmedseqs{$species} = substr ($myseqs{$species}, $start, $subseqlength);
		}
		
		# check for positions where focal sequence has gaps 
		my %gaps; 
		$focalseq = $mytrimmedseqs{$focal};
		$length = length $focalseq;
		@splitseq = split (//, $focalseq);
		for (my $i=0; $i<$length; $i++) {
			if ($splitseq[$i] eq "-") {
				$gaps{$i} = 1;
			}
		} 
		
		# remove positions where the focal sequence has gaps
		my %degappedseqs; 
		foreach $species (sort(keys(%mytrimmedseqs))) {
			my $preseq = $mytrimmedseqs{$species};
			$length = length $preseq;
			@splitseq = split (//, $preseq);
			my $degapped = "";
			for (my $i=0; $i<$length; $i++) {
				if (!exists $gaps{$i}) {
					$degapped .= $splitseq[$i];				
				}
			}
			$degappedseqs{$species} = $degapped;
		}
		
		# remove codons with gaps (i.e. 1-3 '-')
		foreach $species (sort(keys(%degappedseqs))) {
			my $degapped = "";
			$length = length $degappedseqs{$species};
			for (my $i=0; $i<$length-3; $i+=3) {
				my $three = substr($degappedseqs{$species}, $i, 3);
				if ($three !~ /-/) {
					$degapped .= $three;				
				} elsif ($species eq $focal) {
					print "\n$name has focal sequence with gaps!\n"}
			}
			print OUT ">$species\n$degapped\n";
		}
		close OUT;
		print PRANK "perl $tooldir/translatorx_mod.pl -i $degappedaligndir/$name -o $prankoutdir/$name\_out -p P\n";
	}	
}
print "Finito.                                              \n";	
closedir (DIR);
close PRANK;
exit;

###############################################################################
#	Slurps in all data; Receives FILENAME as argument; outputs an array of data
#
sub slurp {
	my ($fileName) = @_;
	my @fileData = ();
	open(GET_DATA, $fileName);
	while(<GET_DATA>) {
		chomp($_);
		push (@fileData, $_);
	}
	close GET_DATA;
	return @fileData;
}	
	

	
