#! /perl/bin/perl

use diagnostics;
use warnings;
use strict;

#############################################################################

my (%speciesshortname) = (
						'PL'	=> 'PL',
						'EH'    => 'EH',
						'GS'	=> 'GS',
						'FG'	=> 'FG',
						'GC'	=> 'GC',
						'NN'    => 'NN',
						'EG'    => 'EG',
						'PC'    => 'PC',
						'HL'    => 'HL',
						'TA'    => 'TA',
						'EM'    => 'EM',
						'CC'    => 'CC',
						'SM'    => 'SM',
						'CB'    => 'CB',
						'AF'    => 'AF',
						'PA'    => 'PA',
						'PP'    => 'PP',
						'AP'    => 'AP',
						'EC'    => 'EC',
						'OH'    => 'OH',
						);
					
my @mycds = ("PL", "EH", "GS", "FG", "GC", "NN", "EG", "PC", "HL", "TA", "EM", "CC", "SM", "CB", "AF", "PA", "PP", "AP", "EC", "OH");

my $minCDSlength = 150;

foreach my $species (@mycds) {
	print "\n";
	my $shortname = $speciesshortname{$species};
	my $speciesdir = "0-$shortname";
	mkdir ($speciesdir, 0777);
	my $filename = "00-all_cds_fasta/$species.fna";
	open (INFILE, $filename) or die "Can not open input $filename!\n";
	open (RM, ">$speciesdir/_$species\_has_been_renamed_$shortname\_");
	close RM;
	open (LOG, ">$speciesdir/$shortname.cds_renaming_log.txt");
	print LOG "original\tnew";
	my %sequence;
	my $c=0;
	my $newseqname = "";
	while (my $line=<INFILE>) {
		chomp($line);
		if ($line =~ /^>/) {
			$c++;
			my $seqname = $line;
			$seqname =~ s/>//;
			$newseqname = "$shortname\_$c";
			print LOG "\n$seqname\t$newseqname";
			print "parsing CDS of $species (->$shortname): $c   \r";
			$sequence{$newseqname} = "";
		} else {
			$sequence{$newseqname} .= "$line";
		}
	}
	close LOG;
	open (TXT, ">$speciesdir/$shortname.txt");
	open (FAS, ">$speciesdir/$shortname.fna");
	foreach my $who (sort (keys(%sequence))) {
		my $seqlength = length $sequence{$who};
		if ($seqlength >= $minCDSlength) {
			print FAS ">$who\n$sequence{$who}\n";
			print TXT "$who\t$sequence{$who}\n";
		}	
	}
	close TXT;
	close FAS;
}	
print "\n\n";
exit;
