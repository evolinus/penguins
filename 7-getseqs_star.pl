#! /perl/bin/perl

use diagnostics;
use warnings;
use strict;
#############################################################################

print "\n\n***********************************************************************\n";

my @myspecies = ("AF", "EH", "GS", "FG", "GC", "NN", "EG", "PC", "HL", "TA", "EM", "CC", "SM", "CB", "PL", "PA", "PP", "AP", "EC", "OH"); # focal species must be first species

my $fractionidentity_threshold 	= 0.7; 		######     0-1; change if necessary
my $fractionaligned_threshold 	= 0.6; 		######     0-1; change if necessary
my $percentidentity_threshold 	= int($fractionidentity_threshold * 100);
my $percentalign_threshold 		= int($fractionaligned_threshold * 100);

#my $aligner = "muscle";
my $aligner = "mafft";

my $howmanymissing = 0;
my $nspecies = scalar (@myspecies);
my $focal = $myspecies[0];

my $blastoutdir = "6-reducedblastout";
my $outdir = "8-seqs";
my $alignoutdir = "10-alignments";

#######################

my %goodspecies;
my $beaststring = "star$nspecies";
foreach my $beast (@myspecies) {
	$beaststring .= "-".$beast;
}


my %mysequences;
my $newnameHit;
takeseqs (@myspecies);
print "\n";

my $focaldir = "$beaststring";
$focaldir = $outdir."/".$focaldir;
mkdir ($focaldir, 0777);
my $repdir = "reports"; 
$repdir = $focaldir."/".$repdir;
mkdir ($repdir, 0777);
my $aligneroutdir = "$beaststring"; 
$aligneroutdir = $alignoutdir."/".$aligneroutdir;
mkdir ($aligneroutdir, 0777);
my $fastadir = "fasta"; 
$fastadir = $focaldir."/".$fastadir;
mkdir ($fastadir, 0777);
my $matchesfile = "_matches.txt";
$matchesfile = $focaldir."/".$matchesfile;
open (MATCHES, ">$matchesfile");	
print MATCHES "$focal\t";
my $newmatchesfile = "_matches_taken.txt";
$newmatchesfile = $focaldir."/".$newmatchesfile;
open (MATCHESTAKEN, ">$newmatchesfile");	
print MATCHESTAKEN "$focal\t";

###################################
###

my @allmatches = processblast (@myspecies);

###
###################################

my $a = $allmatches[0];
my %focalspecieshash = %$a;
print "\n";
my $count = 0;
foreach my $key (sort (keys(%focalspecieshash))) {	
	$count++;
	print MATCHES "\n$key\t";
	print "Writing blast data ... $count                \r";
	for (my $s=0; $s<$nspecies-1; $s++)	{
		if (exists $focalspecieshash{$key}[$s])	{ 
			print MATCHES "$focalspecieshash{$key}[$s]\t";
		} else {  
			print MATCHES "NA\t";
		}
	}
	for (my $second = 1; $second<$nspecies; $second++) {
		my $b = $allmatches[$second];
		my %secondhash = %$b;
		if (exists $focalspecieshash{$key}[$second-1]) {
			if ( exists $secondhash{$focalspecieshash{$key}[$second-1]}[0] ) { 
				print MATCHES "$secondhash{$focalspecieshash{$key}[$second-1]}[0]\t";
			} else { 
				print MATCHES "NA\t";
			}
		} else { 
			print MATCHES "NA\t";
		}
	}
}
close MATCHES;
print "Writing blast data ... done          \n";

my $alignaligner_infile = "9-align_$beaststring.$aligner.sh";
open (ALI, ">$alignaligner_infile");	
my $numberpairs = ($nspecies - 1);
	
print "\n";	
my @blastmatches = slurp ($matchesfile);
shift @blastmatches;
foreach (@blastmatches) {
	my @data = split /\t/, $_;
	my $gene = $data[0];
	print "Processing reciprocal blast results ... $data[0]          \r";
	my $comparisons = scalar (@data) - 1;
	my $na = 0;
	my $blastmismatch = 0;
	my $blastNA = 0;
	for (my $which=0; $which<$nspecies; $which++) {
		$goodspecies{$myspecies[$which]} = 0;
	}			
	foreach my $comp (@data) {
		if ($comp eq "NA") { $na++ }
	}
	if ($na >= 0 ) {
		my @wrong = ();
		my $steps = 0;
		my $blasthit;
		for (my $pair=1; $pair<$numberpairs+1; $pair++) {
			if ( $data[$nspecies+$pair-1] ne $data[0] && $data[$nspecies+$pair-1] ne 'NA') {
				$blastmismatch++;
				$goodspecies{$myspecies[$pair]} = 1;
			}
			if ( $data[$nspecies+$pair-1] ne $data[0] && $data[$nspecies+$pair-1] eq 'NA') {
				$blastNA++;
				$goodspecies{$myspecies[$pair]} = 1;
			}
		}
		if ( ($blastmismatch + $blastNA) <= $howmanymissing) {
		#if ($blastmismatch == 0 && $blastNA <= $howmanymissing) {
		#if ($blastmismatch == 0) {
			my $filename = $fastadir."/".$data[0].".orthologs.fas";
			open (OUT, ">$filename");
			for (my $which=0; $which<$nspecies; $which++) {
				my $this = $data[$which];
				if ($this ne 'NA' && $goodspecies{$myspecies[$which]} == 0) {
				#if ($this ne 'NA') {
					my $newthis = "";
					if($which == 0)	{
						$newthis = $this;
					} else {
						$newthis = $data[0]."_".$this;
					}
					if (exists $mysequences{$newthis}) {
						print OUT ">".$this."\n".$mysequences{$newthis}."\n";
					} else {
						print OUT ">".$this."\n".$mysequences{$this}."\n";
					}
				}	
			}
			close OUT;
			
			my $newline = join ("\t", @data);
			print MATCHESTAKEN "\n$newline";
			
			if ($aligner eq 'mafft') {
				my $outfile = "$aligneroutdir/$data[0].mafft.fasta";
				print ALI "mafft --localpair --maxiterate 1000 --adjustdirection $filename > $outfile\n";
			} elsif ($aligner eq 'muscle') {
				my $outfile = "$aligneroutdir/$data[0].muscle.fasta";
				print ALI "muscle -in $filename -out $outfile\n";
			}
		}	
	}
}
print "Processing reciprocal blast results ... done                   \n";
close ALI;

print "\n***********************************************************************\n\n\n";
exit;



	
###############################################################################
#	take blast output and retrieve putative orthologues pairs
#	input is name of species1 and species2
#	output is a hash with key=species1 and argument=species2 names

sub processblast {	
	my @myspecies = @_;
	my @mymatches;
	my $howmany = scalar (@myspecies);
	my %thegenes = ();
	my $slot = 0;
	
	for (my $i=1; $i<$howmany; $i++) {
		my $species1 = $myspecies[0];
		my $species2 = $myspecies[$i];
		print MATCHES "$species1-$species2\_match\t";
		print MATCHESTAKEN "$species1-$species2\_match\t";
		my $one2two_blastname = $species1."2".$species2."_my.out";
		$one2two_blastname = $blastoutdir."/".$one2two_blastname;
		my $one2two_blastedfilereport = "$species1-$species2.report.txt";
		$one2two_blastedfilereport = $repdir."/".$one2two_blastedfilereport; 
		open (REP, ">$one2two_blastedfilereport");
		print REP "$species1(S1)\t$species2(S2)\tLengthS1\tLengthS2\tLengthAlignment\tFractS1\tFractS2\tOkAlign\tFractIdentity\tOkIdentity\tTaken\n";
		my @one2two_blast = slurp($one2two_blastname);
		shift (@one2two_blast);
		my $nameQuery_bis = "nul";
		my $nameHit_bis = "nul";
		foreach (@one2two_blast) {
			my @blastdata = split /\t/, $_;
			my $nameQuery = $blastdata[0];
			my $nameHit   = $blastdata[1];
			my $start = $blastdata[8];
			my $end   = $blastdata[9];
			print "Processing '$species1' to '$species2' blast ... $nameQuery -> $nameHit\r";
			my $length1 = length($mysequences{$nameQuery});
			my $length2 = length($mysequences{$nameHit});
			my $seq = $mysequences{$nameHit};
			my @sequence = split(//, $seq);	
			if ($end < $start) {
				@sequence = complement($seq); 
			}
			$newnameHit = $nameQuery."_".$nameHit;
			$mysequences{$newnameHit} = join ("", @sequence);
			my $fractionidentity = $blastdata[2];
			my $fractionaligned1 = sprintf "%.2f", $blastdata[3] / $length1;
			my $fractionaligned2 = sprintf "%.2f", $blastdata[3] / $length2;
			my $okfractionidentity = "no";
			my $okfractionaligned = "no";
			my $taken = "no";
			if ($fractionidentity >= $fractionidentity_threshold) {	
				$okfractionidentity = "yes";	
			}
			if ($fractionaligned1 >= $fractionaligned_threshold || $fractionaligned2 >= $fractionaligned_threshold) {
				$okfractionaligned = "yes";
			}
			if ($okfractionaligned eq "yes" && $okfractionidentity eq "yes") {
				$taken = "yes";
			}
			print REP "$nameQuery\t$nameHit\t$length1\t$length2\t$blastdata[3]\t$fractionaligned1\t$fractionaligned2\t$okfractionaligned\t$fractionidentity\t$okfractionidentity\t$taken\n";
			if ($taken eq "yes" && $nameQuery ne $nameQuery_bis) {
				$thegenes{$nameQuery}[$slot] = $nameHit;
			} elsif ($taken eq "no" && $nameQuery ne $nameQuery_bis) {
				$thegenes{$nameQuery}[$slot] = "NA";
			}
			$nameQuery_bis = $nameQuery;
			$nameHit_bis = $nameHit;	
		}
		print "Processing '$species1' to '$species2' blast ... done                                 \n";
		$slot++;
		close REP;	
	}
	$mymatches[0] = \%thegenes;	

	for (my $i=1; $i<$howmany; $i++) {
		my %thegenes2 = ();
		$slot = 0;
		my $species1 = $myspecies[$i];
		my $species2 = $myspecies[0];
		print MATCHES "$species1-$species2\_match\t";
		print MATCHESTAKEN "$species1-$species2\_match\t";
		my $one2two_blastname = $species1."2".$species2."_my.out";
		$one2two_blastname = $blastoutdir."/".$one2two_blastname;
		my $one2two_blastedfilereport = "$species1-$species2.report.txt";
		$one2two_blastedfilereport = $repdir."/".$one2two_blastedfilereport; 
		open (REP, ">$one2two_blastedfilereport");
		print REP "$species1(S1)\t$species2(S2)\tLengthS1\tLengthS2\tLengthAlignment\tFractS1\tFractS2\tOkAlign\tFractIdentity\tOkIdentity\tTaken\n";
		my @one2two_blast = slurp($one2two_blastname);
		shift (@one2two_blast);
		my $nameQuery_bis = "nul";
		my $nameHit_bis = "nul";
		foreach (@one2two_blast) {
			my @blastdata = split /\t/, $_;
			my $nameQuery = $blastdata[0];
			my $nameHit = $blastdata[1];
			my $start = $blastdata[8];
			my $end = $blastdata[9];
			print "Processing '$species1' to '$species2' blast results ... $nameQuery - $nameHit\r";
			my $length1 = length($mysequences{$nameQuery});
			my $length2 = length($mysequences{$nameHit});
			my $fractionidentity = $blastdata[2];
			my $fractionaligned1 = sprintf "%.2f", $blastdata[3] / $length1;
			my $fractionaligned2 = sprintf "%.2f", $blastdata[3] / $length2;
			my $okfractionidentity = "no";
			my $okfractionaligned = "no";
			my $taken = "no";
			if ($fractionidentity >= $fractionidentity_threshold) {	
				$okfractionidentity = "yes";	
			}
			if ($fractionaligned1 >= $fractionaligned_threshold || $fractionaligned2 >= $fractionaligned_threshold) {
				$okfractionaligned = "yes";
			}
			if ($okfractionaligned eq "yes" && $okfractionidentity eq "yes") {
				$taken = "yes";
			}
			print REP "$nameQuery\t$nameHit\t$length1\t$length2\t$blastdata[3]\t$fractionaligned1\t$fractionaligned2\t$okfractionaligned\t$fractionidentity\t$okfractionidentity\t$taken\n";
			if ($taken eq "yes" && $nameQuery ne $nameQuery_bis) {
				$thegenes2{$nameQuery}[$slot] = $nameHit;
			} elsif ($taken eq "no" && $nameQuery ne $nameQuery_bis) {
				$thegenes2{$nameQuery}[$slot] = "NA";
			}
			$nameQuery_bis = $nameQuery;
			$nameHit_bis = $nameHit;	
		}
		print "Processing '$species1' to '$species2' blast results ... done                                 \n";
		$slot++;
		close REP;	
		$mymatches[$i] = \%thegenes2;
	}			
	return @mymatches;
}


###############################################################################
#	import sequences from eaqch species' fasta file 
#	input is name of species
#	output is a hash with all sequences, where key=species and argument=sequence

sub takeseqs {
	my @thespecies = @_;
	foreach my $spec (@thespecies) {
		my $filename = $spec.".txt";
		my $filedir = "0-".$spec;
		$filename = $filedir."/".$filename; ################
		my $counter = 0;
		my @data = slurp($filename);
		foreach (@data) {
			$counter++;
			print "Importing '$spec' sequences ... $counter\r";
			my @theseq = split /\t/, $_;
			$mysequences{$theseq[0]} = $theseq[1];
		}
		print "Importing '$spec' sequences ... done   \n";
	}
}

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
		push (@fileData, $_);
	}
	close GET_DATA;
	return @fileData;
}	
	

###############################################################################
#	user interaction commands
#
sub promptUser	{
	my($promptString,$default) = @_;
	my $input;
	if ($default) {
		print $promptString, "[", $default, "] : ";   	
	} else {
		print $promptString, " ";
	}
	$| = 1;
	chomp($input = <STDIN>);
	if ($default && !$input) {
	  	return $default;
	} else {
		return $input;	
	}
}

###############################################################################
#	make sequence complement
#
sub complement {
	my ($inseq) = @_;
	my $seqlength = length ($inseq);
	my @sequence = split //, $inseq;
	my @complsequence = ();
	my $j = 1;
	for (my $i = 0; $i < $seqlength; $i++) {
		my $done = 0;
		if ($sequence[$seqlength-$j] eq 'A' || $sequence[$seqlength-$j] eq 'a') {	$complsequence[$i] = "T"; $done = 1;	}
		if ($sequence[$seqlength-$j] eq 'C' || $sequence[$seqlength-$j] eq 'c') {	$complsequence[$i] = "G"; $done = 1;	}
		if ($sequence[$seqlength-$j] eq 'G' || $sequence[$seqlength-$j] eq 'g') {	$complsequence[$i] = "C"; $done = 1;	}
		if ($sequence[$seqlength-$j] eq 'T' || $sequence[$seqlength-$j] eq 't') {	$complsequence[$i] = "A"; $done = 1;	}
		if ($sequence[$seqlength-$j] eq 'N' || $sequence[$seqlength-$j] eq 'n') {	$complsequence[$i] = "N"; $done = 1;	}
		if ($sequence[$seqlength-$j] eq '-') {	$complsequence[$i] = "-"; $done = 1;	}
		elsif ($done == 0) {	$complsequence[$i] = $sequence[$seqlength-$j]	}
		$j++;
	}
	return @complsequence;
}
	
	
