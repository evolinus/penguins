#! /perl/bin/perl

use diagnostics;
use warnings;
use strict;
#############################################################################

my $dir = "4-blastout";
my $outdir = "6-reducedblastout";

opendir (DIR, $dir) or die "can't opendir $dir: $!";
while (defined(my $file = readdir(DIR))) {
	if ($file =~ /.out/ && $file !~ /my/) {
		my $outfile = $file;
		$outfile =~ s/\.out/_my.out/;
		open (OUT, ">$outdir/$outfile");
		print OUT "query(gene)\thit(contig)\tav%identity\ttotalignmentlength\ttotmismatch\ttotgapopen\tquery_protostart\tquery_protoend\thit_protostart\thit_protoend\tmyevalue\tavbitscore";
		my $infile = $dir."/".$file;
		my @data = slurp($infile);
		my %count;
		my %data;
		my $takingthis = "";
		my $this = "";
		my ($identity, $alignmentlength, $mismatch, $gapopen, $queryfrom, $queryto, $querylength, $hitfrom, $hitto, $hitlength, $evalue, $bitscore) = 0;
		my ($totidentity, $totalignmentlength, $totmismatch, $totgapopen, $totquerylength, $tothitlength, $myevalue, $totbitscore) = 0;
		foreach my $line (@data) {
			$line =~ s/\*/NA/g;
			my @blastdata = split /\t/, $line;
			my $queryname = $blastdata[0];
			my $hitname = $blastdata[1];
			$this = "$queryname$hitname";
			$identity = $blastdata[2];
			$alignmentlength = $blastdata[3];
			$mismatch = $blastdata[4];
			$gapopen = $blastdata[5];
			$queryfrom = $blastdata[6];
			$queryto = $blastdata[7];
			$querylength = abs($queryto-$queryfrom) +1;
			$hitfrom = $blastdata[8];
			$hitto = $blastdata[9];
			$hitlength = abs($hitto-$hitfrom) +1;
			$evalue = $blastdata[10];
			$evalue =~ s/[0-9.]+e//;
			$evalue =~ s/0\.0/0/;
			$bitscore = $blastdata[11];
			
			print "Processing $file: $queryname                         \r";
			if (exists $count{$this}) {
				$count{$this}++;
				${$data{$this}}[2] += $identity;
				${$data{$this}}[3] += $alignmentlength;
				${$data{$this}}[4] += $mismatch;
				${$data{$this}}[5] += $gapopen;
				${$data{$this}}[6] += $querylength;
				${$data{$this}}[7] += $hitlength;
				if (${$data{$this}}[8] < $evalue) {
					${$data{$this}}[8] = $evalue;
				}
				${$data{$this}}[9] += $bitscore;
			} else {
				$count{$this} = 1;
				${$data{$this}}[0] = $queryname;
				${$data{$this}}[1] = $hitname;
				${$data{$this}}[2] = $identity;
				${$data{$this}}[3] = $alignmentlength;
				${$data{$this}}[4] = $mismatch;
				${$data{$this}}[5] = $gapopen;
				${$data{$this}}[6] = $querylength;
				${$data{$this}}[7] = $hitlength;
				${$data{$this}}[8] = $evalue;
				${$data{$this}}[9] = $bitscore;
			}
		}

		my %newdata;
		foreach my $key (sort(keys(%data))) {
			my $avidentity = ${$data{$key}}[2] / $count{$key};
			${$data{$key}}[2] = $avidentity;
			my $myevalue = "1e".${$data{$key}}[8];
			$myevalue =~ s/^1e0$/0/;
			${$data{$key}}[8] = $myevalue;
			my $avbitscore = ${$data{$key}}[9] / $count{$key};
			${$data{$key}}[9] = $avbitscore;
			$avbitscore = 1000-$avbitscore;
			my $newkey = ${$data{$key}}[0]."_".$avbitscore;
			@{$newdata{$newkey}} = @{$data{$key}};
		}

		print "Processing $file: done                     \n";		
		foreach my $key (sort(keys(%newdata))) {
			#print OUT "\n$key";
			print OUT "\n${$newdata{$key}}[0]";
			print OUT "\t${$newdata{$key}}[1]";
			print OUT "\t${$newdata{$key}}[2]";
			print OUT "\t${$newdata{$key}}[3]";
			print OUT "\t${$newdata{$key}}[4]";
			print OUT "\t${$newdata{$key}}[5]";
			print OUT "\t0\t${$newdata{$key}}[6]";
			print OUT "\t0\t${$newdata{$key}}[7]";
			print OUT "\t${$newdata{$key}}[8]";
			print OUT "\t${$newdata{$key}}[9]";
		}
		close OUT;
	}	
}
closedir (DIR);
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
		push (@fileData, $_);
	}
	close GET_DATA;
	return @fileData;
}
