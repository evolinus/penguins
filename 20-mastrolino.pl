#!/usr/bin/perl 

use diagnostics;
use strict;
use warnings;
######################################################################

my $seed		= 5;	# number of AA in seed (shortest AA stretch in which to analyse divergence); for DNA in seed -> AA seed x 3
my $diffAAseed  = 3;	# max number of different AA  allowed in seed
my $diffDNAseed = 8;	# max number of different DNA allowed in seed

my $diffAAi  = 0.6;		# fraction of differences allowed in AA seq
my $diffDNAi = 0.4;		# fraction of differences allowed in DNA seq

my $minAAlengt = 50;	# min length (in AA) of the final alignment (if <$minAAlengt, alignment is discarded)
my $minCons = 0.9;		# min fraction of the alignment after cleaning (if <$minCons, alignment is discarded)

my $myoutgroup = "OH";	# seq of this species will be ignored in flagging
#my $myoutgroup = "lino";	# seq of this species will be ignored in flagging

my $mydir  = "star20-AF-EH-GS-FG-GC-NN-EG-PC-HL-TA-EM-CC-SM-CB-PL-PA-PP-AP-EC-OH";

my $indir = "19-alignments_nostops/$mydir";
my $outdir = "21-alignments_clean_mastrolino/$mydir";
mkdir ($outdir, 0777);
open (LOG, ">21-alignments_clean_mastrolino/$mydir.log.txt");
print LOG "number of AA in seed (shortest AA stretch in which analyse divergence); for DNA in seed -> AA seed x 3\n";
print LOG "$seed";
print LOG "\nmax number of different AA allowed in seed\n";
print LOG "$diffAAseed";
print LOG "\nmax number of different DNA allowed in seed\n";
print LOG "$diffDNAseed";
print LOG "\nfraction of differences allowed in AA seq\n";
print LOG "$diffAAi";
print LOG "\nfraction of differences allowed in DNA seq\n";
print LOG "$diffDNAi";
print LOG "\nmin length (in AA) of the alignment\n";
print LOG "$minAAlengt";
print LOG "\nfile\toriginallength\tnewlength\tremovedDNA\tremovedAA\t%retained\tretained";

#
my $cnt = 0;
opendir my $dh, $indir or die "Can't open $indir: $!";
while (my $filename = readdir($dh)) {
    $cnt++ if $filename =~ /fas/;
}
my $c = 0;
opendir (DIR, $indir) or die "can't opendir $indir: $!";
ALIGNMENT: while (defined (my $file = readdir(DIR))) {
	if ($file =~ /.fasta$/) {
		$c++;
		print "Processing $file ($c of $cnt)    \r";
		print LOG "\n$file";
		my $infile = $indir."/".$file;	# DNA input sequence file
		open (IN, $infile) or die print "Cannot open file \"$infile\"\n\n";
		$/ = "\n>";
		my $count = 0;
		my $lengthAA = 0;
		my %seqsAA;
		my %seqsDNA;
		my %seqsCODON;
	
		# get sequences and translate
		while (<IN>) {
			$count++;
			(my $species, my @seqseq) = split /\n/, $_;
			chomp $species;
			$species =~ s/>//g;
			my $seqDNA = join('',@seqseq);
			$seqDNA =~ s/ //g;
			my $lengthDNA = length $seqDNA;
			if ($lengthDNA < (3*$minAAlengt) ) {
				next ALIGNMENT;
			}
			my $seqAA = '';
			my $codon;
			for (my $i=0; $i<$lengthDNA-2; $i+=3) {
				$codon = substr($seqDNA,$i,3);
				my $aaaa = codon2aa($codon);
				$seqAA .= codon2aa($codon);
				#print "\n$codon $aaaa";
			}
			my @seqDNAq = $seqDNA =~ /.../g;
			$seqsCODON{$species} = [@seqDNAq];
			@seqDNAq = $seqDNA =~ /./g;
			$seqsDNA{$species} = [@seqDNAq];
			my @seqAAq = split "", $seqAA;
			$seqsAA{$species} = [@seqAAq];
			$lengthAA = length $seqAA;
			#print "\n$species\n$seqAA\n$seqDNA\n";
			#print "$species$lengthDNA\n";
		}		
		close IN;	
		$file =~ s/mafft.fasta/.fasta/;
		$file =~ s/.out.nt_ali.fasta/.fasta/;
		open (OUT, ">$outdir/$file");
		my $lengthDNA = 3*$lengthAA;
		print LOG "\t$lengthDNA";
		
		# make majority consensus for AA and DNA
		#print "\n\n\nlength AA: $lengthAA \t DNA: $lengthDNA\n";
		#print "cons\t";
		my @consensusAA = ();
		for (my $i=0; $i<$lengthAA; $i++) {
			my %countAA;
			foreach my $key (sort(keys(%seqsAA))) {
				my $myAA = ${$seqsAA{$key}}[$i];
				if (exists $countAA{$myAA}) {
					$countAA{$myAA}++;
				} else {
					$countAA{$myAA} = 1;
				}
			}
			my $cons = (sort {$countAA{$b} <=> $countAA{$a}} keys %countAA)[0];
			push (@consensusAA, $cons);
			#print "$cons";
		}
		#print "\n\n\nlength AA: $lengthAA \t DNA: $lengthDNA\n";
		my @consensusDNA = ();
		for (my $i=0; $i<$lengthDNA; $i++) {
			my %countDNA;
			foreach my $key (sort(keys(%seqsDNA))) {
				my $myDNA = ${$seqsDNA{$key}}[$i];
				if ($myDNA ne '-') {
					if (exists $countDNA{$myDNA}) {
						$countDNA{$myDNA}++;
					} else {
						$countDNA{$myDNA} = 1;
					}
				}	
			}
			my $cons = (sort {$countDNA{$b} <=> $countDNA{$a}} keys %countDNA)[0];
			push (@consensusDNA, $cons);
		}
	
		# use sliding window (length = $seed) to flag, in each sequence, 
		# those AA different from consensus (no more than $diffAA and $diffDNA differences are allowed in the window)
		my %masker;
		foreach my $key (sort(keys(%seqsAA))) {
			if ($key ne $myoutgroup) { # masker raccoglie tutte le sequenze, ma non outgroup; usare nome di outgroup farlocco se si vuole includerlo
				my $concatenateseq = join "", @{$seqsAA{$key}};
				#print "\n\n$key\t$concatenateseq";
				my @maskedseq = ();
				for (0..($lengthAA-1)){
					push @maskedseq, 0;
				}
		
				my $diffAA  =  $diffAAi;
				my $diffDNA =  $diffDNAi;			
				my $countdiffAA = 0;
				my $countdiffDNA = 0;

				for (my $i=0; $i<4; $i++) {
					my $myAA = ${$seqsAA{$key}}[$i];
					my $conAA = $consensusAA[$i];
					if ($myAA ne $conAA && $myAA ne '-') {
						$countdiffAA++;
					}
				}	
				for (my $j=0; $j<12; $j++) {
					my $myDNA = ${$seqsDNA{$key}}[$j];
					my $conDNA = $consensusDNA[$j];
					if ($myDNA ne $conDNA && $myDNA ne '-') {
						$countdiffDNA++;
					}	
				}
				if ($countdiffAA >= 3 && $countdiffDNA >= 6) {
					$maskedseq[0] = 1;
					$maskedseq[1] = 1;
					$maskedseq[2] = 1;
				}						

				$countdiffAA = 0;
				$countdiffDNA = 0;
				for (my $i=$lengthAA-3; $i<$lengthAA; $i++) {
					my $myAA = ${$seqsAA{$key}}[$i];
					my $conAA = $consensusAA[$i];
					if ($myAA ne $conAA) {
						$countdiffAA++;
					}
				}	
				for (my $j=$lengthDNA-12; $j<$lengthDNA; $j++) {
					my $myDNA = ${$seqsDNA{$key}}[$j];
					my $conDNA = $consensusDNA[$j];
					if ($myDNA ne $conDNA) {
						$countdiffDNA++;
					}	
				}
				if ($countdiffAA >= 3 && $countdiffDNA >= 6) {
					$maskedseq[$lengthAA-1] = 1;
					$maskedseq[$lengthAA-2] = 1;
					$maskedseq[$lengthAA-3] = 1;
				}						
		
				$countdiffAA = 0;
				$countdiffDNA = 0;
				STEP: for (my $i=3; $i<$lengthAA-$seed; $i++) {
					#print "\nAA $i\t";
					$countdiffAA = 0;
					$countdiffDNA = 0;
					my $myAA = ${$seqsAA{$key}}[$i];
					my $conAA = $consensusAA[$i];
					my $conDNA = "";
			
					if ($myAA ne $conAA) {
						for (my $j=$i; $j<$i+$seed; $j++) {
							$myAA = ${$seqsAA{$key}}[$j];
							$conAA = $consensusAA[$j];
							if ($myAA ne $conAA && $myAA ne '-') {
								$countdiffAA++;
							}
						#print "$myAA.$conAA ";
						}	

						for (my $j=($i*3); $j<($i*3)+($seed*3); $j++) {
							my $myDNA = ${$seqsDNA{$key}}[$j];
							my $conDNA = $consensusDNA[$j];
							if ($myDNA ne $conDNA && $myDNA ne '-') {
								$countdiffDNA++;
							}
						#print "dna $myDNA.$conDNA ";
						}	
					} else {
						next STEP;
					}
					#print "\ndiffAA: $countdiffAA\tdiffDNA: $countdiffDNA";
			
					if ($countdiffAA >= $diffAAseed && $countdiffDNA >= $diffDNAseed) {
						my $from = $i;
						my $to = $i+1;
						#print "\nSliding from $from\n";
						while ($from < $lengthAA) {
							$myAA = ${$seqsAA{$key}}[$to];
							$conAA = $consensusAA[$to];
							if ($myAA ne $conAA && $myAA ne '-') {
								$countdiffAA++;
							}
							#print "$to\t$myAA\t$conAA\n";
							my $nAA = $to-$from+1;
					
							for (my $d=($to*3); $d<($to*3)+3; $d++) {
								my $myDNA = ${$seqsDNA{$key}}[$d];
								my $conDNA = $consensusDNA[$d];
								if ($myDNA ne $conDNA && $myDNA ne '-') {
									$countdiffDNA++;
								}						
							}
							my $nDNA = $nAA*3;
				
							my $fractdiffAA  = $countdiffAA  / $nAA;	
							my $fractdiffDNA = $countdiffDNA / $nDNA;	
				
							my $aa1  = ${$seqsAA{$key}}[$to+1];
							my $con1 = $consensusAA[$to+1];
							my $aa2  = ${$seqsAA{$key}}[$to+2];
							my $con2 = $consensusAA[$to+2];
							my $aa3  = ${$seqsAA{$key}}[$to+3];
							my $con3 = $consensusAA[$to+3];
							if ( ($aa1 eq $con1 && $aa2 eq $con2 && $aa3 eq $con3) 
									|| ($aa1 eq $con1 && $aa2 eq $con2 && $aa3 ne $con3) 
									|| ($fractdiffAA < $diffAA && $fractdiffDNA < $diffDNA) 
									|| $to >= $lengthAA-4) {
								if ($fractdiffAA < $diffAA && $fractdiffDNA < $diffDNA) {
									$to--;
								}
								#print "\nMasking from $from to $to";
								if ($to >= $lengthAA-4) {
									$to+=2;
								}
								for (my $k=$from; $k<($to+1); $k++) {
									$myAA = ${$seqsAA{$key}}[$k];
									$conAA = $consensusAA[$k];
									if ($myAA ne $conAA) {
										$maskedseq[$k] = 1;
									} elsif ($k == $to && $myAA eq $conAA) {
										$maskedseq[$k] = 1; 
									#} elsif ($k == ($to-1) && $myAA eq $conAA && ${$seqsAA{$key}}[$k+1] ne $consensusAA[$k+1]) {
									#	$maskedseq[$k] = 3;
									} else {
										$maskedseq[$k] = 1; # set to 0 if isolated conserved AA are to be retained
									}						
								}
								next STEP;
							}
							$to++;	
						}
					}		
				}
				$masker{$key} = [@maskedseq];
			}
		}	
		#print "\n";	
		
		# make a 0's masker for te outgroup
		my @maskedseq = ();
		for (my $n=0; $n<$lengthAA; $n++) {
			$maskedseq[$n] = 0;
		}
		$masker{$myoutgroup} = [@maskedseq];
			
		# clean sequences based on flagged AA:
		# AA are flagged if AA is flagged in at least one sequence
		# cleaned sequences are back-translated into DNA
		my %cleanedseq;
		my $consensusseq = "";
		for (my $i=0; $i<$lengthAA-1; $i++) {    # $lengthAA-1 is to avoid stop codons
			my $countmasked;
			foreach my $species (sort(keys(%masker))) {
				$countmasked += $masker{$species}[$i];
			}
			if ($countmasked == 0) {
				foreach my $species (sort(keys(%masker))) {
					my $codon = $seqsCODON{$species}[$i];
					$cleanedseq{$species} .= $codon;
				} 
				$consensusseq .= "...";
			} else {
				foreach my $species (sort(keys(%masker))) {
					my $codon = $seqsCODON{$species}[$i];
# keep following line only while testing
					#$cleanedseq{$species} .= $codon; 
				}
				$consensusseq .= "xxx";
			}
		}
		my $cleanedalilength = 0;
		foreach my $species (sort(keys(%cleanedseq))) {
			#my $shortspecies = substr ($species, 0, 3);
			#print OUT ">$shortspecies\n$cleanedseq{$species}\n";
			print OUT ">$species\n$cleanedseq{$species}\n";
			$cleanedalilength = length $cleanedseq{$species};
		}
		print LOG "\t$cleanedalilength";
		my $diffDNA = $lengthDNA - $cleanedalilength;
		my $diffAA = $diffDNA/3;
		print LOG "\t$diffDNA";
		print LOG "\t$diffAA";
		my $retained = $cleanedalilength/$lengthDNA;
		my $percRetained = int($retained*100);
		print LOG "\t$percRetained";
		if ($retained >= $minCons) {
			print LOG "\tyes";
		} else {
			print LOG "\tno";
		}
		

# keep below only if testing

# 		foreach my $species (sort(keys(%masker))) {
# 			my $shortspecies = substr ($species, 0, 3);
# 			print OUT "\n>$shortspecies\n";
# 			foreach my $n (@{$masker{$species}}) {
# 				print OUT "$n$n$n";
# 			}
# 		}
# 		print OUT "\n>con\n$consensusseq";

# keep above only if testing

		close OUT;
		if ($cleanedalilength < (3*$minAAlengt) && $retained < $minCons) {
			system ("rm $outdir/$file");
			next ALIGNMENT;
		}
	}
}	
print "\n";
close LOG;
exit;

########################################################
sub codon2aa {
	my ($codon) = @_;
	$codon = uc $codon;
	my (%g) = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
	
	if (exists $g{$codon})   { return $g{$codon}; }
	elsif ($codon = ~/-/g)   { return '-'; }
	elsif ($codon = ~/GC./i) { return 'A'; }
	elsif ($codon = ~/GG./i) { return 'G'; }
	elsif ($codon = ~/CC./i) { return 'P'; }
	elsif ($codon = ~/AC./i) { return 'T'; }
	elsif ($codon = ~/GT./i) { return 'V'; }
	elsif ($codon = ~/CG./i) { return 'R'; }
	elsif ($codon = ~/TC./i) { return 'S'; }
	else                     { return 'x';
	}
}
