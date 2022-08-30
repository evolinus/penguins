#!/usr/bin/perl

use diagnostics;
use warnings;
use strict;

#############################################################################
my $mydir  = "prankete_clean";

#
my @myforegrounds = ("AF", "AP", "(AF,AP)");
#
####

my $alignmentdir = "21-alignments_clean_mastrolino/$mydir";

my $mytree = "22-mytree.tree";
my $clean = 0;
my $foregroundtree = $mytree;

my $pamldir = "24-pamlout/$mydir\_clade";
mkdir ($pamldir, 0777);
my @codemlfiles = ("2NG.dN", "2NG.dS", "2NG.t", "lnf", "rst", "rst1", "rub");

my $myforeground = "";
my $name = "";
my %list;
my %liks;
my %nparamters;
my ($test, $alt, $null, $altdef, $nulldef, $delta, $df, $chisprob, $deltaround);

#
###########################
# only shared models, for foreground-based models see below
my @models = (	"n_M0",   "a_free");

				
# 	parameters lists       	     (treefile,       model; NSsites; fix_omega; omega; description)
$list{"n_M0"}			= [($mytree,         0,      0,        0,       1,	"Neutral model - M0 (same selective pressure across the phylogeny)")];
$list{"a_free"}			= [($mytree,         1,      0,        0,       1,	"Alternative (free-ratio) model (each branch can have different selective pressure)")];
$list{"Sn_M1a"}			= [($mytree,         0,      1,        0,       1,	"Site nearly neutral model - M1a (sites are evolving nearly neutrally)")];
$list{"Sa_M2a"}			= [($mytree,         0,      2,        0,       1,	"Site alternative model - M2a (some sites are evolving nearly neutrally, others are under positive selection)")];
$list{"Sn_M7"}			= [($mytree,         0,      7,        0,       1,	"Site nearly neutral model beta - M7 (sites are evolving nearly neutrally, omega follows a beta distribution)")];
$list{"Sa_M8"}			= [($mytree,         0,      8,        0,       1,	"Site alternative model - M8 (some sites are evolving nearly neutrally, others are under positive selection, omega follows a beta distribution)")];
$list{"CmC_M2a_rel"}	= [($mytree,         0,     22,        0,       1,	"Clade C null model - M2a_rel (some sites are evolving nearly neutrally, others are under positive selection, modified M2a as in Weadick&Chang 2012 MBE)")];
#
###########################

opendir (DIR, $alignmentdir) or die "can't opendir $alignmentdir: $!";
while (defined(my $myalignment = readdir(DIR))) {
	if ($myalignment =~ /fasta/) {
		my $name = $myalignment;
		$name =~ s/.fasta//;
		my $genedir = $name;
		$genedir = "$pamldir/$genedir";
		mkdir ($genedir, 0777);
		$myalignment = "$alignmentdir/$myalignment";
		print "PAML analyses on '$myalignment'\n";

		my $subdir = "_all";
		my $mydir = "$genedir/$subdir";
		mkdir ($mydir, 0777);
		my $nucdir = "$mydir/nuc";
		mkdir ($nucdir, 0777);
		my $ctldir = "$mydir/ctl";
		mkdir ($ctldir, 0777);
		my $outdir = "$mydir/out";
		mkdir ($outdir, 0777);
		my $codemldir = "$mydir/codemlfiles";
		mkdir ($codemldir, 0777);

		my %mysequences;
		%mysequences = fa_reader("$myalignment");

		my $howmany = 0;
		my $alignmentlength;
		for my $species (sort(keys(%mysequences))) {
			$alignmentlength = length (${$mysequences{$species}}[1]);
			$howmany++;
		}

		##################
		# preparing and running PAML
		#
		# prepare nuc file
		my $Nuc = "$name.nuc";
		$Nuc = "$nucdir/$Nuc";
		open (NUC, ">$Nuc");
		print NUC "\t$howmany\t$alignmentlength\n";
		for my $species (sort (keys(%mysequences))) {
			my $seq = ${$mysequences{$species}}[1];
			my $seqlength = length ($seq);
			my $newseq = "";
			for (my $i=0; $i<$seqlength; $i+=3) {
				my $three = substr($seq, $i, 3);
				if ($three eq 'TAG' || $three eq 'TAA' || $three eq 'TGA') {
						$newseq .= "---";
				} else {
					$newseq .= "$three";
				}
			}
			print NUC  "${$mysequences{$species}}[0]\n$newseq\n";
		}
		close NUC;

		# run codeml
		foreach my $model (@models) {
			# make ctl file
			my $CtlFile = "$ctldir/$model.ctl";
			my $OutFile = "$outdir/$model.out";
			my $tree = ${$list{$model}}[0];
			my $mod  = ${$list{$model}}[1];
			my $NSs  = ${$list{$model}}[2];
			my $fixw = ${$list{$model}}[3];
			my $w    = ${$list{$model}}[4];
			open (CTL, ">$CtlFile");
			print CTL "seqfile = $Nuc\ntreefile = $tree\noutfile = $OutFile\n";
			print CTL "noisy = 1\nverbose = 0\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\nclock = 0\naaDist = 0\naaRatefile = dat/lg.dat\n";
			print CTL "model = $mod\nNSsites = $NSs\n";
			print CTL "icode = 0\nMgene = 0\nfix_kappa = 0\nkappa = 2\n";
			print CTL "fix_omega = $fixw\nomega = $w\n";
			print CTL "fix_alpha = 1\nalpha = 0\nMalpha = 0\nncatG = 10\ngetSE = 0\nRateAncestor = 1\nSmall_Diff = .5e-6\ncleandata = $clean";
			close CTL;
			# run codeml
			print "\n#############################################\nRunning codeml $CtlFile\n";
			system("/giorgio/data3/pirri/paml4.8/bin/codeml  $CtlFile") == 0;# or die "codeml failed";
			print "PAML codeml '$model' has finished!        \n#############################################\n\n";
			foreach my $ff (@codemlfiles) {
				system("mv $ff $codemldir/$model.$clean.$ff.txt");
			}

		}

		my $results1file = "$genedir/$name.results.txt";
		open (OUT, ">$results1file");			
		##############
		## dN, dS, omega
		print OUT "\n\nomegafix";
		my $BnOutFile = "$genedir/_all/out/n_M0.out";
		my @data = slurp ($BnOutFile);
		if (!exists $data[1]) {
			print OUT "\nNA";
		} else {
			foreach my $linedata (@data) {
				$linedata =~ s/^\s+//;
				my @line = split (/\s+/, $linedata);
				my $siz =  scalar (@line);
				if ($siz > 1) {
					if ($line[0] =~ /^omega/ ) {
						print OUT "\n$line[3]";
					}
				}	
			}
		}
		print OUT "\n\nbranch\tt\tN\tS\tdN/dS\tdN\tdS\tN*dN\tS*dS";
		my $BaOutFile = "$genedir/_all/out/a_free.out";
		@data = slurp ($BaOutFile);
		if (!exists $data[1]) {
			print OUT "\nNA";
			print OUT "\tNA" x 8;
		} else {
			foreach my $linedata (@data) {
				$linedata =~ s/^\s+//;
				my @line = split (/\s+/, $linedata);
				my $siz = scalar (@line);
				if ($siz > 1) {
					if ($siz > 4) {
						if ($line[3] =~ /^ls/ && $line[5] == 0) {
							print OUT "\nNA";
							print OUT "\tNA" x 8;
						}	
					}
					if ( ($line[0] =~ /[56789]{1}\.\./ || $line[0] =~ /10\.\./) && $line[1] =~ /^[0-9]+\.[0-9]+/ ) {
						print OUT "\n$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]";
					}
				}
			}
		}
		close OUT;

		foreach $myforeground (@myforegrounds) {

			$myforeground =~ s/[\(\)\,]+/ /g;
			$myforeground =~ s/[ ]+/-/g;
			$myforeground =~ s/^-//;
			$myforeground =~ s/-$//;

			###########################
			#
			my @modelsF = (	"Bn.$myforeground",  "Ba.$myforeground",);
				
			# 	parameters lists       	     (treefile,       model; NSsites; fix_omega; omega; description)
			$list{"Bn.$myforeground"} 	   = [($foregroundtree, 2,      0,        1,       1,	"Branch neutral model (same selective pressure across the phylogeny)")];
			$list{"Ba.$myforeground"} 	   = [($foregroundtree, 2,      0,        0,       1,	"Branch alternative model (foreground branch under different selective pressure than the rest of the phylogeny)")];
			$list{"BSn.$myforeground"}     = [($foregroundtree, 2,      2,        1,       1,	"Branch-site neutral model (same selective pressure across sites and across the phylogeny)")];
			$list{"BSa.$myforeground"}     = [($foregroundtree, 2,      2,        0,       1,	"Branch-site alternative model (some sites in the foreground branch are under different (positive) selective pressure than the rest of the phylogeny)")];
			$list{"CmC.$myforeground"} 	   = [($foregroundtree, 3,      2,        0,       1,	"Clade model C (clade has different selective pressure than the rest of the phylogeny)")];
			#
			###########################

			$subdir = "$myforeground";
			$mydir = "$genedir/$subdir";
			mkdir ($mydir, 0777);
			$ctldir = "$mydir/ctl";
			mkdir ($ctldir, 0777);
			$outdir = "$mydir/out";
			mkdir ($outdir, 0777);
			$codemldir = "$mydir/codemlfiles";
			mkdir ($codemldir, 0777);
			my $foregroundtreedir = "$mydir/tree";
			mkdir ($foregroundtreedir, 0777);

			$foregroundtree = "$foregroundtreedir/$myforeground.tree";
			print "\nForeground species '$myforeground'.\n";
			my $outfile4tree = "$genedir/_all/out/n_M0.out";
			taketree ($outfile4tree, $myforeground);

			# run codeml
			foreach my $model (@modelsF) {
				# make ctl file
				my $CtlFile = "$ctldir/$model.ctl";
				my $OutFile = "$outdir/$model.out";
				my $tree = $foregroundtree;
				my $mod  = ${$list{$model}}[1];
				my $NSs  = ${$list{$model}}[2];
				my $fixw = ${$list{$model}}[3];
				my $w    = ${$list{$model}}[4];
				open (CTL, ">$CtlFile");
				print CTL "seqfile = $Nuc\ntreefile = $tree\noutfile = $OutFile\n";
				print CTL "noisy = 1\nverbose = 0\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\nclock = 0\naaDist = 0\naaRatefile = dat/lg.dat\n";
				print CTL "model = $mod\nNSsites = $NSs\n";
				print CTL "icode = 0\nMgene = 0\nfix_kappa = 0\nkappa = 2\n";
				print CTL "fix_omega = $fixw\nomega = $w\n";
				print CTL "fix_alpha = 1\nalpha = 0\nMalpha = 0\nncatG = 10\ngetSE = 0\nRateAncestor = 1\nSmall_Diff = .5e-6\ncleandata = $clean";
				close CTL;
				# run codeml
				print "\n#############################################\nRunning codeml $CtlFile\n";
				system("/giorgio/data3/pirri/paml4.8/bin/codeml  $CtlFile") == 0;# or die "codeml failed";
				print "PAML codeml '$model' has finished!        \n#############################################\n\n";
				foreach my $ff (@codemlfiles) {
					system("mv $ff $codemldir/$model.$clean.$ff.txt");
				}
			}	

			##################
			# write results
			#
			my $resultsfile = "$genedir/$name.$myforeground.results.txt";
			open (OUT, ">$resultsfile");
			print OUT "test\talt\tnull\taltLik\tnullLik\tparam";

			## likelihoods
			foreach my $model (@models) {
				my $what = ${$list{$model}}[5];
				my $tree = ${$list{$model}}[0];
				my $mod  = ${$list{$model}}[1];
				my $NSs  = ${$list{$model}}[2];
				my $fixw = ${$list{$model}}[3];
				my $w    = ${$list{$model}}[4];
				my $OutFile = "$genedir/_all/out/$model.out";
				my @data = slurp ($OutFile);
				if (!exists $data[1]) {
					$liks{$model} = "NA";
					next;
				}
				foreach my $linedata (@data) {
					$linedata =~ s/[\)\(:;,'#]//g;
					$linedata =~ s/^\s+//;
					my @line = split (/\s+/, $linedata);
					my $siz = scalar (@line);
					if ($siz > 1) {
						if ($line[0] =~ /^lnL/) {
							$liks{$model} = $line[4];
							$nparamters{$model} = $line[3];
						}
					}
				}
			}
			foreach my $model (@modelsF) {
				my $what = ${$list{$model}}[5];
				my $tree = $foregroundtree;
				my $mod  = ${$list{$model}}[1];
				my $NSs  = ${$list{$model}}[2];
				my $fixw = ${$list{$model}}[3];
				my $w    = ${$list{$model}}[4];
				my $OutFile = "$genedir/$myforeground/out/$model.out";
				my @data = slurp ($OutFile);
				if (!exists $data[1]) {
					$liks{$model} = "NA";
					next;
				}
				foreach my $linedata (@data) {
					$linedata =~ s/[\)\(:;,'#]//g;
					$linedata =~ s/^\s+//;
					my @line = split (/\s+/, $linedata);
					my $siz = scalar (@line);
					if ($siz > 1) {
						if ($line[0] =~ /^lnL/) {
							$liks{$model} = $line[4];
							$nparamters{$model} = $line[3];
						}
					}
				}
			}

			#########
			## likelihood ratio tests			
			my $testtype = "";
			my $nullmodel = "";
			my $altmodel = "";
			#########
			## site tests
#			$testtype  = "SITE-TEST 12";
#			$altmodel  = "Sa_M2a";
#			$nullmodel = "Sn_M1a";
#			lrt_test ($testtype, $altmodel, $nullmodel);
#			$testtype  = "SITE-TEST 78";
#			$altmodel  = "Sa_M8";
#			$nullmodel = "Sn_M7";
#			lrt_test ($testtype, $altmodel, $nullmodel);

			#########
			## branch tests
 			$testtype  = "BRANCH-TEST";
 			$altmodel  = "Ba.$myforeground";
 			$nullmodel = "n_M0";
			lrt_test ($testtype, $altmodel, $nullmodel);
#			$testtype  = "BRANCH-TEST 2";
#			$altmodel  = "Ba.$myforeground";
#			$nullmodel = "Bn.$myforeground";
#			lrt_test ($testtype, $altmodel, $nullmodel);

			#########
			## clade tests
#			$testtype  = "CLADE-TEST";
# 			$altmodel  = "CmC.$myforeground";
#			$nullmodel = "Sn_M1a";
# 			lrt_test ($testtype, $altmodel, $nullmodel);
# 			$testtype  = "CLADE-TEST 2";
# 			$altmodel  = "CmC.$myforeground";
# 			$nullmodel = "CmC_M2a_rel";
# 			lrt_test ($testtype, $altmodel, $nullmodel);

			#########
			## branch-site tests
# 			$testtype  = "BRANCHSITE-TEST";
# 			$altmodel  = "BSa.$myforeground";
# 			$nullmodel = "BSn.$myforeground";
# 			lrt_test ($testtype, $altmodel, $nullmodel);
			
			close OUT;
			
			my $compressdir = "$genedir/$myforeground";
			system ("tar -zcvf $compressdir.tar.gz $compressdir/");
			system ("rm -r $compressdir/");
		}
		my $compressdir = "$genedir/_all";
		system ("tar -zcvf $compressdir.tar.gz $compressdir/");
		system ("rm -r $compressdir/");
	}
}
##################

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

###############################################################################
#	read fa file: each seq must be in a single line
#
sub fa_reader  {
	my ($infile) = @_;
	my @data = slurp ($infile);
	my $n = scalar @data;
	my $seqinalign = "";

	my %seqs = ();
	for (my $i=0; $i<$n; $i+=2) {
		my $name = $data[$i];
		$name =~ s/>//;
		my $seq = $data[$i+1];
		${$seqs{$name}}[0] = $name;
		${$seqs{$name}}[1] = $seq;
	}
	return %seqs;
}

###############################################################################
# making tree with foreground/background
#
sub taketree {
	my ($infile, $myforeground) = @_;
	open (TREE, ">$foregroundtree");
	my @data = slurp ($infile);
	my $ok = 0;
	if ($myforeground !~ /-/) {
		foreach my $linedata (@data) {
			$linedata =~ s/^\s+//;
			my @line = split (/\s+/, $linedata);
			my $siz = scalar (@line);
			if ($siz > 1) {
				if ($line[0] =~ /^tree/) {
					$ok = 1;
				}
				if ($linedata =~ /$myforeground/ && $ok == 1) {
					$linedata =~ s/$myforeground: ([0-9.]+)/$myforeground: $1 #1/;
					print TREE "$linedata";
					last;
				}
			}
		}
	} else {
		my @leaves = split /-/, $myforeground;
		my $first = $leaves[0];
		my $last = $leaves[-1];
		foreach my $linedata (@data) {
			$linedata =~ s/^\s+//;
			my @line = split (/\s+/, $linedata);
			my $siz = scalar (@line);
			if ($siz > 1) {
				if ($line[0] =~ /^tree/) {
					$ok = 1;
				}
				if ($linedata =~ /$first/ && $linedata =~ /$last/ && $ok == 1) {
					$linedata =~ s/($first.*$last: [0-9.]+[\)]+: [0-9.]+)/$1 #1/; # ((AF: 0.000004, AP: 0.002465 #1): 0.008224,
					print TREE "$linedata";
					last;
				}
			}
		}
	}
	close TREE;
}

###############################################################################
# likelihood ratio test: usage -> lrt_test (name, alt, null);
#
sub lrt_test {
	my ($test, $alt, $null) = @_;
	if ($liks{$alt} ne 'NA' && $liks{$null} ne 'NA') {
		$df = $nparamters{$alt} - $nparamters{$null};
	} else {
		$df = "NA";
	}
	print OUT "\n$test\t$alt\t$null\t$liks{$alt}\t$liks{$null}\t$df";
 }
