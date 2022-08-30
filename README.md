# penguins
Scripts to replicate analyses done in Pirri et al. 2202


Pipeline to obtain orthogroups using a Reciprocal Best Hit (RBH) approach, align them, clean the alignment and perform PAML tests
lino.ometto@unipv.it

All scripts, whether in perl (.pl) or bash (.sh), must be run within the folder where this _README_ is located. They should be run in the order given by the number preceding their name (e.g. 1-..., 2-..., etc).

The pipeline has 'variables' (folder names, species names, BLAST parameters, etc) set for the analysis presented in the associated manuscript.
* Be careful to change the 'variables' in the various scripts to suit your needs *


1) Put the .fasta files containing the CDS of the different species in the folder '00-all_cds_fasta'. 


2) Use the '1-parsefasta.pl' script ==>
* IN: the .fasta files contained in '00-all_cds_fasta'.
* OUT: folders (one per species) with fasta and txt files containing CDS and a log file with old and new CDS names
Only CDS files in the variable list 'my @mycds = (......)' will be taken into account.
Output files and folders are created using the matches given in the variable '%speciesshortname':
e.g.: 'mySpecies' => 'MYSPEC' indicates that CDS contained in the file '00-all_cds_fasta/mySpecies.fasta' will be rewritten in the file 'O-MYSPEC/MYSPEC.fasta'	
	0-MYSPEC
		|__  _mySpecies_has_been_renamed_MYSPEC_
		|__  MYSPEC.cds_renaming_log.txt
		|__  MYSPEC.fasta
		|__  MYSPEC.txt
		
Individual CDSs will be renamed as '>MYSPEC_X', where X is a number ranging from 1 to N, where N=number of CDSs, as indicated in the log file 'MYSPEC.cds_renaming_log.txt'.
You may decide to filter CDSs by length, e.g. so as not to use fragments or CDSs that are unlikely to be used; if so, change the variable '$minCDSlength' (default is 150 bp)


3) Use the script '2-make_blast_database.sh' to create BLAST databases (you must have BLAST+ installed on your computer) ==>
* IN: the .fasta files contained in the 0-.... folders indicated in the variable 'SPECIES=(.....)'.
* OUT: all database files needed to use BLAST
	

4) Make reciprocal BLASTn between all species with the '3-BLASTn.sh' script ==>
* IN: the .fasta files contained in the 0-... folders indicated in the variable 'SPECIES=(.....)'.
* OUT: tab-delimited BLASTn outputs in the '4-blastout' folder, e.g. 'MYSPEC2YOURSPEC.out' contains the BLAST hits with MYSPEC CDS as query and YOURSPEC as database
If you want you can also use tBLASTx (in case modify the script) and modify the arguments used by BLAST


5) Use the '5-reduceblastout.pl' script ==>
* IN: BLAST results contained in '4-blastout'. 
* OUT: tables with BLAST hits where redundant results (same hits) have been removed
Output files are saved in '6-reducedblastout', e.g. 'MYSPEC2YOURSPEC_my.out'.


6) Use the script '7-getseqs_star.pl' or '7-getseqs_web.pl' to identify sets of sequences that should correspond to orthologues ==>
These scripts identify sequences that are reciprocal best hits (RBH) between species. The difference between the two is that whereas with the script '7-getseqs_star.pl' the putative orthologues will be considered those CDSs that are RBHs between a focal species and a second species, with '7-getseqs_web.pl' only CDSs where each CDS is the RBH between all pairwise comparisons are considered orthologues.
Specify which species are to be analysed in 'my @myspecies = ("MYSPEC ", " YOURSPEC ", " THEIRSPEC"); # focal species must be first species'.
In case you wish to modify search parameters: 
my $fractionidentity_threshold = 0.7;
my $fractionaligned_threshold = 0.6;
they specify the minimum values of identity and fraction of aligned CDS to define a hit as a potential orthologue (as from BLAST output)
* IN: BLAST outputs in '6-reducedblastout' and CDS sequences in folders 0-...
* OUT: a folder in '8-seqs' with a .fasta file for each orthologue group (the name of each file is equal to the name of the CDS of the focal species) and other intermediate files/with info
e.g.:	
	8-seqs
		|__  web3-MYSPEC-YOURSPEC-THEIRSPEC
					|__  _matches_taken.txt
					|__  _matches.txt
					|__  fasta
					|		|__ MYSPEC_974.orthologs.fas ...
					|
					|__  reports
							|__  YOURSPEC-MYSPEC.report.txt ...
											
In addition, a file is written with commands to perform alignments with mafft or muscle (see variable 'my $aligner = "mafft";'; you must have mafft and/or muscle installed on your computer):
e.g. 9-align_web3-MYSPEC-YOURSPEC-THEIRSPEC.mafft.sh
and a folder is created where the alignments are saved in '10-alignments', in this example '10-alignments/web3-MYSPEC-YOURSPEC-THEIRSPEC'


7) Use the script '9-align_.... .sh' to do the alignments of the orthologues ==>
* IN: e.g. sequences in '8-seqs/web3-MYSPEC-YOURSPEC/'
* OUT: e.g. alignments saved in '10-alignments/web3-MYSPEC-YOURSPEC-THEIRSPEC/'


8) Parse the alignments with '11-parse_alignments.pl' (the sequences are put in one line for convenience and the names are re-set if they have been changed by the alignment program) ==>
Specify in the script the name of the folder in which the alignments are contained (see e.g. my $mydir = "web3-MYSPEC-YOURSPEC-THEIRSPEC";)
* IN: e.g. alignments saved in '10-alignments/web3-MYSPEC-YOURSPEC-THEIRSPEC/...mafft.fasta'.
* OUT: e.g. alignments saved in '12-alignments_parsed/web3-MYSPEC-YOURSPEC-THEIRSPEC/...mafft.fasta'
a file is also written, e.g. '12-alignments_parsed/web3-MYSPEC-YOURSPEC-THEIRSPEC.ortho_names.txt', with the new and original names of the orthologues
	

9) with the script '13-trim_and_degap.pl'  ==>
the alignments are trimmed and cleaned up according to the sequence of the focal species (which is then assumed to be in phase and of length multiple of three, hence made from intact codons). In addition, the portions of the alignment where there are gaps are removed (this is because there may be the unlikely possibility that certain sequences contain partial codons due to errors in sequencing and thus gaps of length 1 or 2 would occur in the alignment). Note that all PAML analyses are normally only done on the portions of the alignment where there are no gaps ('clean = 1' option), so no information is lost for PAML analyses. In case modify the script to change this behaviour.
Specify in the script the name of the folder in which the alignments are contained (see e.g. my $mydir = "web3-MYSPEC-YOURSPEC-THEIRSPEC";)
* IN: e.g. alignments saved in '12-alignments_parsed/web3-MYSPEC-YOURSPEC-THEIRSPEC'.
* OUT: a) e.g. alignments saved in '14-translatorx_in/web3-MYSPEC-YOURSPEC-THEIRSPEC'
	b) script file to make alignments with 'translatorX', e.g. 15-prankcommands_web3-MYSPEC-YOURSPEC-THEIRSPEC.sh	


10) To align with TranslatorX the orthologues use the script 15-prankcommands_....sh ==>
Specify in the script the name of the folder in which the alignments are contained (see e.g. my $mydir = "web3-MYSPEC-YOURSPEC-THEIRSPEC";)
By default, alignments are made using PRANK and are based on the amino acid sequence corresponding to the nucleotide sequence.
This script makes use of the perl script that can be downloaded at: http://translatorx.co.uk/
You must have prank installed on your computer: http://wasabiapp.org/software/prank/
If you want to use another alignment program, you have to modify the script '13-trim_and_degap.pl' on line 112:
print PRANK "perl $tooldir/translatorx_mod.pl -i $degappedaligndir/$name -o $prankoutdir/$name_out -p P\n";
where 'P' in the parameter "-p P" indicates that prank is used; if "-p M" is put, muscle is used, if "-p F" is put, mafft is used
* IN: e.g. alignments saved in '14-translatorx_in/web3-MYSPEC-YOURSPEC-THEIRSPEC'.
* OUT: e.g. alignments saved in '16-translatorx_out/web3-MYSPEC-YOURSPEC-THEIRSPEC'
	 

11) TranslatorX produces many intermediate files, if you want to remove them use the script '17-cleantxout.sh'.  ==>
Remember to change the name of the folder where the alignments are, e.g. 'folder=web3-MYSPEC-YOURSPEC-THEIRSPEC'
	 

12) In some CDS, internal stop codons may be present, either due to the sequence being a pseudogene, an error in the sequencing, or due to an error in the alignments+trimming mentioned above. 
To be on the safe side, and because PAML approach assumes the sequences represent functional coding sequences, alignments containing such sequences can be removed from the subsequent analysis, but in case they can be used by replacing the TAA/TGA/TAG codon with NNN.
Use the script '18-check4stops.pl'  ==>
to delete or simply return alignments with internal stop codons to the terminal
In the script change the variable '$wantstops' to delete ('my $wantstops = "no";') or not ('my $wantstops = "yes";') these alignments.
Specify in the script the name of the folder in which the alignments are contained (see e.g. my $mydir = 'web3-MYSPEC-YOURSPEC-THEIRSPEC';)
* IN: e.g. alignments saved in '16-translatorx_out/web3-MYSPEC-YOURSPEC-THEIRSPEC'.
* OUT: e.g. alignments saved in '19-alignments_nostops/web3-MYSPEC-YOURSPEC-THEIRSPEC'.


13) to eliminate problematic regions of the alignments (i.e. containing too many differences to be considered "true" orthologous positions) use the script '20-mastrolino.pl'  ==>
This is to prevent PAML from identifying as fast evolving sites positions that in reality are only the result of a poorly done alignment (which could also be the result of incorrect sequences, so the starting dataset should be the best possible to avoid bias in the analysis)
This script analyses the alignment and eliminates the portions that do not meet the following conditions (one can always change these values):
	my $seed = 5; # number of AA in seed (shortest AA stretch in which to analyse divergence); for DNA in seed -> AA seed x 3
	my $diffAAseed = 3; # max number of different AA allowed in seed
	my $diffDNAseed = 8; # max number of different DNA allowed in seed
	my $diffAAi = 0.6; # fraction of differences allowed in AA seq
	my $diffDNAi = 0.4; # fraction of differences allowed in DNA seq
	my $minAAlengt = 50; # min length (in AA) of the alignment
	my $minCons = 0.9; # min fraction of the alignment after cleaning
	my $myoutgroup = "MYSPEC"; # seq of this species will be ignored in flagging
The $seed is the initial length of the portion of the alignment where the number of differences between sequences is to be seen (see parameters $diffDNAseed and $diffAAseed)
If the minimum conditions for considering this portion of the alignment as "OK" are not met, then the length is increased by 1 aa, and the parameters ($diffDNAi and $diffAAi) are tested further.
to prevent a highly divergent sequence (e.g. a distant outgroup) from giving a bias, you can choose to ignore it by entering the name of the outgroup (if you put a home name the outgroup will automatically be included in the analysis)
* IN: e.g. alignments saved in '19-alignments_nostops/web3-MYSPEC-YOURSPEC-THEIRSPEC'
* OUT: e.g. alignments saved in '21-alignments_clean_mastops/web3-MYSPEC-YOURSPEC-THEIRSPEC'
A file is also written, e.g. '21-alignments_clean_mastrolino/web3-MYSPEC-YOURSPEC-THEIRSPEC.log.txt' with the alignment lengths before and after cleaning


14) If you wish to run PAML without removing gapped regions of the alignment whenever they do not constitute a huge portion of the alignment itself (as in our analyses; PAML option 'clean = 0'), use script '21b-cleanalignments.pl' ==>
in case change the following variables:
my $maxgapped = 0.20;  # 0-1; maximum fraction of sequences that can have a gap at a certain position
my @musthave = ('AP', 'AF');  # name of sequences that must always have no gap
my $minlength = 150;  # minimum length of alignment after degapping (should be multiple of 3)
* IN: e.g. alignments saved in '19-alignments_nostops/web3-MYSPEC-YOURSPEC-THEIRSPEC'
* OUT: e.g. alignments saved in '21-alignments_clean_mastrolino/web3-MYSPEC-YOURSPEC-THEIRSPEC_clean'


15) to run PAML analyses use '23-runpaml_clade.pl' and '23-runpaml.pl' ==>
The scripts will produce input files using the alignments 
* IN: e.g. alignments saved in '21-alignments_clean_mastrolino/web3-MYSPEC-YOURSPEC-THEIRSPEC_clean'
* OUT: e.g. PAML outputs saved in '24-pamlout'
PAML analyses will need a phylogenetic tree, as the example given in '22-mytree.tree'


16) the two scripts '25-pamlHarvester.dnds.v2.pl' and '26-pamlHarvester.tests.pl' will take all the PAML output results for each single gene and merge them into text tables
* IN: e.g. PAML outputs saved in '24-pamlout'
* OUT: e.g. tables saved in '27-paml_results/codeml_output'


