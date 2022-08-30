
SPECIES=('PL' 'EH' 'GS' 'FG' 'GC' 'NN' 'EG' 'PC' 'HL' 'TA' 'EM' 'CC' 'SM' 'CB' 'AF' 'PA' 'PP' 'AP' 'EC' 'OH');

NSPECIES=${#SPECIES[@]}
for (( i=0;i<$NSPECIES;i++)); do
	./ncbi-blast-2.9.0+/bin/makeblastdb -in 0-${SPECIES[${i}]}/${SPECIES[${i}]}.fna -dbtype nucl -parse_seqids -out 0-${SPECIES[${i}]}/${SPECIES[${i}]} -title "${SPECIES[${i}]} transcripts"
done
