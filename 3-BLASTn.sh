
SPECIES=('PL' 'EH' 'GS' 'FG' 'GC' 'NN' 'EG' 'PC' 'HL' 'TA' 'EM' 'CC' 'SM' 'CB' 'AF' 'PA' 'PP' 'AP' 'EC' 'OH');

NSPECIES=${#SPECIES[@]}
for (( i=0;i<$NSPECIES;i++)); do
	j=i+1
	for (( j=i+1;j<$NSPECIES;j++)); do	
		echo "blasting ${SPECIES[${i}]} vs ${SPECIES[${j}]}"
		./ncbi-blast-2.9.0+/bin/blastn -db 0-${SPECIES[${j}]}/${SPECIES[${j}]} -query 0-${SPECIES[${i}]}/${SPECIES[${i}]}.fna -out 4-blastout/${SPECIES[${i}]}2${SPECIES[${j}]}.out -soft_masking F -evalue 1e-15 -max_target_seqs 10 -best_hit_overhang 0.18 -outfmt 6
		echo "blasting ${SPECIES[${j}]} vs ${SPECIES[${i}]}"
		./ncbi-blast-2.9.0+/bin/blastn -db 0-${SPECIES[${i}]}/${SPECIES[${i}]} -query 0-${SPECIES[${j}]}/${SPECIES[${j}]}.fna -out 4-blastout/${SPECIES[${j}]}2${SPECIES[${i}]}.out -soft_masking F -evalue 1e-15 -max_target_seqs 10 -best_hit_overhang 0.18 -outfmt 6
	done
done
