#!/bin/sh -i
#
module load StdEnv/2020
module load gcc/9.3.0
module load blast+/2.12.0
module load python/3.7
module load feht
pip install biopython
pip install pandas
pip install openpyxl
pip install os
#
#1: Assemblies.fasta, 2: ReferenceGene.fasta, 3: metadata.xsls
#
#Run local BLAST of reference gene on multifasta assembly
#Save results as TSV. Results will include sequence aligned hits. These are gene variants
makeblastdb -in $1 -dbtype nucl
for file in *.fna
do
	gene="${file%.*}"
	mkdir $gene
#makeblastdb -in $1 -dbtype nucl #
	blastn -query $file -db $1 -out ./$gene/BLASTresults.tsv -outfmt "7 sseqid length sstart send pident qcovs sseq" -max_target_seqs 300000
#
#Filter hits where %identity is at least 80 and coverage is at leat 95
#Write new multifasta file with these filtered hits. This fasta has all variants of gene
	awk '$5 > 80 && $6 >95 {print ">" $1 "\n" $7}' ./$gene/BLASTresults.tsv | tail -n +5 > ./gene_variants.fasta
#
#Save the sequence ids of hits (i.e, contigs with refgene) from blast results as txt file. This txt file is used to assess PAV
	grep ">" $1 | cut -d '_' -f1 | cut -d '>' -f2 | sort | uniq > ./present_and_absent.txt
	grep ">" ./gene_variants.fasta | cut -d '_' -f1 | cut -d '>' -f2 | sort | uniq > ./present.txt #
#
#cp $2 $gene/$2
#cd $gene/
	python3 binmat_builder.py $2 #present_and_absent.txt present.txt gene_variants.fasta
#
	mv *.txt ./$gene/
	mv gene_variants.fasta ./$gene/
#rm $2
#find . -type f ! -name '*.txt' -delete
#rm present.txt
#rm present_and_absent.txt
#cd ../
done
mkdir RefSeqs
mv *.fna RefSeqs/
#feht -i metadata-PAV.txt -d pav-data.txt > $gene/fehtPAV.txt
#feht -i metadata-SAV.txt -d sav-data.txt > $gene/fehtSAV.txt
#feht -i metadata-SAV.txt -d allelic-data.txt > $gene/fehtAllelic.txt
