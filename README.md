# KN9204 Annotation

#### Input files
````
# TriAnnot file.
id_KN9204_annotation_only.bed
# RNA-seq file.
id_KN9204_rnaseq_only.bed
# ISO-seq file.
KN9204_isoseq_only.bed
# Transposable element file.(chr1-chr7, subgenome A, B, D and Un)
TE_KN9204.*.txt
# reference genome. (chr1-chr7, subgenome A, B, D and Un)
KN9204.*.fasta
# CDS blast database of Chinese Spring
CS_CDS_db
# protein blast database of Chinese Spring
CS_PROTEIN_db
````

#### Prerequisite
- python-2.7  (https://www.python.org/)
- TAMA   (https://github.com/GenomeRIK/tama)
- ncbi-blast-2.7.1 (https://ftp.ncbi.nlm.nih.gov/blast/executables/)
- genometools-1.5.9 (http://genometools.org/)

#### Usage
```shell script
export workspace_zhu=/home/sczhuhd/workspace_zhu
cd $workspace_zhu/workspace_tama/merge/
# generate filelist.txt
echo -e "id_KN9204_annotation_only.bed  capped  3,3,3   annotation\nKN9204_isoseq_only.bed capped  1,1,1   iso_seq\nid_KN9204_rnaseq_only.bed  no_cap  2,2,2   rna_seq" > filelist.txt
# merge annotation files
python $workspace/workspace_tama/tama-master/tama_merge.py -a 100 -z 100 -m 50 -f filelist.txt -p KN_9204 -d merge_dup
cd $workspace/workspace_final_rna/remove_TE/
# reformat the TE file.
cat $workspace/workspace_tama/merge/KN_9204.bed |awk '{print $1,$2,$3}' > $workspace/workspace_final_rna/remove_TE/KN_9204.bed
# export TE coverage of each gene.
gcc search.c -o searchTE
$workspace/workspace_final_rna/remove_TE/searchTE
cd $workspace/workspace_bed2gff/origin_KN9204/
cp $workspace/workspace_final_rna/remove_TE/result.txt ./
cp $workspace/workspace_tama/merge/KN_9204.bed ./
cp $workspace/workspace_tama/merge/KN_9204_gene_report.txt ./
cp $workspace/workspace_tama/merge/KN_9204_trans_report.txt ./
# convert bed to gff
$workspace/workspace_genometools/genometools-1.5.9/bin/gt bed_to_gff3 -featuretype gene -blocktype exon $workspace/workspace_tama/merge/KN_9204.bed > KN_9204.bed.gff3
# reformat the gff file
python gff2newgff.py KN_9204.bed.gff3
# translate all genes
python bed2gff.py $workspace/workspace_tama/merge/KN_9204.bed
mv IGDB_v1.0_KN_9204_201811_cds.fasta cds_origin_KN9204.fasta
mv IGDB_v1.0_KN_9204_201811_pep.fasta prot_origin_KN9204.fasta
# blast CS_CDS_db
/home/sczhuhd/software/mysoftware/ncbi-blast-2.7.1+/bin/blastn -db $workspace/workspace_final_rna/blast/CS_CDS_db -query cds_origin_KN9204.fasta -out cds_origin_KN9204.fasta.outfmt6 -evalue 0.00001 -max_target_seqs 5 -num_threads 20 -outfmt 6
# format blast result.
python parse_prot.py prot_origin_KN9204.fasta;
# blast CS_PROTEIN_db
/home/sczhuhd/software/mysoftware/ncbi-blast-2.7.1+/bin/blastp -db $workspace/workspace_final_rna/blast/CS_PROTEIN_db -query match_prot_origin_KN9204.fasta -out match_prot_origin_KN9204.fasta.out -evalue 0.00001 -max_target_seqs 5 -num_threads 20 -outfmt 6
# export hign confidence and low confidence genes.
python final_report.py
cp KN9204_LC.bed $workspace/workspace_bed2gff/
cp KN9204_HC.bed $workspace/workspace_bed2gff/
cd $workspace/workspace_bed2gff/
# translate all high confidence genes.
python bed2gff.py KN9204_HC.bed
# translate all low confidence genes.
python bed2gff.py KN9204_LC.bed
```

#### Output files
- IGDB_v1.0_KN9204_HC_201811.gff3
- IGDB_v1.0_KN9204_HC_201811.gtf
- IGDB_v1.0_KN9204_HC_201811_cds.fasta
- IGDB_v1.0_KN9204_HC_201811_pep.fasta
- IGDB_v1.0_KN9204_HC_201811_transcripts.fasta
- IGDB_v1.0_KN9204_LC_201811.gff3
- IGDB_v1.0_KN9204_LC_201811.gtf
- IGDB_v1.0_KN9204_LC_201811_cds.fasta
- IGDB_v1.0_KN9204_LC_201811_pep.fasta
- IGDB_v1.0_KN9204_LC_201811_transcripts.fasta
