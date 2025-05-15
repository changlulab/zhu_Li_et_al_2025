#!/bin/bash

INDEX=INDEXGOESHERE
Group=GROUPGOESHERE
Subgroup=SUBgroupGOESHERE

cd /projects/lu_lab/Gaoshan/MIA/RNA_seq/$Group/$Subgroup

mkdir fastqc
mkdir trimmed

cd Raw_Data/

gunzip *.gz
ls *.fastq | while read id; do (trim_galore -o /projects/lu_lab/Gaoshan/MIA/RNA_seq/$Group/$Subgroup/trimmed $id);done

cd ..
cd trimmed

fastqc -o /projects/lu_lab/Gaoshan/MIA/RNA_seq/$Group/$Subgroup/fastqc *.fq
cd ..

FILES=$PWD/trimmed/*.fq
SAM=.sam
LOG=.log
BAM=.bam
TXT=.txt
GTF=.ensGene.gtf
FLAG=.flagstat

mkdir Aligned_SAM
mkdir Aligned_results
mkdir Aligned_BAM
mkdir Count

for fn in $FILES
do
sample=$(basename $fn)
f=${sample%%_*}
echo $f

hisat2 -p 16 -x /home/gaoshanli/Data/mm10/Sequence/hisat2Index/$INDEX/genome -U $PWD/trimmed/$sample -S $PWD/Aligned_SAM/$f$SAM 2> $PWD/Aligned_results/$f$LOG

samtools flagstat $PWD/Aligned_SAM/$f$SAM > $PWD/Aligned_results/$f$FLAG

samtools view -bS $PWD/Aligned_SAM/$f$SAM | samtools sort -o $PWD/Aligned_BAM/$f$BAM

samtools index $PWD/Aligned_BAM/$f$BAM

featureCounts -t exon -g gene_name -a /home/gaoshanli/Data/mm10/gtf/$INDEX/$INDEX$GTF -o $PWD/Count/$f$TXT $PWD/Aligned_BAM/$f$BAM 2> $PWD/Count/$f$LOG

done



