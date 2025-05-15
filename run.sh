#!/bin/bash

echo "choose the Index you want to use and press [ENTER] :"
ls /home/gaoshanli/Data/mm10/gtf/
read INDEX
cd /home/gaoshanli/Data/mm10/gtf/$INDEX 2> /dev/null
until [ "$PWD" == "/home/gaoshanli/Data/mm10/gtf/$INDEX" ]; do
echo "Please choose a real index"
ls /home/gaoshanli/Data/mm10/gtf/
read INDEX
cd /home/gaoshanli/Data/mm10/gtf/$INDEX 2> /dev/null;
done

echo "choose the group you want to process and press [ENTER] :"
ls /projects/lu_lab/Gaoshan/MIA/RNA_seq/
read GROUP
cd /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP 2> /dev/null
until ["$PWD" == "/projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP"]; do
echo "Please choose a real group"
ls /projects/lu_lab/Gaoshan/MIA/RNA_seq/
read GROUP
cd /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP 2> /dev/null;

echo "choose the subgroup you want to process and press [ENTER] :"
ls /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP
read SUBGROUP
cd /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP/$SUBGROUP 2> /dev/null
until ["$PWD" == "/projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP/$SUBGROUP"]; do
echo "Please choose a real subgroup"
ls /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP
read SUBGROUP
cd /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP/$SUBGROUP 2> /dev/null;

sed "s/GROUPGOESHERE/"$GROUP"/" ./base_RNAseqanalysis.sh > /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP/$SUBGROUP/RNAseqanalysis_1.sh

sed "s/SUBgroupGOESHERE/"$SUBGROUP"/" /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP/$SUBGROUP/RNAseqanalysis_1.sh > /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP/$SUBGROUP/RNAseqanalysis_2.sh

sed "s/INDEXGOESHERE/"$INDEX"/" /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP/$SUBGROUP/RNAseqanalysis_2.sh > /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP/$SUBGROUP/RNAseqanalysis.sh

cd /projects/lu_lab/Gaoshan/MIA/RNA_seq/$GROUP/$SUBGROUP

rm ./RNAseqanalysis_1.sh
rm ./RNAseqanalysis_2.sh
chmod u+x ./RNAseqanalysis.sh

sbatch -A chipseq --mem-per-cpu=20G -t 2-00:00:00 -p normal_q ./RNAseqanalysis.sh
