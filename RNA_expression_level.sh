#!/bin/bash

dir=/projects/lu_lab/Gaoshan/MIA/Taiji/RNA_data
dir1=/projects/lu_lab/Gaoshan/MIA/Taiji/rna_data
dir2=/projects/lu_lab/Gaoshan/MIA/Taiji/rna_data_ready

FILES=${dir}/*.txt

TXT=.txt


for fn in $FILES
do
echo `basename "$fn"`
f=`basename "${fn%.*}"`
echo $f

sed -i '1d' ${dir}/$f$TXT

sed -i '1d' ${dir}/$f$TXT

awk '{ print $1, $7 }' ${dir}/$f$TXT >${dir1}/$f$TXT

tr ' ' \\t <${dir1}/$f$TXT > ${dir2}/$f$TXT

done
