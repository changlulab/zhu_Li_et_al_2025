#!/bin/bash

echo "choose the group you want to process and press [ENTER] :"
ls /projects/lu_lab/Gaoshan/MIA/ChIP_seq/
read GROUP
cd /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP 2> /dev/null

##############################################################################
until [ "$PWD" == "/projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP" ]; do
#the space in the bracket (at the begining and at the end) is very important##

echo "Please choose a real group"
ls /projects/lu_lab/Gaoshan/MIA/ChIP_seq/
read GROUP
cd /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP 2> /dev/null;
done

echo "choose the hisotne you want to process and press [ENTER] :"
ls /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP
read HISTONE
cd /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE 2> /dev/null
until [ "$PWD" == "/projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE" ]; do
echo "Please choose a real histone mark"
ls /projects/lu_lab/Gaoshan/MIA/ChIP_seq/
read HISTONE
cd /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE 2> /dev/null;
done

echo "choose the subgroup you want to process and press [ENTER] :"
ls /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE
read SUBGROUP
cd /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP 2> /dev/null
until [ "$PWD" == "/projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP" ]; do
echo "Please choose a real subgroup"
ls /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE
read SUBGROUP
cd /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP 2> /dev/null;
done

cd /projects/lu_lab/Gaoshan/MIA/ChIP_seq/
sed "s/GROUPGOESHERE/"$GROUP"/" ./base_precorr_enhancer.sh > /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr_1.sh

sed "s/SUBgroupGOESHERE/"$SUBGROUP"/" /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr_1.sh > /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr_2.sh

sed "s/HISTONEGOESHERE/"$HISTONE"/" /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr_2.sh > /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr.sh

cd /projects/lu_lab/Gaoshan/MIA/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP

rm ./precorr_1.sh ./precorr_2.sh
chmod u+x ./precorr.sh

sbatch -p normal_q --account=chipseq -t 3-00:00:00 --mem=20GB ./precorr.sh
