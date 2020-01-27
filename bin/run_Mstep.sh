#!/usr/bin/bash


rv=$1
i=$2
wkDir=$3
hypcurrent=$4
Annofile=$5

# echo Curerent working directory
# pwd

echo Run Mstep with $rv, $i step, $hypcurrent, $wkDir

Rscript --vanilla /home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/bin/Mstep.r $rv $i $hypcurrent $Annofile $wkDir>> $wkDir/Rout.txt

exit
