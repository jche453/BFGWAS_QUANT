#!/usr/bin/bash

wkdir=$1
mkfile=$2
j=$3

cd $wkdir

echo Run make clean
make -f ${mkfile} clean


echo Run make with $mkdir and j=$j parallel jobs

make -k -C ${wkdir} -f ${mkfile} -j ${j} > ${wkdir}/make.output 2> ${wkdir}/make.err

exit


