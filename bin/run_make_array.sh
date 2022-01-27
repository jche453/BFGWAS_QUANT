#!/usr/bin/bash

j=$1
i=$SGE_TASK_ID
wkdir=/mnt/YangFSS/data/jchen/simulation_30casual_h2_25_a23switch/wkdir_$SGE_TASK_ID
mkfile=${wkdir}/_BFGWAS.mk
cd $wkdir

echo Run make clean
make -f ${mkfile} clean


echo Run make with $mkdir and j=$j parallel jobs

make -k -C ${wkdir} -f ${mkfile} -j ${j} > ${wkdir}/make.output 2> ${wkdir}/make.err

exit
