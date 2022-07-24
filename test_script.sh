BFGWAS_dir="/home/jyang/GIT/BFGWAS_QUANT"
geno_dir=/mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/LDdetect_SegmentedVCF
wkdir="/home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir"
cd $wkdir

############################################################
############ Generate LDcorr matrix and GWAS Zscore for simulation data
pheno=/home/jyang/GIT/BFGWAS_QUANT/Example/ExData/sim_pheno.txt
filehead=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/sim_20blocks_filehead.txt
# filehead=/home/jyang/GIT/BFGWAS_QUANT/Example/ExData/filehead_4block.txt

Zscore_dir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Sim_Zscore
anno_dir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Anno
LDdir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD

cat $filehead | while read line ; do
	rsync ${LDdir}/${line}.LDcorr.txt.gz* /home/jyang/GIT/BFGWAS_QUANT/Example/ExData/RefLD/
	rsync ${Zscore_dir}/${line}.Zscore.txt.gz* /home/jyang/GIT/BFGWAS_QUANT/Example/ExData/Zscore/
	rsync ${anno_dir}/Anno_${line}.txt.gz /home/jyang/GIT/BFGWAS_QUANT/Example/ExData/Anno/
done

# LDwindow=1000000
## Run all blocks sequencially
${BFGWAS_dir}/bin/GetRefLD.sh --wkdir ${wkdir} \
--toolE ${BFGWAS_dir}/bin/Estep_mcmc \
--geno_dir ${geno_dir} --filehead ${filehead} --pheno ${pheno} \
--GTfield GT --genofile_type vcf --maf 0.001 --LDwindow ${LDwindow}  \
--Zscore_dir ${Zscore_dir} --LDdir ${LDdir}

qsub -q b.q -j y -pe smp 5 -wd ${wkdir} -N GetRefLD -t 1-20 -tc 14 /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Scripts/ShellScripts/run_GetRefLD_array.sh

# line=WGS_1898_samples_CHR_19_11284028_13471127
/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Scripts/ShellScripts/GetRefLD_per_block.sh ${line}

############################################################
############# Test Junyu Simulation Data
# simudir="/mnt/YangFSS/data/jchen/BFGWAS_simulation/01-SimulationOld/simulation_10casual_h2_50"
# filehead=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/sim_20blocks_filehead.txt

hfile=/home/jyang/GIT/BFGWAS_QUANT/Example/ExData/hypval_4anno.txt
Nsample=1893 # sample size
Anum=4 # 4 annotations
em=3 # EM steps
Nburnin=10000  # Burn-in iterations in MCMC
Nmcmc=10000  # MCMC iteration number
mkfile=${wkdir}/simu_BFGWAS.mk

${BFGWAS_dir}/bin/gen_mkf.pl \
--wkdir ${wkdir} --tool_dir ${BFGWAS_dir} \
--filehead ${filehead} --LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--anno_dir ${anno_dir} --AnnoNumber ${Anum} --hfile ${hfile} \
--maf 0.01 --Nsample ${Nsample} --a_gamma 1.001 --b_gamma 1 \
--Nburnin ${Nburnin} --Nmcmc ${Nmcmc} --NmcmcLast ${Nmcmc} \
--em ${em} --mf ${mkfile}

j=5 # Number of cores to run the job in parallel
# make -f ${mkfile} clean
qsub -q b.q -j y -pe smp ${j} -wd ${wkdir} -N BFGWAS ${BFGWAS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}

###### Test single block and Mstep
/home/jyang/GIT/BFGWAS_QUANT/bin/Estep_mcmc -inputSS -Zscore /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Sim_Zscore/WGS_1898_samples_CHR_19_3019660_4348967.Zscore.txt.gz -LDcorr /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD/WGS_1898_samples_CHR_19_3019660_4348967.LDcorr.txt.gz -a /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Anno/Anno_WGS_1898_samples_CHR_19_3019660_4348967.txt.gz -hfile /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/hypval.current -n 1893 -maf 0.01 -bvsrm -smin 0 -smax 5 -win 100 -o WGS_1898_samples_CHR_19_3019660_4348967 -w 10 -s 10 -AnnoNumber 4 -seed 2022

Rscript --vanilla /home/jyang/GIT/BFGWAS_QUANT/bin/Mstep.r /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/hyptemp0.txt 0 1.001 1 1893 /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/hypval.current /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/paramtemp0.txt /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/EM_result.txt


