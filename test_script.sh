BFGWAS_dir="/home/jyang/GIT/BFGWAS_QUANT"
LDdir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD
geno_dir=/mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/LDdetect_SegmentedVCF
wkdir="/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir"
cd $wkdir

############################################################
############ Generate LDcorr matrix and GWAS Zscore for simulation data
pheno=/mnt/YangFSS/data/jchen/BFGWAS_simulation/01-SimulationOld/simulation_10casual_h2_50/wkdir_1/pheno.txt
filehead=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/sim_20blocks_filehead.txt
Zscore_dir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Sim_Zscore
# LDwindow=1000000

## Run all blocks sequencially
${BFGWAS_dir}/bin/GetRefLD.sh --wkdir ${wkdir} \
--toolE ${BFGWAS_dir}/bin/Estep_mcmc \
--geno_dir ${geno_dir} --filehead ${filehead} --pheno ${pheno} \
--GTfield GT --genofile_type vcf --maf 0.001 --LDwindow ${LDwindow}  \
--Zscore_dir ${Zscore_dir} --LDdir ${LDdir}

# line=WGS_1898_samples_CHR_19_11284028_13471127
/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Scripts/ShellScripts/GetRefLD_per_block.sh ${line}

############################################################
############# Test Junyu Simulation Data
# simudir="/mnt/YangFSS/data/jchen/BFGWAS_simulation/01-SimulationOld/simulation_10casual_h2_50"
filehead=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/sim_20blocks_filehead.txt
anno_dir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Anno
hfile=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/hypval_4anno.txt
Nsample=1893 # sample size
Anum=4 # 4 annotations
maf=0.01
em=2 # EM steps
Nburnin=10000  # Burn-in iterations in MCMC
Nmcmc=10000  # MCMC iteration number
mkfile=${wkdir}/simu_BFGWAS.mk

${BFGWAS_dir}/bin/gen_mkf.pl \
--wkdir ${wkdir} --tool_dir ${BFGWAS_dir} \
--filehead ${filehead} --LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--anno_dir ${anno_dir} --AnnoNumber ${Anum} --hfile ${hfile} \
--maf ${maf} --Nsample ${Nsample} \
--Nburnin ${Nburnin} --Nmcmc ${Nmcmc} --NmcmcLast ${Nmcmc} \
--em ${em} --mf ${mkfile}

j=20 # Number of cores to run the job in parallel
# make -f ${mkfile} clean
qsub -q b.q -j y -pe smp ${j} -wd ${wkdir} -N BFGWAS ${BFGWAS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}

###### Test single block and Mstep
/home/jyang/GIT/BFGWAS_QUANT/bin/Estep_mcmc -inputSS -Zscore /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Sim_Zscore/WGS_1898_samples_CHR_19_23467746_28557893.Zscore.txt.gz -LDcorr /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD/WGS_1898_samples_CHR_19_23467746_28557893.LDcorr.txt.gz -a /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Anno/Anno_WGS_1898_samples_CHR_19_23467746_28557893.txt.gz -hfile /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir/hypval.current -maf 0.01 -n 1893 -bvsrm -smin 0 -smax 5 -win 100 -o WGS_1898_samples_CHR_19_23467746_28557893 -w 10 -s 10 -AnnoNumber 4 -seed 2022

Rscript --vanilla /home/jyang/GIT/BFGWAS_QUANT/bin/Mstep.r /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir/Eoutput/hyptemp0.txt 0 1.0001 1 1893 /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir/hypval.current /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir/Eoutput/paramtemp0.txt /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir/Eoutput/EM_result.txt


