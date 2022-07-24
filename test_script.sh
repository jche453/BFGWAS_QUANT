BFGWAS_dir="/home/jyang/GIT/BFGWAS_QUANT"
geno_dir=/mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/LDdetect_SegmentedVCF
wkdir="/home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir"
cd $wkdir

############################################################
############ Set input argument values
Nsample=1893 # sample size
Anum=4 # 4 annotations
em=3 # EM steps
Nburnin=10000  # Burn-in iterations in MCMC
Nmcmc=10000  # MCMC iteration number
hfile=${BFGWAS_dir}/Example/ExData/hypval_4anno.txt #  Initial prior parameter values
mkfile=${wkdir}/simu_BFGWAS.mk ## Make file directory

############# Test Simulation Data with 4 blocks
##### Set directories for filehead, GWAS summary data, Annotation Data, and Reference LD files
filehead=${BFGWAS_dir}/Example/ExData/filehead_4block.txt
Zscore_dir=${BFGWAS_dir}/Example/ExData/Zscore
LDdir=${BFGWAS_dir}/Example/ExData/RefLD
anno_dir=${BFGWAS_dir}/Example/ExData/Anno

# filehead=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/sim_20blocks_filehead.txt
cat $filehead | while read line ; do
echo $line
#rsync /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD/${line}.LDcorr.txt.gz*  /mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/RefLD/
rsync /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD/${line}.LDcorr.txt.gz* /home/jyang/GIT/BFGWAS_QUANT/Example/ExData/RefLD/
done

############# Test Simulation Data with  20 blocks
# filehead=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/sim_20blocks_filehead.txt
# Zscore_dir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Sim_Zscore
# LDdir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD
# anno_dir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Anno

########### Generate make file with all jobs
${BFGWAS_dir}/bin/gen_mkf.pl \
--wkdir ${wkdir} --tool_dir ${BFGWAS_dir} \
--filehead ${filehead} --LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--anno_dir ${anno_dir} --AnnoNumber ${Anum} --hfile ${hfile} \
--maf 0.01 --Nsample ${Nsample} --a_gamma 1.001 --b_gamma 1 \
--Nburnin ${Nburnin} --Nmcmc ${Nmcmc} --NmcmcLast ${Nmcmc} \
--em ${em} --mf ${mkfile}

########## Submit a job to run the Makefile
j=5 # Number of cores to run jobs in parallel
# make -f ${mkfile} clean
qsub -q b.q -j y -pe smp ${j} -wd ${wkdir} -N BFGWAS ${BFGWAS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}

###### Test E-step of a single block
/home/jyang/GIT/BFGWAS_QUANT/bin/Estep_mcmc -inputSS -Zscore /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Sim_Zscore/WGS_1898_samples_CHR_19_3019660_4348967.Zscore.txt.gz -LDcorr /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD/WGS_1898_samples_CHR_19_3019660_4348967.LDcorr.txt.gz -a /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Anno/Anno_WGS_1898_samples_CHR_19_3019660_4348967.txt.gz -hfile /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/hypval.current -n 1893 -maf 0.01 -bvsrm -smin 0 -smax 5 -win 100 -o WGS_1898_samples_CHR_19_3019660_4348967 -w 10 -s 10 -AnnoNumber 4 -seed 2022

###### Test M-step
Rscript --vanilla /home/jyang/GIT/BFGWAS_QUANT/bin/Mstep.r /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/hyptemp0.txt 0 1.001 1 1893 /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/hypval.current /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/paramtemp0.txt /home/jyang/GIT/BFGWAS_QUANT/Example/Test_wkdir/Eoutput/EM_result.txt


############ Generate LDcorr matrix and GWAS Zscore for simulation data (20 blocks)
pheno=${BFGWAS_dir}/Example/ExData/sim_pheno.txt
# filehead=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/sim_20blocks_filehead.txt
#Zscore_dir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/Sim_Zscore
#LDdir=/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/RefLD
# LDwindow=1000000

## Run all blocks sequencially
## >32GB memory might be needed
${BFGWAS_dir}/bin/GetRefLD.sh --wkdir ${wkdir} \
--toolE ${BFGWAS_dir}/bin/Estep_mcmc \
--geno_dir ${geno_dir} --filehead ${filehead} --pheno ${pheno} \
--GTfield GT --genofile_type vcf --maf 0.001 --LDwindow ${LDwindow}  \
--Zscore_dir ${Zscore_dir} --LDdir ${LDdir}

qsub -q b.q -j y -pe smp 4 -wd ${wkdir} -N GetRefLD -t 1-20 -tc 14 /home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Scripts/ShellScripts/run_GetRefLD_array.sh