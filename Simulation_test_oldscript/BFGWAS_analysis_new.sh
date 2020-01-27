######## Please use our data disk (/mnt/YangFSS/data) for input and output of all BFGWAS jobs
######## Example scripts are under /home/jyang/ResearchProjects/ROSMAP_GWAS/BFGWAS/Scripts
######### Tool directory
bfGWAS_SS_dir="/home/jyang/GIT/bfGWAS_SS"

######## Specify file directories and variable values for BFGWAS
chr=2
study=simulation_test # project name
N=1893  # sample size
number=$(($chr + 1))
pheno_var=$(sed -n ${number}p /mnt/YangFSS/data/jchen/simulation_20blocks/simulation_info.txt | awk '{print $5}')
# genome-block names
filehead=/mnt/YangFSS/data/jchen/simulation_20blocks/FileHead_1703_anno_chr19.bed
# FileHead=/home/jyang/Collaborations/IrwinSAGE/SummaryData/FileHead_AA_2581.txt

# phenotype directory
pheno=/mnt/YangFSS/data/jchen/simulation_20blocks/wkdir_${chr}/pheno.txt
# genotype directory
geno_dir=/mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/LDdetect_SegmentedVCF

######### Please first generate summary statistics (LD and score statistics files) files for best computation speed
# LD coefficient directory
LD_dir=/mnt/YangFSS/data/jchen/simulation_20blocks/LD_statistic
# Score statistics directory
Score_dir=/mnt/YangFSS/data/jchen/simulation_20blocks/wkdir_${chr}/Score
LDwindow=1

# Submit Array jobs with GetRefLD.sh
#qsub -j y -wd ${Score_dir} -N GetRefLD_${study} -t 1-1703 -tc 100 -pe smp 4 /mnt/icebreaker/data2/home/jchen/bfGWAS/GetRefLD.sh ${geno_dir} ${pheno} ${filehead} ${LD_dir} ${Score_dir} ${LDwindow}

####### Command to run BFGWAS for genome-wide blocks
##  run_BFGWAS.sh will call the perl script gen_mkf.pl to generate a makefile with all BFGWAS_SS jobs run by run_make.sh
## The executible file /home/jyang/GIT/bfGWAS_SS/bin/Estep_mcmc was run by run_Estep.sh, which is called by make
## Then clean all previous outputs and re-run the generated Makefile
## Please test with ~10 genome-blocks instead of running the whole genome-wide data

# Hyper parameter file
hfile=/mnt/YangFSS/data/jchen/simulation_20blocks/old_hypval.current

# Annotation file directory and Annotation code for categories (optional)
annoDir=/mnt/YangFSS/data/jchen/simulation_20blocks/simulation_test_oldscript/Anno
annoCode=/mnt/YangFSS/data/jchen/simulation_20blocks/simulation_test_oldscript/AnnocodeROSMAP.txt

# Set up Working directory
wkdir=/mnt/YangFSS/data/jchen/simulation_20blocks/simulation_test_oldscript
mkdir -p $wkdir
cd $wkdir

###### Call gen_mkf.pl to generate Make file
# directory for run_Estep.sh and run_Mstep.sh
EMdir=/home/jyang/ResearchProjects/ROSMAP_GWAS/BFGWAS/Scripts
# Specify computation specs
em=5 # EM steps
burnin=10000 # Burn-in iterations in MCMC
Nmcmc=10000 # MCMC iteration number
mkfile=${wkdir}/${study}_BFGWAS.mk

/home/jyang/ResearchProjects/ROSMAP_GWAS/BFGWAS/Scripts/gen_mkf.pl \
--EMdir ${EMdir} --hyp ${hfile} \
-n ${N} --pv ${pheno_var} \
-w ${wkdir} --geno sumstat \
-f ${filehead} -l local \
--ad ${annoDir} --ac ${annoCode} \
--LDdir ${LD_dir} --Scoredir ${Score_dir} \
-j BFGWAS_${study} --em ${em} -b ${burnin} -N ${Nmcmc} \
--mf ${mkfile}


######### Submit the job for running the makefile (run_Estep.sh and run_Mstep.sh were called by make)
j=25 # Number of cores to run the job in parallel
qsub -q b.q -j y -pe smp ${j} -wd ${wkdir} -N BFGWAS_${study} /home/jyang/ResearchProjects/ROSMAP_GWAS/BFGWAS/Scripts/run_make.sh ${wkdir} ${mkfile} ${j}
