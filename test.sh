bfGWAS_SS_dir="/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation"

######## Specify file directories and variable values for BFGWAS
study=simulation_geno # project name
N=1893  # sample size
chr=10
number=$(($chr + 1))
pheno_var=$(sed -n ${number}p /mnt/YangFSS/data/jchen/20block_simulation/simulation_info.txt | awk '{print $5}')

# genome-block names
filehead=/mnt/YangFSS/data/jchen/20block_simulation/FileHead.txt

# phenotype directory
pheno=/mnt/YangFSS/data/jchen/20block_simulation/wkdir_10/pheno.txt
# genotype directory
geno_dir=/mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/LDdetect_SegmentedVCF
Score_dir=/mnt/YangFSS/data/jchen/20block_simulation/wkdir_10/Score
LD_dir=/mnt/YangFSS/data/jchen/20block_simulation/LDstatistic
LDwindow=1
AnnoNumber=4

hfile=/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/simulation/hypval.current

annoDir=/mnt/YangFSS/data/jchen/20block_simulation/Anno

# Set up Working directory
wkdir=/mnt/YangFSS/data/jchen/20block_simulation/wkdir_10
cd $wkdir

###### Call gen_mkf.pl to generate Make file
# directory for run_Estep.sh and run_Mstep.sh
EMdir=/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/bin
# Specify computation specs
em=5 # EM steps
burnin=100000 # Burn-in iterations in MCMC
Nmcmc=100000 # MCMC iteration number
mkfile=${wkdir}/${study}_BFGWAS.mk
line=WGS_1898_samples_CHR_19_11284028_13471127
/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/bin/Estep_mcmc -inputSS \
-score ${Score_dir}/${line}.score.txt.gz \
-LDcorr ${LD_dir}/${line}.LDcorr.txt.gz \
-a ${annoDir}/Anno_${line}.gz \
-hfile ${hfile} \
-n ${N} -pv ${pheno_var} -maf 0.01 -r2 0.001 -bvsrm -smin 0 -smax 10 -win 100 \
-o ${line} -w ${burnin} -s ${Nmcmc} -initype 3 -AnnoNumber ${AnnoNumber}\
> ${wkdir}/OUT/${line}.output.txt
