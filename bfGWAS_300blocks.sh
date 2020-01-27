######## Please use our data disk (/mnt/YangFSS/data) for input and output of all BFGWAS jobs
######## Example scripts are under /home/jyang/ResearchProjects/ROSMAP_GWAS/BFGWAS/Scripts
######### Tool directory

for chr in {2..300};
do
  bfGWAS_SS_dir="/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation"
  ######## Specify file directories and variable values for BFGWAS
  N=1893  # sample size
  number=$(($chr + 1))
  pheno_var=$(sed -n ${number}p /mnt/YangFSS/data/jchen/simulation_20blocks/simulation_info.txt | awk '{print $5}')
  # genome-block names
  filehead=/mnt/YangFSS/data/jchen/simulation_20blocks/FileHead_1703_anno_chr19.bed

  # phenotype directory
  pheno=/mnt/YangFSS/data/jchen/simulation_20blocks/wkdir_${chr}/pheno.txt
  # genotype directory
  geno_dir=/mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/LDdetect_SegmentedVCF

  ######### Please first generate summary statistics (LD and score statistics files) files for best computation speed
  # LD coefficient directory
  LD_dir=/mnt/YangFSS/data/jchen/simulation_20blocks/LD_statistic
  # Score statistics directory
  Score_dir=/mnt/YangFSS/data/jchen/simulation_20blocks/wkdir_${chr}/Score
  if [[ ! -e $Score_dir ]]; then
    mkdir $Score_dir
  elif [[ ! -d $Score_dir ]]; then
    echo "$Score_dir already exists" 1>&2
  fi
  LDwindow=1

  # Submit Array jobs with GetRefLD.sh
  qsub -j y -wd ${Score_dir} -N Score_${chr} -t 1-20 -tc 100 -pe smp 5 /mnt/icebreaker/data2/home/jchen/bfGWAS/GetRefLD.sh ${geno_dir} ${pheno} ${filehead} ${LD_dir} ${Score_dir} ${LDwindow}
done


####### Command to run BFGWAS for genome-wide blocks
##  run_BFGWAS.sh will call the perl script gen_mkf.pl to generate a makefile with all BFGWAS_SS jobs run by run_make.sh
## The executible file /home/jyang/GIT/bfGWAS_SS/bin/Estep_mcmc was run by run_Estep.sh, which is called by make
## Then clean all previous outputs and re-run the generated Makefile
## Please test with ~10 genome-blocks instead of running the whole genome-wide data

for chr in {2..2};
do
  bfGWAS_SS_dir="/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation"
  ######## Specify file directories and variable values for BFGWAS
  N=1893  # sample size
  number=$(($chr + 1))
  pheno_var=$(sed -n ${number}p /mnt/YangFSS/data/jchen/simulation_20blocks/simulation_info.txt | awk '{print $5}')
  # genome-block names
  filehead=/mnt/YangFSS/data/jchen/simulation_20blocks/FileHead_1703_anno_chr19.bed
  #filehead=/mnt/YangFSS/data/jchen/simulation_20blocks/FileHead_1703_anno_chr19_new.bed

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
  Anum=4
  hfile=/mnt/YangFSS/data/jchen/simulation_20blocks/hypval.current

  annoDir=/mnt/YangFSS/data/jchen/simulation_20blocks/Anno

  # Set up Working directory
  #wkdir=/mnt/YangFSS/data/jchen/simulation_20blocks/wkdir_4/Theta_10-6
  wkdir=/mnt/YangFSS/data/jchen/simulation_20blocks/wkdir_${chr}
  cd $wkdir

  ###### Call gen_mkf.pl to generate Make file
  # directory for run_Estep.sh and run_Mstep.sh
  EMdir=/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/bin
  # Specify computation specs
  em=5 # EM steps
  burnin=10000  # Burn-in iterations in MCMC
  Nmcmc=10000  # MCMC iteration number
  mkfile=${wkdir}/${study}_BFGWAS.mk

  /home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/bin/gen_mkf.pl \
  --EMdir ${EMdir} --hyp ${hfile} \
  -n ${N} --pv ${pheno_var} \
  -w ${wkdir} --geno sumstat \
  -f ${filehead} -l local \
  --ad ${annoDir} \
  --LDdir ${LD_dir} --Scoredir ${Score_dir} \
  -j BFGWAS_${chr} --em ${em} -b ${burnin} -N ${Nmcmc} \
  --mf ${mkfile} --AnnoNumber ${Anum}

  j=20 # Number of cores to run the job in parallel
  qsub -q b.q -j y -pe smp ${j} -wd ${wkdir} -N BFGWAS_${chr} /home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/bin/run_make.sh ${wkdir} ${mkfile} ${j}

done
