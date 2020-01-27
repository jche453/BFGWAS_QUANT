qlogin
module load tabix/0.2.6

#Take 5000 SNPs to initial simulation
cd /home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/simulation

zcat /mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_19.recalibrated_variants.vcf.gz | head -n 300 | grep CHROM > chr19.vcf
zcat /mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_19.recalibrated_variants.vcf.gz | sed -n '301,5301p'  >> chr19.vcf

/home/jyang/GIT/bfGWAS_SS/bin/Estep_mcmc -vcf chr19.vcf \
          -p /home/jchen/bfGWAS/Phenotype/cogdx_wgs.txt \
          -o simulation -GTfield GT -saveGeno -maf 0

sed -i 's/#//g' ./output/simulation.geno

#run simulation_Jingjing.R

#Extract SNP with MAF larger 0.05
cd /home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/simulation
zcat /mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_19.recalibrated_variants.vcf.gz | head -n 300 | grep CHROM > chr19_maf05.vcf
cat simulation_geno.txt| while read line ; do
	chr=$(echo $line | awk -F[:_/] '{print $1}' )
	start=$(echo $line | awk -F[:_/] '{print $2}' )
	ref=$(echo $line | awk -F[:_/] '{print $3}' )
	alt=$(echo $line | awk -F[:_/] '{print $4}' )
	echo ${chr}:${start}-${start}
	tabix /mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_${chr}.recalibrated_variants.vcf.gz ${chr}:${start}-${start} | awk -v ref=${ref} '$4==ref{print }' >> chr19_maf05.vcf
done


/home/jyang/GIT/bfGWAS_SS/bin/Estep_mcmc -vcf chr19_maf05.vcf \
          -p /home/jchen/bfGWAS/Phenotype/cogdx_wgs.txt \
          -o simulation_maf05 -GTfield GT -saveGeno -maf 0

sed -i 's/#//g' ./output/simulation_maf05.geno

#gzip the annotation file
gzip Anno_chr19_maf05.txt
gzip chr19_maf05.vcf
