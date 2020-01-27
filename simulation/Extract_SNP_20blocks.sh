#!/usr/bin/bash

## Segment VCF files into 2.5MB blocks
bed_dir=/mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/LDdetect_SegmentedVCF_CADDsubset
bed_file=${bed_dir}/FileHead_1703_anno_chr19.bed
data_dir=/mnt/YangFSS/data/WingoLab-Data/WGS_JointCall/LDdetect_SegmentedVCF
out_dir=/mnt/YangFSS/data/jchen/20block_simulation/original_genotype
cd ${out_dir}

#module load vcftools/0.1.13
#vcf-query -l ${data_dir}/WGS_1898_samples_CHR_10_100241302_100668400.vcf.gz > ${out_dir}/sampleID.txt

cat ${bed_file} | while read line ; do

	seg_name=$(echo ${line} | awk '{print $4}')

	echo $seg_name

	/home/jyang/GIT/bfGWAS_SS/bin/Estep_mcmc -vcf ${data_dir}/${seg_name}.vcf.gz -p ${out_dir}/sampleID.txt -o ${seg_name} -GTfield GT -saveGeno -maf 0.01
	sed -i 's/#//g' ./output/${seg_name}.geno
done
