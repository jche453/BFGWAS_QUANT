Bayesian Functional Genome-wide Association Study (BFGWAS)
=======
**BFGWAS_QUANT** is developed based on our previous tool [**BFGWAS**](https://github.com/yjingj/BFGWAS_SS) tool as proposed  by [Yang J. et. al, AJHG, 2017](https://doi.org/10.1016/j.ajhg.2017.08.002). 

- With individual-level GWAS data, this **BFGWAS_QUANT** version first compute single variant Z-score statistics (GWAS summary statistics) and LD correlation files from individual-level GWAS data, and then load these two summary statistics files to run MCMC.
- With summary-level GWAS data, this **BFGWAS_QUANT** version  can read Z-score statistics file and LD correlation files generated from a reference panel (*need to be generated by BFGWAS_QUANT*) to run MCMC.
- **BFGWAS_QUANT** runs MCMC in parallel for multiple genome-blocks, which is handled through a `Makefile` generated by a PERL script `./bin/gen_mkf.pl`. EM-MCMC algorithm is handled based on dependency jobs in the generated `Makefile`. Users will only need to submit the generated `Makefile` as a job for the analysis results.
- **BFGWAS_QUANT** is based on a multivariable Bayesian variable selection model, which accounts for multivariate functional annotations and LD to fine-map and prioritize GWAS hits.

---

- [Software Installation](#software-installation)
- [Input Files](#input-files)
	- [Individual GWAS Input Files](#input-file-1)
		- [1. Genotype VCF Files](#1-genotype-vcf-files)
		- [2. Phenotype File](#2-phenotype-file)
	- [Summary GWAS Input Files](#input-file-2)
		- [1. GWAS Zscore File](#1-gwas-zscore-file)
		- [2. Reference LD File](#2-reference-LD-file)
	- [Other Input Files](#input-files)
		- [1. Genome Block Prefix File](#1-prefix-file)
		- [2. Prior Parameter File](#2-prior-file)
		- [3. Annotation Code File](#3-anno-code-file)
		- [4. Annotation File](#4-anno-file)
- [Example Usage](#example-usage)
	- [Step 1. Obtain Summary Statistics](#step-1-obtain-summary-statistics)
	- [Step 2. Generate Make File](#step-2-create-makefile)
	- [Step 3. Run Make File](#step-3-run-makefile)
- [Output Files](#output-files)
- [Analyze Results](#analyze-results)
- [Flags](#flags)

---


## Software Installation

### 1. Compile **Estep_mcmc**
* Install required C++ libraries C++ libraries **zlib**, **gsl**, **eigen3**, **lapack**, **atlas**, **blas** are used to develop this tool. Please install these libraries to your system and include the library path `-I[path to libraries]` accordingly in the C++ compilation command line in the `Makefile`.

* Compile C++ library `./libStatGen/libStatGen.a` under your system by using the following commands:

```
cd BFGWAS_SS/libStatGen/;
make clean;
make
```

* Compile C++ source code for the executible file *./bin/Estep_mcmc* that will be used to run the Estep MCMC algorithm, by using the following commands under `BFGWAS_SS/` directory:

```
cd BFGWAS_SS/;
make clean ;
make
```

* Even though a compiled executible file *./bin/Estep_mcmc* from our cluster is provided on GITHUB, please still compile one for your own system. The `BFGWAS_SS/Makefile` might need to be adapted for one's own system.

### 2. Additional Requirements
* Tool [**TABIX**](https://www.htslib.org/doc/tabix.html)
* R library [**data.table**](https://github.com/Rdatatable/data.table/wiki/Installation), [**tidyverse**](https://www.tidyverse.org/), [**optimx**](https://cran.r-project.org/web/packages/optimx/index.html)

## Input Files
**Example data files are provided under `./Example/ExData/`.**

### Individual GWAS Input Files
#### 1. Genotype VCF Files

- **[VCF Genotype files](http://samtools.github.io/hts-specs/VCFv4.1.pdf)** are required for using individual-level GWAS data and for generating LD correlation files. The VCF genotype files should be one per genome block (variants of the same chromosome should be in the same block), sorted by position, and then zipped by `bgzip` (file names are of `[filehead].vcf.gz`). 
- All VCF genotype files should be put under the same parent directory.
- Genotype files are supposed to be segmented based on LD. Genome blocks are expected to be approximately independent with ~5K-10K SNPs per block.
- Example segmentation information derived by [LDetect](https://bitbucket.org/nygcresearch/ldetect/src/master/) are provided in `./Example/ExData/Segmentations/lddetect_1KG_*_hg19.bed ` for EUR, AFR, and ASN populations.

#### 2. Phenotype File

### Summary GWAS Input Files
#### 1. GWAS Zscore File
#### 2. Reference LD File

### Other Input Files 
#### 1. Genome Block Prefix File

* **Genome Block Prefix File** of VCF genotype files as in `./Example/ExData/filehead_4block.txt` is required. Each row of the list file is the file head of the VCF file of one genome block as in `[filehead].vcf.gz`. Note that the VCF file extension suffix `.vcf.gz` should not be included.

#### 2. Prior Parameter File


#### 3. Annotation File


## Example Usage (`./test_script.sh`)

#### Step 1. Obtain Summary Statistics

#### Step 2. Generate Make File

#### Step 3. Run Make File


## Output Files


## Analyze Results 
- `/Example/AnalyzeResults/Analysis.r`

## Flags
