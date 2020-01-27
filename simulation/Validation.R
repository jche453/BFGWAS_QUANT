library(statip)
library(wtest)

setwd("/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/test_output")
load("/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/simulation/parameter.rda")
geno = read.table("/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/simulation/output/simulation.geno",stringsAsFactors = F, header = T)
geno_table = as.data.frame(t(geno[,-c(1:5)]))
geno_table = round(geno_table,0)
geno$maf = apply(geno_table,2,maf)
geno = geno [geno$maf >= 0.01,]
snp_n = nrow(geno)
para = data.frame(geno$ID, beta, gamma, pi, sigma_2)

Annodata= read.table("/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/test_output/Eoutput/AnnoScore5.txt")
Anum = ncol(Annodata) - 3
temp_col_names <- c("SNP","gamma", "beta")
for(i in 1:Anum){
  temp_col_names <- c(temp_col_names, 
                      paste("Annotation", i, sep = "_"))
}
colnames(Annodata) <- temp_col_names

prehyp <- read.table("/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/test_output/hypval.current", header=F)
a = prehyp[,1]
b = prehyp[,2]

log_pi = a[1]
for (i in 1:(Anum)) {
  log_pi = log_pi + a[i + 1] * Annodata[,i+3]
}
Annodata$pi = exp(log_pi)/(1 + exp(log_pi))

log_sigma_2 = b[1]
for (i in 1:(Anum)) {
  log_sigma_2 = log_sigma_2 + b[i + 1] * Annodata[,i+3]
}
Annodata$sigma_2 = exp(log_sigma_2)

combine = merge(para, Annodata[c(1:3,7,8)], by.x = "geno.ID", by.y = "SNP")
cor(combine$beta.x,combine$beta.y) ^2
#0.5308783
cor(combine$gamma.x,combine$gamma.y) ^2
#0.762831

# of casual snps 
nrow(combine[combine$gamma.x > 0 & combine$gamma.y > 0.1 , ])
#5
nrow(combine[combine$gamma.y > 0.1, ])


sum(combine$gamma.x)
#5
sum(combine$gamma.y)
#9.9762

#Pheno R2
#x: true
pheno_x = read.table("/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/simulation/pheno.txt")
#y: model
pheno_geno = merge(geno, combine,by.x = "ID", by.y = "geno.ID")
pheno_y = t(pheno_geno[,-c(1:5,1108: length(pheno_geno))]) %*% pheno_geno$beta.y
cor(pheno_x[,2], pheno_y) ^ 2
#0.1844257

#Geno R2
geno_x = pheno_geno[pheno_geno$gamma.x == 1,]
geno_y = pheno_geno[pheno_geno$gamma.y > 0.1,]
LD = cor(t(geno_x[,-c(1:5,1108: length(pheno_geno))]), t(geno_y[,-c(1:5,1108: length(pheno_geno))])) ^ 2

param = read.table("/mnt/icebreaker/data2/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/test_output/output/chr19_maf05.paramtemp", header = T)