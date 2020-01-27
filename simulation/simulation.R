setwd("/home/jchen/bfGWAS/bfGWAS_QuantitativeAnnotation/Test/simulation")
library(statip)
library(wtest)

#Set parameters for simulation
geno = read.table("output/simulation.geno",stringsAsFactors = F, header = T)

Anum = 3
a = c(-12, 4, 0.8, 0)
b = c(0.1, 0.11, 3, 0)
h2 = 0.25

##### Check MAF of genotypes
geno_table = as.data.frame(t(geno[,-c(1:5)]))
geno_table = round(geno_table,0)
geno$maf = apply(geno_table,2,maf)
geno = geno [geno$maf >= 0.01,]
snp_n = nrow(geno)
#write.table(geno$ID, "simulation_geno.txt", sep = "\t", quote = F, row.names = F , col.names = F)
#Run Extract SNPs with maf larger than 0.01
#geno = read.table("output/simulation_maf05.geno",stringsAsFactors = F, header = T)

#Annotation file
A = geno[,1:5]
set.seed(2019)
annotation = round(rnorm(snp_n, mean = 0, sd = 1),4)
annotation_matrix = annotation
for (i in 1:(Anum-1)) {
  annotation_new = round(rnorm(snp_n, mean = 0, sd = 1),4)
  annotation = paste(annotation,annotation_new,sep = ",")
  annotation_matrix =cbind(annotation_matrix,annotation_new)
}
A = data.frame(A,annotation)
write.table(A,"Annotation.txt",col.names = F,quote = F,row.names = F, sep = "\t")

##### Set target proportion of SNPs, average of pi
log_pi = a[1]
for (i in 1:(Anum)) {
  log_pi = log_pi + a[i + 1] * annotation_matrix[,i]
}
pi = exp(log_pi)/(1 + exp(log_pi))
png("pi.png")
plot(pi)
dev.off()
mean(pi)

log_sigma_2 = b[1]
for (i in 1:(Anum)) {
  log_sigma_2 = log_sigma_2 + b[i + 1] * annotation_matrix[,i]
}
sigma_2 = exp(log_sigma_2)
mean(sigma_2)
png("sigma_2.png")
hist(sigma_2)
dev.off()

generate_gamma = function(x) {rbern(1, x)}
gamma = sapply(pi,generate_gamma)
mean(gamma)
length(gamma[gamma>0])
#Number of casual snp is 5

### Generate beta 
generate_beta = function(x) {rnorm(1,0,sqrt(x))}
beta = sapply(sigma_2, generate_beta) * gamma

pheno =  t(geno[,-c(1:5,ncol(geno))]) %*% beta  # without error term
var1 = var(pheno)
var_error = var1 * (1/h2 - 1)
var_error / (var1 + var_error)
sd = sqrt(var_error)
error = rnorm(1102, mean = 0, sd = sd)

#Save(error, beta, pi, sigma2)
save(error, beta, gamma, sigma_2, pi, file = "parameter.rda")

pheno =  pheno  + error
var1 / var(pheno)
var(pheno)
#3.081258
png("pheno.png")
hist(pheno)
dev.off()
pheno_1 = data.frame(rownames(pheno), pheno[,1])
pheno_1[,1] = gsub("[/.]", "-", pheno_1[,1])
write.table(pheno_1, "pheno.txt", col.names = F, row.names = F, quote = F)

#Annodata = data.frame(A$ID, gamma, beta, annotation_matrix)
#write.table(Annodata,"Annodata.txt",col.names = F,quote = F,row.names = F)
