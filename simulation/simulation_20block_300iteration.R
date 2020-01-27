#!/bin/env Rscript
library(statip)
library(wtest)

Anum = 4
a = c(-9.5, 4, 0.8, 2.3, 0)
b = c(0.04, 0.11, 1.6, 0.07, 0)
h2 = 0.25

load("/mnt/YangFSS/data/jchen/20block_simulation/Annotation_final.rda")
load("/mnt/YangFSS/data/jchen/20block_simulation/geno_final.rda")
#94342

##### Set target proportion of SNPs, average of pi
log_pi = a[1]
for (i in 1:(Anum)) {
  log_pi = log_pi + a[i + 1] * A_final[, i + 5]
}
pi = exp(log_pi)/(1 + exp(log_pi))
png("/mnt/YangFSS/data/jchen/20block_simulation/pi.png")
plot(pi)
dev.off()
mean(pi)

log_sigma_2 = b[1]
for (i in 1:(Anum)) {
  log_sigma_2 = log_sigma_2 + b[i + 1] * A_final[, i + 5]
}
sigma_2 = exp(log_sigma_2)
mean(sigma_2)
png("/mnt/YangFSS/data/jchen/20block_simulation/sigma_2.png")
hist(sigma_2)
dev.off()

generate_gamma = function(x) {rbern(1, x)}
gamma = sapply(pi,generate_gamma)
mean(gamma)
length(gamma[gamma>0])
#Number of casual snp is 18

### Generate beta 
generate_beta = function(x) {rnorm(1,0,sqrt(x))}
beta = sapply(sigma_2, generate_beta) * gamma
range(beta)

pheno =  t(geno_final[,-c(1:5,ncol(geno_final))]) %*% beta  # without error term
var1 = var(pheno)
var_error = var1 * (1/h2 - 1)
var_error / (var1 + var_error)
sd = sqrt(var_error)
error = rnorm(length(pheno), mean = 0, sd = sd)

#Save(error, beta, pi, sigma2)
save(error, beta, gamma, sigma_2, pi, file = "/mnt/YangFSS/data/jchen/20block_simulation/parameter.rda")

pheno =  pheno  + error
var1 / var(pheno)
var(pheno)
#58.12276
png("/mnt/YangFSS/data/jchen/20block_simulation/pheno.png")
hist(pheno)
dev.off()
pheno_1 = data.frame(rownames(pheno), pheno[,1])
pheno_1[,1] = gsub("[/.]", "-", pheno_1[,1])
startchange = function(x){
  if (startsWith(x, "X")) {
    newname = substr(x, 2, nchar(x))
  } else {
    newname = x
  }
  return(newname)
}
pheno_1[,1] = sapply(pheno_1[,1], startchange)

write.table(pheno_1, "/mnt/YangFSS/data/jchen/20block_simulation/pheno.txt", col.names = F, row.names = F, quote = F)

#Annodata = data.frame(A$ID, gamma, beta, annotation_matrix)
#write.table(Annodata,"Annodata.txt",col.names = F,quote = F,row.names = F)
