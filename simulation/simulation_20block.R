library(statip)
library(wtest)

Anum = 4
a = c(-9.5, 4, 0.8, 2.3, 0)
b = c(0.04, 0.11, 1.6, 0.07, 0)
h2 = 0.25

DataDir = "/mnt/YangFSS/data/jchen/simulation_20blocks/original_genotype/output"
BedPattern = "*.geno"
myFiles <- list.files(path = DataDir, pattern = BedPattern)

#Assemble genotype and annotation data
A_final = NULL
geno_final = NULL
for (file in myFiles) {
  geno = read.table(paste0(DataDir, "/", file), stringsAsFactors = F, header = T)
  rownames(geno) = geno$ID
  samplename = substr(file, 1, nchar(file) - (nchar(BedPattern) - 1))
  
  #Annotation file
  #1.95CredibleSet
  CredS = read.table(paste0("/mnt/YangFSS/data/jchen/Brain_Cortex_annotation/95CredibleSet/", samplename, ".anno")
                     , stringsAsFactors = F, header = F)
  colnames(CredS) = c("chr", "bp", "rs", "cm", "CredS")
  A = merge(geno[,1:5], CredS[,c("bp", "CredS")], by.x = "POS", by.y = "bp", all.x = T)
  
  AllciseQTL = read.table(paste0("/mnt/YangFSS/data/jchen/Brain_Cortex_annotation/AllciseQTL/", samplename, ".anno")
                          , stringsAsFactors = F, header = F)
  colnames(AllciseQTL) = c("chr", "bp", "rs", "cm", "AllciseQTL")
  A = merge(A, AllciseQTL[,c("bp", "AllciseQTL")], by.x = "POS", by.y = "bp", all.x = T)
  
  MaxCPP = read.table(paste0("/mnt/YangFSS/data/jchen/Brain_Cortex_annotation/MaxCPP/", samplename, ".anno")
                      , stringsAsFactors = F, header = F)
  colnames(MaxCPP) = c("chr", "bp", "rs", "cm", "MaxCPP")
  A = merge(A, MaxCPP[,c("bp", "MaxCPP")], by.x = "POS", by.y = "bp", all.x = T)
  
  A[is.na(A)] <- 0
  A$zeroAnno =  round(rnorm(nrow(A), mean = 0, sd = 1),4)
  A$annotation = paste0(A$CredS, ",", A$AllciseQTL, ",", A$MaxCPP, ",", A$zeroAnno)
  write.table(A[,c(2, 1, 3:5, 10)],paste0("/mnt/YangFSS/data/jchen/simulation_20blocks/Anno/Anno_", 
                                   samplename),col.names = F,quote = F,row.names = F, sep = "\t")
  geno = geno[A$ID,]
  A_final = rbind(A_final, A)
  geno_final = rbind(geno_final, geno)
}

save(A_final, file = "/mnt/YangFSS/data/jchen/simulation_20blocks/Annotation_final.rda")
save(geno_final, file = "/mnt/YangFSS/data/jchen/simulation_20blocks/geno_final.rda")
#122745


##### Set target proportion of SNPs, average of pi
log_pi = a[1]
for (i in 1:(Anum)) {
  log_pi = log_pi + a[i + 1] * A_final[, i + 5]
}
pi = exp(log_pi)/(1 + exp(log_pi))
png("/mnt/YangFSS/data/jchen/simulation_20blocks/test/pi.png")
plot(pi)
dev.off()
mean(pi)

log_sigma_2 = b[1]
for (i in 1:(Anum)) {
  log_sigma_2 = log_sigma_2 + b[i + 1] * A_final[, i + 5]
}
sigma_2 = exp(log_sigma_2)
mean(sigma_2)
png("/mnt/YangFSS/data/jchen/simulation_20blocks/test/sigma_2.png")
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
save(error, beta, gamma, sigma_2, pi, file = "/mnt/YangFSS/data/jchen/simulation_20blocks/test/parameter.rda")

pheno =  pheno  + error
var1 / var(pheno)
var(pheno)
#58.12276
png("/mnt/YangFSS/data/jchen/simulation_20blocks/pheno.png")
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

write.table(pheno_1, "/mnt/YangFSS/data/jchen/simulation_20blocks/pheno.txt", col.names = F, row.names = F, quote = F)

#Annodata = data.frame(A$ID, gamma, beta, annotation_matrix)
#write.table(Annodata,"Annodata.txt",col.names = F,quote = F,row.names = F)
