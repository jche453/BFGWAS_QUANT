library(statip)
library(wtest)

setwd("/mnt/YangFSS/data/jchen/simulation_20blocks/wkdir_1")
load("parameter.rda")
load("/mnt/YangFSS/data/jchen/simulation_20blocks/geno_final.rda")
geno = geno_final
#geno_final = geno_final[,1:5]
geno_final$beta = beta
geno_final$gamma = gamma
geno_final$sigma_2 = sigma_2
geno_final$pi = pi

Annodata= read.table("Eoutput/AnnoScore5.txt")
Anum = ncol(Annodata) - 3
temp_col_names <- c("SNP","gamma", "beta")
for(i in 1:Anum){
  temp_col_names <- c(temp_col_names, 
                      paste("Annotation", i, sep = "_"))
}
colnames(Annodata) <- temp_col_names

prehyp <- read.table("/mnt/YangFSS/data/jchen/simulation_20blocks/wkdir_1/hypval.current-5", header=F)
a = prehyp[,1]
b = prehyp[,2]

log_pi = a[1]
for (i in 1:(Anum)) {
  log_pi = log_pi + a[i + 1] * Annodata[,i+3]
}
Annodata$pi = exp(log_pi)/(1 + exp(log_pi))
png("/mnt/YangFSS/data/jchen/simulation_20blocks/wkdir_1/pi_pos.png")
plot(Annodata$pi)
dev.off()
sum(Annodata$pi)
#266.3301
dim(Annodata[Annodata$gamma > 0.1, ])
#269

log_sigma_2 = b[1]
for (i in 1:(Anum)) {
  log_sigma_2 = log_sigma_2 + b[i + 1] * Annodata[,i+3]
}
Annodata$sigma_2 = exp(log_sigma_2)

combine = merge(geno_final, Annodata[c(1:3,8:9)], by.x = "ID", by.y = "SNP")

# of casual snps 
nrow(combine[combine$gamma.x > 0 & combine$gamma.y > 0.1 , ])
#2
nrow(Annodata[Annodata$gamma >0.1, ])
#269

sum(combine$gamma.x)
#20
sum(Annodata$gamma)
#199.9728

#Pheno R2
#x: true
pheno_x = read.table("pheno.txt", stringsAsFactors = F)
pheno_x$V1[1:698] = paste0("X", pheno_x$V1[1:698])

#y: model
combine =combine[combine$beta.y > 0 ,]
pheno_y = t(combine[,(6: 1899)]) %*% combine$beta.y
pheno_y = data.frame(pheno_y)
pheno_y$V1 = rownames(pheno_y)
pheno_y$V1 = gsub("[/.]", "-", pheno_y$V1)

pheno = merge(pheno_y, pheno_x, by = "V1")
cor(pheno[,2], pheno[,3]) ^ 2
#0.100503

cc = cor(t(combine[combine$gamma.x > 0, 6:1899]), t(combine[combine$gamma.y > 0, 6:1899])) ^2
cc[cc < 0.3] =0