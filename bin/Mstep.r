library(optimx)
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=F)

ptm <- proc.time()

######## Need to pass args(hypfile, paramfile, k, hypcurrent_file) from bash
args <- commandArgs(TRUE)
#hyptemp#.txt
#hypfile = args[[1]]
rv = args[[1]]
#the number of iterition 
k = as.numeric(args[[2]])
#hyper.current
hypcurrent_file = args[[3]]
#Annotation file <- [SNPs,c(SNPname, gamma, beta, annotation1, annotation2, annotation3...)]
Annofile = args[[4]]
print(hypcurrent_file)

wkDir = args[[5]]
########## Load beta and gamma 
Annodata = read.table(Annofile, sep="\t", header=FALSE)
#Number of annotation files
Anum = ncol(Annodata) - 3
temp_col_names <- c("SNP","gamma", "beta")
for(i in 1:Anum){
	temp_col_names <- c(temp_col_names, 
	paste("Annotation", i, sep = "_"))
}
colnames(Annodata) <- temp_col_names

############hyptemp.txt
#hypdata = read.table(hypfile, sep="\t", header=FALSE)
#temp_col_names_hyp <- c("block", "loglike", "GV", "rv")
#colnames(hypdata) <-  temp_col_names_hyp

########### Update hyper parameter values
#rv is phenotype variance
#rv = mean(hypdata[, "rv"])
tau = 1.0 / as.numeric(rv)
#pve = sum(hypdata[, "GV"])

prehyp <- read.table(hypcurrent_file, header=F)
print("hyper parameter values before MCMC: ")
print(prehyp)
#a_old = prehyp[,1]
#b_old = prehyp[,2]
a_old = upper = c(-20, rep(0,Anum))
b_old = upper = c(1, rep(0,Anum))
a_variance = prehyp[,3]
b_variance = prehyp[,4]
inverse_a_Var_matrix = diag(1/a_variance)
inverse_b_Var_matrix = diag(1/b_variance)

#### updating hyper a vector and b vector
#hypcurrent <- NULL
hypmat <- NULL

Annodata = Annodata[Annodata$gamma > 0,]
beta_temp = Annodata[, "beta"]
gamma_temp = Annodata[,"gamma"]
A_temp = cbind(rep(1,nrow(Annodata)),Annodata[,-c(1:3)])
A_temp = as.matrix(A_temp)

####################################################
a_fn <- function(a) {
  asum = sum(gamma_temp * (A_temp %*% a) - log(1 + exp(A_temp %*% a))) - sum(1/2 * t(a) %*% (inverse_a_Var_matrix * a)) 
  asum = as.vector(asum)
  return(-asum)
}

a_gr <- function(a) {
  a_gr_sum = apply(gamma_temp * A_temp - as.vector(1/(1 + exp(-A_temp %*% a))) * A_temp,2,sum) - apply(inverse_a_Var_matrix * a, 2, sum)
  a_gr_sum = as.vector(a_gr_sum)
  return(-a_gr_sum)
}

a_Hess <- function(a) {
  a_Hess_sum = matrix(0,nrow(inverse_a_Var_matrix),ncol(inverse_a_Var_matrix))
  for (i in 1:nrow(A_temp)) {
    a_Hess_value = as.numeric(-exp(-A_temp[i,] %*% a)/(1 + exp(-A_temp[i,] %*% a)) ^ 2) * outer(A_temp[i,],A_temp[i,]) 
    a_Hess_sum = a_Hess_sum + a_Hess_value
  }
  a_Hess_sum = a_Hess_sum - inverse_a_Var_matrix
  a_Hess_sum = as.matrix(a_Hess_sum)
  return(-(a_Hess_sum))
}

a_temp = optimx(a_old, a_fn, method='L-BFGS-B', gr = a_gr, hess = a_Hess, upper = c(-14, rep(10,Anum)))[1:nrow(inverse_a_Var_matrix)]
a_new_Variance = 1/diag(a_Hess(a_old))
####################################################
b_fn <- function(b){
  bsum = sum(gamma_temp/2 *(- (A_temp %*% b) - tau * beta_temp ^ 2 * (1/exp(A_temp %*% b)))) - sum(1/2 * t(b) %*% (inverse_b_Var_matrix * b)) 
  bsum = as.vector(bsum)
  return(-bsum)
}

b_gr <- function(b){
  b_gr_sum = apply(gamma_temp/2 *(- A_temp + as.vector(tau * beta_temp ^ 2 * (1/exp(A_temp %*% b))) * A_temp), 2, sum) - apply(inverse_b_Var_matrix * b, 2, sum) 
  b_gr_sum = as.vector(b_gr_sum)
  return(-b_gr_sum)
}

b_Hess <- function(b) {
  b_Hess_sum = matrix(0,nrow(inverse_b_Var_matrix),ncol(inverse_b_Var_matrix))
  for (i in 1:nrow(A_temp)) {
    b_Hess_value = as.numeric(- gamma_temp[i]/2 * tau * beta_temp[i] ^ 2 * exp(-A_temp[i,] %*% b)) * outer(A_temp[i,],A_temp[i,])
    b_Hess_sum = b_Hess_sum + b_Hess_value
  }
  b_Hess_sum = b_Hess_sum - inverse_b_Var_matrix
  b_Hess_sum = as.matrix(b_Hess_sum)
  return(-b_Hess_sum)
}

b_temp = optimx(b_old, b_fn, method='L-BFGS-B', gr = b_gr, hess = b_Hess)[1:nrow(inverse_b_Var_matrix)]
b_new_Variance = 1/diag(b_Hess(b_old))
#####################################################
#hypcurrent <- c(a_temp, b_temp)
hypmat <- data.frame(as.vector(t(a_temp)), as.vector(t(b_temp)), a_new_Variance, b_new_Variance)

########## Write out updated hyper parameter values
colnames(hypmat) <- c("#a", "b", "a_variance", "b_variance")
print("hyper parameter values updates after MCMC: ")
print(hypmat)
write.table(format(hypmat, scientific=TRUE), 
            file=hypcurrent_file, 
            quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)

write.table(Annodata,
            file=paste0(wkDir , "/Anno_data", k, ".txt"), 
            quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)

write.table(format(hypmat, scientific=TRUE), 
            file=paste0(hypcurrent_file, "-", k), 
            quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)


#### Summarize log-likelihood
#loglike_total = sum(hypdata$loglike)

########## Write out updated hyper parameter values and se to EM_result_file
#hypcurrent = c(pve, loglike_total, hypcurrent)
#hypcurrent <- format(hypcurrent, scientific = TRUE)
#print("write to hypcurrent file with hyper parameter values after MCMC: ")
#print(c(k, hypcurrent))
#write.table(matrix(c(k, hypcurrent), nrow=1), file = EM_result_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)

print("EM step time cost (in minutes) : ")
print((proc.time() - ptm)/60)








