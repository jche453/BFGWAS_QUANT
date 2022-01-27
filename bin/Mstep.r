library(optimx)
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=F)

ptm <- proc.time()

######## Need to pass args(hypfile, paramfile, k, hypcurrent_file) from bash
args <- commandArgs(TRUE)
rv = args[[1]]
#the number of iterition 
k = as.numeric(args[[2]])
#hyper.current
hypcurrent_file = args[[3]]
#Annotation file <- [SNPs,c(SNPname, gamma, beta, annotation1, annotation2, annotation3...)]
Annofile = args[[4]]
print(hypcurrent_file)
wkDir = args[[5]]
abgamma = as.numeric(args[[6]])
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
a_gamma = abgamma
b_gamma = abgamma

########### Update hyper parameter values
#rv is phenotype variance, rv = mean(hypdata[, "rv"])
tau = 1.0 / as.numeric(rv)
prehyp <- read.table(hypcurrent_file, header=F)
print("hyper parameter values before MCMC: ")
print(prehyp)

a_old = prehyp[,1]
a_variance = prehyp[,3]
inverse_a_Var_matrix = diag(1/a_variance)

#### updating hyper a vector and b vector
hypmat <- NULL

if (nrow(Annodata[Annodata$gamma > 0,]) > 0){
  Annodata = Annodata[Annodata$gamma > 0,]
}
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

a_new_Variance = diag(solve(a_Hess(a_old)))
#e = try( optimx(a_old, a_fn, method='L-BFGS-B', gr = a_gr, hess = a_Hess, upper = c(-12, rep(10,Anum)))[1:nrow(inverse_a_Var_matrix)], TRUE)
#if (class(e) != "try-error"){
#  a_temp = optimx(a_old, a_fn, method='L-BFGS-B', gr = a_gr, hess = a_Hess, upper = c(-12, rep(10,Anum)))[1:nrow(inverse_a_Var_matrix)]
#}else{
  a_old = upper = c(-40, rep(0,Anum))
  a_temp = optimx(a_old, a_fn, method='L-BFGS-B', gr = a_gr, hess = a_Hess, upper = c(-12, rep(10,Anum)))[1:nrow(inverse_a_Var_matrix)]
#}
####################################################
Est_sigma2 <- function(sigma2, m, tau, a, b){
  sigma2_hat = (sum(sigma2 * m) * tau + 2 * b) / (sum(m) + 2 * (a + 1))
  return(sigma2_hat)
}

CI_fish_sigma2 <- function(sigma2, m, tau, a, b){
  sigma2_hat = Est_sigma2(sigma2, m, tau, a, b)
  if( (sum(m) * tau - sum(m)/2 - (a+1) + 2*b/sigma2_hat) < 0){
    se_sigma2=0
  }else{
    se_sigma2 = sigma2_hat * sqrt(1/(sum(m) * tau - sum(m)/2 - (a+1) + 2*b/sigma2_hat))
  }
  return(c(sigma2_hat, se_sigma2))
}

sigma2_temp = CI_fish_sigma2(gamma_temp, beta_temp, tau, a_gamma, b_gamma)

b_temp = rep(sigma2_temp[1], length(a_temp))
b_new_Variance = rep(sigma2_temp[2], length(a_temp))
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








