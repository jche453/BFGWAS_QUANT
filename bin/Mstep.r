library(optimx)
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=F)

# Package names
packages <- c("data.table", "tidyverse", "optimx")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#################################
ptm <- proc.time()

######## Read input arguments
args <- commandArgs(TRUE)
hypfile=args[[1]]
k = as.numeric(args[[2]]) #the number of iterition
a_gamma = as.numeric(args[[3]])
b_gamma = as.numeric(args[[4]])
gwas_n = as.numeric(args[[5]]) # GWAS sample size
hypcurrent_file = args[[6]] # hyper.current
# sum(gamma) and sum(gamma * beta^2) are saved in hypfile
param_file = args[[7]] # Annotation data and gamma
EM_result_file = args[[8]] # Save EM_result_file

########### Read prior hyper parameter values
#rv is phenotype variance, rv = mean(hypdata[, "rv"])
# hypcurrent_file = "/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_Data/hypval_10anno.txt"
prehyp <- read.table(hypcurrent_file, header=F)
print("hyper parameter values before MCMC: ")
print(prehyp)
avec_old = str_split(prehyp[1, 2], ",") %>% unlist() %>% as.numeric()
tau_beta_old = as.numeric(prehyp[2, 2])

########## Read hyptemp file
# hypfile = "/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir/output/WGS_1898_samples_CHR_19_43862455_44744108.hyptemp"
hyp_dt = fread(hypfile, sep = "\t", header = TRUE)
sum_gamma = sum(hyp_dt$sum_gamma)
sum_beta2 = sum(hyp_dt$sum_beta2)

########## Read annotation data
# param_file="/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir/output/WGS_1898_samples_CHR_19_43862455_44744108.paramtemp"
param_dt = fread(param_file, sep="\t", header=TRUE)
gamma_temp = param_dt$Pi
split_anno <- function(x){
  xvec = str_split(x, ",") %>% unlist()%>% t() %>% data.frame()
  return(xvec)
}
Anno_df = map_dfr(param_dt$Anno, split_anno)
Anum = ncol(Anno_df)
A_temp = data.frame(intercept = 1, Anno_df) %>% unlist() %>% as.numeric() %>% matrix(nrow = nrow(param_dt))
gamma_A = gamma_temp * A_temp

####################################################
####### Solve for a_vec
a_fn <- function(a) {
  Ata = A_temp %*% a
  asum = sum(gamma_temp * Ata - log(1 + exp(Ata))) - 0.5 * sum(a * a)
  return(-asum)
}

a_gr <- function(a) {
  Ata = A_temp %*% a
  a_gr_sum = apply( (gamma_A - as.vector(1/(1 + exp(-Ata))) * A_temp) , 2, sum) - a
  return( as.vector(-a_gr_sum) )
}

a_Hess <- function(a) {
  a_num = length(a)
  exp_mAta = exp(-A_temp %*% a)
  a_Hess_sum = matrix(0, a_num, a_num)
  for (i in 1:nrow(A_temp)) {
    a_Hess_sum = a_Hess_sum + (exp_mAta[i] / ((1 + exp_mAta[i])^2) ) * outer(A_temp[i,], A_temp[i,])
  }
  a_Hess_sum = a_Hess_sum + diag(rep(1, a_num))
  return( as.matrix(a_Hess_sum))
}

a_temp = optimx(avec_old, fn = a_fn, method='L-BFGS-B', gr = a_gr, hess = a_Hess,
                upper = c(avec_old[1], rep(10, Anum)), lower = c(avec_old[1], rep(-10, Anum)) ) [1:length(avec_old)]

####################################################
# gwas_n = 2453; a_gamma = 1.01; b_gamma = 1
Est_tau_beta <- function(sum_gamma, sum_beta2, gwas_n, a, b){
  tau_beta_hat =  ( sum_gamma + 2 * (a - 1)) / ( gwas_n * sum_beta2 + 2 * b)
  return(tau_beta_hat)
}

CI_fish_tau_beta <- function(sum_gamma, sum_beta2, gwas_n, a, b){
  tau_beta_hat =  ( sum_gamma + 2 * (a - 1)) / ( gwas_n * sum_beta2 + 2 * b)
  temp = 0.5 * sum_gamma + (a - 1)
  if( temp  < 0){
    se_tau_beta = 0
  }else{
    se_tau_beta = tau_beta_hat / sqrt(temp)
  }
  return(c(tau_beta_hat, se_tau_beta))
}

tau_beta_temp = CI_fish_tau_beta(sum_gamma, sum_beta2, gwas_n, a_gamma, b_gamma)
print("Estimates for tau_beta, and tau_beta_SE: ")
print(tau_beta_temp)

#####################################################
hypmat <- data.table(`#hyper_parameter` = c("a", "tau_beta"), value = c(paste(a_temp, collapse = ",") , tau_beta_temp[1]))

########## Write out updated hyper parameter values
print("hyper parameter values updates after MCMC: ")
print(hypmat)
# hypcurrent_file="/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir/hyper.current.txt"
write.table(format(hypmat, scientific=TRUE), 
            file=hypcurrent_file,
            quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)

########## Write out updated hyper parameter values and se to EM_result_file
R2 = sum(hyp_dt$r2)
loglike = sum(hyp_dt$log_post_likelihood)
print(paste("Sum PIP = ", sum_gamma))
print(paste("Regression R2 = ", R2))
print(paste("Posterior log likelihood = ",  loglike))

# EM_result_file="/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir/EM_result.txt"
if(k==0){
  write.table(data.frame(EM_iteration = k, R2 = R2, Loglike = loglike,
                        tau_beta = tau_beta_temp[1], avec = paste(a_temp, collapse = ",")),
              file = EM_result_file,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append=FALSE)
}else{
  write.table(data.frame(EM_iteration = k, R2 = R2, Loglike = loglike,
                        tau_beta = tau_beta_temp[1], avec = paste(a_temp, collapse = ",")),
              file = EM_result_file,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append=TRUE)
}

print("EM step time cost (in minutes) : ")
print((proc.time() - ptm)/60)








