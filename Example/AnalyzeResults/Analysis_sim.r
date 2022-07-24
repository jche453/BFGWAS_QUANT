# rm(list=ls(all=TRUE))
library(data.table)
library(tidyverse)

####### Source Util Functions and set data directories
source("/home/jyang/GIT/BFGWAS_QUANT/bin/R_funcs.r")
setwd("/home/jyang/GIT/BFGWAS_QUANT/Example/")

######## Compare results
paramdata_bfgwas = LoadEMdata(filename="./Test_wkdir/Eoutput/paramtemp3.txt", header = TRUE)
#head(paramdata_bfgwas)
sum(paramdata_bfgwas$Pi)

## Manhantton plot
paramdata_bfgwas_sig <- filter(paramdata_bfgwas, Pi > 0.1)
dim(paramdata_bfgwas_sig)
paramdata_bfgwas_sig

###### Manhantton Plot
ggplot(paramdata_bfgwas, aes(x=POS, y = -log10(Pval), color = Pi)) +
	geom_point() +scale_color_gradient(low="blue", high="red") +
	facet_grid(cols = vars(CHR), scales = "free") +
	# geom_point(paramdata_bfgwas_sig, aes(x=POS, y = -log10(Pval), color = Pi)) +
	geom_hline(yintercept=-log10(5e-8))
ggsave("./AnalyzeResults/mp.pdf")

###### Effect size plot: Bayesian estimates vs. marginal estimates
ggplot(paramdata_bfgwas[paramdata_bfgwas$Beta>0, ], aes(x = mBeta, y = Beta, col = Pi)) +
	geom_point() + geom_abline(intercept=0, slope = 1) +
	scale_color_gradient(low="blue", high="red")
ggsave("./AnalyzeResults/beta.pdf")

######## Results of hyper parameter estimates #######
test_hyp <- LoadEMhyp(filename = "./Test_wkdir/Eoutput/EM_result.txt", header = TRUE)
print(test_hyp)

######## END ################










