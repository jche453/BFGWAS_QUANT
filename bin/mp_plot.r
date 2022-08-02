# Rscript to generate Manhattan plots highlighed by posterior probabilities PI

Sys.setlocale('LC_ALL','C') 
options(stringsAsFactors=F)

library(tidyverse) # filter data
library(data.table) # load gwas data
library(autoimage) # legend plot
library(plotrix) # axis.break()

#### read input arguments #######
args=(commandArgs(TRUE))
if(length(args)==0) {
    stop("Error: No arguments supplied. You need to provide at least param.txt in the command line")
} else {
    for (i in 1:length(args)) eval(parse(text=args[[i]]))
}

##############
sig_level = as.numeric(sig_level)
yRedLine <- -log10(sig_level)
break.top=20;
pngwidth <- 12800
pngheight <- 6800
plotRedLine <- TRUE;
top.size = 0.3
#paramFile="/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/Test_wkdir_AD/Eoutput/paramtemp3.txt"
#out_prefix="/home/jyang/ResearchProjects/BFGWAS_QUANT_Test/AnalyzeResults/AD"

############# Load gwas data 
if(!file.exists(paramFile)) stop("Error: Input file not found")
file_out <- paste0(out_prefix, "_mp.png")
print("Load BFGWAS Results from:")
print(paramFile)
print("Make Manhanttan Plot : ")
print(file_out)

gwas <- fread(paramFile, sep = "\t", header = TRUE)
setnames(gwas, c("CHR", "POS", "ID", "REF", "ALT", "MAF", "Pi", "Beta", "mBeta", "Chisq", "Pval", "Rank", "Anno"))
gwas <- gwas[!is.na(gwas$Pval), ]

#### Calculate -log10_pval
pval <- gwas$Pval
mlog10_pval <- -log10(pval)
mlog10_pval[mlog10_pval == Inf] <- max(mlog10_pval[mlog10_pval < Inf]) + runif(sum(mlog10_pval == Inf))

# set chromosome colors / prepare data frame and reorder variants
# gwas is a data.table here
chrs <- sort(unique(gwas$CHR))
p = dim(gwas)[1]
pcol = rep("grey40", length = p)
pcol[gwas$CHR %% 2 == 0] <- "grey60"

gwas_man <- data.frame(chr=gwas$CHR, bp=gwas$POS, mlog10_pval = mlog10_pval, pi = gwas$Pi, check.names=F)
gwas_man <- data.frame(gwas_man, plotPos=NA, highlightColor=NA, pch=20,
    pipColor=NA, pcol = pcol, check.names=F)
gwas_man <- dplyr::arrange(gwas_man, chr, bp)
rm(gwas)

rbPal <- colorRampPalette(c('lightpink','darkred'))
if(sum(gwas_man$pi>0.1068)>0){
    gwas_man$pipColor[gwas_man$pi>0.1068] <- rbPal(10)[as.numeric(cut(gwas_man$pi[gwas_man$pi>0.1068],breaks = 10))]
}

# determine gap between chromosomes
chrGAP <- 1.5E9/(pngwidth /16)
endPos <- 0
chrs = sort(unique(gwas_man$chr))
plotPos <- numeric(0); chrLab <- numeric(0); chrEnd <- numeric(0);
for(chr_num in chrs){
    chrTemp <- dplyr::filter(gwas_man, chr == chr_num)
    chrPOS <- chrTemp$bp - min(chrTemp$bp, na.rm=TRUE) + endPos + 1
    chrLab <- c(chrLab, mean(chrPOS, na.rm=T))
    endPos <- max(chrPOS, na.rm=T) + chrGAP
    plotPos <- c(plotPos, chrPOS)
    yline_pos <- max(chrPOS, na.rm=T) + chrGAP/2
    chrEnd <- c(chrEnd, yline_pos)
}
gwas_man$plotPos <- plotPos
print(head(gwas_man))

###Testing with 1000 points 
#gwas_man_temp = gwas_man
#gwas_man <- gwas_man_temp[sort(sample(1:nrow(gwas_man_temp), 1000)), ]

## Function to make legend
legend.col <- function(col, lev){
    opar <- par
    n <- length(col)
    bx <- par("usr")
    box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
    bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3]) / n
    xx <- rep(box.cx, each = 2)
    par(xpd = TRUE)
    for(i in 1:n){
        yy <- c(box.cy[1] + (box.sy * (i - 1)),
        box.cy[1] + (box.sy * (i)),
        box.cy[1] + (box.sy * (i)),
        box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
    }
    par(new = TRUE)
    plot(0, 0, type = "n",
    ylim = c(min(lev), max(lev)),
    yaxt = "n", ylab = "",
    xaxt = "n", xlab = "",
    frame.plot = FALSE)
    axis(side = 4, las = 2, tick = FALSE, line = .25, cex.axis=1.2, cex.lab=1.2)
    par <- opar
}


##### Making plots
png(filename = file_out, width = pngwidth, height = pngheight, pointsize = 16, res=600)
par(mar=c(5.1, 5.1, 2.1, 3.5), mgp = c(3, 0.7, 0), las=1)
    x = gwas_man$plotPos
    y = gwas_man$mlog10_pval

    if(max(y,na.rm=T) > break.top){

        # Manhattan plot with two different y axis scales
        # set axis labels of both scales
        lab1 <- pretty(c(0,break.top), n=ceiling(12 * (1-top.size)))
        lab1 <- c(lab1[lab1 < break.top], break.top)
        lab2 <- pretty(c(break.top,max(y,na.rm=T)), n=max(3,floor(12 * top.size)))
        lab2 <- lab2[lab2 > max(lab1)]
   
        # resulting range of top scale in bottom scale units
        top.range = break.top/(1 - top.size) - break.top
        top.data = max(lab2) - break.top
        # function to rescale the top part
        rescale = function(y) { break.top+(y-break.top)/(top.data/top.range)}

       # png(filename = file_out, width = pngwidth, height = pngheight, pointsize = 16, res=600)
        # plot bottom part / rescaled
        plot(x[y<break.top], y[y<break.top], ylim = c(0, break.top+top.range),
            axes=FALSE,
            pch=20, cex=1.2, cex.lab=1.2, cex.axis=1.2, xaxt="n",
            col=gwas_man$pcol[y<break.top], ylab=expression(-log[10]*italic(P)), xlab="Chromosome", bty="L", main="")
        # plot top part
        points(x[y>break.top],rescale(y[y>break.top]),pch=20,
            col=gwas_man$pcol[y>break.top], cex=1.2)

        # SNPs with PIP > 0.1068
        gwas_sub1_top = gwas_man[gwas_man$pi > 0.1068 & gwas_man$pi <= 0.5 & gwas_man$mlog10_pval>= break.top, ]
        gwas_sub1_bottom = gwas_man[gwas_man$pi > 0.1068 & gwas_man$pi <= 0.5 & gwas_man$mlog10_pval < break.top, ]
        points(gwas_sub1_top$plotPos, rescale(gwas_sub1_top$mlog10_pval),
                    pch=gwas_sub1_top$pch, col=gwas_sub1_top$pipColor, cex=1.2)
        points(gwas_sub1_bottom$plotPos, gwas_sub1_bottom$mlog10_pval,
                    pch=gwas_sub1_bottom$pch, col=gwas_sub1_bottom$pipColor, cex=1.2)
        # SNPs with PIP > 0.5
        gwas_sub1_top = gwas_man[gwas_man$pi > 0.5 & gwas_man$mlog10_pval>= break.top, ]
        gwas_sub1_bottom = gwas_man[gwas_man$pi > 0.5 & gwas_man$mlog10_pval < break.top, ]
        points(gwas_sub1_top$plotPos, rescale(gwas_sub1_top$mlog10_pval),
                    pch=17, col=gwas_sub1_top$pipColor, cex=1.2)
        points(gwas_sub1_bottom$plotPos, gwas_sub1_bottom$mlog10_pval,
                    pch=17, col=gwas_sub1_bottom$pipColor, cex=1.2)

        # add axes and axis labels
        if(length(chrs) > 1) {
                axis(1, at=chrLab[seq(1,length(chrLab),by=2)],
                    labels=chrs[1:length(chrLab)][seq(1,length(chrLab), by=2)],
                    las=1, tick=F, cex.axis=1.2, cex.lab=1.2, line=1)
                axis(1,at=chrLab[seq(2,length(chrLab),by=2)],
                    labels=chrs[1:length(chrLab)][seq(2,length(chrLab), by=2)],
                    las=1,tick=F, cex.axis=1.2,cex.lab=1.2,line=0)
            } else {
                axis(1, at=chrLab[1], labels=chrs[1], las=1,tick=F,cex.axis=1.2,cex.lab=1.2,line=2)
            }
        axis(side=2, at=lab1, cex.axis=1.2, cex.lab=1.2)
        axis(side=2, at=rescale(lab2), labels=lab2, cex.axis=1.2, cex.lab=1.2)

        # plot axis breaks and indicate line of axis break
        box()
        axis.break(axis=2, breakpos=break.top, style="zigzag", brw=0.02)
        axis.break(axis=4, breakpos=break.top, style="zigzag", brw=0.02)
        abline(h=break.top, lwd=1.5, lty=2, col="grey")
        if(plotRedLine) {
            if(yRedLine <= break.top) {
                abline(h=yRedLine, lwd=1.5, col="red", lty=2)
            } else {
                abline(h=rescale(yRedLine), lwd=1.5, col="magenta", lty=2)
            }
        }
    } else {
        # Standard Manhattan plot if no association signal above break
        # png(filename = file_out, width = pngwidth, height = pngheight, pointsize = 16, res=600)
        par(mar=c(5.1, 5.1, 2.1, 3.5), mgp = c(3, 0.7, 0), las=1)
        plot(x, y, xaxt="n", pch=20, cex=1.2, cex.lab=1.2,
            cex.axis=1.2, xaxt="n",
            col=gwas_man$pcol, ylab=expression(-log[10]*italic(P)),
            xlab="", bty="L", main="", ylim = c(0, max(c(yRedLine*1.1, y))) )
        if(length(chrs) > 1) {
            axis(1,at=chrLab[seq(1,length(chrLab),by=2)],
                labels=chrs[1:length(chrLab)][seq(1,length(chrLab),by=2)],
                las=1,tick=F,cex.axis=1.2,cex.lab=1.2,line=1)
            axis(1,at=chrLab[seq(2,length(chrLab),by=2)],
                labels=chrs[1:length(chrLab)][seq(2,length(chrLab),by=2)],
                las=1,tick=F,cex.axis=1.2,cex.lab=1.2,line=0)
            } else {
            axis(1,at=chrLab[1], labels=chrs[1], las=1,tick=F,cex.axis=1.2,cex.lab=1.2,line=2)
                }
        gwas_sub1 = gwas_man[gwas_man$pi > 0.1068 & gwas_man$pi <= 0.5, ]
        points(gwas_sub1$plotPos, gwas_sub1$mlog10_pval,
                    pch=gwas_sub1$pch, col=gwas_sub1$pipColor, cex=1.2)
        gwas_sub1 = gwas_man[gwas_man$pi > 0.5, ]
        points(gwas_sub1$plotPos, gwas_sub1$mlog10_pval,
                    pch=17, col=gwas_sub1$pipColor, cex=1.2)
        if(plotRedLine) abline(h=yRedLine,lwd=1.5,col="magenta",lty=2)
    }

legend("topleft", legend = "CPP>0.5", pch = 17, col = rbPal(10)[6], cex = 1.2, pt.cex = 1.2)

cpp = round(seq(0.1068, 1, length.out = 10), 2)
par(mar=c(5.1, 5.1, 2.1, 4.1), mgp = c(2, 2, 2), las=1)
legend.col(col = rbPal(10), lev = cpp)
dev.off()


