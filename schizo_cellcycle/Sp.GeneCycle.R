library("GeneCycle")
library("qvalue")
library("longitudinal")

Spf <-read.csv("Sp.genes_in.tsv",sep="\t",header=T,row.names=1)
summary(Spf)
#head(Spf)
#t(Spf)
Sp <- as.longitudinal(t(Spf),time=c(0,20,40,60,80,100,120,140))
dim(Sp)
get.time.repeats(Sp)

pdf("Spombe_plots.pdf")
plot(Sp,1:20)
plot(Sp,21:40)


cons = is.constant(Sp)
sum(cons) 

find.periodic.genes <- function(dataset)
{
  cat("Computing p-values ...\n")
  pval = fisher.g.test(dataset)   
  n1 = sum(qvalue(pval, lambda=0)$qvalues < 0.05)
  n2 = sum(qvalue(pval)$qvalues < 0.05)
  fdr.out <- fdrtool(pval, statistic="pvalue", plot=FALSE)
  n3 <- sum(fdr.out$qval < 0.05) 
  n4 <- sum(fdr.out$lfdr < 0.2) 
  cat("Conservative estimate (Wichert et al.) =", n1, "\n")
  cat("Less conservative estimate =", n2, "\n")
  cat("Semiparametric Fdr < 0.05 (fdrtool) =", n3, "\n")
  cat("Semiparametric fdr < 0.2 (fdrtool) =", n4, "\n")
}
#find.periodic.genes(Sp)

pval = fisher.g.test(Sp)   
fdr.out <- fdrtool(pval, statistic="pvalue", plot=FALSE)
fd_filt <- (fdr.out$lfdr < 0.2)
fd_filt
#Sp_filter <- Sp[fd_filt==TRUE]
#Sp_filter

#pval.estimate.eta0(pval.Sp, method="conservative")
#pval.estimate.eta0(pval.Sp, method="adaptive")
#pval.estimate.eta0(pval.Sp, method="bootstrap")
#pval.estimate.eta0(pval.Sp, method="smoother")

#avgp.Sp <- avgp(Sp,"S.pombe")

#dominant.freqs(Sp, 3)
#dominant.freqs(Sp,6)


#periodogram(Sp)

#spe5 <- robust.spectrum(Sp)

#g.statistic(spe5)
#pval = robust.g.test(spe5)

#pval
#pval = robust.g.test(spe5, 4)
#pval
