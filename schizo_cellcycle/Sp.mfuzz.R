library("longitudinal")
library("Biobase")
library("Mfuzz")
exprs <- as.matrix(read.csv("Sp.genes_in.tsv", header=TRUE, sep="\t", row.names=1, as.is=TRUE))

class(exprs)
dim(exprs)
colnames(exprs) <- c(0,20,40,60,80,100,120,140)	

Spexpr <- ExpressionSet(assayData=exprs)
Sp.s <- standardise(Spexpr)

m1 <- mestimate(Sp.s)
m1

#cselection(Sp.s,m1,crange=seq(4,32,4),repeats=5,visu=TRUE)
#Dmin(Sp.s,m1,crange=seq(4,40,4),repeats=3,visu=TRUE)


cl <- mfuzz(Sp.s,c=8,m=m1)
acore.list <- acore(Sp.s,cl,min.acore=0.5)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("CLUSTER")
write.csv(file="Sp.mfuzz_k8.clusters.csv",dat,row.names=TRUE)

pdf("Sp.mfuzz_k8.pdf")
mfuzz.plot(Sp.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)


cl <- mfuzz(Sp.s,c=12,m=m1)
pdf("Sp.mfuzz_k12.pdf")
mfuzz.plot(Sp.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("CLUSTER")
write.csv(file="Sp.mfuzz_k12.clusters.csv",dat,row.names=TRUE)

cl <- mfuzz(Sp.s,c=16,m=m1)
pdf("Sp.mfuzz_k16.pdf")
mfuzz.plot(Sp.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("CLUSTER")
write.csv(file="Sp.mfuzz_k16.clusters.csv",dat,row.names=TRUE)

cl <- mfuzz(Sp.s,c=20,m=m1)
pdf("Sp.mfuzz_k20.pdf")
mfuzz.plot(Sp.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("CLUSTER")
write.csv(file="Sp.mfuzz_k20.clusters.csv",dat,row.names=TRUE)


cl <- mfuzz(Sp.s,c=24,m=m1)
pdf("Sp.mfuzz_k24.pdf")
mfuzz.plot(Sp.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("GENE")
write.csv(file="Sp.mfuzz_k24.clusters.csv",dat,row.names=TRUE)
