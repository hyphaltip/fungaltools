library("Biobase")
library("Mfuzz")
exprs <- as.matrix(read.csv("Sj.genes_in.tsv", header=TRUE, sep="\t", row.names=1, as.is=TRUE))

class(exprs)
dim(exprs)
colnames(exprs) <- c(0,20,40,60,80,100,120,140)	

Sjexpr <- ExpressionSet(assayData=exprs)
Sj.s <- standardise(Sjexpr)

m1 <- mestimate(Sj.s)
m1

#cselection(Sj.s,m1,crange=seq(4,32,4),repeats=5,visu=TRUE)
#Dmin(Sj.s,m1,crange=seq(4,40,4),repeats=3,visu=TRUE)


cl <- mfuzz(Sj.s,c=8,m=m1)
acore.list <- acore(Sj.s,cl,min.acore=0.5)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("CLUSTER")
write.csv(file="Sj.mfuzz_k8.clusters.csv",dat,row.names=TRUE)

pdf("Sj.mfuzz_k8.pdf")
mfuzz.plot(Sj.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)


cl <- mfuzz(Sj.s,c=12,m=m1)
pdf("Sj.mfuzz_k12.pdf")
mfuzz.plot(Sj.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("CLUSTER")
write.csv(file="Sj.mfuzz_k12.clusters.csv",dat,row.names=TRUE)

cl <- mfuzz(Sj.s,c=16,m=m1)
pdf("Sj.mfuzz_k16.pdf")
mfuzz.plot(Sj.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("CLUSTER")
write.csv(file="Sj.mfuzz_k16.clusters.csv",dat,row.names=TRUE)

cl <- mfuzz(Sj.s,c=20,m=m1)
pdf("Sj.mfuzz_k20.pdf")
mfuzz.plot(Sj.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("CLUSTER")
write.csv(file="Sj.mfuzz_k20.clusters.csv",dat,row.names=TRUE)


cl <- mfuzz(Sj.s,c=24,m=m1)
pdf("Sj.mfuzz_k24.pdf")
mfuzz.plot(Sj.s,cl=cl,mfrow=c(3,2),time.labels=colnames(exprs),new.window = FALSE)
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
dat <- data.frame(cl$cluster)
colnames(dat)<-c("GENE")
write.csv(file="Sj.mfuzz_k24.clusters.csv",dat,row.names=TRUE)
