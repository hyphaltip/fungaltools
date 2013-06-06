#   EXAMPLE FOR THE USE OF fdr.cycle
#
#   Copyright (C) Matthias E. Futschik,
#   Institute of Theoretical Biology,
#   Charite, Humboldt-University
#   Berlin, Germany
#   (matthias.futschik@charite.de;m.futschik@staff.hu-berlin.de)
#
#
#   This code is free for academic research use.
#   For all other use, please contact Mattthias Futschik.
#   If you employing this code, please cite
#   Matthias E. Futschik & Hanspeter Herzel, Bioinformatics, 2008
#

library(Biobase) # Necessary bioconductor package 
library(Mfuzz)   # Bioconductor package - to obtain the yeast cell cycle dataset for the example
                 # and for data preprocessing 

source("d:/Matthias/R/fdr.cycle.R") # MODIFY TO THE DIRECTORY WHERE fdr.cycle.R is saved or
                                    # copy and paste fdr.cycle.R in the workspace  


data(yeast) # loading the reduced CDC28 yeast set 
ow <- options("warn")
options(warn=-1) # to switch off possible warnings caused by old exprs object for yeast data 

# Data preprocessing 
yeast <- filter.NA(yeast)
yeast  <- fill.NA(yeast) # for illustration only; rather use knn method
yeast <- standardise(yeast)
# 
T.yeast <- 85   # cell cycle period (t=85min)
times.yeast <-  pData(yeast)$time  # time of measurements
#
yeast.test <- yeast[1:600,] # To speed up the example
#

NN <- 20 # number of generated background models
         # Here, a rather small number was chosen for demonstration purpose.
 
# Calculation of FDRs
# i) based on random permutation as background model
fdr.rr <- fdr.cycle(eset=yeast.test,T=T.yeast,
                    times=times.yeast,background.model="rr",N=NN)
# ii) based on Gaussian distribution 
fdr.g <- fdr.cycle(eset=yeast.test,T=T.yeast,
                   times=times.yeast,background.model="gauss",N=NN)
# iii) based on AR(1) models as background
fdr.ar1 <- fdr.cycle(eset=yeast.test,T=T.yeast,
                     times=times.yeast,background.model="ar1",N=NN)


# Number of significant genes based on diff. background models
sum(fdr.rr$fdr < 0.1) 
sum(fdr.g$fdr < 0.1)
sum(fdr.ar1$fdr < 0.1)

# Plot top scoring gene
plot(times.yeast,exprs(yeast.test)[order(fdr.ar1$fdr)[1],],type="o",
     xlab="Time",ylab="Expression",
     main=paste(geneNames(yeast.test)[order(fdr.ar1$fdr)[1]],"-- FDR:",
       fdr.ar1$fdr[order(fdr.ar1$fdr)[1]]))

# List significant genes
fdr.ar1$fdr[which(fdr.ar1$fdr < 0.1)]

options(ow) # returning to the original warning options
