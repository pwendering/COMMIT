#!/usr/bin/Rscript

# CSS normalization of OTU abundances
# adapted from Ruben Garrido-Oter

options(warn=-1)

library(metagenomeSeq, quietly = T, warn.conflicts = F)
library(biomformat, quietly = T, warn.conflicts = F)

args <- commandArgs(TRUE)

inFile <- args[1]
outFile <- args[2]

# read OTU table form file and convert to MRexperiment and biomformat dgCMatrix
b <- read_biom(inFile)
b_exp <- biom2MRexperiment(b)
b <- biom_data(b)

# for the calculation of the scaling factors,
# columns with none or only one entry have to be removed form the matrix
b <- matrix(b, ncol = ncol(b), nrow = nrow(b))
singletons = apply(X = b, MARGIN = 2, FUN = function(x) {sum(x!=0)<=1})
b <- b[,!singletons]
p <- cumNormStatFast(b)

# normalize the OTU table
b <- cumNorm(b_exp, p=p)

# convert the result back to biom format and write the file
write_biom(MRexperiment2biom(b, norm=T, log=T), outFile)

