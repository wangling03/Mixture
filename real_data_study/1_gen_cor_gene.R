rm(list=ls())

library(HiClimR)


data_chb=as.matrix(read.delim("/datafile/CHB_unrelated_norm_march2007.txt", row.names = 1, header= TRUE))
data_jpt=as.matrix(read.delim("/datafile/JPT_unrelated_norm_march2007.txt", row.names = 1, header= TRUE))

data_chb=as.data.frame(data_chb)
data_jpt=as.data.frame(data_jpt)

data_asian = merge(data_chb, data_jpt, by="row.names", all=TRUE)

rownames(data_asian)=data_asian$Row.names
data_asian=data_asian[2:length(data_asian)]

## remove CHRNA6 (i.e. GI_21361147-S)

row.names.remove = c("GI_21361147-S")

data_chrna6 = data_asian[(row.names(data_asian) %in% row.names.remove), ]
data_asian2 = data_asian[!(row.names(data_asian) %in% row.names.remove), ]

## Use fastCor function to compute the correlation matrix
t0 <- proc.time() 
genecor <- fastCor(t(data_asian2),nSplit = 20, upperTri = TRUE, optBLAS = TRUE, verbose = TRUE)
proc.time() - t0


####

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

genecor=makeSymm(genecor)


################

gene_hc_count06=apply(genecor, 2, function(x) length(x[x>0.6]))

write.table(gene_hc_count06, "/datafile/gene_hc_count06.txt", sep="\t")

