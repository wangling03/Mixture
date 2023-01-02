rm(list=ls())

library(HiClimR)


#data_chb=as.matrix(read.delim("E:\\document related\\Ling MSU\\papers\\Genetic data\\HapMap\\CHB_unrelated_norm_march2007\\CHB_unrelated_norm_march2007.txt", row.names = 1, header= TRUE))
#data_jpt=as.matrix(read.delim("E:\\document related\\Ling MSU\\papers\\Genetic data\\HapMap\\JPT_unrelated_norm_march2007\\JPT_unrelated_norm_march2007.txt", row.names = 1, header= TRUE))

data_chb=as.matrix(read.delim("../real_data_study2/datafile/CHB_unrelated_norm_march2007.txt", row.names = 1, header= TRUE))
data_jpt=as.matrix(read.delim("../real_data_study2/datafile/JPT_unrelated_norm_march2007.txt", row.names = 1, header= TRUE))

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

genecor[1:5,1:4]

####

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

genecor=makeSymm(genecor)

#write.csv(genecor, "/mnt/scratch/wangli35/genecor.csv")

################

gene_hc_count05=apply(genecor, 2, function(x) length(x[x>0.5]))

write.table(gene_hc_count05, "./gene_hc_count05.txt", sep="\t")

################

gene_hc_count06=apply(genecor, 2, function(x) length(x[x>0.6]))

write.table(gene_hc_count06, "./gene_hc_count06.txt", sep="\t")

################

gene_hc_count07=apply(genecor, 2, function(x) length(x[x>0.7]))

write.table(gene_hc_count07, "./gene_hc_count07.txt", sep="\t")

################

gene_hc_count08=apply(genecor, 2, function(x) length(x[x>0.8]))

write.table(gene_hc_count08, "./gene_hc_count08.txt", sep="\t")

################

gene_hc_count09=apply(genecor, 2, function(x) length(x[x>0.9]))

write.table(gene_hc_count09, "./gene_hc_count09.txt", sep="\t")

################

gene_hc_count2=apply(genecor, 2, function(x) length(x[x>=0.5]))

write.table(gene_hc_count2, "./gene_hc_count2.txt", sep="\t")

#######
length(gene_hc_count05[gene_hc_count05>2200])
##6629
length(gene_hc_count05[gene_hc_count05>=2200])
##6634

length(gene_hc_count06[gene_hc_count06>2200])
#2269
length(gene_hc_count06[gene_hc_count06>=2200])
#2270

length(gene_hc_count07[gene_hc_count07>2200])
#8
length(gene_hc_count07[gene_hc_count07>=2200])
#8