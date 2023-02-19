rm(list=ls())

library(dplyr)
library(HiClimR)


corcnt06_imp=as.matrix(read.delim("./datafile/gene_hc_count06.txt", row.names = 1, header= TRUE))


hist(corcnt06_imp)
quantile(corcnt06_imp, probs = seq(0, 1, 0.25))


corcnt_keep=subset(corcnt06_imp,  corcnt06_imp > 2200)

num_gene=length(corcnt_keep)

print(num_gene)

genename_keep=row.names(corcnt_keep)


data_chb=as.matrix(read.delim("datafile/CHB_unrelated_norm_march2007.txt", row.names = 1, header= TRUE))
data_jpt=as.matrix(read.delim("datafile/JPT_unrelated_norm_march2007.txt", row.names = 1, header= TRUE))


data_chb=as.data.frame(data_chb)
data_jpt=as.data.frame(data_jpt)

data_asian = merge(data_chb, data_jpt, by="row.names", all=TRUE)

rownames(data_asian)=data_asian$Row.names
data_asian=data_asian[2:length(data_asian)]

## remove CHRNA6 (i.e. GI_21361147-S)

row.names.remove = c("GI_21361147-S")

data_chrna6 = data_asian[(row.names(data_asian) %in% row.names.remove), ]
data_asian2 = data_asian[!(row.names(data_asian) %in% row.names.remove), ]

data_hc=data_asian2[(row.names(data_asian2) %in% genename_keep),]

data_hc_all = dplyr::bind_rows(data_chrna6, data_hc) 

write.table(data_hc_all, "./datafile/data_hc_06.txt", sep="\t")
write.csv(data_hc_all, file = "./datafile/data_hc_06.csv")


############################################################
## Use fastCor function to compute the correlation matrix
t0 <- proc.time() 
genecor <- fastCor(t(data_hc),nSplit = 1, upperTri = TRUE, optBLAS = TRUE, verbose = TRUE)
proc.time() - t0

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

genecor=makeSymm(genecor)

diagv= rep(1,num_gene)

diag(genecor) <- diagv

setEPS()


eps_graph_name="./datafile/corcnt06.eps"
pdf_graph_name="./datafile/corcnt06.pdf"
png_graph_name="./datafile/corcnt06.png"

#print(eps_graph_name)

#postscript(eps_graph_name)
#pdf(file=pdf_graph_name)
png(png_graph_name)

Sigma_m=abs(as.matrix(genecor))
filled.contour(x=seq(from=1,to=num_gene,length=num_gene),
               y=seq(from=1,to=num_gene,length=num_gene),
               z=Sigma_m,
               zlim=seq(0, 1),
               nlevels=15,
               axes=TRUE,
               #col=hsv(h=seq(from=-1,to=1,length=1000),s=1,v=1),
               color.palette=topo.colors,
               key.axes = axis(4, cex.axis=1.8),
               plot.axes={
                 axis(1,cex.axis=1.5)
                 axis(2,cex.axis=1.5)
                 
               },

)

dev.off()

############################################################

