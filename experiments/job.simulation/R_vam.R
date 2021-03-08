#!/usr/bin/env Rscript

start_time =proc.time()

library(VAM)
library(Seurat)

# Load parameters 
args = commandArgs(trailingOnly=TRUE)
SC_FILE=args[1]
GS_FILE=args[2]
OUT_FILE=args[3]

print(paste0('SC_FILE: ', SC_FILE))
print(paste0('GS_FILE: ', GS_FILE))
print(paste0('OUT_FILE: ', OUT_FILE))

# Load single-cell data
raw_counts<-read.table(file=SC_FILE, row.names = 1, header=TRUE, sep='\t')
mydata <- CreateSeuratObject(counts = raw_counts, min.cells = 3, min.genes = 200, project = "mydata_scRNAseq")
mydata = NormalizeData(mydata)

# Load .gs data
df_gs = read.table(file=GS_FILE, header=TRUE, sep='\t')
gene.set.id.list = list()
for (i in 1:dim(df_gs)[1]){
    gene.set.id.list[[i]] = strsplit(df_gs[i,2], ",")[[1]]
    names(gene.set.id.list)[i] = df_gs[i,1]
}

# Run VAM and write results
mydata.vam = vamForSeurat(seurat.data=mydata,
                          gene.set.collection=gene.set.id.list,
                          center=F, gamma=T, sample.cov=F, return.dist=T)
df_pval = mydata.vam@assays$VAMcdf[,]
df_pval = data.frame(t(1-df_pval))
write.table(df_pval, file=OUT_FILE, sep='\t')

print('Finished')
print(proc.time() - start_time)