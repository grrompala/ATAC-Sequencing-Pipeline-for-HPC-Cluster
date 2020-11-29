library(genomation)
library(stringr)
library(DESeq2)
library(biomaRt)
library(variancePartition)
library(doParallel)

setwd("/sc/arion/projects/DADisorders/greg.working/ATAC.pipeline")

# Metadata file needed
meta <- read.csv("counts.meta.csv",header=T,row.names = 1)
# Filtering steps for grabbing GLU celltype
working.meta <- meta[meta$celltype=="GLU",]
rownames(working.meta) <- str_replace_all(rownames(working.meta),"-",".")

filt.count <- read.table("Counts/Filtered.counts.txt",header=T)
norm.counts <- read.table("Counts/Normed.counts.txt",header=T)

# Differential expression
filt.count <- filt.count[,rownames(working.meta)]
colnames(filt.count)==rownames(working.meta)
matrix <- DESeqDataSetFromMatrix(filt.count,working.meta,design=~SEX+AGE+batch+GROUP)
analysis <- DESeq(matrix)
results <- results(analysis)
# Annotate the results
anno.file <- read.table("Counts/Peak.counts.annnotated.txt",header=T)
anno <- anno.file[,c(1,67:76)]
anno <- anno[!duplicated(anno$rownames.feat.counts.),]
rownames(anno) <- anno$rownames.feat.counts.
anno <- anno[rownames(results),]
final <- data.frame(cbind(results,anno))
dir.create("Diff.Expression")
write.table(final,file="Diff.Expression/DEtable.txt",col.names =TRUE,row.names=TRUE)

#Variance Partition
norm.counts <- norm.counts[,rownames(working.meta)]
working.meta$batch <- as.factor(working.meta$batch)                           
form <- ~ AGE + (1|batch) + (1|GROUP) + (1|SEX)
cl <- makeCluster(3)
registerDoParallel(cl)
varPart <- fitExtractVarPartModel(norm.counts,form,working.meta )
dir.create("Variance.Partition")
vp <- sortCols(varPart)
write.table(vp,file="Variance.Partition/varPart.txt")

