library(genomation)
library(stringr)
library(DESeq2)
library(biomaRt)
library(variancePartition)


# Metadata file needed
meta <- read.csv("counts.meta.csv",header=T,row.names = 1)

# Filtering steps for grabbing GLU celltype
working.meta <- meta[meta$celltype=="GLU",]

working.meta$SOURCE <- ifelse(grepl(pattern="S",rownames(working.meta)),"Bank","Hurd")

colnames(working.meta)[2:3] <- c("BATCH","OPIOID")

# Raw Counts file needed -- columns need to be handled! See #Clean Up Column Names Below
feat.counts <- read.table("Peak.counts.txt",header=T,row.names=1)

load("ensembl.RData") # load in ensembl mart object



## hg19 GTF to hg38 BED12
#ml ucsc-utils
#awk '{if($3 != "gene") print $0}' hg19.gtf | grep -v "^#" | gtfToGenePred /dev/stdin /dev/stdout | genePredToBed stdin hg19.bed



# Clean up column names 
colnames(feat.counts) <- str_replace_all(colnames(feat.counts),c("X|.GLU.bam|GLU.|.bam"),c(""))
colnames(feat.counts)[6:length(colnames(feat.counts))] <- paste("GLU-",colnames(feat.counts)[6:length(colnames(feat.counts))],sep="")

gene.annotation <- readTranscriptFeatures("hg19.bed", remove.unusual = TRUE, up.flank = 1000,
  down.flank = 1000, unique.prom = TRUE)
  
feat.counts$Chr <- str_replace_all(feat.counts$Chr,"chr","")

anno <- annotateWithGeneParts(as(feat.counts[,c(1:3)],"GRanges"),gene.annotation)
gene.part <- getMembers(anno)
tss.dis <-  getAssociationWithTSS(anno) 
feat.annotated <- cbind(feat.counts,gene.part,tss.dis,rownames(feat.counts))

values <- unique(feat.annotated$feature.name)
data <- getBM(attributes=c("hgnc_symbol","description","gene_biotype","ensembl_transcript_id"), filters = "ensembl_transcript_id", values = values, mart= ensembl)
feat.annotated <- merge(feat.annotated,data,by.x="feature.name",by.y="ensembl_transcript_id",all.x=T)



# Filtered and Normalized counts
feat.counts <- feat.counts[rowMeans(feat.counts[,c(6:length(colnames(feat.counts)))])>10,]
filt.count <- feat.counts[,6:length(colnames(feat.counts))]
normed.counts <- varianceStabilizingTransformation(as.matrix(filt.count))
dir.create("Counts")
write.table(feat.annotated,file="Counts/Peak.counts.annnotated.txt",col.names=T,row.names=T)
write.table(filt.count,file="Counts/Filtered.counts.txt",col.names=T,row.names=T)
write.table(normed.counts,file="Counts/Normed.counts.txt",col.names=T,row.names=T)

# PCA calculations
suppressPackageStartupMessages(library(ggfortify))
PCA <- prcomp(t(normed.counts))
importance <- summary(PCA)
sample.PCAs <- PCA$x
dir.create("PCA.output")
write.table(importance$importance,file="PCA.output/PCA.importance.txt",col.names=T,row.names=T)
write.table(sample.PCAs,file="PCA.output/PCA.each.sample.txt",col.names=T,row.names=T) 

# Grafting metadata with outcomes
  #PCA
Glu.Meta.Stats <- write.table(cbind(working.meta,sample.PCAs[,1:2]),file="PCA.output/Shiny_PCA_1and2.txt",col.names=T,row.names=T)
  #Peak number
                      

# Differential expression
#filt.count <- filt.count[,rownames(working.meta)]
#colnames(filt.count)==rownames(working.meta)
#matrix <- DESeqDataSetFromMatrix(filt.count,working.meta,design=~SEX+AGE+BATCH+OPIOID)
#analysis <- DESeq(matrix)
#results <- results(analysis)
#dir.create("Diff.Expression")
#write.csv(results,file="DEtable.csv",header=T,row.names=1)

# Variance Partition
# norm.counts <- norm.counts[,rownames(working.meta)
# working.meta$batch <- as.factor(working.meta$batch)                           
# form <- ~ AGE + (1|BATCH) + (1|OPIOID) + (1|SEX)
# varPart <- fitExtractVarPartModel(norm.counts,form,working.meta )




