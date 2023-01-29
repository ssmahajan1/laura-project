## R script for tximport of STAR/salmon counts, consensus clustering
# Jeremiah Holt and Sidharth Mahajan
# For Laura's BXD project
library(tidyverse)
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(EnhancedVolcano)
library(circlize)
library(ComplexHeatmap)
library(ConsensusClusterPlus)

source("C:/Users/ssmah/OneDrive - University of Tennessee/Sid/laura-project/R_functions.R")
#C:\Users\ssmah\OneDrive - University of Tennessee\Sid\laura-project
setwd("C:/Users/ssmah/OneDrive - University of Tennessee/Laura/Sid/laura-project")

# get file names and extract samples from each
files = list.files(path="./", pattern="*\\.sf", recursive=T)
samples = gsub("_.*","", files)
samples = gsub(".*/", "", samples)
names(files) = samples
length(files) #107

# extract condition from sample-conditions file
conditionTable = read_tsv("C:/Users/ssmah/OneDrive - University of Tennessee/Medical School/Research/Genomics_Project/Laura/laura_sample_conditions_n107.txt") %>%
  #filter(remove.sample != "remove") %>%
  arrange(SampleID) %>%
  mutate(SampleID.updated = if_else(is.na(SampleID.new), SampleID, SampleID.new)) %>%
  mutate(SampleID.updated = if_else(is.na(SampleID.new.JH), SampleID.updated, SampleID.new.JH))

dim(conditionTable) #107 7

# missing samples
all(conditionTable$SampleID == names(files)) #TRUE

#missing1 = conditionTable %>% filter(!SampleID %in% names(files))
#missing2 = files[which(!names(files) %in% conditionTable$SampleID)]

# filter conditionTable 
#conditionTable = conditionTable %>% 
#  filter(SampleID %in% names(files)) %>%
#  arrange(SampleID)
#dim(conditionTable) #104 7

#all(conditionTable$SampleID == names(files)) #TRUE

# change sample name on file list
names(files) = conditionTable$SampleID.updated

# remove samples to filter
conditionTable = conditionTable %>% 
  filter(is.na(remove.sample)) %>% # removes specified samples 
  filter(tissue == "Tumor") %>%
  arrange(SampleID.updated)
dim(conditionTable) #95 7 

files = files[conditionTable$SampleID.updated]
length(files) #95

all(conditionTable$SampleID.updated == names(files)) #TRUE

# make colData table for DESeq2
sampleTable = data.frame(sample=names(files),
                         #sample.updated = conditionTable$SampleID.updated,
                         #remove = conditionTable$remove.sample,
                         file=files,
                         tissue=factor(conditionTable$tissue),
                         group=factor(conditionTable$group))

table(sampleTable$tissue, sampleTable$group)

#####################################
### Make TxDb object for tximport ###
#####################################

# TxDb object
txdb = makeTxDbFromGFF("C:/Users/ssmah/OneDrive - University of Tennessee/Medical School/Research/Genomics_Project/Laura/ensembl_mm10.gtf", dataSource = "ensembl?", organism = "Mus musculus")
k = AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
#tx2gene = tx2gene %>% mutate(TXNAME = gsub("\\..*", "", TXNAME),
#                             GENEID = gsub("\\..*", "", GENEID)) # remove version suffix from transcripts and genes

#################################
### Read files using tximport ###
#################################

all(file.exists(files))
txi = tximport(files, 
               type="salmon", 
               tx2gene=tx2gene,
               ignoreTxVersion=T)
names(txi)

########################
### Plug into DESeq2 ###
########################

ddsTxi = DESeqDataSetFromTximport(txi,
                                  colData = sampleTable,
                                  design = ~group) 
dim(ddsTxi) #47606    95
head(counts(ddsTxi))

# pre filter for genes that dont have at least 10 total counts
keep = rowSums(counts(ddsTxi)) >= 10
dds = ddsTxi[keep,]
dim(dds) #28582    95

# more prefiltering - protein coding? 
# read in protein coding genes
IDsToKeep = read_tsv("C:/Users/ssmah/OneDrive - University of Tennessee/Medical School/Research/Genomics_Project/Laura/protein-coding-genes-laura-BXD.txt") %>%
  filter(ensembl_gene_id %in% rownames(dds))
dim(IDsToKeep) #18523     2

# filter genes and remove controls 
dds = dds[IDsToKeep$ensembl_gene_id,]
dim(dds) #18523    95

# check meta data
colData(dds)

# run DESeq
dds = DESeq(dds) 

# write counts
norm.df = as.data.frame(counts(dds, normalized=T)) %>% 
  mutate(symbol = toupper(mapIds(org.Mm.eg.db,
                                 keys=row.names(counts(dds, normalized=T)),
                                 column="SYMBOL",
                                 keytype="ENSEMBL",
                                 multiVals="first")), 
         .before=1) %>%
  mutate(ID = row.names(counts(dds, normalized=T)), .before=1)
dim(norm.df) #18523    97

write_tsv(norm.df, "salmon_mm10_counts_library_size_normalized_BXD_project_DESeq2_protein_coding_18523_n97.txt")

# write pseudocounts
as.data.frame(log2(counts(dds, normalized=T)+1)) %>% 
  mutate(symbol = toupper(mapIds(org.Mm.eg.db,
                                 keys=row.names(counts(dds, normalized=T)),
                                 column="SYMBOL",
                                 keytype="ENSEMBL",
                                 multiVals="first")), 
         .before=1) %>%
  mutate(ID = row.names(counts(dds, normalized=T)), .before=1) %>%
  write_tsv("salmon_mm10_pseudocounts_library_size_normalized_BXD_project_DESeq2_protein_coding_18523_n97.txt")

################################################################################
# normalize for downstream analysis, e.g. clustering (not DEG)                 #
# use VST over rlog for sample n > 30                                          #
################################################################################
vsd = vst(dds, blind=T) # Variance stabilizing transformation (VST)  
head(assay(vsd)) #average ~ (+6)
res = as.data.frame(assay(vsd))
dim(res) #18523    95

res %>% 
  mutate(GeneSymbol = toupper(mapIds(org.Mm.eg.db,
                                     keys=row.names(counts(dds, normalized=T)),
                                     column="SYMBOL",
                                     keytype="ENSEMBL",
                                     multiVals="first")), .before=1) %>%
  mutate(GeneID = rownames(res), .before=1) %>%
  write_tsv("salmon_mm10_counts_VST_normalized_log2_BXD_project_DESeq2_protein_coding_18523_n95.txt")


# plot and compare
boxplot(log2(counts(dds, normalized=F)+1), xlab ="Sample", ylab="Raw pseudocounts (log2+1)")
boxplot(res[,-c(1,2)], xlab ="Sample", ylab="Log2 normalized counts (VST)") 
boxplot(log2(counts(dds, normalized=T)+1), xlab ="Sample", ylab="Library sized normalized pseudocounts (log2+1)")

# call plotting functions from 2_R_functions.R
p = plotDistribution(dds, 
                     plot_type="boxplot",
                     norm=F)
print(p + ggtitle("Psuedocounts: log2(TPM+1)"))

p = plotDistribution(dds, 
                     plot_type="boxplot",
                     norm=T)
print(p + ggtitle("Library Size Normalized Psuedocounts: log2(TPM+1)"))

# density
p = plotDistribution(dds, 
                     plot_type="density",
                     norm=F)
print(p + ggtitle("Psuedocounts: log2(TPM+1)"))

p = plotDistribution(dds, 
                     plot_type="density",
                     norm=T)
print(p + ggtitle("Library Size Normalized Psuedocounts: log2(TPM+1)"))

# normalized
p = plotDistributionNorm(vsd, plot_type="boxplot")
print(p + ggtitle("VST Normalized Expression"))

p = plotDistributionNorm(vsd, plot_type="density")
print(p + ggtitle("VST Normalized Expression"))


##############################
### PCA plots and clusters ###
##############################

# built in DESeq2 functions

plotPCA(vsd, intgroup=c("group"), ntop=1000) + geom_text(aes(label=name),vjust=2)

CA(vsd, intgroup=c("tissue"), ntop=1000) + geom_text(aes(label=name),vjust=2)

########################################
### Get quick heatmaps for MAD genes ###
########################################

#simple heatmap function
mad.heatmap = function(matrix, sampleTable, name){
  
  # filtered MAD
  mad = apply(matrix,1,mad)
  mad.cutoff= quantile(mad, 0.75) 
  mad.subset = matrix[which(mad>mad.cutoff),]
  dim(mad.subset) #4427   25
  
  mad.matrix.scaled = t(scale(t(data.matrix(mad.subset))))
  rownames(mad.matrix.scaled) = rownames(mad.subset)
  head(mad.matrix.scaled)
  dim(mad.matrix.scaled) #4427   25
  
  
  #heatmap color function
  col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
  
  all(sampleTable$sample == colnames(mad.matrix.scaled))
  
  # heatmaps
  ha = HeatmapAnnotation(
    group = sampleTable$group,
    tissue = sampleTable$tissue,
    simple_anno_size = unit(0.75, "cm")
  )
  
  hm = Heatmap(mad.matrix.scaled,
               column_title = paste("All samples n=", ncol(mad.matrix.scaled), sep=""),
               row_title = paste("Top 25% most variable genes n=", nrow(mad.matrix.scaled), sep=""),
               show_row_names = F,
               cluster_rows = T,
               row_km = 3, 
               row_km_repeats = 3,
               row_gap = unit(1, "mm"),
               # row_split = dif$UP.DOWN,
               #column_split=sampleTable$group,
               cluster_columns =T,
               column_km = 3,
               column_km_repeats = 100,
               column_gap = unit(1, "mm"),
               top_annotation = ha
  )
  
  
  pdf(paste(name,".pdf", sep=""), width = 14, height = 14)
  draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
}
# convert to matrix
results.matrix = data.matrix(res)
rownames(results.matrix) = rownames(res)
dim(results.matrix)
#############
###filtering
############
mad = apply(results.matrix,1,mad)
mad.cutoff= quantile(mad, 0.75) 
results.matrix = results.matrix[which(mad>mad.cutoff),]
dim(results.matrix) #4427   25



results.matrix.centered = results.matrix - apply(results.matrix,1, median)
#results.matrix = t(scale(t(results.matrix)))
head(results.matrix)
dim(results.matrix) #18620    98

# call function
mad.heatmap(results.matrix, sampleTable, "BXD_top_25_MAD_heatmap_cluster_both_protein_coding_18523_n95_km3")

##################################
###### ConsensusClusterPlus ######
##################################

ccp.res <- ConsensusClusterPlus(results.matrix.centered,maxK=9,
                                reps=1000,
                                pItem=0.8,
                                pFeature=1,
                                clusterAlg="hc", 
                                title="cluster_test_centered.1000.filter",
                                #innerLinkage="average", 
                                #finalLinkage="average", 
                                distance= "pearson",
                                #ml=NULL,
                                #tmyPal=NULL,
                                seed = 12345,
                                plot="png",
                                #weightsItem=NULL,
                                #weightsFeature=NULL,
                                #verbose=F
)
ccp.res[[4]]$consensusClass
icl.test <- calcICL(ccp.res,title="cluster_test", plot=NULL)


##################
###### DEGs ######
##################

plotCounts(dds, "ENSMUSG00000059552", intgroup=c("tissue"))
plotCounts(dds, "ENSMUSG00000022105", intgroup=c("tissue"))

############################
library(ensembldb)

# get Biotype information for each gene using biomaRt
sensembl = useMart("ensembl",
                   dataset="mmusculus_gene_ensembl",
                   host="useast.ensembl.org")

goids.mm10 = getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), 
                   filters = 'ensembl_gene_id', 
                   values = rownames(dds), 
                   mart = sensembl)

table(goids.mm10$gene_biotype)
detach("package:ensembldb")

IDsToKeep = goids.mm10 %>%
  dplyr::filter(gene_biotype == "protein_coding") # keep protein coding 
dim(IDsToKeep) #18620     2

write_tsv(IDsToKeep,"protein-coding-genes-laura-BXD.txt")

####################################




