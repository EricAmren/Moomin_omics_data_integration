library(limma)
library(gplots)
library(dplyr)
setwd("Documents/Moomin_omics_data_integration/microarray_VS_RNAseq/Microarray/")

source("cyber-t_R_scripts/bayesreg.R")

## Reading raw data

#Raw data consist of 6 TSV files : 2 conditions and 3 replicates each.

# Conditions are :
  
# * Leaf lysate
# * Minimal media

# Limma functions are used here to extract relevant field.
setwd("microarray_data/raw_data/")
summaryTable<-read.delim("summaryTable.txt",check.names=FALSE,stringsAsFactors=FALSE)
microarrayData <-read.maimages(summaryTable[,"FileName"],source="agilent",green.only=TRUE)

## Probes filtering : keeping only probes that are 10% brighter than negative probes.

neg95 <- apply(microarrayData$E[microarrayData$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95)) #95 percentile of negative probes
cutoff <- matrix(1.1*neg95,nrow(microarrayData),ncol(microarrayData),byrow=TRUE) # keep probes that are 10% brighter than the negative probes
filtered <- rowSums(microarrayData$E > cutoff) >= 3
noNegData <- microarrayData[microarrayData$genes$ControlType==0 & filtered,]

## Data processing : rearranging data in a way cyber-t can use.

setwd("../../..")
df <- as.data.frame(noNegData$E)
df <- cbind(noNegData$genes$SystematicName,df)
df <- aggregate(df[,2:7], list(df$`noNegData$genes$SystematicName`), median) # Merging duplicates and keeping median values
colnames(df) <- c("gene_id", "C1", "C2", "C3", "E1", "E2", "E3" )
df <- df[c("gene_id", "E1", "E2", "E3", "C1", "C2", "C3")] # Rearranging df so our values are in the same order than RNAseq data
df <- format(df,  scientific = TRUE)
write.table(df, file="Microarray/microarray_data/processed_data/processed_MA_results.txt", col.names=FALSE, row.names=FALSE, sep =",", quote = FALSE)

## Call of Cyber-t "bayesT" function

# Cyber-t is a method that compute a PPDE and LogFC for probes that got change in intensity (AKA differentially expressed genes).
# Due to the large number of gene studied, we use 'winSize = 101' following cyber-t manual. 'bayes' parameter is set to '1' as we want to use a bayesian method.

dataFile<-read.table("Microarray/microarray_data/processed_data/processed_MA_results.txt",sep=",", header = TRUE, row.names = 1)
cyber_t_results <- bayesT(dataFile,numC=3,numE=3,ppde=1, bayes=1, winSize=101, conf=4)
logFC <- log(cyber_t_results$meanE / cyber_t_results$meanC)

## Building reference table for later (a table that matches systematic gene names and probe names.)
# Re-using a table from older study to associate gene names with probe names. (gene names <-> probe names <-> ECs name)
blast <- read.table("misc/SakaiRNASeq/final_blast_Table.txt", sep='\t', header=TRUE)
corresponding_names <- merge(dplyr::select(blast, query.name, subject), dplyr::select(noNegData$genes, ProbeName, SystematicName), by.x="query.name", by.y="ProbeName")
corresponding_names <- unique(corresponding_names)
write.table(corresponding_names, file="Comparison/reference_table.txt")

## Writing cyber-t results
# Output genes matched with their logFC and PPDE.

cyber_t_output <- data.frame(rownames(cyber_t_results), cyber_t_results$ppde.p, logFC)
cyber_t_output <- cyber_t_output[order(cyber_t_output$cyber_t_results.ppde.p, decreasing=TRUE), ]
write.table(cyber_t_output, file="Microarray/results/cyber_t_results.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

colnames(cyber_t_output) <- c("SystematicName", "PPDE", "logFC")
references <- read.table("Comparison/reference_table.txt")
cyber_t_output2 <- merge(cyber_t_output, references, by.x="SystematicName", by.y="SystematicName")

cyber_t_output3 <- dplyr::select(cyber_t_output2, subject, PPDE, logFC)
cyber_t_output3 <- aggregate(cyber_t_output3[,2:3], list(cyber_t_output3$subject), median)
colnames(cyber_t_output3) <- c("subject", "PPDE", "logFC")
write.table(cyber_t_output3, file="Microarray/results/cyber_t_final_results.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Differential expression analysis with RNAseq data and EBseq
library(memisc)

## Import Cyber-t results for microarray

cyber_t_results <- read.table("Microarray/results/cyber_t_final_results.txt")
colnames(cyber_t_results) <- c("gene_id", "cyber_t_PPDE", "cyber_t_logFC")
head(cyber_t_results)
cyber_t_results <- format(cyber_t_results, scientific = FALSE)

## Repair name discrepancies between RNAseq data and microarray data
references <- read.table("Comparison/reference_table.txt")
EBseq_corresponding_genes <- merge(small_EBseq_results, references, by.x="gene_id", by.y="subject")
colnames(EBseq_corresponding_genes) <- c("ECs","EBseq_PPDE","EBseq_logFC","query.name", "SystematicName")

## Merging all results in one big df
res_summary <- merge(cyber_t_results, small_EBseq_results, by="gene_id")
res_summary <- merge(res_summary, EdgeR_corresponding_genes, by="gene_id") 
res_summary <- merge(res_summary, small_limma_results, by="gene_id")
res_summary <- res_summary[c("gene_id", "cyber_t_logFC", "EBseq_logFC", "limma_logFC", "EdgeR_logFC", "cyber_t_PPDE", "EBseq_PPDE", "limma_Pval", "EdgeR_Pval")]

bayes_threshold = 0.97
freq_threshold = 0.01

full_res_summary <- merge(cyber_t_results, small_EBseq_results, by="gene_id", all=TRUE)
full_res_summary <- merge(full_res_summary, EdgeR_corresponding_genes, by="gene_id", all=TRUE) 
full_res_summary <- merge(full_res_summary, small_limma_results, by="gene_id", all=TRUE)
full_res_summary <- full_res_summary[c("gene_id", "cyber_t_logFC", "EBseq_logFC", "limma_logFC", "EdgeR_logFC", "cyber_t_PPDE", "EBseq_PPDE", "limma_Pval", "EdgeR_Pval")]

full_only_significant_values<-full_res_summary[(full_res_summary$cyber_t_PPDE > bayes_threshold | full_res_summary$EBseq_PPDE > bayes_threshold | full_res_summary$limma_Pval < freq_threshold | full_res_summary$EdgeR_Pval < freq_threshold),]
full_only_significant_values <- full_only_significant_values[!(rowSums(is.na(full_only_significant_values))==NCOL(full_only_significant_values)),]

## Filtering out unsignificant values
only_significant_values<-res_summary[(res_summary$cyber_t_PPDE > bayes_threshold & res_summary$EBseq_PPDE > bayes_threshold & res_summary$limma_Pval < freq_threshold & res_summary$EdgeR_Pval < freq_threshold),]




# Format and export data for Moomin input

cyber_t_exp_mat <- data.frame(only_significant_values$gene_id, only_significant_values$cyber_t_PPDE, only_significant_values$cyber_t_logFC)
colnames(cyber_t_exp_mat) <- c("GeneID", "PPDE", "FC")
write.table(cyber_t_exp_mat, file = "Moomin_input/cyber_t_res.tsv", quote=F, sep="\t", row.names = F)

