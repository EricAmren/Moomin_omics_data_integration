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
summaryTable<-read.delim("summaryTable.txt",check.names=FALSE,stringsAsFactors=FALSE) # Read data of files listed in the summaryTable.txt
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

## Import Cyber-t results for microarray

cyber_t_results <- read.table("Microarray/results/cyber_t_final_results.txt")
colnames(cyber_t_results) <- c("GeneID", "PPDE", "FC")
cyber_t_results <- format(cyber_t_results, scientific = FALSE)

## Filtering out genes with ppde below a threshold of 97%
bayes_threshold = 0.97
filtered_cyber_t_results <- cyber_t_results[(cyber_t_results$PPDE > bayes_threshold),]

# Format and export data for Moomin input

write.table(filtered_cyber_t_results, file = "Moomin_input/cyber_t_res2.tsv", quote=F, sep="\t", row.names = F)
