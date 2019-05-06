library(limma)
library(gplots)
library(dplyr)
setwd("/home/ericcumunel/Documents/internship/rnaseq_MA_comparison/SakaiRNASeq/1.Microarray_Data_Analysis")
source("bayesreg.R")

#read data in 
summaryTable<-read.delim("summaryTable.txt",check.names=FALSE,stringsAsFactors=FALSE)
microarrayData <-read.maimages(summaryTable[,"FileName"],source="agilent",green.only=TRUE)

#normalisation : looks like a bad idea ! TODO : Reverting afterward 
correctedMicroarrayData <-backgroundCorrect(microarrayData,method="normexp")
correctedMicroarrayData <-normalizeBetweenArrays(correctedMicroarrayData,method="quantile")
correctedMicroarrayData <- microarrayData

#filtering of probes
neg95 <- apply(correctedMicroarrayData$E[correctedMicroarrayData$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95)) #95 percentile of negative probes
cutoff <- matrix(1.1*neg95,nrow(correctedMicroarrayData),ncol(correctedMicroarrayData),byrow=TRUE) # keep probes that are 10% brighter than the negative probes
filtered <- rowSums(correctedMicroarrayData$E > cutoff) >= 3

#removing of negative probes
noNegData <- correctedMicroarrayData[correctedMicroarrayData$genes$ControlType==0 & filtered,]

#Normalized dataframe
df <- as.data.frame(noNegData$E)
df <- cbind(noNegData$genes$SystematicName,df)
normalized_df <- aggregate(df[,2:7], list(df$`noNegData$genes$SystematicName`), median)
colnames(normalized_df) <- c("gene_id", "C1", "C2", "C3", "E1", "E2", "E3" )
normalized_df <- normalized_df[c("gene_id", "E1", "E2", "E3", "C1", "C2", "C3")]
normalized_df <- format(normalized_df,  scientific = TRUE)
write.table(normalized_df, file="normalized_results.txt", col.names=FALSE, row.names=FALSE, sep =",", quote = FALSE)
dataFile<-read.table("normalized_results.txt",sep=",", header = TRUE, row.names = 1)
norm_results <- bayesT(dataFile,numC=3,numE=3,ppde=1, bayes=1, winSize=101, conf=4)
# norm_FC <- (norm_results$meanE / norm_results$meanC) - 1
logFC <- log(norm_results$meanE / norm_results$meanC)


cybert_norm_results <- data.frame(rownames(norm_results), norm_results$ppde.p, norm_results$fold)
# cybert_norm_results <- data.frame(rownames(norm_results), norm_results$ppde.p, norm_FC)
cybert_norm_results <- data.frame(rownames(norm_results), norm_results$ppde.p, logFC)
cybert_norm_results <- cybert_norm_results[order(cybert_norm_results$norm_results.ppde.p, decreasing=TRUE), ]
write.table(cybert_norm_results, file="cybert_norm_results.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Not Normalized dataframe
o_neg95 <- apply(microarrayData$E[microarrayData$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95)) #95 percentile of negative probes
o_cutoff <- matrix(1.1*o_neg95,nrow(microarrayData),ncol(microarrayData),byrow=TRUE) # keep probes that are 10% brighter than the negative probes
o_filtered <- rowSums(microarrayData$E > o_cutoff) >= 3

#removing of negative probes
o_noNegData <- microarrayData[microarrayData$genes$ControlType==0 & o_filtered,]

df2 <- as.data.frame(o_noNegData$E)
df2 <- cbind(o_noNegData$genes$SystematicName,df2)
not_normalized_df <- aggregate(df2[,2:7], list(df2$`o_noNegData$genes$SystematicName`), median)

write.table(not_normalized_df, file="results.txt", col.names=FALSE, row.names=FALSE, sep =",", quote = FALSE )
o_dataFile<-read.table("results.txt",sep=",", header = TRUE, row.names = 1)
no_norm_results <-bayesT(o_dataFile,numC=3,numE=3,ppde=1,bayes=1,winSize = 101)
cybert_no_norm_results <- data.frame(rownames(no_norm_results), no_norm_results$ppde.p, no_norm_results$fold)
cybert_no_norm_results <- cybert_no_norm_results[order(cybert_no_norm_results$no_norm_results.ppde.p, decreasing=TRUE), ]
write.table(cybert_no_norm_results, file="cybert_no_norm_results.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## gene names <-> probe names <-> ECs name
noNegData$genes
blast <- read.table("/home/ericcumunel/Documents/internship/rnaseq_MA_comparison/SakaiRNASeq/3.Microarray_Vs_RNASeq/3.table merge/final_blast_Table.txt", sep='\t', header=TRUE)
corresponding_names <- merge(dplyr::select(blast, query.name, subject), dplyr::select(noNegData$genes, ProbeName, SystematicName), by.x="query.name", by.y="ProbeName")
corresponding_names <- unique(corresponding_names)
write.table(corresponding_names, file="reference_table.txt")
