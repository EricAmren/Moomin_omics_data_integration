# setwd("/home/ericcumunel/Documents/internship/rnaseq_MA_comparison/SakaiRNASeq/1.Microarray_Data_Analysis")
# norm_results<-read.table("cyber-t_results/norm_results.txt",sep="\t", header = TRUE)
# norm_svn_results<-read.table("cyber-t_results/norm_svn_results.txt",sep="\t", header = TRUE)
# norm_log_results<-read.table("cyber-t_results/norm_log_results.txt",sep="\t", header = TRUE)
# 
# no_norm_results<-read.table("cyber-t_results/nonorm_results.txt",sep="\t", header = TRUE)
# no_norm_svn_results<-read.table("cyber-t_results/nonorm_svn_results.txt",sep="\t", header = TRUE)
# no_norm_log_results<-read.table("cyber-t_results/nonorm_log.results.txt",sep="\t", header = TRUE)
# 
# cybert_results<-read.table("/home/ericcumunel/Documents/internship/rnaseq_MA_comparison/RNAseq/top1000.tsv", sep="\t", header = TRUE)
# 
# limma_results<-read.table("microarrayTopTable.txt",sep='\t', header = TRUE)
# limma_genes <- limma_results$SystematicName
# ranking_df = data.frame(limma=limma_results$SystematicName[1:50], norm=norm_results$Lab_0[1:50], norm_svn=norm_svn_results$Lab_0[1:50], norm_log=norm_log_results$Lab_0[1:50])
# 
# length(intersect(limma_results$SystematicName[1:100],norm_results$Lab_0[1:100]))
# length(intersect(limma_results$SystematicName[1:100],norm_svn_results$Lab_0[1:100]))
# length(intersect(limma_results$SystematicName[1:100],norm_log_results$Lab_0[1:100]))
# 
# length(intersect(limma_results$SystematicName[1:100],no_norm_results$Lab_0[1:100]))
# length(intersect(limma_results$SystematicName[1:100],no_norm_svn_results$Lab_0[1:100]))
# length(intersect(limma_results$SystematicName[1:100],no_norm_log_results$Lab_0[1:100]))
# 


## Compare countings
setwd("/home/ericcumunel/Documents/internship/rnaseq_MA_comparison/comparison")

# Results from the paper for microarray (using limma)
limma_results<-read.table("microarrayTopTable.txt",sep='\t', header = TRUE)

# Results of cyber-t for microarray
cybert_results <- read.table("/home/ericcumunel/Documents/internship/rnaseq_MA_comparison/SakaiRNASeq/1.Microarray_Data_Analysis/cybert_norm_results.txt")
colnames(cybert_results) <- c("gene_id", "cybert_PPDE", "cybert_logFC")
cybert_results <- format(cybert_results, scientific = TRUE)

# Results of EBseq for RNAseq
EBseq_results <- read.table("/home/ericcumunel/Documents/internship/rnaseq_MA_comparison/RNAseq/final_EBseq_results.tsv", header=TRUE)
head(EBseq_results)
small_EBseq_results <- data.frame(gsub("ECs_", "ECs", EBseq_results$names), EBseq_results$PPDE, EBseq_results$logFC)

colnames(small_EBseq_results) <- c("gene_id", "PPDE", "logFC")
references <- read.table("reference_table.txt")
# summary <- data.frame(cybert_results)

EB_with_real_names <- merge(small_EBseq_results, references, by.x="gene_id", by.y="subject")
head(EB_with_real_names)
colnames(EB_with_real_names) <- c("ECs","EBseq_PPDE","EB_logFC","query.name", "SystematicName")

# Results from the paper for RNAseq (featureCounts)
rnaseq_results <- read.table("rnaseqTopTable.txt", skip=1)
colnames(rnaseq_results) <- c("ECs", "logFC", "t", "P.Value", "adj.P.Val", "b")
rnaseq_with_real_names <- merge(rnaseq_results, references, by.x="ECs", by.y="subject")
rnaseq_with_real_names <- data.frame(rnaseq_with_real_names$SystematicName, rnaseq_with_real_names$P.Value,rnaseq_with_real_names$logFC)
colnames(rnaseq_with_real_names) <- c("gene_id", "rna_Pval", "rna_logFC")

# Merging all results in one big df
ppde_summary <- merge(cybert_results, EB_with_real_names, by.x="gene_id", by.y="SystematicName")
ppde_summary <- merge(ppde_summary, rnaseq_with_real_names, by="gene_id")
head(ppde_summary)
# Filtering out missing values
ppde_summary[complete.cases(ppde_summary), ]


final_summary <- merge(ppde_summary, limma_results, by.x="gene_id", by.y="SystematicName")
final_summary <- final_summary[,c("gene_id","cybert_PPDE", "cybert_logFC", "EBseq_PPDE", "EB_logFC", "rna_Pval", "rna_logFC", "logFC", "P.Value")]
colnames(final_summary) <- c("gene_id","cybert_PPDE", "cybert_logFC", "EBseq_PPDE", "EB_logFC", "rna_Pval", "rna_logFC", "limma_logFC", "limma_Pval")
final_summary <- final_summary[c("gene_id","cybert_PPDE", "cybert_logFC", "EBseq_PPDE", "EB_logFC", "rna_Pval", "rna_logFC", "limma_Pval", "limma_logFC")]
final_summary <- final_summary[complete.cases(final_summary),]


## Statistical Tests :

# Limma Vs FeatureCounts
plot(final_summary$limma_logFC, final_summary$rna_logFC)
rho_limma_vs_FC <- cor.test(~ as.numeric(final_summary$limma_logFC) + as.numeric(final_summary$rna_logFC), data=final_summary, method = "spearman", continuity = TRUE)$estimate
r_squared_limma_vs_FC <- summary(lm(final_summary$limma_logFC ~ final_summary$rna_logFC, data=final_summary))$adj.r.squared


# EBseq Vs FeatureCounts
plot(final_summary$EB_logFC, final_summary$rna_logFC)
rho_EB_vs_FC <- cor.test(~ as.numeric(final_summary$EB_logFC) + as.numeric(final_summary$rna_logFC), data=final_summary, method = "spearman", continuity = TRUE)$estimate
r_squared_EB_vs_FC <- summary(lm(final_summary$EB_logFC ~ final_summary$rna_logFC, data=final_summary))$adj.r.squared

# Cyber-t Vs Limma
plot(final_summary$cybert_logFC, final_summary$limma_logFC)
rho_cybert_VS_limma <- cor.test(~ as.numeric(final_summary$cybert_logFC) + as.numeric(final_summary$limma_logFC), data=final_summary, method = "spearman", continuity = TRUE)$estimate
r_squared_cybert_vs_limma <- summary(lm(final_summary$cybert_logFC ~ final_summary$limma_logFC, data=final_summary))$adj.r.squared

# Cyber-t Vs EBseq
plot(final_summary$cybert_logFC, final_summary$EB_logFC)
rho_cybert_vs_EB <- cor.test(~ as.numeric(final_summary$cybert_logFC) + as.numeric(final_summary$EB_logFC), data=final_summary, method = "spearman", continuity = TRUE)$estimate
r_squared_cybert_vs_EB <- summary(lm(final_summary$cybert_logFC ~ final_summary$EB_logFC, data=final_summary))$adj.r.squared


