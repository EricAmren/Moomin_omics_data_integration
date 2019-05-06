# install.packages("memisc")
library(memisc)
setwd("/home/ericcumunel/Documents/internship/rnaseq_MA_comparison/RNAseq/")

par(mfrow=c(2,3))
files <- list.files(path="/home/ericcumunel/Documents/internship/rnaseq_MA_comparison/RNAseq/new_featureCounts_results", pattern="*.tsv$", full.names=TRUE, recursive=FALSE)

get_title <- function(x) {
  title=cases(
    "sample 1 : Control"=grepl("ERR2955747", x),
    "sample 2 : Control"=grepl("ERR2955748", x),
    "sample 3 : Control"=grepl("ERR2955749", x),
    "sample 4 : with spinach leaves"=grepl("ERR2955750", x),
    "sample 5 : with spinach leaves"=grepl("ERR2955751", x),
    "sample 6 : with spinach leaves"=grepl("ERR2955752", x),
    "wtfisthat"=TRUE
  )
}

data = read.csv(file = "new_featureCounts_results/trimmed_ERR2955747.tsv", header = TRUE, sep = "\t", skip=1)
gene_id <- data[,1]
counts_list <- list()
counts_list[["gene_id"]] <- gene_id
listplot <- function(x) {
  data = read.csv(file = x, header = TRUE, sep = "\t", skip=1)
  # plot(data[,length(data)], type = "p",
  #      main = get_title(x),
  #      xlab = "genes",
  #      ylab = "counts"
  # )
  counts <- data[,length(data)]
  sample_name <<- sub('\\.tsv$', '',strsplit(x, "_")[[1]][6])
  counts_list[[sample_name]] <<- counts
}
lapply(files,listplot)

output_table <- data.frame(counts_list)
output_table <- output_table[c("gene_id", "ERR2955750", "ERR2955751", "ERR2955752", "ERR2955747", "ERR2955748", "ERR2955749" )]

write.table(output_table, file="counts.txt", row.names=FALSE,sep ="\t", quote = FALSE )


library(EBSeq)
counts = read.csv("counts.txt", header = TRUE, sep = "\t")
data=as.matrix(counts)
rnames=counts[,1]
rownames(data)=rnames
data=data[,-1]
class(data)="numeric"

conditions = as.factor(c("C1","C1","C1","C2","C2","C2"))
Sizes=MedianNorm(data)
EBOut=EBTest(Data=data,Conditions=conditions,sizeFactors=Sizes, maxround=5)
EBDERes=GetDEResults(EBOut, FDR=0.05)

EBseq_results <-data.frame(EBDERes$PPMat, counts_list[2:7])
EBseq_results <- EBseq_results[order(EBseq_results$PPDE, decreasing=TRUE),]

#check convergence
rounds = 5
diff(EBOut$Alpha[(rounds-2):rounds,])
diff(EBOut$Beta[(rounds-2):rounds,])
diff(EBOut$P[(rounds-2):rounds,])

#check model fit
par(mfrow=c(1,1))
QQP(EBOut)
DenNHist(EBOut)

EBFC <- PostFC(EBOut)
df_EBFC <- as.data.frame(EBFC["PostFC"])
df_EBFC <- data.frame(ECs=rownames(df_EBFC), logFC=log(df_EBFC[,1]))
df1 <- data.frame(EBDERes$PPMat)
head(df1)
df1$names <- rownames(EBDERes$PPMat)
# colnames(df1) <- c("ECs", "PPEE", "PPDE")
# 
# head(df1)
# head(df_EBFC)
head(EBDERes$PPMat)
final_EBseq_results <- merge(df1, df_EBFC, by.x="names", by.y="ECs")
write.table(final_EBseq_results, file="final_EBseq_results.tsv", sep="\t", row.names = FALSE, quote = FALSE)

