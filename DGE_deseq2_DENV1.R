library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(IHW)

setwd("D://Projects/kenezam/gene_counts/GSE69602/")


df1 <- read.csv("GSE69602_ribo_mock_novirus_notime_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df2 <- read.csv("GSE69602_ribo_mock_novirus_notime_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df3 <- read.csv("GSE69602_ribo_virus_DENV1_48h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df4 <- read.csv("GSE69602_ribo_virus_DENV1_48h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df5 <- read.csv("GSE69602_ribo_virus_DENV1_72h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df6 <- read.csv("GSE69602_ribo_virus_DENV1_72h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

dfr1 <- read.csv("GSE69602_rna_mock_novirus_notime_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr2 <- read.csv("GSE69602_rna_mock_novirus_notime_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr3 <- read.csv("GSE69602_rna_virus_DENV1_48h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr4 <- read.csv("GSE69602_rna_virus_DENV1_48h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr5 <- read.csv("GSE69602_rna_virus_DENV1_72h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr6 <- read.csv("GSE69602_rna_virus_DENV1_72h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

cts <- purrr::reduce(
  list(df1, df2, df3, df4, df5, df6, dfr1, dfr2, dfr3, dfr4, dfr5, dfr6),
  full_join,
  by = "N_unmapped")

id = gsub("\\..*","",cts$N_unmapped)
cts <- cbind(id, cts)

library("AnnotationDbi")
library("org.Hs.eg.db")
cts$id <- AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                keys = cts$id,
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "list")

colnames(cts) <- c("GENES_SYMBOL",
                   "ENSEMBL_ID",
                   "ribo_mock_novirus_notime_1",
                   "ribo_mock_novirus_notime_2",
                   "ribo_virus_DENV1_48h_1",
                   "ribo_virus_DENV1_48h_2",
                   "ribo_virus_DENV1_72h_1",
                   "ribo_virus_DENV1_72h_2",
                   "rna_mock_novirus_notime_1",
                   "rna_mock_novirus_notime_2",
                   "rna_virus_DENV1_48h_1",
                   "rna_virus_DENV1_48h_2",
                   "rna_virus_DENV1_72h_1",
                   "rna_virus_DENV1_72h_2")

#w
d <- cts %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_DENV1 <- sum(d[1,])

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

type <- c('mock', 'mock', 'd1_48', 'd1_48', 'd1_72', 'd1_72', 'mock', 'mock', 'd1_48', 'd1_48', 'd1_72', 'd1_72')
type <- factor(type)
reg <- c('ribo', 'ribo','ribo', 'ribo', 'ribo', 'ribo', 'rna','rna','rna','rna','rna','rna')
reg <- factor(reg)
group  = reg:type

coldata <- data.frame(type=type,reg=reg)

#good filtering
keep <- filterByExpr(cts[,3:ncol(cts)], group = group)
cts_f <- cts[keep,]

x <- cts_f[,c(3:ncol(cts_f))]

rownames(coldata) <- colnames(x)

rownames(x) <- paste(cts_f$GENES_SYMBOL, cts_f$ENSEMBL_ID, sep="__")

#48h
dds <- DESeq2::DESeqDataSetFromMatrix(countData = x[type != "d1_72"],
                                      colData = coldata[type != "d1_72",], 
                                      design = ~reg + type + reg:type)

dds <- DESeq(dds, test = "LRT", reduced = ~reg + type)

dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

resultsNames(dds) # lists the coefficients
res <- results(dds, filterFun=ihw)
res_la <- results(dds, name="regrna.typemock", test="Wald", lfcThreshold=1, altHypothesis="lessAbs", filterFun=ihw)
res_l <- results(dds, name="regrna.typemock", test="Wald", lfcThreshold=0, altHypothesis="less", filterFun=ihw)

head(res[order(res$pvalue), ], 10)
summary(res)

## Checking the normalization
epsilon <- 1
col.strain <- c("mock"="blue","DENV1"="orange")

#RIBO
colnames(dds.norm) <- c("ribo_mock_novirus_notime_1", "ribo_mock_novirus_notime_2", "ribo_virus_DENV1_48h_1", "ribo_virus_DENV1_48h_2", "rna_mock_novirus_notime_1", "rna_mock_novirus_notime_2", "rna_virus_DENV1_48h_1", "rna_virus_DENV1_48h_2")
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)[,1:4]+epsilon),  col=col.strain[as.vector(coldata$type[1:4])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)[!is.na(res@listData$padj),1:4]+epsilon),  col=col.strain[as.vector(coldata$type[1:4])], cex.axis=0.7, 
        las=1, xlab="log2(filtered & normalized counts)", horizontal=TRUE, main="filtered & normalized counts") 
plotDensity(log2(counts(dds.norm)+epsilon)[,1:4],  col=col.strain[as.vector(coldata$type[1:4])], lwd=2,
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)[!is.na(res@listData$padj),1:4]+epsilon), col=col.strain[as.vector(coldata$type[1:4])], lwd=2,
            xlab="log2(filtered & normalized counts)", cex.lab=0.7, panel.first=grid())

#RNA
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)[,5:8]+epsilon),  col=col.strain[as.vector(coldata$type[5:8])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),5:8]+epsilon),  col=col.strain[as.vector(coldata$type[5:8])], cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds.norm)[,5:8]+epsilon),  col=col.strain[as.vector(coldata$type[5:8])], lwd=2,
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),5:8]+epsilon), col=col.strain[as.vector(coldata$type[5:8])], lwd=2,
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())

res_RO_DENV1_48 <- as.data.frame(res)
res_RO_DENV1_48_la <- as.data.frame(res_la)
res_RO_DENV1_48_l <- as.data.frame(res_l)

write.table(res_RO_DENV1_48, str_c("./",'res_RO_ds2_48h.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_DENV1_48_la, str_c("./",'res_RO_ds2_48h_la.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_DENV1_48_l, str_c("./",'res_RO_ds2_48h_l.tsv'), sep="\t", quote = F, row.names = T)

par(mfrow=c(2,2),cex.lab=0.7)
hist(res_RO_DENV1_48_la$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald LessAbs")
hist(res_RO_DENV1_48_l$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald 1-sided")
hist(res_RO_DENV1_48$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="LRT 2-sided")

#72h
dds <- DESeq2::DESeqDataSetFromMatrix(countData = x[type != "d1_48"],
                                      colData = coldata[type != "d1_48",], 
                                      design = ~reg + type + reg:type)
dds <- DESeq(dds, test = "LRT", reduced = ~reg + type)
resultsNames(dds) # lists the coefficients

res <- results(dds, filterFun=ihw)
res_la <- results(dds, name="regrna.typemock", test="Wald", lfcThreshold=1, altHypothesis="lessAbs", filterFun=ihw)
res_l <- results(dds, name="regrna.typemock", test="Wald", lfcThreshold=0, altHypothesis="less", filterFun=ihw)

head(res[order(res$pvalue), ], 10)
summary(res)

## Checking the normalization
epsilon <- 1
col.strain <- c("mock"="blue","DENV1"="orange")

#RIBO
colnames(dds.norm) <- c("ribo_mock_novirus_notime_1", "ribo_mock_novirus_notime_2", "ribo_virus_DENV1_72h_1", "ribo_virus_DENV1_72h_2", "rna_mock_novirus_notime_1", "rna_mock_novirus_notime_2", "rna_virus_DENV1_72h_1", "rna_virus_DENV1_72h_2")
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)[,1:4]+epsilon),  col=col.strain[as.vector(coldata$type[1:4])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)[!is.na(res@listData$padj),1:4]+epsilon),  col=col.strain[as.vector(coldata$type[1:4])], cex.axis=0.7, 
        las=1, xlab="log2(filtered & normalized counts)", horizontal=TRUE, main="filtered & normalized counts") 
plotDensity(log2(counts(dds.norm)+epsilon)[,1:4],  col=col.strain[as.vector(coldata$type[1:4])], lwd=2,
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)[!is.na(res@listData$padj),1:4]+epsilon), col=col.strain[as.vector(coldata$type[1:4])], lwd=2,
            xlab="log2(filtered & normalized counts)", cex.lab=0.7, panel.first=grid())

#RNA
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)[,5:8]+epsilon),  col=col.strain[as.vector(coldata$type[5:8])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),5:8]+epsilon),  col=col.strain[as.vector(coldata$type[5:8])], cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds.norm)[,5:8]+epsilon),  col=col.strain[as.vector(coldata$type[5:8])], lwd=2,
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),5:8]+epsilon), col=col.strain[as.vector(coldata$type[5:8])], lwd=2,
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())

res_RO_DENV1_72 <- as.data.frame(res)
res_RO_DENV1_72_la <- as.data.frame(res_la)
res_RO_DENV1_72_l <- as.data.frame(res_l)

write.table(res_RO_DENV1_72, str_c("./",'res_RO_ds2_72h.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_DENV1_72_la, str_c("./",'res_RO_ds2_72h_la.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_DENV1_72_l, str_c("./",'res_RO_ds2_72h_l.tsv'), sep="\t", quote = F, row.names = T)

par(mfrow=c(2,2),cex.lab=0.7)
hist(res_RO_DENV1_72_la$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald LessAbs")
hist(res_RO_DENV1_72_l$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald 1-sided")
hist(res_RO_DENV1_72$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="LRT 2-sided")