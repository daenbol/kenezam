#!!!!  MOUSE



library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(IHW)
library(affy)

setwd("D://Projects/kenezam/gene_counts/GSE137757")

df1 <- read.csv("GSE137757_ribo_virus_reovirus_18h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df2 <- read.csv("GSE137757_ribo_virus_reovirus_18h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df3 <- read.csv("GSE137757_ribo_mock_novirus_18h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df4 <- read.csv("GSE137757_ribo_mock_novirus_18h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

dfr1 <- read.csv("GSE137757_rna_virus_reovirus_18h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr2 <- read.csv("GSE137757_rna_virus_reovirus_18h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr3 <- read.csv("GSE137757_rna_mock_novirus_18h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr4 <- read.csv("GSE137757_rna_mock_novirus_18h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

cts <- purrr::reduce(
  list(df1, df2, df3, df4, dfr1, dfr2, dfr3, dfr4),
  full_join,
  by = "N_unmapped")

id = gsub("\\..*","",cts$N_unmapped)
cts <- cbind(id, cts)

#good filtering
keep <- filterByExpr(cts[,3:ncol(cts)], group = group)
cts_f <- cts[keep,]

library("AnnotationDbi")
library("org.Hs.eg.db")
cts$id <- AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                keys = cts$id,
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "list")

colnames(cts) <- c("GENES_SYMBOL",
                   "ENSEMBL_ID",
                   "ribo_virus_reovirus_18h_1",
                   "ribo_virus_reovirus_18h_2",
                   "ribo_mock_novirus_18h_1",
                   "ribo_mock_novirus_18h_2",
                   "rna_virus_reovirus_18h_1",
                   "rna_virus_reovirus_18h_2",
                   "rna_mock_novirus_18h_1",
                   "rna_mock_novirus_18h_2")

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

type <- c('reovirus', 'reovirus', 'mock', 'mock', 'reovirus', 'reovirus', 'mock', 'mock')
type <- factor(type)
reg <- c('ribo', 'ribo','ribo', 'ribo', 'rna','rna','rna','rna')
rep <- c(rep(c('1', '2'),4))
rep <- factor(rep)
reg <- factor(reg)
group  = reg:type

coldata <- data.frame(type=type,reg=reg,rep=rep)

x <- cts_f[,c(3:ncol(cts_f))]

rownames(coldata) <- colnames(x)

rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                      colData = coldata, 
                                      design = ~reg + type + reg:type)
dds <- DESeq(dds, test = "LRT", reduced = ~reg + type)
dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)


## Checking the normalization
epsilon <- 1
col.strain <- c("mock"="blue","reovirus"="orange")

#RIBO
colnames(dds.norm) <- c("ribo_reo_1", "ribo_reo_2", "ribo_mock_1", "ribo_mock_2", "rna_reo_1", "rna_reo_2", "rna_mock_1", "rna_mock_2")
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

## Computing mean and variance
norm.counts <- counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),5:8]
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)

mean.var.col <- densCols(x=log2(mean.counts), y=log2(variance.counts))
plot(x=log2(mean.counts), y=log2(variance.counts), pch=16, cex=0.5, 
     col=mean.var.col, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")

resultsNames(dds) # lists the coefficients
res <- results(dds)
res_la <- results(dds, name="regrna.typereovirus", test="Wald", lfcThreshold=1, altHypothesis="lessAbs", filterFun=ihw)
res_l <- results(dds, name="regrna.typereovirus", test="Wald", lfcThreshold=0, altHypothesis="less", filterFun=ihw)

#????????????
filterGenes <- res$baseMean > 0
res$pvalue[ !filterGenes ] <- NA
res$padj <- p.adjust( res$pvalue, method="BH")

#ok
myCPM <- cpm(cts[,3:ncol(cts)])
plot(myCPM[,3],cts[,3],ylim=c(0,50),xlim=c(0,0.001))
abline(h = 10, col = "red")
abline(v = 0.15, col = "blue")

head(res[order(res$pvalue), ], 10)
summary(res)

res_RO_reo_la <- as.data.frame(res_la)
res_RO_reo_l <- as.data.frame(res_l)

write.table(res_RO_reo_la[order(res_RO_reo_la$padj),], str_c("./",'res_RO_la.tsv'), sep="\t", quote = F, row.names = T, )
write.table(res_RO_reo_l, str_c("./",'res_RO_l.tsv'), sep="\t", quote = F, row.names = T)

write.table(cts, str_c("./",'sum_reo_cts.tsv'), sep="\t", quote = F, row.names = T)

tst <- read.csv("res_RO_la.tsv", sep="\t", header = T)

par(mfrow=c(2,2),cex.lab=0.7)
hist(res_RO_reo_la$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald LessAbs")
hist(res_RO_reo_l$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald 1-sided")
hist(res$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="LRT 2-sided")