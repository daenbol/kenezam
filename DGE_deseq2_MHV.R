#Unable to remove hidden batches

library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(IHW)
library(IHW)
library(sva)
library(ggfortify)
library(affy)

setwd("D://Projects/kenezam/gene_counts/MHV/")


df1 <- read.csv("EMTAB8650_ribo_mock_novirus_5h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df2 <- read.csv("EMTAB8650_ribo_mock_novirus_5h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df3 <- read.csv("EMTAB8650_ribo_cntrl_novirus_6h_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df4 <- read.csv("EMTAB8650_ribo_mock_novirus_8h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df5 <- read.csv("EMTAB8650_ribo_virus_MHVA59_5h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df6 <- read.csv("EMTAB8650_ribo_virus_MHVA59_5h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df7 <- read.csv("EMTAB8650_ribo_virus_MHVA59_8h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

dfr1 <- read.csv("EMTAB8650_rna_mock_novirus_5h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr2 <- read.csv("EMTAB8650_rna_mock_novirus_5h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr3 <- read.csv("EMTAB8650_rna_mock_novirus_8h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr4 <- read.csv("EMTAB8650_rna_virus_MHVA59_5h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr5 <- read.csv("EMTAB8650_rna_virus_MHVA59_5h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr6 <- read.csv("EMTAB8650_rna_virus_MHVA59_8h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

cts <- purrr::reduce(
  list(df1, df2, df3, df4, df5, df6, df7, dfr1, dfr2, dfr3, dfr4, dfr5, dfr6),
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
                   "ribo_mock_novirus_5h_1",
                   "ribo_mock_novirus_5h_2",
                   "ribo_cntrl_novirus_6h",
                   "ribo_mock_novirus_8h_1",
                   "ribo_virus_MHVA59_5h_1",
                   "ribo_virus_MHVA59_5h_2",
                   "ribo_virus_MHVA59_8h_1",
                   "rna_mock_novirus_5h_1",
                   "rna_mock_novirus_5h_2",
                   "rna_mock_novirus_8h_1",
                   "rna_virus_MHVA59_5h_1",
                   "rna_virus_MHVA59_5h_2",
                   "rna_virus_MHVA59_8h_1")

#w
d <- cts[,-c(5,6,9,12,15)] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_MHV <- sum(d[1,])

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

#Filtering using filterByExpr
type <- c('mock', 'mock', 'ctrl', 'mock', 'MHVA59', 'MHVA59', 'MHVA59', 'mock', 'mock', 'mock', 'MHVA59', 'MHVA59', 'MHVA59')
type <- factor(type)
reg <- c('ribo', 'ribo','ribo', 'ribo', 'ribo', 'ribo', 'ribo', 'rna','rna','rna','rna','rna','rna')
reg <- factor(reg)
rep <- c('1', '2', '1', '1', '1', '2', '1', '1', '2', '1', '1', '2', '1')
rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type

coldata <- data.frame(type=type,reg=reg,rep=rep)

#good filtering
keep <- filterByExpr(cts[,3:ncol(cts)], group = group[-c(3,4,7,10,13)])
cts_f <- cts[keep,]

### COMMON PART
x <- cts_f[,c(3:ncol(cts_f))]

rownames(coldata) <- colnames(x)

rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")[keep]

xm <- x[,-c(3,4,7,10,13)]
cold <- coldata[-c(3,4,7,10,13),]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = xm,
                                      colData = cold,
                                      design = ~reg + type + reg:type)

dds <- DESeq(dds, test = "LRT", reduced = ~reg + type)

dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

resultsNames(dds) # lists the coefficients
res <- results(dds, filterFun=ihw)
res_la <- results(dds, name="regribo.typeMHVA59", test="Wald", lfcThreshold=1, altHypothesis="lessAbs", filterFun=ihw)
res_l <- results(dds, name="regribo.typeMHVA59", test="Wald", lfcThreshold=0, altHypothesis="less", filterFun=ihw)

head(res[order(res$pvalue), ], 20)
summary(res)

## Checking the normalization
epsilon <- 1
col.strain <- c("mock"="blue","MHVA59"="orange")

#RIBO
colnames(dds.norm) <- c("ribo_MHV_1", "ribo_MHV_2", "ribo_mock_1", "ribo_mock_2", "rna_MHV_1", "rna_MHV_2", "rna_mock_1", "rna_mock_2")
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)[,1:4]+epsilon),  col=col.strain[as.vector(cold$type[1:4])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),1:4]+epsilon),  col=col.strain[as.vector(cold$type[1:4])], cex.axis=0.7, 
        las=1, xlab="log2(filtered & normalized counts)", horizontal=TRUE, main="filtered & normalized counts") 
plotDensity(log2(counts(dds.norm)+epsilon)[,1:4],  col=col.strain[as.vector(cold$type[1:4])], lwd=2,
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),1:4]+epsilon), col=col.strain[as.vector(cold$type[1:4])], lwd=2,
            xlab="log2(filtered & normalized counts)", cex.lab=0.7, panel.first=grid())

#RNA
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)[,5:8]+epsilon),  col=col.strain[as.vector(cold$type[5:8])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),5:8]+epsilon),  col=col.strain[as.vector(cold$type[5:8])], cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds.norm)[,5:8]+epsilon),  col=col.strain[as.vector(cold$type[5:8])], lwd=2,
            xlab="log2(counts+1)", cex.lab=0.7, panel.first=grid())
plotDensity(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),5:8]+epsilon), col=col.strain[as.vector(cold$type[5:8])], lwd=2,
            xlab="log2(normalized counts+1)", cex.lab=0.7, panel.first=grid())

res_RO_MHV <- as.data.frame(res)
res_RO_MHV_la <- as.data.frame(res_la)
res_RO_MHV_l <- as.data.frame(res_l)

write.table(res_RO_MHV, str_c("./",'RO_MHV_ds2.tsv'), sep="\t", quote = F)
write.table(res_RO_MHV_la, str_c("./",'RO_MHV_ds2_la.tsv'), sep="\t", quote = F)
write.table(res_RO_MHV_l, str_c("./",'RO_MHV_ds2_l.tsv'), sep="\t", quote = F)

par(mfrow=c(2,2),cex.lab=0.7)
hist(res_RO_MHV_la$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald LessAbs")
hist(res_RO_MHV_l$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald 1-sided")
hist(res$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="LRT 2-sided")