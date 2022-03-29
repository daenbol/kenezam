#Ok rep
remove.packages("rlang")
library(purrr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(IHW)
library(sva)
library(ggplot2)
library(rlang)
library(ggfortify)

setwd("D://Projects/kenezam/gene_counts/GSE101760/")


df1 <- read.csv("GSE101760_ribo_ctrl_novirus_12h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df2 <- read.csv("GSE101760_ribo_ctrl_novirus_12h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df3 <- read.csv("GSE101760_ribo_mock_novirus_12h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df4 <- read.csv("GSE101760_ribo_mock_novirus_12h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df5 <- read.csv("GSE101760_ribo_virus_IAV_12h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df6 <- read.csv("GSE101760_ribo_virus_IAV_12h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

dfr1 <- read.csv("GSE101760_rna_ctrl_novirus_12h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr2 <- read.csv("GSE101760_rna_ctrl_novirus_12h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr3 <- read.csv("GSE101760_rna_mock_novirus_12h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr4 <- read.csv("GSE101760_rna_mock_novirus_12h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr5 <- read.csv("GSE101760_rna_virus_IAV_12h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr6 <- read.csv("GSE101760_rna_virus_IAV_12h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

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
                   "ribo_ctrl_novirus_12h_1",
                   "ribo_ctrl_novirus_12h_2",
                   "ribo_mock_novirus_12h_1",
                   "ribo_mock_novirus_12h_2",
                   "ribo_virus_IAV_12h_1",
                   "ribo_virus_IAV_12h_2",
                   "rna_ctrl_novirus_12h_1",
                   "rna_ctrl_novirus_12h_2",
                   "rna_mock_novirus_12h_1",
                   "rna_mock_novirus_12h_2",
                   "rna_virus_IAV_12h_1",
                   "rna_virus_IAV_12h_2")

#w
#d <- cts %>%
#  summarize_if(is.numeric, sum, na.rm=TRUE)
#w_IAV <- sum(d[1,])

#w
d <- cts[,-c(3,4,9,10)] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_IAV <- sum(d[1,])

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

#Filtering using filterByExpr
#type <- c('ctrl', 'ctrl', 'mock', 'mock', 'IAV', 'IAV', 'ctrl', 'ctrl', 'mock', 'mock', 'IAV', 'IAV')
type <- c('mock', 'mock', 'mock', 'mock', 'IAV', 'IAV', 'mock', 'mock', 'mock', 'mock', 'IAV', 'IAV')
type <- factor(type)
reg <- c('ribo', 'ribo','ribo', 'ribo', 'ribo', 'ribo', 'rna','rna','rna','rna','rna','rna')
reg <- factor(reg)
rep <- c(rep(c('1', '2'),6))
rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type

coldata <- data.frame(type=type,reg=reg,rep=rep)

#good filtering
keep <- filterByExpr(cts[,3:ncol(cts)], group = group)
cts_f <- cts[keep,]

x <- cts_f[,c(3:ncol(cts_f))]

rownames(coldata) <- colnames(x)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                      colData = coldata, 
                                      design = ~reg + type + rep + reg:type)
dds <- DESeq(dds, test = "LRT", reduced = ~reg + type + rep)

dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

#SVA

#ALL
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ group, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
#x <- tibble::column_to_rownames(cts[2:ncol(x)], 'ENSEMBL_ID')
#xx <- log2(x + 2)
n.sv = num.sv(as.matrix(dat),mod,method="leek")
n.sv
#1
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

t <- svseq$sv
t <- cbind(t, as.factor(type), as.factor(reg))
colnames(t) <- c("SV1", "SV2", "type", "reg")
ggp <- ggplot(t, aes(x=SV1, y=SV2)) + 
  geom_point(aes(shape=as.factor(type), color=as.factor(reg)))
ggp

#2
svseq <- svaseq(dat, mod, mod0, n.sv = 7)

sv.pca <- prcomp(svseq$sv,
                 center = TRUE,
                 scale. = TRUE)

k2 <- kmeans(svseq$sv, centers = 2, nstart = 25)
str(k2)

library(factoextra)
fviz_cluster(k2, data = svseq$sv)

library(ggfortify)
sv.pca.plot <- autoplot(sv.pca,
                        data = svseq$sv)

sv.pca.plot

#RIBO
dds <- DESeq2::DESeqDataSetFromMatrix(countData = x[reg != "rna"],
                                      colData = coldata[reg != "rna",], 
                                      design = ~type)
dds <- DESeq(dds)

dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ group, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
#x <- tibble::column_to_rownames(cts[2:ncol(x)], 'ENSEMBL_ID')
#xx <- log2(x + 2)
n.sv = num.sv(as.matrix(dat),mod,method="leek")
n.sv
#1
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

t <- svseq$sv
t <- cbind(t, as.factor(type), as.factor(reg))
colnames(t) <- c("SV1", "SV2", "type", "reg")
ggp <- ggplot(t, aes(x=SV1, y=SV2)) + 
  geom_point(aes(shape=as.factor(type), color=as.factor(reg)))
ggp

#2
svseq <- svaseq(dat, mod, mod0, n.sv = 7)

sv.pca <- prcomp(svseq$sv,
                 center = TRUE,
                 scale. = TRUE)

k2 <- kmeans(svseq$sv, centers = 2, nstart = 25)
str(k2)

library(factoextra)
fviz_cluster(k2, data = svseq$sv)

library(ggfortify)
sv.pca.plot <- autoplot(sv.pca,
                        data = svseq$sv)

sv.pca.plot

dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

resultsNames(dds) # lists the coefficients
res <- results(dds, filterFun=ihw)
res_la <- results(dds, name="regribo.typeIAV", test="Wald", lfcThreshold=1, altHypothesis="lessAbs", filterFun=ihw)
res_l <- results(dds, name="regribo.typeIAV", test="Wald", lfcThreshold=0, altHypothesis="less", filterFun=ihw)

head(res[order(res$pvalue), ], 20)
summary(res)

## Checking the normalization
epsilon <- 1
col.strain <- c("mock"="blue","IAV"="orange")

#RIBO
colnames(dds.norm) <- c("ribo_ctrl_1", "ribo_ctrl_2", "ribo_mock_1", "ribo_mock_2", "ribo_IAV_1", "ribo_IAV_2", "rna_ctrl_1", "rna_ctrl_2", "rna_mock_1", "rna_mock_2", "rna_IAV_1", "rna_IAV_2")
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)[,1:6]+epsilon),  col=col.strain[as.vector(coldata$type[1:6])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)[!is.na(res@listData$padj),1:6]+epsilon),  col=col.strain[as.vector(coldata$type[1:6])], cex.axis=0.7, 
        las=1, xlab="log2(filtered & normalized counts)", horizontal=TRUE, main="filtered & normalized counts") 
plotDensity(log2(counts(dds.norm)+epsilon)[,1:6],  col=col.strain[as.vector(coldata$type[1:6])], lwd=2,
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)[!is.na(res@listData$padj),1:6]+epsilon), col=col.strain[as.vector(coldata$type[1:6])], lwd=2,
            xlab="log2(filtered & normalized counts)", cex.lab=0.7, panel.first=grid())

#RNA
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)[,7:12]+epsilon),  col=col.strain[as.vector(coldata$type[7:12])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),7:12]+epsilon),  col=col.strain[as.vector(coldata$type[7:12])], cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds.norm)[,7:12]+epsilon),  col=col.strain[as.vector(coldata$type[7:12])], lwd=2,
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),7:12]+epsilon), col=col.strain[as.vector(coldata$type[7:12])], lwd=2,
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())

res_RO_IAV <- as.data.frame(res)
res_RO_IAV_la <- as.data.frame(res_la)
res_RO_IAV_l <- as.data.frame(res_l)

write.table(res_RO_IAV, str_c("./",'res_RO_ds2_b.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_IAV_la, str_c("./",'res_RO_ds2_b_la.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_IAV_l, str_c("./",'res_RO_ds2_b_l.tsv'), sep="\t", quote = F, row.names = T)

par(mfrow=c(2,2),cex.lab=0.7)
hist(res_RO_IAV_la$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald LessAbs")
hist(res_RO_IAV_l$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald 1-sided")
hist(res_RO_IAV$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="LRT 2-sided")