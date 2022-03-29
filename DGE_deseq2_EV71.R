#Ok

library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(IHW)

setwd("D://Projects/kenezam/gene_counts/GSE103308/")

df1 <- read.csv("GSE103308_ribo_mock_novirus_1.75h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df2 <- read.csv("GSE103308_ribo_mock_novirus_3.5h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df3 <- read.csv("GSE103308_ribo_mock_novirus_3.5h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df4 <- read.csv("GSE103308_ribo_mock_novirus_6.25h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df5 <- read.csv("GSE103308_ribo_mock_novirus_6.25h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df6 <- read.csv("GSE103308_ribo_virus_EV71_1.75h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df7 <- read.csv("GSE103308_ribo_virus_EV71_3.5h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df8 <- read.csv("GSE103308_ribo_virus_EV71_3.5h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df9 <- read.csv("GSE103308_ribo_virus_EV71_6.25h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df10 <- read.csv("GSE103308_ribo_virus_EV71_6.25h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

dfr1 <- read.csv("GSE103308_rna_mock_novirus_1.75h_1_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr2 <- read.csv("GSE103308_rna_mock_novirus_3.5h_1_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr3 <- read.csv("GSE103308_rna_mock_novirus_3.5h_2_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr4 <- read.csv("GSE103308_rna_mock_novirus_6.25h_1_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr5 <- read.csv("GSE103308_rna_mock_novirus_6.25h_2_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr6 <- read.csv("GSE103308_rna_virus_EV71_1.75h_1_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr7 <- read.csv("GSE103308_rna_virus_EV71_3.5h_1_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr8 <- read.csv("GSE103308_rna_virus_EV71_3.5h_2_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr9 <- read.csv("GSE103308_rna_virus_EV71_6.25h_1_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr10 <- read.csv("GSE103308_rna_virus_EV71_6.25h_2_4s.fastq/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

cts <- purrr::reduce(
  #list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, dfr1, dfr2, dfr3, dfr4, dfr5, dfr6, dfr7, dfr8, dfr9, dfr10),
  list(df2, df3, df4, df5, df7, df8, df9, df10, dfr2, dfr3, dfr4, dfr5, dfr7, dfr8, dfr9, dfr10),
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
                   "ribo_mock_novirus_3.5h_1",
                   "ribo_mock_novirus_3.5h_2",
                   "ribo_mock_novirus_6.25h_1",
                   "ribo_mock_novirus_6.25h_2",
                   "ribo_virus_EV71_3.5h_1",
                   "ribo_virus_EV71_3.5h_2",
                   "ribo_virus_EV71_6.25h_1",
                   "ribo_virus_EV71_6.25h_2",
                   "rna_mock_novirus_3.5h_1",
                   "rna_mock_novirus_3.5h_2",
                   "rna_mock_novirus_6.25h_1",
                   "rna_mock_novirus_6.25h_2",
                   "rna_virus_EV71_3.5h_1",
                   "rna_virus_EV71_3.5h_2",
                   "rna_virus_EV71_6.25h_1",
                   "rna_virus_EV71_6.25h_2")

#w
d <- cts %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_EV71 <- sum(d[1,])

#w
d <- cts[,c(3,4,7,8,11,12,15,16)] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_EV71_35h <- sum(d[1,])

#w
d <- cts[,-c(1,2,3,4,7,8,11,12,15,16)] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_EV71_625h <- sum(d[1,])

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

time <- c(rep(c(3.5, 3.5, 6.25, 6.25),4))
time <- factor(time)
type <- c(rep('mock', 4), rep('EV71', 4), rep('mock', 4), rep('EV71', 4))
type <- factor(type)
reg <- c(rep('ribo', 8), rep('rna', 8))
reg <- factor(reg)
rep <- c(rep(c('1', '2'),8))
rep <- factor(rep)
#time <- relevel(time, ref = 3.5)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type
#ann <- data.frame(type=type,reg=reg,rep=rep,time=time)
#write.table(ann, str_c("./",'ann.tsv'), sep="\t", quote = F, row.names = F)

coldata <- data.frame(type=type,reg=reg,time=time)
colda <- data.frame(group=group,time=time)

#good filtering
keep <- filterByExpr(cts[,3:ncol(cts)], group = group)
cts_f <- cts[keep,]

x <- cts_f[,c(3:ncol(cts_f))]

rownames(coldata) <- colnames(x)

rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")

#1. ???????????????????? ???????????????? ?????????????????? ??????????
dds35 <- DESeq2::DESeqDataSetFromMatrix(countData = dplyr::select(x, contains("3.5")),
                                      colData = coldata[coldata$time != 6.25,c(1,2)],
                                      design = ~reg + type + reg:type)

dds625 <- DESeq2::DESeqDataSetFromMatrix(countData = dplyr::select(x, contains("6.25")),
                                        colData = coldata[coldata$time != 3.5,c(1,2)],
                                        design = ~reg + type + reg:type)

#dds <- dds625
dds <- dds35

dds <- DESeq(dds, test = "LRT", reduced = ~reg + type)
dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

resultsNames(dds) # lists the coefficients
res <- results(dds, filterFun=ihw)
res_la <- results(dds, name="regribo.typeEV71", test="Wald", lfcThreshold=1, altHypothesis="lessAbs", filterFun=ihw)
res_l <- results(dds, name="regribo.typeEV71", test="Wald", lfcThreshold=0, altHypothesis="less", filterFun=ihw)
head(res[order(res$pvalue), ], 10)
summary(res)

## Checking the normalization
epsilon <- 1
col.strain <- c("mock"="blue","EV71"="orange")

#RIBO
colnames(dds.norm) <- c("ribo_mock_1", "ribo_mock_2", "ribo_EV71_1", "ribo_EV71_2", "rna_mock_1", "rna_mock_2", "rna_EV71_1", "rna_EV71_2")
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

res_RO_EV71 <- as.data.frame(res)
res_RO_EV71_la <- as.data.frame(res_la)
res_RO_EV71_l <- as.data.frame(res_l)

write.table(res_RO_EV71, str_c("./",'res_RO_ds2_35.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_EV71_la, str_c("./",'res_RO_ds2_35_la.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_EV71_l, str_c("./",'res_RO_ds2_35_l.tsv'), sep="\t", quote = F, row.names = T)

#2. ???? ???????????????? ???????????????? ????-?????????????? ?? ??????????????
dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                      colData = coldata,
                                      #colData = colda,
                                      design = ~reg + type + reg:type + reg:type:time)
                                      #design = ~group + group:time)

dds <- DESeq(dds, test = "LRT", reduced = ~group)

dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

resultsNames(dds) # lists the coefficients
res <- results(dds)
res_la <- results(dds, name="regribo.typesarscov2", test="Wald", lfcThreshold=1, altHypothesis="lessAbs", filterFun=ihw)
res_l <- results(dds, name="regribo.typesarscov2", test="Wald", lfcThreshold=0, altHypothesis="less", filterFun=ihw)
head(res[order(res$pvalue), ], 10)
summary(res)

## Checking the normalization
epsilon <- 1
col.strain <- c("mock"="blue","EV71"="orange")

#RIBO
colnames(dds.norm) <- c("ribo_mock_1", "ribo_mock_2", "ribo_EV71_1", "ribo_EV71_2", "rna_mock_1", "rna_mock_2", "rna_EV71_1", "rna_EV71_2")par(mfrow=c(2,2),cex.lab=0.7)
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