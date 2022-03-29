#LAST due to deseq design problems
#Hundreds of thousends of ribo counts

library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(IHW)
library(sva)

remove.packages("DESeq2")
library(devtools)
install_github("daenbol/DESeq2")

setwd("D://Projects/kenezam/gene_counts/SARSCOV2/")

df1 <- read.csv("PRJNA704763_ribo_mock_novirus_24h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df2 <- read.csv("PRJNA704763_ribo_mock_novirus_24h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df3 <- read.csv("PRJNA704763_ribo_mock_novirus_24h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df4 <- read.csv("PRJNA704763_ribo_virus_SARSCOV2_6h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df5 <- read.csv("PRJNA704763_ribo_virus_SARSCOV2_6h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df6 <- read.csv("PRJNA704763_ribo_virus_SARSCOV2_6h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df7 <- read.csv("PRJNA704763_ribo_virus_SARSCOV2_24h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df8 <- read.csv("PRJNA704763_ribo_virus_SARSCOV2_24h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df9 <- read.csv("PRJNA704763_ribo_virus_SARSCOV2_24h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

dfr1 <- read.csv("PRJNA704763_rna_mock_novirus_24h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr2 <- read.csv("PRJNA704763_rna_mock_novirus_24h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr3 <- read.csv("PRJNA704763_rna_mock_novirus_24h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr4 <- read.csv("PRJNA704763_rna_virus_SARSCOV2_6h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr5 <- read.csv("PRJNA704763_rna_virus_SARSCOV2_6h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr6 <- read.csv("PRJNA704763_rna_virus_SARSCOV2_6h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr7 <- read.csv("PRJNA704763_rna_virus_SARSCOV2_24h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr8 <- read.csv("PRJNA704763_rna_virus_SARSCOV2_24h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr9 <- read.csv("PRJNA704763_rna_virus_SARSCOV2_24h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

cts <- purrr::reduce(
  list(df1, df2, df3, df4, df5, df6, df7, df8, df9, dfr1, dfr2, dfr3, dfr4, dfr5, dfr6, dfr7, dfr8, dfr9),
  full_join,
  by = "N_unmapped")

#Not actual while using filterByExpr
#cts <- cts[rowSums(dplyr::select(cts,where(is.numeric)))>0,]
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
                   "ribo_mock_novirus_24h_1",
                   "ribo_mock_novirus_24h_2",
                   "ribo_mock_novirus_24h_3",
                   "ribo_virus_SARS-CoV-2_6h_1",
                   "ribo_virus_SARS-CoV-2_6h_2",
                   "ribo_virus_SARS-CoV-2_6h_3",
                   "ribo_virus_SARS-CoV-2_24h_1",
                   "ribo_virus_SARS-CoV-2_24h_2",
                   "ribo_virus_SARS-CoV-2_24h_3",
                   
                   "rna_mock_novirus_24h_1",
                   "rna_mock_novirus_24h_2",
                   "rna_mock_novirus_24h_3",
                   "rna_virus_SARS-CoV-2_6h_1",
                   "rna_virus_SARS-CoV-2_6h_2",
                   "rna_virus_SARS-CoV-2_6h_3",
                   "rna_virus_SARS-CoV-2_24h_1",
                   "rna_virus_SARS-CoV-2_24h_2",
                   "rna_virus_SARS-CoV-2_24h_3")

#w_24
d <- cts[,-c(6,7,8,15,16,17)] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_SC2_m_24h <- sum(d[1,])

#w_6
d <- cts[,c(c(3:8),c(12:17))] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_SC2_m_6h <- sum(d[1,])

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

time <- c(rep(c(rep(24, 3), rep(6, 3), rep(24, 3)),2))
time <- factor(time)
type <- c(rep('mock', 3), rep('sarscov2', 6), rep('mock', 3), rep('sarscov2', 6))
type <- factor(type)
reg <- c(rep('ribo', 9), rep('rna', 9))
reg <- factor(reg)
rep <- c(rep(c('1', '2', '3'), 6))
rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type

coldata <- data.frame(type=type,reg=reg,rep=rep,time=time)

#good filtering
keep <- filterByExpr(cts[,3:ncol(cts)], group = group[time != 6])
cts_f <- cts[keep,]

### COMMON PART
x <- cts_f[,c(3:ncol(cts_f))]

rownames(coldata) <- colnames(x)

rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")
#a <- ~reg + type + reg:type + reg:type:time
dds <- DESeq2::DESeqDataSetFromMatrix(countData = x[time != 6],
                                      colData = coldata[time != 6,], 
                                      design = ~reg + type + reg:type)

dds <- DESeq(dds, test = "LRT", reduced = ~reg + type)

dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)
#SVA

type <- c(rep('mock', 3), rep('sarscov2', 3), rep('mock', 3), rep('sarscov2', 3))
type <- factor(type)
reg <- c(rep('ribo', 6), rep('rna', 6))
reg <- factor(reg)
rep <- c(rep(c('1', '2', '3'), 4))
rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type

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
#no clusters

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

type <- c(rep('mock', 3), rep('sarscov2', 3))
type <- factor(type)
reg <- c(rep('ribo', 6))
reg <- factor(reg)
rep <- c(rep(c('1', '2', '3'), 2))
rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type

dds <- DESeq(dds)

dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ type, colData(dds))
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
  geom_point(aes(color=as.factor(type)))
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

resultsNames(dds) # lists the coefficients
res <- results(dds, filterFun=ihw)
res_la <- results(dds, name="regribo.typesarscov2", test="Wald", lfcThreshold=1, altHypothesis="lessAbs", filterFun=ihw)
res_l <- results(dds, name="regribo.typesarscov2", test="Wald", lfcThreshold=0, altHypothesis="less", filterFun=ihw)

head(res[order(res$pvalue), ], 10)
summary(res)

## Checking the normalization
epsilon <- 1
col.strain <- c("mock"="blue","sarscov2"="orange")

#RIBO
#colnames(dds.norm) <- c("ribo_mock_novirus_notime_1", "ribo_mock_novirus_notime_2", "ribo_virus_DENV1_48h_1", "ribo_virus_DENV1_48h_2", "rna_mock_novirus_notime_1", "rna_mock_novirus_notime_2", "rna_virus_DENV1_48h_1", "rna_virus_DENV1_48h_2")
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(cts[,3:ncol(cts)][time != 6][,1:6]+epsilon),  col=col.strain[as.vector(coldata$type[1:6])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),1:6]+epsilon),  col=col.strain[as.vector(coldata$type[1:6])], cex.axis=0.7, 
        las=1, xlab="log2(filtered & normalized counts)", horizontal=TRUE, main="filtered & normalized counts") 
plotDensity(log2(cts[,3:ncol(cts)][time != 6]+epsilon)[,1:6],  col=col.strain[as.vector(coldata$type[1:6])], lwd=2,
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)[!is.na(res_l@listData$padj),1:6]+epsilon), col=col.strain[as.vector(coldata$type[1:6])], lwd=2,
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

res_RO_SC2_m_24h <- as.data.frame(res)
res_RO_SC2_m_24h_la <- as.data.frame(res_la)
res_RO_SC2_m_24h_l <- as.data.frame(res_l)

write.table(res_RO_SC2_m_24h, str_c("./",'res_RO_ds2_24h.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_m_24h_la, str_c("./",'res_RO_ds2_24h_la.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_m_24h_l, str_c("./",'res_RO_ds2_24h_l.tsv'), sep="\t", quote = F, row.names = T)

par(mfrow=c(2,2),cex.lab=0.7)
hist(res_RO_DENV1_48_la$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald LessAbs")
hist(res_RO_DENV1_48_l$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald 1-sided")
hist(res_RO_DENV1_48$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="LRT 2-sided")

#6h
time <- c(rep(c(rep(24, 3), rep(6, 3), rep(24, 3)),2))
time <- factor(time)
type <- c(rep('mock', 3), rep('sarscov2', 6), rep('mock', 3), rep('sarscov2', 6))
type <- factor(type)
reg <- c(rep('ribo', 9), rep('rna', 9))
reg <- factor(reg)
rep <- c(rep(c('1', '2', '3'), 6))
rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type

dds <- DESeq2::DESeqDataSetFromMatrix(countData = x[(time != 24) | (type != "sarscov2")],
                                      colData = coldata[(time != 24) | (type != "sarscov2"),], 
                                      design = ~reg + type + reg:type)
dds <- DESeq(dds, test = "LRT", reduced = ~reg + type)
dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

resultsNames(dds) # lists the coefficients

res <- results(dds, filterFun=ihw)
res_la <- results(dds, name="regribo.typesarscov2", test="Wald", lfcThreshold=1, altHypothesis="lessAbs", filterFun=ihw)
res_l <- results(dds, name="regribo.typesarscov2", test="Wald", lfcThreshold=0, altHypothesis="less", filterFun=ihw)

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

res_RO_SC2_6 <- as.data.frame(res)
res_RO_SC2_6_la <- as.data.frame(res_la)
res_RO_SC2_6_l <- as.data.frame(res_l)

write.table(res_RO_SC2_6, str_c("./",'res_RO_ds2_6h.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_6_la, str_c("./",'res_RO_ds2_6h_la.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_6_l, str_c("./",'res_RO_ds2_6h_l.tsv'), sep="\t", quote = F, row.names = T)

par(mfrow=c(2,2),cex.lab=0.7)
hist(res_RO_DENV1_72_la$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald LessAbs")
hist(res_RO_DENV1_72_l$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald 1-sided")
hist(res_RO_DENV1_72$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="LRT 2-sided")