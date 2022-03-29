#Ok. Authors

library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(IHW)
library(tidyr)

#1. Calculated counts and extended counts
setwd("D://Projects/kenezam/gene_counts/GSE162323/")

df1 <- read.csv("GSE162323_ribo_virus_SARS-CoV-2_8h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df2 <- read.csv("GSE162323_ribo_virus_SARS-CoV-2_8h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df3 <- read.csv("GSE162323_ribo_mock_novirus_notime_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df4 <- read.csv("GSE162323_ribo_mock_novirus_notime_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

dfm1 <- read.csv("GSE162323_ribo_virus_SARS-CoV-2_3h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfm2 <- read.csv("GSE162323_ribo_virus_SARS-CoV-2_5h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]


dfr1 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_8h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr2 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_8h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr3 <- read.csv("GSE162323_rna_mock_novirus_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr4 <- read.csv("GSE162323_rna_mock_novirus_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

dfrm1 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_3h_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfrm2 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_3h_1_total_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfrm3 <- read.csv("GSE162323_rna_mock_novirus_1_total_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfrm4 <- read.csv("GSE162323_rna_mock_novirus_2_total_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfrm5 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_5h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfrm6 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_5h_1_total_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfrm7 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_8h_1_total_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfrm8 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_8h_2_total_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

cts_ex <- purrr::reduce(
  list(df1, df2, df3, df4, dfr1, dfr2, dfr3, dfr4, dfm1, dfm2, dfrm1, dfrm2, dfrm3, dfrm4, dfrm5, dfrm6, dfrm7, dfrm8),
  full_join,
  by = "N_unmapped")

cts <- purrr::reduce(
  list(df1, df2, df3, df4, dfr1, dfr2, dfr3, dfr4),
  full_join,
  by = "N_unmapped")

id = gsub("\\..*","",cts$N_unmapped)
cts <- cbind(id, cts)
cts_ex <- cbind(id, cts_ex)

library("AnnotationDbi")
library("org.Hs.eg.db")
cts$id <- AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                keys = cts$id,
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "list")
cts_ex$id <- AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                keys = cts_ex$id,
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "list")

colnames(cts) <- c("GENES_SYMBOL",
                   "ENSEMBL_ID",
                   "ribo_virus_sarscov2_8h_1",
                   "ribo_virus_sarscov2_8h_2",
                   "ribo_mock_novirus_8h_1",
                   "ribo_mock_novirus_8h_2",
                   "rna_virus_sarscov2_8h_1",
                   "rna_virus_sarscov2_8h_2",
                   "rna_mock_novirus_8h_1",
                   "rna_mock_novirus_8h_2")

colnames(cts_ex) <- c("GENES_SYMBOL",
                   "ENSEMBL_ID",
                   "ribo_virus_sarscov2_8h_1",
                   "ribo_virus_sarscov2_8h_2",
                   "ribo_mock_novirus_8h_1",
                   "ribo_mock_novirus_8h_2",
                   "ribo_virus_sarscov2_3h_1",
                   "ribo_virus_sarscov2_5h_2",
                   "rna_virus_sarscov2_8h_1",
                   "rna_virus_sarscov2_8h_2",
                   "rna_mock_novirus_8h_1",
                   "rna_mock_novirus_8h_2",
                   "rna_virus_sarscov2_3h_1",
                   "rna_virus_sarscov2_3h_t",
                   "rna_mock_novirus_notime_1_t",
                   "rna_mock_novirus_notime_2_t",
                   "rna_virus_sarscov2_5h",
                   "rna_virus_sarscov2_5h_t",
                   "rna_virus_sarscov2_8h_1_t",
                   "rna_virus_sarscov2_8h_2_t")

cts_ex$GENES_SYMBOL <- gsub('[c()"",]', "", cts_ex$GENES_SYMBOL)
cts_ex_rs = separate_rows(cts_ex,1,sep = " ")
cts_ex_rs_b <- cts_ex_rs

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

#2. paper counts
setwd("D://Projects/kenezam/gene_counts/GSE162323/Authors/")

cts_a <- read.csv("table.tsv", sep="\t")

library("AnnotationDbi")
library("org.Hs.eg.db")

id = gsub("\\..*","",cts_a$ucsc...ensembl.ID)

cts_a$id <- AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                  keys = id,
                                  column = "ENSEMBL",
                                  keytype = "UCSCKG",
                                  multiVals = "list")

#2.1
cts_a_ro <- cts_a[,c(1,2,19,20,15,16,17,18,7,8,3,4,5,11,9,10,6,12,13,14)]

colnames(cts_a_ro) <- c("GENES_SYMBOL",
                        "ENSEMBL_ID",
                        "ribo_virus_sarscov2_8h_1",
                        "ribo_virus_sarscov2_8h_2",
                        "ribo_mock_novirus_8h_1",
                        "ribo_mock_novirus_8h_2",
                        "ribo_virus_sarscov2_3h_1",
                        "ribo_virus_sarscov2_5h_2",
                        "rna_virus_sarscov2_8h_1",
                        "rna_virus_sarscov2_8h_2",
                        "rna_mock_novirus_8h_1",
                        "rna_mock_novirus_8h_2",
                        "rna_virus_sarscov2_3h_1",
                        "rna_virus_sarscov2_3h_t",
                        "rna_mock_novirus_notime_1_t",
                        "rna_mock_novirus_notime_2_t",
                        "rna_virus_sarscov2_5h",
                        "rna_virus_sarscov2_5h_t",
                        "rna_virus_sarscov2_8h_1_t",
                        "rna_virus_sarscov2_8h_2_t")

type <- c('sarscov2', 'sarscov2', 'mock', 'mock', 'sarscov2', 'sarscov2', 'sarscov2', 'sarscov2', 'mock', 'mock', 'sarscov2', 'sarscov2', 'mock', 'mock', 'sarscov2', 'sarscov2', 'sarscov2', 'sarscov2')
type <- factor(type)
reg <- c(rep('ribo', 6), rep('rna',12))
reg <- factor(reg)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group  = reg:type

cts <- cts_a_ro

#2.2
cts_a_ro_c <- cts_a_ro[,-c(c(7,8),c(13:20))]
cts_a_ro_c[is.na(cts_a_ro_c)] = 0
#OK
lapply(cts_a_ro_c[,3:ncol(cts_a_ro_c)],sum)

type <- c(rep(c('sarscov2', 'sarscov2', 'mock', 'mock'),2))
type <- factor(type)
reg <- c(rep('ribo', 4), rep('rna',4))
reg <- factor(reg)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group  = reg:type

cts <- cts_a_ro_c

#
d <- cts %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_SC2 <- sum(d[1,])

#good filtering
keep <- filterByExpr(cts[,3:ncol(cts)], group = group)
cts_f <- cts[keep,]

library("AnnotationDbi")
library("org.Hs.eg.db")
cts_f$id <- AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                keys = cts_f$id,
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "list")

coldata <- data.frame(type=type,reg=reg)

x <- cts_f[,c(3:ncol(cts_f))]

rownames(coldata) <- colnames(x)

rownames(x) <- paste(cts_f$GENES_SYMBOL, cts_f$ENSEMBL_ID, sep="__")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
      colData = coldata, 
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
col.strain <- c("mock"="blue","sarscov2"="orange")

#RIBO
colnames(dds.norm) <- c("ribo_sc2_1", "ribo_sc2_2", "ribo_mock_1", "ribo_mock_2", "rna_sc2_1", "rna_sc2_2", "rna_mock_1", "rna_mock_2")
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

res_RO_SC2 <- as.data.frame(res)
res_RO_SC2_la <- as.data.frame(res_la)
res_RO_SC2_l <- as.data.frame(res_l)

write.table(res_RO_SC2, str_c("./",'res_RO_ds2.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_la, str_c("./",'res_RO_ds2_la.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_l, str_c("./",'res_RO_ds2_l.tsv'), sep="\t", quote = F, row.names = T)

par(mfrow=c(2,2),cex.lab=0.7)
hist(res_RO_SC2_la$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald LessAbs")
hist(res_RO_SC2_l$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="Wald 1-sided")
hist(res$pvalue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="LRT 2-sided")