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

setwd("D://Projects/kenezam/gene_counts/GSE158930/")

df1 <- read.csv("GSE158930_ribo_mock_novirus_4h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df2 <- read.csv("GSE158930_ribo_mock_novirus_4h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df3 <- read.csv("GSE158930_ribo_mock_novirus_4h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df4 <- read.csv("GSE158930_ribo_mock_novirus_96h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df5 <- read.csv("GSE158930_ribo_mock_novirus_96h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df6 <- read.csv("GSE158930_ribo_mock_novirus_96h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df7 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_4h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df8 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_4h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df9 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_4h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df10 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_24h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df11 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_24h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df12 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_24h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df13 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_24h_4_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df14 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_48h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df15 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_48h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df16 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_48h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df17 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_48h_4_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df18 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_72h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df19 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_72h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df20 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_72h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df21 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_72h_4_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df22 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_96h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df23 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_96h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df24 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_96h_3_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df25 <- read.csv("GSE158930_ribo_virus_SARS-CoV-2_96h_4_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

dfr1 <- read.csv("GSE158930_rna_mock_novirus_4h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr2 <- read.csv("GSE158930_rna_mock_novirus_4h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr3 <- read.csv("GSE158930_rna_mock_novirus_4h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr4 <- read.csv("GSE158930_rna_mock_novirus_96h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr5 <- read.csv("GSE158930_rna_mock_novirus_96h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr6 <- read.csv("GSE158930_rna_mock_novirus_96h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr7 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_4h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr8 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_4h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr9 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_4h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr10 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_24h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr11 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_24h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr12 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_24h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr13 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_24h_4_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr14 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_48h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr15 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_48h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr16 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_48h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr17 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_48h_4_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr18 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_72h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr19 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_72h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr20 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_72h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr21 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_72h_4_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr22 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_96h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr23 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_96h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr24 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_96h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr25 <- read.csv("GSE158930_rna_virus_SARS-CoV-2_96h_4_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

cts <- purrr::reduce(
  list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18, df19, df20, df21, df22, df23, df24, df25, dfr1, dfr2, dfr3, dfr4, dfr5, dfr6, dfr7, dfr8, dfr9, dfr10, dfr11, dfr12, dfr13, dfr14, dfr15, dfr16, dfr17, dfr18, dfr19, dfr20, dfr21, dfr22, dfr23, dfr24, dfr25),
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
                   "ribo_mock_novirus_4h_1",
                   "ribo_mock_novirus_4h_2",
                   "ribo_mock_novirus_4h_3",
                   "ribo_mock_novirus_96h_1",
                   "ribo_mock_novirus_96h_2",
                   "ribo_mock_novirus_96h_3",
                   "ribo_virus_SARS-CoV-2_4h_1",
                   "ribo_virus_SARS-CoV-2_4h_2",
                   "ribo_virus_SARS-CoV-2_4h_3",
                   "ribo_virus_SARS-CoV-2_24h_1",
                   "ribo_virus_SARS-CoV-2_24h_2",
                   "ribo_virus_SARS-CoV-2_24h_3",
                   "ribo_virus_SARS-CoV-2_24h_4",
                   "ribo_virus_SARS-CoV-2_48h_1",
                   "ribo_virus_SARS-CoV-2_48h_2",
                   "ribo_virus_SARS-CoV-2_48h_3",
                   "ribo_virus_SARS-CoV-2_48h_4",
                   "ribo_virus_SARS-CoV-2_72h_1",
                   "ribo_virus_SARS-CoV-2_72h_2",
                   "ribo_virus_SARS-CoV-2_72h_3",
                   "ribo_virus_SARS-CoV-2_72h_4",
                   "ribo_virus_SARS-CoV-2_96h_1",
                   "ribo_virus_SARS-CoV-2_96h_2",
                   "ribo_virus_SARS-CoV-2_96h_3",
                   "ribo_virus_SARS-CoV-2_96h_4",
                   "rna_mock_novirus_4h_1",
                   "rna_mock_novirus_4h_2",
                   "rna_mock_novirus_4h_3",
                   "rna_mock_novirus_96h_1",
                   "rna_mock_novirus_96h_2",
                   "rna_mock_novirus_96h_3",
                   "rna_virus_SARS-CoV-2_4h_1",
                   "rna_virus_SARS-CoV-2_4h_2",
                   "rna_virus_SARS-CoV-2_4h_3",
                   "rna_virus_SARS-CoV-2_24h_1",
                   "rna_virus_SARS-CoV-2_24h_2",
                   "rna_virus_SARS-CoV-2_24h_3",
                   "rna_virus_SARS-CoV-2_24h_4",
                   "rna_virus_SARS-CoV-2_48h_1",
                   "rna_virus_SARS-CoV-2_48h_2",
                   "rna_virus_SARS-CoV-2_48h_3",
                   "rna_virus_SARS-CoV-2_48h_4",
                   "rna_virus_SARS-CoV-2_72h_1",
                   "rna_virus_SARS-CoV-2_72h_2",
                   "rna_virus_SARS-CoV-2_72h_3",
                   "rna_virus_SARS-CoV-2_72h_4",
                   "rna_virus_SARS-CoV-2_96h_1",
                   "rna_virus_SARS-CoV-2_96h_2",
                   "rna_virus_SARS-CoV-2_96h_3",
                   "rna_virus_SARS-CoV-2_96h_4")

#w
d <- cts %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_SC2_b <- sum(d[1,])

#w
d <- cts[,c(c(3:5),c(9:11),c(28:30),c(34:36))] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_SC2_b_4h <- sum(d[1,])

#w
d <- cts[,c(c(6:8),c(24:27),c(31:33),c((ncol(cts)-3):ncol(cts)))] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_SC2_b_96h <- sum(d[1,])

#w
d <- cts[,c(c(3:5),c(12:15),c(28:30),c(37:40))] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_SC2_b_24h <- sum(d[1,])

#w
d <- cts[,c(c(3:5),c(16:19),c(28:30),c(41:44))] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_SC2_b_48h <- sum(d[1,])

#w
d <- cts[,c(c(6:8),c(20:23),c(31:33),c(45:48))] %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
w_SC2_b_96h <- sum(d[1,])

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

#1st attempt
#time <- c(rep(c(1.75, rep(c(3.5, 6.25),2)), 4))
time <- c(rep(c(rep(4, 3), rep(96, 3), rep(4, 3), rep(24, 4), rep(48, 4), rep(72, 4), rep(96, 4)),2))
time <- factor(time)
#type <- c(rep('mock', 5), rep('EV71', 5), rep('mock', 5), rep('EV71', 5))
type <- c(rep('mock', 6), rep('sarscov2', 19), rep('mock', 6), rep('sarscov2', 19))
type <- factor(type)
reg <- c(rep('ribo', 25), rep('rna', 25))
reg <- factor(reg)
rep <- c(rep(c(rep(c('1', '2', '3'), 3), rep(c('1', '2', '3', '4'), 4)),2))
rep <- factor(rep)
#time <- relevel(time, ref = 3.5)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
batch <- c(1,0)[1 + (rep %in% c(1,2))]
batch <- factor(batch)
batch <- relevel(batch, ref = '0')
group = reg:type
coldata <- data.frame(type=type,reg=reg,time=time)
#coldata_c <- coldata[-c(c(10:21), c(35:46)),]
coldata_c <- coldata[time == 4,]
coldata_a <- data.frame(group=group,time=time)

### COMMON PART
#good filtering
keep <- filterByExpr(cts[,3:ncol(cts)], group = group)
cts_f <- cts[keep,]

x <- cts_f[,c(3:ncol(cts_f))]

rownames(coldata) <- colnames(x)

rownames(x) <- paste(cts$GENES_SYMBOL[keep], cts$ENSEMBL_ID[keep], sep="__")
#a <- ~reg + type + reg:type + reg:type:time
a <- factor(c("ribo:mock:4", "ribo:mock:4", "ribo:mock:4", "ribo:mock:96", "ribo:mock:96", "ribo:mock:96", "ribo:sarscov2:4", "ribo:sarscov2:4", "ribo:sarscov2:4", "ribo:sarscov2:96", "ribo:sarscov2:96", "ribo:sarscov2:96", "ribo:sarscov2:96", "rna:mock:4", "rna:mock:4", "rna:mock:4", "rna:mock:96", "rna:mock:96", "rna:mock:4", "rna:mock:4", "rna:mock:4", "rna:mock:96", "rna:mock:96", "rna:mock:96", "rna:sarscov2:4", "rna:sarscov2:4", "rna:sarscov2:4", "rna:sarscov2:96", "rna:sarscov2:96", "rna:sarscov2:96", "rna:sarscov2:96", "ribo:mock:24", "ribo:mock:24", "ribo:mock:24", "ribo:mock:48", "ribo:mock:48", "ribo:mock:48", "ribo:mock:72", "ribo:mock:72", "ribo:mock:72", "rna:mock:24 rna:mock:24", "rna:mock:24", "rna:mock:48", "rna:mock:48", "rna:mock:48", "rna:mock:72", "rna:mock:72", "rna:mock:72"))
aa <- reg[-c(c(10:21), c(35:46))]
bb <- type[-c(c(10:21), c(35:46))]
cc <- (reg:type)[-c(c(10:21), c(35:46))]
dd <- (reg:type:time)[-c(c(10:21), c(35:46))]
trace(DESeqDataSet, edit = T)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = x[,(time == 4 | time == 96)],
                                      colData = coldata[(time == 4 | time == 96),c(1,2)], 
                                      design = ~group + group:time)

#z 4h
reg <- reg[time == 4]
type <- type[time == 4]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = x[,(time == 4)],
                                      colData = coldata[(time == 4),c(1,2)], 
                                      design = ~reg + type + reg:type)

dds <- DESeq(dds, test = "LRT", reduced = ~reg + type)
dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

#SVA
dat  <- counts(dds.norm, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat_m  <- dat[idx, ]
mod  <- model.matrix(~ as.factor(reg+type), colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
#x <- tibble::column_to_rownames(cts[2:ncol(x)], 'ENSEMBL_ID')
#xx <- log2(x + 2)
apply(dat_m, 2, function(x) any(is.na(x)))
n.sv = num.sv(as.matrix(dat_m),mod,method="leek")
n.sv
svseq <- svaseq(dat_m, mod, mod0, n.sv = 2)

sv.pca <- prcomp(svseq$sv,
                 center = TRUE,
                 scale. = TRUE)
library(ggfortify)
a <- cbind(svseq$sv, as.character(reg:type))
colnames(a) <- c("SV1", "SV2", "group")
sv.pca.plot <- autoplot(sv.pca,
                        data = a,
                        colour = "group")

sv.pca.plot

#z all
time <- c(rep(c(rep(4, 3), rep(96, 3), rep(4, 3), rep(24, 4), rep(48, 4), rep(72, 4), rep(96, 4)),2))
time <- factor(time)
#type <- c(rep('mock', 5), rep('EV71', 5), rep('mock', 5), rep('EV71', 5))
type <- c(rep('mock', 6), rep('sarscov2', 19), rep('mock', 6), rep('sarscov2', 19))
type <- factor(type)
reg <- c(rep('ribo', 25), rep('rna', 25))
reg <- factor(reg)
rep <- c(rep(c(rep(c('1', '2', '3'), 3), rep(c('1', '2', '3', '4'), 4)),2))
rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
batch <- c(1,0)[1 + (rep %in% c(1,2))]
batch <- factor(batch)
batch <- relevel(batch, ref = '0')
group = reg:type
coldata <- data.frame(type=type,reg=reg,time=time)
coldata_c <- coldata[time == 4,]
coldata_a <- data.frame(group=group,time=time)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                      colData = coldata[,c(1,2)], 
                                      design = ~reg + type + reg:type)

dds <- DESeq(dds, test = "LRT", reduced = ~reg + type)
dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

#SVA
dat  <- counts(dds.norm, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat_m  <- dat[idx, ]
mod  <- model.matrix(~ reg+type+time, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
apply(dat_m, 2, function(x) any(is.na(x)))
n.sv = num.sv(as.matrix(dat_m),mod,method="leek")
n.sv
svseq <- svaseq(dat_m, mod, mod0, n.sv = 1)

sv.pca <- prcomp(svseq$sv,
                 center = TRUE,
                 scale. = TRUE)
library(ggfortify)
a <- cbind(svseq$sv, as.character(rep), as.character(rep))
colnames(a) <- c("SV1", "SV2", "group", "rep")
sv.pca.plot <- autoplot(sv.pca,
                        data = a,
                        colour = "group")

sv.pca.plot

coldata <- data.frame(type=type,reg=reg,sv1=a[,1],sv2=a[,2])

sv <- as.numeric(a[,1])
coldata <- data.frame(type=type,reg=reg,sv=sv)

sv1 <- as.numeric(coldata$sv1)
sv2 <- as.numeric(coldata$sv2)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                      colData = coldata, 
                                      design = ~reg + type + sv + reg:type)

dds <- DESeq(dds, test = "LRT", reduced = ~reg + type + sv)

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

res_RO_SC2_b_4h <- as.data.frame(res)
res_RO_SC2_b_4h_la <- as.data.frame(res_la)
res_RO_SC2_b_4h_l <- as.data.frame(res_l)

write.table(res_RO_SC2_b_4h, str_c("./",'res_RO_ds2_4h.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_b_4h_la, str_c("./",'res_RO_ds2_4h_la.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_b_4h_l, str_c("./",'res_RO_ds2_4h_l.tsv'), sep="\t", quote = F, row.names = T)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = x[,time == 96],
                                      colData = coldata[time == 96,c(1,2,5)], 
                                      design = ~reg + type + batch + reg:type)

dds <- DESeq(dds, test = "LRT", reduced = ~reg + type + batch)

#SVA
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ reg:type, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
#x <- tibble::column_to_rownames(cts[2:ncol(x)], 'ENSEMBL_ID')
#xx <- log2(x + 2)
n.sv = num.sv(as.matrix(dat),mod,method="leek")
n.sv
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

sv.pca <- prcomp(svseq$sv,
                 center = TRUE,
                 scale. = TRUE)
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

res_RO_SC2_b_96h <- as.data.frame(res)
res_RO_SC2_b_96h_la <- as.data.frame(res_la)
res_RO_SC2_b_96h_l <- as.data.frame(res_l)

write.table(res_RO_SC2_b_96h, str_c("./",'res_RO_ds2_96h.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_b_96h_la, str_c("./",'res_RO_ds2_96h_la.tsv'), sep="\t", quote = F, row.names = T)
write.table(res_RO_SC2_b_96h_l, str_c("./",'res_RO_ds2_96h_l.tsv'), sep="\t", quote = F, row.names = T)
