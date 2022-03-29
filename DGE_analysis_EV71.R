library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)

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

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

#time <- c(rep(c(1.75, rep(c(3.5, 6.25),2)), 4))
time <- c(rep(c(3.5, 3.5, 6.25, 6.25),4))
time <- factor(time)
#type <- c(rep('mock', 5), rep('EV71', 5), rep('mock', 5), rep('EV71', 5))
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
ann <- data.frame(type=type,reg=reg,rep=rep,time=time)
write.table(ann, str_c("./",'ann.tsv'), sep="\t", quote = F, row.names = F)

x <- cts[,c(3:ncol(cts))]
rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")
y <- DGEList(counts=x,group=group)

keep <- filterByExpr(y, group=group)
y <- y[keep,,keep.lib.sizes=FALSE]

#glMDSPlot(y[,1:6], groups=group, folder = "glimma-plots/")
glimmaMDS(y[,1:10], groups=ann[1:10,], gene.selection = "common", folder="mds")
#glMDSPlot(y[,7:12], groups=group, folder = "glimma-plots/")
glimmaMDS(y[,11:20], groups=ann[11:20,], gene.selection = "common", folder="mds")
y <- calcNormFactors(y)

plot(density(log2(y$counts[,1]*1000000/(y$samples$lib.size[1]+1))), from = -2, to = 20, lwd = 3)
for (k in 2:6){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$lib.size[k]+1))),col=k)
}

plot(density(log2(y$counts[,7]*1000000/(y$samples$norm.factors[1]*y$samples$lib.size[7]+1))), from = -2, to = 20, lwd = 3)
for (k in 8:12){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$norm.factors[k]*y$samples$lib.size[k]+1))),col=k)
}

design <- model.matrix(~0 + group + rep + group:time)
ncol(design) <= qr(design)$rank
y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)
#grouprna:wt grouprna:ko groupribo:wt groupribo:ko
#1. 35
qlf_RO_EV71_35 <- glmQLFTest(fit, contrast = c(1, -1, -1, 1, rep(0,5)))
qlf_rna_EV71_35 <- glmQLFTest(fit, contrast = c(c(-1, 1), rep(0, 7)))
qlf_ribo_EV71_35 <- glmQLFTest(fit, contrast = c(rep(0,2), c(-1, 1), rep(0,5)))

res_RO_EV71_35 <- topTags(qlf_RO_EV71_35, n = Inf)[[1]]
res_rna_EV71_35 <- topTags(qlf_rna_EV71_35, n =Inf)[[1]]
res_ribo_EV71_35 <- topTags(qlf_ribo_EV71_35, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_EV71_35$row_names <- rownames(res_RO_EV71_35)
res_RO_EV71_35 <- left_join(res_RO_EV71_35, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_EV71_35$PValue, covariates = rowMeans(dplyr::select(res_RO_EV71_35, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_EV71_35$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_EV71_35$PValue, method="BH")
sum(FDR < 0.1)

res_ribo_48$row_names <- rownames(res_ribo_48)
res_ribo_48 <- left_join(res_ribo_48, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_ribo_48$PValue, covariates = rowMeans(dplyr::select(res_ribo_48, contains("ribo_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_ribo_48$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_ribo_48$PValue, method="BH")
sum(FDR < 0.1)

#2. 625
qlf_RO_EV71_625 <- glmQLFTest(fit, contrast = c(1, -1, -1, 1, 0, 1, -1, -1, 1))
qlf_rna_EV71_625 <- glmQLFTest(fit, contrast = c(c(-1, 1), rep(0, 3), c(-1, 1), c(0, 0)))
qlf_ribo_EV71_625 <- glmQLFTest(fit, contrast = c(rep(0,2), c(-1, 1), rep(0,3), c(-1, 1)))

res_RO_EV71_625 <- topTags(qlf_RO_EV71_625, n = Inf)[[1]]
res_rna_EV71_625 <- topTags(qlf_rna_EV71_625, n =Inf)[[1]]
res_ribo_EV71_625 <- topTags(qlf_ribo_EV71_625, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_EV71_625$row_names <- rownames(res_RO_EV71_625)
res_RO_EV71_625 <- left_join(res_RO_EV71_625, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_EV71_625$PValue, covariates = rowMeans(dplyr::select(res_RO_EV71_625, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_EV71_625$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_EV71_625$PValue, method="BH")
sum(FDR < 0.1)

res_ribo_48$row_names <- rownames(res_ribo_48)
res_ribo_48 <- left_join(res_ribo_48, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_ribo_48$PValue, covariates = rowMeans(dplyr::select(res_ribo_48, contains("ribo_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_ribo_48$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_ribo_48$PValue, method="BH")
sum(FDR < 0.1)

res_rna_48$row_names <- rownames(res_rna_48)
res_rna_48 <- left_join(res_ribo_48, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_rna_48$PValue, covariates = rowMeans(dplyr::select(res_rna_48, contains("rna_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_rna_48$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_rna_48$PValue, method="BH")
sum(FDR < 0.1)

write.table(res_rna_48, str_c("./",'RNA_48.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_ribo_48, str_c("./",'Ribo_48.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_RO_48, str_c("./",'RO_48.tsv'), sep="\t", quote = F, row.names = F)

#2. 72
qlf_RO_72 <- glmQLFTest(fit, contrast = c(1, 0, -1, -1, 0, 1))
qlf_rna_72 <- glmQLFTest(fit, contrast = c(-1, 0, 1, 0, 0, 0))
qlf_ribo_72 <- glmQLFTest(fit, contrast = c(0, 0, 0, -1, 0, 1))

res_RO_72 <-  topTags(qlf_RO_72, n = Inf)[[1]]
res_rna_72 <-  topTags(qlf_rna_72, n =Inf)[[1]]
res_ribo_72 <-  topTags(qlf_ribo_72, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_72$row_names <- rownames(res_RO_72)
res_RO_72 <- left_join(res_RO_72, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_72$PValue, covariates = rowMeans(dplyr::select(res_RO_72, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_72$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_72$PValue, method="BH")
sum(FDR < 0.1)

res_ribo_72$row_names <- rownames(res_ribo_72)
res_ribo_72 <- left_join(res_ribo_72, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_ribo_72$PValue, covariates = rowMeans(dplyr::select(res_ribo_72, contains("ribo_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_ribo_72$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_ribo_72$PValue, method="BH")
sum(FDR < 0.1)

res_rna_72$row_names <- rownames(res_rna_72)
res_rna_72 <- left_join(res_ribo_72, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_rna_72$PValue, covariates = rowMeans(dplyr::select(res_rna_72, contains("rna_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_rna_72$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_rna_72$PValue, method="BH")
sum(FDR < 0.1)

write.table(res_rna_72, str_c("./",'RNA_72.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_ribo_72, str_c("./",'Ribo_72.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_RO_72, str_c("./",'RO_72.tsv'), sep="\t", quote = F, row.names = F)