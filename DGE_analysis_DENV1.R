library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)

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

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

#Filtering using filterByExpr
type <- c('mock', 'mock', 'd1_48', 'd1_48', 'd1_72', 'd1_72', 'mock', 'mock', 'd1_48', 'd1_48', 'd1_72', 'd1_72')
type <- factor(type)
reg <- c('ribo', 'ribo','ribo', 'ribo', 'ribo', 'ribo', 'rna','rna','rna','rna','rna','rna')
reg <- factor(reg)
rep <- c(rep(c('1', '2'),6))
rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type

x <- cts[,c(3:ncol(cts))]
rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")
y <- DGEList(counts=x,group=group)

keep <- filterByExpr(y, group=group)
y <- y[keep,,keep.lib.sizes=FALSE]

glMDSPlot(y[,1:6], groups=group, folder = "glimma-plots/")
glMDSPlot(y[,7:12], groups=group, folder = "glimma-plots/")

y <- calcNormFactors(y)

plot(density(log2(y$counts[,1]*1000000/(y$samples$lib.size[1]+1))), from = -2, to = 20, lwd = 3)
for (k in 2:6){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$lib.size[k]+1))),col=k)
}

plot(density(log2(y$counts[,7]*1000000/(y$samples$norm.factors[1]*y$samples$lib.size[7]+1))), from = -2, to = 20, lwd = 3)
for (k in 8:12){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$norm.factors[k]*y$samples$lib.size[k]+1))),col=k)
}

design <- model.matrix(~0 + group)
y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)

#1. 48
qlf_RO_48 <- glmQLFTest(fit, contrast = c(1, -1, 0, -1, 1, 0))
qlf_rna_48 <- glmQLFTest(fit, contrast = c(c(-1, 1), rep(0, 4)))
qlf_ribo <- glmQLFTest(fit, contrast = c(rep(0,3), c(-1, 1), c(0)))

res_RO_48 <-  topTags(qlf_RO_48, n = Inf)[[1]]
res_rna_48 <-  topTags(qlf_rna_48, n =Inf)[[1]]
res_ribo_48 <-  topTags(qlf_ribo, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_48$row_names <- rownames(res_RO_48)
res_RO_48 <- left_join(res_RO_48, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_48$PValue, covariates = rowMeans(dplyr::select(res_RO_48, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_48$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_48$PValue, method="BH")
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