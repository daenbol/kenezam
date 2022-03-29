library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(sva)

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

ann <- data.frame(type=type,reg=reg,rep=rep)

x <- cts[,c(3:ncol(cts))]
rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")
y <- DGEList(counts=x,group=group)

keep <- filterByExpr(y, group=group)
y <- y[keep,,keep.lib.sizes=FALSE]

#glMDSPlot(y[,1:6], groups=group, folder = "glimma-plots/")
glimmaMDS(y[,1:7], groups=ann[1:7,], gene.selection = "common", folder="mds")
#glMDSPlot(y[,7:12], groups=group, folder = "glimma-plots/")
glimmaMDS(y[,8:13], groups=ann[8:13,], gene.selection = "common", folder="mds")
y <- calcNormFactors(y)

plot(density(log2(y$counts[,1]*1000000/(y$samples$lib.size[1]+1))), from = -2, to = 20, lwd = 3)
for (k in 2:7){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$lib.size[k]+1))),col=k)
}

plot(density(log2(y$counts[,8]*1000000/(y$samples$norm.factors[1]*y$samples$lib.size[7]+1))), from = -2, to = 20, lwd = 3)
for (k in 9:13){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$norm.factors[k]*y$samples$lib.size[k]+1))),col=k)
}

design <- model.matrix(~0 + group)
ncol(design) <= qr(design)$rank

#redesign
type <- c('mock', 'mock', 'MHVA59', 'MHVA59', 'mock', 'mock', 'MHVA59', 'MHVA59')
type <- factor(type)
reg <- c('ribo', 'ribo','ribo', 'ribo', 'rna','rna','rna','rna')
reg <- factor(reg)
rep <- c('1', '2', '1', '1', '1', '2', '1', '1', '2', '1', '1', '2', '1')
rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type

ym <- y[,-c(3,4,7,10,13)]
design <- model.matrix(~0 + group)
ncol(design) <= qr(design)$rank

ym <- estimateDisp(ym,design)

fit <- glmQLFit(ym,design)
#grouprna:wt grouprna:ko groupribo:wt groupribo:ko
#1. 48
qlf_RO_MHV <- glmQLFTest(fit, contrast = c(1, -1, -1, 1))
qlf_rna_MHV <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0))
qlf_ribo_MHV <- glmQLFTest(fit, contrast = c(0, 0, -1, 1))

res_RO_MHV <-  topTags(qlf_RO_MHV, n = Inf)[[1]]
res_rna_MHV <-  topTags(qlf_rna_MHV, n =Inf)[[1]]
res_ribo_MHV <-  topTags(qlf_ribo_MHV, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_MHV$row_names <- rownames(res_RO_MHV)
res_RO_MHV <- left_join(res_RO_MHV, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_MHV$PValue, covariates = rowMeans(dplyr::select(res_RO_MHV, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_MHV$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_MHV$PValue, method="BH")
sum(FDR < 0.1)

res_ribo_MHV$row_names <- rownames(res_ribo_MHV)
res_ribo_MHV <- left_join(res_ribo_MHV, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_ribo_MHV$PValue, covariates = rowMeans(dplyr::select(res_ribo_MHV, contains("ribo_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_ribo_MHV$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_ribo_MHV$PValue, method="BH")
sum(FDR < 0.1)

res_rna_MHV$row_names <- rownames(res_rna_MHV)
res_rna_MHV <- left_join(res_ribo_MHV, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_rna_MHV$PValue, covariates = rowMeans(dplyr::select(res_rna_MHV, contains("rna_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_rna_MHV$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_rna_MHV$PValue, method="BH")
sum(FDR < 0.1)

write.table(res_rna_MHV, str_c("./",'RNA_MHV.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_ribo_MHV, str_c("./",'Ribo_MHV.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_RO_MHV, str_c("./",'RO_MHV.tsv'), sep="\t", quote = F, row.names = F)