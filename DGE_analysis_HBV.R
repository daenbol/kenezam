library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(sva)

setwd("D://Projects/kenezam/gene_counts/GSE135860")


df1 <- read.csv("GSE135860_ribo_mock_novirus_72h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df2 <- read.csv("GSE135860_ribo_mock_novirus_72h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df3 <- read.csv("GSE135860_ribo_mock_novirus_72h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df4 <- read.csv("GSE135860_ribo_virus_HBV_72h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df5 <- read.csv("GSE135860_ribo_virus_HBV_72h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df6 <- read.csv("GSE135860_ribo_virus_HBV_72h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

dfr1 <- read.csv("GSE135860_rna_mock_novirus_72h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr2 <- read.csv("GSE135860_rna_mock_novirus_72h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr3 <- read.csv("GSE135860_rna_mock_novirus_72h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr4 <- read.csv("GSE135860_rna_virus_HBV_72h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr5 <- read.csv("GSE135860_rna_virus_HBV_72h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr6 <- read.csv("GSE135860_rna_virus_HBV_72h_3_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

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
                   "ribo_mock_novirus_72h_1",
                   "ribo_mock_novirus_72h_2",
                   "ribo_mock_novirus_72h_3",
                   "ribo_virus_HBV_72h_1",
                   "ribo_virus_HBV_72h_2",
                   "ribo_virus_HBV_72h_3",
                   "rna_mock_novirus_72h_1",
                   "rna_mock_novirus_72h_2",
                   "rna_mock_novirus_72h_3",
                   "rna_virus_HBV_72h_1",
                   "rna_virus_HBV_72h_2",
                   "rna_virus_HBV_72h_3")

#cts_t <- cts
#cts_t$GENES_SYMBOL <- as.factor(as.character(cts_t$GENES_SYMBOL))
  
#n_occur <- data.frame(base::table(cts_t$GENES_SYMBOL))
#n_occur[n_occur$Freq > 1,]

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

#Filtering using filterByExpr
type <- c('mock', 'mock', 'mock', 'HBV', 'HBV', 'HBV', 'mock', 'mock', 'mock', 'HBV', 'HBV', 'HBV')
#type <- c('mock', 'mock', 'mock', 'HBV', 'HBV', 'HBV', 'mock', 'mock', 'HBV', 'HBV', 'HBV')

type <- factor(type)
reg <- c('ribo', 'ribo','ribo', 'ribo', 'ribo', 'ribo', 'rna','rna','rna','rna','rna','rna')
#reg <- c('ribo', 'ribo','ribo', 'ribo', 'ribo', 'ribo', 'rna','rna','rna','rna','rna')

reg <- factor(reg)
rep <- c(rep(c('1', '2', '3'),4))
#rep <- c(rep(c('1', '2', '3'),2),c('1', '3', '1', '2', '3'))

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
glimmaMDS(y[,1:6], groups=ann[1:6,], gene.selection = "common", folder="mds")
#glMDSPlot(y[,7:12], groups=group, folder = "glimma-plots/")
glimmaMDS(y[,7:12], groups=ann[7:12,], gene.selection = "common", folder="mds")
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
ncol(design) <= qr(design)$rank
y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)
#grouprna:wt grouprna:ko groupribo:wt groupribo:ko
#1. HBV
qlf_RO_HBV <- glmQLFTest(fit, contrast = c(1, -1, -1, 1))
qlf_rna_HBV <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0))
qlf_ribo_HBV <- glmQLFTest(fit, contrast = c(0, 0, -1, 1))

res_RO_HBV <-  topTags(qlf_RO_HBV, n = Inf)[[1]]
res_rna_HBV <-  topTags(qlf_rna_HBV, n =Inf)[[1]]
res_ribo_HBV <-  topTags(qlf_ribo_HBV, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_HBV$row_names <- rownames(res_RO_HBV)
res_RO_HBV <- left_join(res_RO_HBV, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_HBV$PValue, covariates = rowMeans(dplyr::select(res_RO_HBV, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_HBV$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_HBV$PValue, method="BH")
sum(FDR < 0.1)

res_ribo_HBV$row_names <- rownames(res_ribo_HBV)
res_ribo_HBV <- left_join(res_ribo_HBV, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_ribo_HBV$PValue, covariates = rowMeans(dplyr::select(res_ribo_HBV, contains("ribo_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_ribo_HBV$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_ribo_HBV$PValue, method="BH")
sum(FDR < 0.1)

res_rna_HBV$row_names <- rownames(res_rna_HBV)
res_rna_HBV <- left_join(res_ribo_HBV, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_rna_HBV$PValue, covariates = rowMeans(dplyr::select(res_rna_HBV, contains("rna_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_rna_HBV$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_rna_HBV$PValue, method="BH")
sum(FDR < 0.1)

write.table(res_rna_HBV, str_c("./",'RNA_HBV.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_ribo_HBV, str_c("./",'Ribo_HBV.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_RO_HBV, str_c("./",'RO_HBV.tsv'), sep="\t", quote = F, row.names = F)

#2. HBV wo distinct sample
type <- c('mock', 'mock', 'mock', 'HBV', 'HBV', 'HBV', 'mock', 'mock', 'HBV', 'HBV', 'HBV')

type <- factor(type)

reg <- c('ribo', 'ribo','ribo', 'ribo', 'ribo', 'ribo', 'rna','rna','rna','rna','rna')

reg <- factor(reg)

rep <- c(rep(c('1', '2', '3'),2),c('1', '3', '1', '2', '3'))

rep <- factor(rep)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group = reg:type
ann <- data.frame(type=type,reg=reg,rep=rep)

x <- cts[,c(3:ncol(cts))][,-c(8)]
rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")
y <- DGEList(counts=x,group=group)

keep <- filterByExpr(y, group=group)
y <- y[keep,,keep.lib.sizes=FALSE]


qlf_RO_HBV_wo <- glmQLFTest(fit, contrast = c(1, -1, -1, 1))
qlf_rna_HBV_wo <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0))
qlf_ribo_HBV_wo <- glmQLFTest(fit, contrast = c(0, 0, -1, 1))

res_RO_HBV_wo <-  topTags(qlf_RO_HBV_wo, n = Inf)[[1]]
res_rna_HBV_wo <-  topTags(qlf_rna_HBV_wo, n =Inf)[[1]]
res_ribo_HBV_wo <-  topTags(qlf_ribo_HBV_wo, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_HBV_wo$row_names <- rownames(res_RO_HBV_wo)
res_RO_HBV_wo <- left_join(res_RO_HBV_wo, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_HBV_wo$PValue, covariates = rowMeans(dplyr::select(res_RO_HBV_wo, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_HBV_wo$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_HBV_wo$PValue, method="BH")
sum(FDR < 0.1)

res_ribo_HBV_wo$row_names <- rownames(res_ribo_HBV_wo)
res_ribo_HBV_wo <- left_join(res_ribo_HBV_wo, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_ribo_HBV_wo$PValue, covariates = rowMeans(dplyr::select(res_ribo_HBV_wo, contains("ribo_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_ribo_HBV_wo$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_ribo_HBV_wo$PValue, method="BH")
sum(FDR < 0.1)

res_rna_HBV_wo$row_names <- rownames(res_rna_HBV_wo)
res_rna_HBV_wo <- left_join(res_ribo_HBV_wo, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_rna_HBV_wo$PValue, covariates = rowMeans(dplyr::select(res_rna_HBV_wo, contains("rna_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_rna_HBV_wo$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_rna_HBV_wo$PValue, method="BH")
sum(FDR < 0.1)

write.table(res_rna_HBV_wo, str_c("./",'RNA_HBV_wo.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_ribo_HBV_wo, str_c("./",'Ribo_HBV_wo.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_RO_HBV_wo, str_c("./",'RO_HBV_wo.tsv'), sep="\t", quote = F, row.names = F)
