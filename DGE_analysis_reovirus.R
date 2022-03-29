library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)

library(stringr)
library(Glimma)

setwd("D://Projects/kenezam/gene_counts/GSE137757")

df1 <- read.csv("GSE137757_ribo_virus_reovirus_18h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df2 <- read.csv("GSE137757_ribo_virus_reovirus_18h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df3 <- read.csv("GSE137757_ribo_mock_novirus_18h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
df4 <- read.csv("GSE137757_ribo_mock_novirus_18h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

dfr1 <- read.csv("GSE137757_rna_virus_reovirus_18h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr2 <- read.csv("GSE137757_rna_virus_reovirus_18h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr3 <- read.csv("GSE137757_rna_mock_novirus_18h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]
dfr4 <- read.csv("GSE137757_rna_mock_novirus_18h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),1:2]

cts <- purrr::reduce(
  list(df1, df2, df3, df4, dfr1, dfr2, dfr3, dfr4),
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
                   "ribo_virus_reovirus_18h_1",
                   "ribo_virus_reovirus_18h_2",
                   "ribo_mock_novirus_18h_1",
                   "ribo_mock_novirus_18h_2",
                   "rna_virus_reovirus_18h_1",
                   "rna_virus_reovirus_18h_2",
                   "rna_mock_novirus_18h_1",
                   "rna_mock_novirus_18h_2")

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

#Filtering using filterByExpr
type <- c('reovirus', 'reovirus', 'mock', 'mock', 'reovirus', 'reovirus', 'mock', 'mock')
type <- factor(type)
reg <- c('ribo', 'ribo','ribo', 'ribo', 'rna','rna','rna','rna')
reg <- factor(reg)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group  = reg:type

x <- cts[,c(3:ncol(cts))]
rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")
y <- DGEList(counts=x,group=group)

keep <- filterByExpr(y, group=group)
y <- y[keep,,keep.lib.sizes=FALSE]

glMDSPlot(y[,1:4], groups=group, folder = "glimma-plots/after_calc_norm_f", gene.selection = "common")

glMDSPlot(y[,5:8], groups=group, folder = "glimma-plots/after_calc_norm_f", gene.selection = "common")

y <- calcNormFactors(y)

plot(density(log2(y$counts[,1]*1000000/(y$samples$norm.factors[1]*y$samples$lib.size[1]+1))), from = -2, to = 20, lwd = 3)
for (k in 1:4){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$norm.factors[k]*y$samples$lib.size[k]+1))),col=k)
}

design <- model.matrix(~0 + group)
y <- estimateDisp(y,design)

y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)

qlf_RO <- glmQLFTest(fit, contrast = c(1, -1, -1, 1))
qlf_rna <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0))
qlf_ribo <- glmQLFTest(fit, contrast = c(0, 0, -1, 1))

res_RO <-  topTags(qlf_RO, n = Inf)[[1]]
res_rna <-  topTags(qlf_rna, n =Inf)[[1]]
res_ribo <-  topTags(qlf_ribo, n =Inf)[[1]]

res_RO_b <- res_RO
res_rna_b <- res_rna
res_ribo_b <- res_ribo

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO$row_names <- rownames(res_RO)
res_RO <- left_join(res_RO, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO$PValue, covariates = rowMeans(dplyr::select(res_RO, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO$PValue, method="BH")
sum(FDR < 0.1)

res_ribo$row_names <- rownames(res_ribo)
res_ribo <- left_join(res_ribo, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_ribo$PValue, covariates = rowMeans(dplyr::select(res_ribo, contains("ribo_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_ribo$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_ribo$PValue, method="BH")
sum(FDR < 0.1)

res_rna$row_names <- rownames(res_rna)
res_rna <- left_join(res_ribo, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_rna$PValue, covariates = rowMeans(dplyr::select(res_rna, contains("rna_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_rna$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_rna$PValue, method="BH")
sum(FDR < 0.1)

write.table(res_rna, str_c("./",'RNA.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_ribo, str_c("./",'Ribo.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_RO, str_c("./",'RO.tsv'), sep="\t", quote = F, row.names = F)

res_rna_reo <- res_rna
res_ribo_reo <- res_ribo
res_RO_reo <- res_RO

hist(res_RO$PValue, freq = FALSE, col="gray", breaks = 50, horizontal=TRUE, main="EdgeR")
