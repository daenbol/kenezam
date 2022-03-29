#Hundreds of thousands of ribo counts

library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(sva)

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
  #list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, dfr1, dfr2, dfr3, dfr4, dfr5, dfr6, dfr7, dfr8, dfr9, dfr10),
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
ann <- data.frame(type=type,reg=reg,rep=rep,time=time)

x <- cts[,c(3:ncol(cts))]
rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")
y <- DGEList(counts=x,group=group)

keep <- filterByExpr(y, group=group)
y <- y[keep,,keep.lib.sizes=FALSE]

#glMDSPlot(y[,1:6], groups=group, folder = "glimma-plots/")
glimmaMDS(y[,1:9], groups=ann[1:9,], gene.selection = "common", folder="mds")
#glMDSPlot(y[,7:12], groups=group, folder = "glimma-plots/")
glimmaMDS(y[,10:18], groups=ann[10:18,], gene.selection = "common", folder="mds")
y <- calcNormFactors(y)

plot(density(log2(y$counts[,1]*1000000/(y$samples$lib.size[1]+1))), from = -2, to = 20, lwd = 3)
for (k in 2:9){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$lib.size[k]+1))),col=k)
}

plot(density(log2(y$counts[,7]*1000000/(y$samples$norm.factors[1]*y$samples$lib.size[7]+1))), from = -2, to = 20, lwd = 3)
for (k in 10:18){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$norm.factors[k]*y$samples$lib.size[k]+1))),col=k)
}

design <- model.matrix(~0 + group + group:time)
design_c <- design[,-c(5,7,9,11,13,15)]
ncol(design_c) <= qr(design_c)$rank
y <- estimateDisp(y,design_c)

fit <- glmQLFit(y,design_c)
#grouprna:wt grouprna:ko groupribo:wt groupribo:ko
#1. 4
qlf_RO_sarscov2_15_4 <- glmQLFTest(fit, contrast = c(1, -1, -1, 1, rep(0,10)))
qlf_rna_sarscov2_15_4 <- glmQLFTest(fit, contrast = c(c(-1, 1), rep(0, 12)))
qlf_ribo_sarscov2_15_4 <- glmQLFTest(fit, contrast = c(rep(0,2), c(-1, 1), rep(0,10)))

res_RO_sarscov2_15_4 <- topTags(qlf_RO_sarscov2_4, n = Inf)[[1]]
res_rna_sarscov2_15_4 <- topTags(qlf_rna_sarscov2_4, n =Inf)[[1]]
res_ribo_sarscov2_15_4 <- topTags(qlf_ribo_sarscov2_4, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_sarscov2_4$row_names <- rownames(res_RO_sarscov2_4)
res_RO_sarscov2_4 <- left_join(res_RO_sarscov2_4, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_sarscov2_4$PValue, covariates = rowMeans(dplyr::select(res_RO_sarscov2_4, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_sarscov2_4$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_sarscov2_4$PValue, method="BH")
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

#2. 24
qlf_RO_sarscov2_24 <- glmQLFTest(fit, contrast = c(1, -1, -1, 1, -1, 1, rep(0, 4), 0.217, -0.217, -0.217, 0.217))
qlf_rna_sarscov2_24 <- glmQLFTest(fit, contrast = c(c(-1, 1), rep(0, 2), 1, rep(0, 5), c(-0.217, 0.217), rep(0,2)))
qlf_ribo_sarscov2_24 <- glmQLFTest(fit, contrast = c(rep(0, 2), c(-1, 1), 0, 1, rep(0, 4), rep(0,2), c(-0.217, 0.217)))

res_RO_sarscov2_24 <- topTags(qlf_RO_sarscov2_24, n = Inf)[[1]]
res_rna_sarscov2_24 <- topTags(qlf_rna_sarscov2_24, n =Inf)[[1]]
res_ribo_sarscov2_24 <- topTags(qlf_ribo_sarscov2_24, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_sarscov2_24$row_names <- rownames(res_RO_sarscov2_24)
res_RO_sarscov2_24 <- left_join(res_RO_sarscov2_24, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_sarscov2_24$PValue, covariates = rowMeans(dplyr::select(res_RO_sarscov2_24, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_sarscov2_24$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_sarscov2_24$PValue, method="BH")
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

#3. 48
qlf_RO_sarscov2_48 <- glmQLFTest(fit, contrast = c(1, 0, -1, -1, 0, 1))
qlf_rna_sarscov2_48 <- glmQLFTest(fit, contrast = c(-1, 0, 1, 0, 0, 0))
qlf_ribo_sarscov2_48 <- glmQLFTest(fit, contrast = c(0, 0, 0, -1, 0, 1))

res_RO_sarscov2_48 <-  topTags(qlf_RO_sarscov2_48, n = Inf)[[1]]
res_rna_sarscov2_48 <-  topTags(qlf_rna_sarscov2_48, n =Inf)[[1]]
res_ribo_sarscov2_48 <-  topTags(qlf_ribo_sarscov2_48, n =Inf)[[1]]

#Another join
x$row_names <- rownames(x)

library("IHW")
res_RO_sarscov2_48$row_names <- rownames(res_RO_sarscov2_48)
res_RO_sarscov2_48 <- left_join(res_RO_sarscov2_48, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_RO_sarscov2_48$PValue, covariates = rowMeans(dplyr::select(res_RO_sarscov2_48, contains("virus"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_RO_sarscov2_48$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_RO_sarscov2_48$PValue, method="BH")
sum(FDR < 0.1)

res_ribo_sarscov2_48$row_names <- rownames(res_ribo_sarscov2_48)
res_ribo_sarscov2_48 <- left_join(res_ribo_sarscov2_48, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_ribo_sarscov2_48$PValue, covariates = rowMeans(dplyr::select(res_ribo_sarscov2_48, contains("ribo_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_ribo_sarscov2_48$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_ribo_sarscov2_48$PValue, method="BH")
sum(FDR < 0.1)

res_rna_sarscov2_48$row_names <- rownames(res_rna_sarscov2_48)
res_rna_sarscov2_48 <- left_join(res_ribo_sarscov2_48, x, by="row_names")
#IHW
ihwRes <- ihw(pvalues = res_rna_sarscov2_48$PValue, covariates = rowMeans(dplyr::select(res_rna_sarscov2_48, contains("rna_"))), alpha = 0.1)
adj_pvalues(ihwRes)
res_rna_sarscov2_48$adj_palues <- adj_pvalues(ihwRes)
count(adj_pvalues(ihwRes) < 0.05)
#BH
FDR <- p.adjust(res_rna_sarscov2_48$PValue, method="BH")
sum(FDR < 0.1)

write.table(res_rna_72, str_c("./",'RNA_72.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_ribo_72, str_c("./",'Ribo_72.tsv'), sep="\t", quote = F, row.names = F)
write.table(res_RO_72, str_c("./",'RO_72.tsv'), sep="\t", quote = F, row.names = F)