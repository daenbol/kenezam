#Ok. Authors

library(purrr)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)
library(stringr)
library(Glimma)
library(sva)

#1. very low ribo. I will use paper counts
setwd("D://Projects/kenezam/gene_counts/GSE162323/")
df1 <- read.csv("GSE162323_ribo_virus_SARS-CoV-2_8h_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df2 <- read.csv("GSE162323_ribo_virus_SARS-CoV-2_8h_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df3 <- read.csv("GSE162323_ribo_mock_novirus_notime_1_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
df4 <- read.csv("GSE162323_ribo_mock_novirus_notime_2_5s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

dfr1 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_8h_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr2 <- read.csv("GSE162323_rna_virus_SARS-CoV-2_8h_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr3 <- read.csv("GSE162323_rna_mock_novirus_1_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]
dfr4 <- read.csv("GSE162323_rna_mock_novirus_2_4s/ReadsPerGene.out.tab", sep="\t")[-c(1:3),c(1,3)]

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
                   "ribo_virus_sarscov2_18h_1",
                   "ribo_virus_sarscov2_18h_2",
                   "ribo_mock_novirus_18h_1",
                   "ribo_mock_novirus_18h_2",
                   "rna_virus_sarscov2_18h_1",
                   "rna_virus_sarscov2_18h_2",
                   "rna_mock_novirus_18h_1",
                   "rna_mock_novirus_18h_2")

cts$GENES_SYMBOL <- gsub('[c()"",]', "", cts$GENES_SYMBOL)

cts_rna <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("rna_")))
cts_ribo <- cbind(cts[,c(1:2)],dplyr::select(cts, contains("ribo_")))

#Filtering using filterByExpr
type <- c('sarscov2', 'sarscov2', 'mock', 'mock', 'sarscov2', 'sarscov2', 'mock', 'mock')

type <- factor(type)
reg <- c('ribo', 'ribo','ribo', 'ribo', 'rna','rna','rna','rna')
reg <- factor(reg)
type <- relevel(type, ref = 'mock')
reg <- relevel(reg, ref = 'rna')
group  = reg:type

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

x <- cts[,c(3:ncol(cts))]
rownames(x) <- paste(cts$GENES_SYMBOL, cts$ENSEMBL_ID, sep="__")
y <- DGEList(counts=x,group=group)

keep <- filterByExpr(y, group=group)
y <- y[keep,,keep.lib.sizes=FALSE]

glMDSPlot(y, groups=group, folder = "glimma-plots/")

y <- calcNormFactors(y)

#ribo
plot(density(log2(y$counts[,1]*1000000/(y$samples$lib.size[1]+1))), from = -2, to = 20, lwd = 3)
for (k in 2:4){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$lib.size[k]+1))),col=k)
}

#rna
plot(density(log2(y$counts[,5]*1000000/(y$samples$lib.size[5]+1))), from = -2, to = 20, lwd = 3)
for (k in 6:48){
  lines(density(log2(y$counts[,k]*1000000/(y$samples$lib.size[k]+1))),col=k)
}

design <- model.matrix(~0 + group)
y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)

qlf_RO <- glmQLFTest(fit, contrast = c(1, -1, -1, 1))
qlf_rna <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0))
qlf_ribo <- glmQLFTest(fit, contrast = c(0, 0, -1, 1))

res_RO <-  topTags(qlf_RO, n = Inf)[[1]]
res_rna <-  topTags(qlf_rna, n =Inf)[[1]]
res_ribo <-  topTags(qlf_ribo, n =Inf)[[1]]

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

res_rna_SC2 <- res_rna
res_ribo_SC2 <- res_ribo
res_RO_SC2 <- res_RO