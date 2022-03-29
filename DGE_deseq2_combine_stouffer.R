library(metap)

setwd("D://Projects/kenezam/gene_counts/")

#MHV: Unable to remove hidden batches
df11 <- read.csv("MHV/RO_MHV_ds2.tsv", sep="\t")
df12 <- read.csv("MHV/RO_MHV_ds2_l.tsv", sep="\t")
df13 <- read.csv("MHV/RO_MHV_ds2_la.tsv", sep="\t")
#SARS-CoV-2_big (GSE158930): Hundreds of thousands of ribo counts
df211 <- read.csv("GSE158930/res_RO_ds2_4h.tsv", sep="\t")
df212 <- read.csv("GSE158930/res_RO_ds2_4h_l.tsv", sep="\t")
df213 <- read.csv("GSE158930/res_RO_ds2_4h_la.tsv", sep="\t")
df221 <- read.csv("GSE158930/res_RO_ds2_96h.tsv", sep="\t")
df222 <- read.csv("GSE158930/res_RO_ds2_96h_l.tsv", sep="\t")
df223 <- read.csv("GSE158930/res_RO_ds2_96h_la.tsv", sep="\t")
#SARS-CoV-2 (GSE162323): Unable to remove hidden batches, low correlation with paper, counts from paper were used
df31 <- read.csv("MHV/RO_MHV_ds2.tsv", sep="\t")
df32 <- read.csv("MHV/RO_MHV_ds2_l.tsv", sep="\t")
df33 <- read.csv("MHV/RO_MHV_ds2_la.tsv", sep="\t")
#SARS-CoV-2_med: Outlier in ribo and strong batch in rna, nice authors
df411 <- read.csv("SARSCOV2/res_RO_ds2_6h.tsv", sep="\t")
df412 <- read.csv("SARSCOV2/res_RO_ds2_6h_l.tsv", sep="\t")
df413 <- read.csv("SARSCOV2/res_RO_ds2_6h_la.tsv", sep="\t")
df421 <- read.csv("SARSCOV2/res_RO_ds2_24h.tsv", sep="\t")
df422 <- read.csv("SARSCOV2/res_RO_ds2_24h_l.tsv", sep="\t")
df423 <- read.csv("SARSCOV2/res_RO_ds2_24h_la.tsv", sep="\t")
#HBV (GSE135860): Length distribution, many SVs but no distinct clusters on SVA-plot, outlier removed
df51 <- read.csv("GSE135860/RO_HBV_ds2_wo.tsv", sep="\t")
df52 <- read.csv("GSE135860/RO_HBV_ds2_wo_l.tsv", sep="\t")
df53 <- read.csv("GSE135860/RO_HBV_ds2_wo_la.tsv", sep="\t")
#IAV (GSE101760): Rep batch
df61 <- read.csv("GSE101760/res_RO_ds2_b.tsv", sep="\t")
df62 <- read.csv("GSE101760/res_RO_ds2_b_l.tsv", sep="\t")
df63 <- read.csv("GSE101760/res_RO_ds2_b_la.tsv", sep="\t")
#EV71 (GSE103308): 2 time points
df711 <- read.csv("GSE103308/res_RO_ds2_35.tsv", sep="\t")
df712 <- read.csv("GSE103308/res_RO_ds2_35_l.tsv", sep="\t")
df713 <- read.csv("GSE103308/res_RO_ds2_35_la.tsv", sep="\t")
df721 <- read.csv("GSE103308/res_RO_ds2_625.tsv", sep="\t")
df722 <- read.csv("GSE103308/res_RO_ds2_625_l.tsv", sep="\t")
df723 <- read.csv("GSE103308/res_RO_ds2_625_la.tsv", sep="\t")
#DENV1 (GSE69602): Hundreds of thousands of ribo counts (close to 1M) + low on rna (3-5M) + problems with author
df811 <- read.csv("GSE69602/res_RO_ds2_48h.tsv", sep="\t")
df812 <- read.csv("GSE69602/res_RO_ds2_48h_l.tsv", sep="\t")
df813 <- read.csv("GSE69602/res_RO_ds2_48h_la.tsv", sep="\t")
df821 <- read.csv("GSE69602/res_RO_ds2_72h.tsv", sep="\t")
df822 <- read.csv("GSE69602/res_RO_ds2_72h_l.tsv", sep="\t")
df823 <- read.csv("GSE69602/res_RO_ds2_72h_la.tsv", sep="\t")

metap::sumz(??(0.01, 0.1, 0.5), c(1,1,1))
metap::sumz(c(0.01, 0.1, 0.5), c(1,0.5,2))

dim(df12)
dim(df212)
dim(df222)
dim(df32)
dim(df412)
dim(df422)
dim(df52)

#1st approach (NAs)
cts <- purrr::reduce(
  list(df12, df212, df222, df32, df412, df422, df52, df62, df712, df722, df812, df822),
  full_join,
  by = "N_unmapped")

#2nd approach (p-value = 1)
df <- cbind(df12, df212, df222, df32, df412, df422, df52, df62, df712, df722, df812, df822)

dfm <- as.data.frame(as_tibble(df, .name_repair = "universal"))
#rownames(dfm) <- rownames(df)
dfm <- dfm %>% rowwise() %>% mutate(Fisher_PValue = c(sumlog_c(c_across(starts_with("pv")))[3][[1]]))
#rownames(dfm) <- rownames(df)
dfm_c <- dfm[!is.na(dfm$Fisher_PValue),]
