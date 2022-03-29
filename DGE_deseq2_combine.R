library(metap)
library(enrichR)


sumlog_c <-
  function(p, log.p = FALSE) {
    keep <- (na.exclude(p) > 0) & (na.exclude(p) <= 1)
    #invalid <- sum(1L * keep) == 1
    if(sum(1L * keep) == 1) {
      warning("Must have at least two valid p values")
      res <- list(chisq = NA_real_, df = NA_integer_,
                  p = na.exclude(p)[1], validp = p[keep])
    } else if (sum(1L * keep) == 0) {
        res <- list(chisq = NA_real_, df = NA_integer_, p = NA_integer_, validp = p[keep])
    } else {
      lnp <- log(na.exclude(p)[keep])
      chisq <- (-2) * sum(lnp)
      df <- 2 * length(lnp)
      if(length(lnp) != length(na.exclude(p))) {
        warning("Some studies omitted")
      }
      res <- list(chisq = chisq, df = df,
                  p = pchisq(chisq, df, lower.tail = FALSE,
                             log.p = log.p), validp = na.exclude(p)[keep])
    }
    class(res) <- c("sumlog", "metap")
    res
  }
print.sumlog_c <- function(x, ...) {
  cat("chisq = ", x$chisq, " with df = ", x$df, " p = ", x$p, "\n")
  invisible(x)
}


#Join attempts
#1
#inner
#DLEU1 DLEU7-AS1__ENSG00000176124.15
#216

#2
#inner
#DLEU1 ENSG00000176124.15
#216
dim(res_RO_reo)
dim(res_RO_SC2)
dim(res_RO_DENV1_72)

df <- cbind(res_RO_reo, res_RO_SC2, res_RO_DENV1_72, res_RO_IAV)

dfm <- as.data.frame(as_tibble(df, .name_repair = "universal"))
#rownames(dfm) <- rownames(df)
dfm <- dfm %>% rowwise() %>% mutate(Fisher_PValue = c(sumlog_c(c_across(starts_with("pv")))[3][[1]]))
#rownames(dfm) <- rownames(df)
dfm_c <- dfm[!is.na(dfm$Fisher_PValue),]
dfm_c$PValue_adjusted <- p.adjust(dfm_c$Fisher_PValue, "BH")
rownames(dfm_c) <- rownames(df)[!is.na(dfm$Fisher_PValue)]
dfm_c <- tibble::rownames_to_column(dfm_c, "Gene")

genid <- AnnotationDbi::mapIds(x = org.Hs.eg.db,
                               keys = res_vl_fpkmRvsNR_up$row_names,
                               column = "ENTREZID",
                               keytype = "SYMBOL",
                               multiVals = "first")

setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()s
df_e <- dfm_c[dfm_c$PValue_adjusted<0.2,]
df_e <- within(df_e, gene_names<-data.frame(do.call('rbind', strsplit(as.character(Gene), '__', fixed=TRUE))))
if (is.null(dbs)) websiteLive <- FALSE

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "Virus_Perturbations_from_GEO_down", "Virus_Perturbations_from_GEO_up", "Pfam_InterPro_Domains", "Pfam_Domains_2019")
if (websiteLive) {
  enriched <- enrichr(c(df_e$gene_names$X1), dbs)
}

#GO Molecular Function 2021

plotEnrich(enriched[["GO_Molecular_Function_2021"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

#GO Cellular Component 2021

plotEnrich(enriched[["GO_Cellular_Component_2021"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

#GO Biological Process 2021

plotEnrich(enriched[["GO_Biological_Process_2021"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

#GO Biological Process 2021

plotEnrich(enriched[["Virus_Perturbations_from_GEO_down"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

#GO Biological Process 2021

plotEnrich(enriched[["Virus_Perturbations_from_GEO_up"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

#GO Biological Process 2021

plotEnrich(enriched[["Pfam_InterPro_Domains"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

#GO Biological Process 2021

plotEnrich(enriched[["Pfam_Domains_2019"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

write.csv2(as.data.frame(dfm_d), str_c("./",'down_fishers_method.csv'), quote = F, row.names = F)
