library(metap)

#DF split and gene names extraction
#reo
res_ribo_reo_up <- res_ribo_reo %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_ribo_reo %>% filter(logFC > 0))$row_names))
res_ribo_reo_down <- res_ribo_reo %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_ribo_reo %>% filter(logFC < 0))$row_names))
res_rna_reo_up <- res_rna_reo %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_rna_reo %>% filter(logFC > 0))$row_names))
res_rna_reo_down <- res_rna_reo %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_rna_reo %>% filter(logFC < 0))$row_names))
res_RO_reo_up <- res_RO_reo %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_RO_reo %>% filter(logFC > 0))$row_names))
res_RO_reo_down <- res_RO_reo %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_RO_reo %>% filter(logFC < 0))$row_names))

res_ribo_reo_up <- res_ribo_reo %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_ribo_reo %>% filter(logFC > 0))$row_names))
res_ribo_reo_down <- res_ribo_reo %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_ribo_reo %>% filter(logFC < 0))$row_names))
res_rna_reo_up <- res_rna_reo %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_rna_reo %>% filter(logFC > 0))$row_names))
res_rna_reo_down <- res_rna_reo %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_rna_reo %>% filter(logFC < 0))$row_names))
res_RO_reo_up <- res_RO_reo %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_RO_reo %>% filter(logFC > 0))$row_names))
res_RO_reo_down <- res_RO_reo %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_RO_reo %>% filter(logFC < 0))$row_names))


res_ribo_SC2_up <- res_ribo_SC2 %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_ribo_SC2 %>% filter(logFC > 0))$row_names))
res_ribo_SC2_down <- res_ribo_SC2 %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_ribo_SC2 %>% filter(logFC < 0))$row_names))
res_rna_SC2_up <- res_rna_SC2 %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_rna_SC2 %>% filter(logFC > 0))$row_names))
res_rna_SC2_down <- res_rna_SC2 %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_rna_SC2 %>% filter(logFC < 0))$row_names))
res_RO_SC2_up <- res_RO_SC2 %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_RO_SC2 %>% filter(logFC > 0))$row_names))
res_RO_SC2_down <- res_RO_SC2 %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_RO_SC2 %>% filter(logFC < 0))$row_names))

res_ribo_SC2_up <- res_ribo_SC2 %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_ribo_SC2 %>% filter(logFC > 0))$row_names))
res_ribo_SC2_down <- res_ribo_SC2 %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_ribo_SC2 %>% filter(logFC < 0))$row_names))
res_rna_SC2_up <- res_rna_SC2 %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_rna_SC2 %>% filter(logFC > 0))$row_names))
res_rna_SC2_down <- res_rna_SC2 %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_rna_SC2 %>% filter(logFC < 0))$row_names))
res_RO_SC2_up <- res_RO_SC2 %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_RO_SC2 %>% filter(logFC > 0))$row_names))
res_RO_SC2_down <- res_RO_SC2 %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_RO_SC2 %>% filter(logFC < 0))$row_names))


res_ribo_48_up <- res_ribo_48 %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_ribo_48 %>% filter(logFC > 0))$row_names))
res_ribo_48_down <- res_ribo_48 %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_ribo_48 %>% filter(logFC < 0))$row_names))
res_rna_48_up <- res_rna_48 %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_rna_48 %>% filter(logFC > 0))$row_names))
res_rna_48_down <- res_rna_48 %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_rna_48 %>% filter(logFC < 0))$row_names))
res_RO_48_up <- res_RO_48 %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_RO_48 %>% filter(logFC > 0))$row_names))
res_RO_48_down <- res_RO_48 %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_RO_48 %>% filter(logFC < 0))$row_names))

res_ribo_48_up <- res_ribo_48 %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_ribo_48 %>% filter(logFC > 0))$row_names))
res_ribo_48_down <- res_ribo_48 %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_ribo_48 %>% filter(logFC < 0))$row_names))
res_rna_48_up <- res_rna_48 %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_rna_48 %>% filter(logFC > 0))$row_names))
res_rna_48_down <- res_rna_48 %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_rna_48 %>% filter(logFC < 0))$row_names))
res_RO_48_up <- res_RO_48 %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_RO_48 %>% filter(logFC > 0))$row_names))
res_RO_48_down <- res_RO_48 %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_RO_48 %>% filter(logFC < 0))$row_names))


res_ribo_72_up <- res_ribo_72 %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_ribo_72 %>% filter(logFC > 0))$row_names))
res_ribo_72_down <- res_ribo_72 %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_ribo_72 %>% filter(logFC < 0))$row_names))
res_rna_72_up <- res_rna_72 %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_rna_72 %>% filter(logFC > 0))$row_names))
res_rna_72_down <- res_rna_72 %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_rna_72 %>% filter(logFC < 0))$row_names))
res_RO_72_up <- res_RO_72 %>% filter(logFC > 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_RO_72 %>% filter(logFC > 0))$row_names))
res_RO_72_down <- res_RO_72 %>% filter(logFC < 0) %>% mutate(Gene = sub(pattern="\\_.*", "", (res_RO_72 %>% filter(logFC < 0))$row_names))

res_ribo_72_up <- res_ribo_72 %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_ribo_72 %>% filter(logFC > 0))$row_names))
res_ribo_72_down <- res_ribo_72 %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_ribo_72 %>% filter(logFC < 0))$row_names))
res_rna_72_up <- res_rna_72 %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_rna_72 %>% filter(logFC > 0))$row_names))
res_rna_72_down <- res_rna_72 %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_rna_72 %>% filter(logFC < 0))$row_names))
res_RO_72_up <- res_RO_72 %>% filter(logFC > 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_RO_72 %>% filter(logFC > 0))$row_names))
res_RO_72_down <- res_RO_72 %>% filter(logFC < 0) %>% mutate(ENSG = sub(pattern=".*\\_", "", (res_RO_72 %>% filter(logFC < 0))$row_names))

#Join attempts
#1
#inner
#DLEU1 DLEU7-AS1__ENSG00000176124.15
#216

#2
#inner
#DLEU1 ENSG00000176124.15
#216

df_u <- purrr::reduce(
  list(res_RO_reo_up, res_RO_SC2_up, res_RO_48_up),
  inner_join,
  by = "ENSG"
)

df_u <- df_u %>% rowwise() %>% mutate(Fisher_PValue = c(sumlog(c_across(starts_with("PV")))[3][[1]]))

dfm_u <- df_u[,c(1:6,17:20,32:35,51)]
dfm_u$PValue_adjusted <- p.adjust(dfm_u$Fisher_PValue, "BH")

dfm_u <- within(dfm_u, row_names.x<-data.frame(do.call('rbind', strsplit(as.character(row_names.x), '__', fixed=TRUE))))
write.csv2(dfm_u, str_c("./",'up_fishers_method.csv'), quote = F, row.names = F)






df_d <- purrr::reduce(
  list(res_RO_reo_down, res_RO_SC2_down, res_RO_48_down),
  inner_join,
  by = "ENSG"
)

df_d <- df_d %>% rowwise() %>% mutate(Fisher_PValue = c(sumlog(c_across(starts_with("PV")))[3][[1]]))

dfm_d <- df_d[,c(1:6,17:20,32:35,51)]
dfm_d$PValue_adjusted <- p.adjust(dfm_d$Fisher_PValue, "BH")

dfm_d <- within(dfm_d, row_names.x<-data.frame(do.call('rbind', strsplit(as.character(row_names.x), '__', fixed=TRUE))))
write.csv2(as.data.frame(dfm_d), str_c("./",'down_fishers_method.csv'), quote = F, row.names = F)
