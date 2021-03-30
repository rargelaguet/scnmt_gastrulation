############################################################
## Regress out technical covariates in the RNA expression ##
############################################################

# Number of expressed genes
foo <- rna_dt[,.(covariate=sum(expr>0)), by=c("id_rna")]
rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
  .[,expr:=lm(formula=expr~covariate)[["residuals"]]+mean(expr), by=c("gene","stage")] %>%
  .[,covariate:=NULL]

# Regress out batch effects in the RNA expression at E7.5
foo <- sample_metadata[stage=="E7.5",c("id_rna","plate")] %>%
  .[,plate:=as.factor(grepl("PS_VE",plate))]

rna_dt <- rbind(
  rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
    .[,expr:=lm(formula=expr~plate)[["residuals"]]+mean(expr), by="gene"] %>% .[,plate:=NULL],
  rna_dt[!id_rna%in%foo$id_rna]
)

# Regress out batch effects in the RNA expression at E6.5
foo <- sample_metadata[stage=="E6.5",c("id_rna","plate")] %>%
  .[,plate:=as.factor(grepl("E6.5_late",plate))]

rna_dt <- rbind(
  rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
    .[,expr:=lm(formula=expr~plate)[["residuals"]]+mean(expr), by="gene"] %>% .[,plate:=NULL],
  rna_dt[!id_rna%in%foo$id_rna]
)

# Mithocondrial content
foo <- rna_dt[grepl("mt-",gene)] %>% .[,.(mt=sum(expr)), by="id_rna"]
rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
  .[,expr:=lm(formula=expr~mt)[["residuals"]]+mean(expr), by=c("gene","stage")] %>%
  .[,mt:=NULL]
