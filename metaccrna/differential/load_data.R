
# Option 1: (lineage_A) vs (lineage_B,lineage_CC)
# Option 2: (lineage_A vs lineage_B) AND (lineageA vs lineage_C)

##############################################
## Load differential RNA expression results ##
##############################################

if (!is.null(io$diff.rna)) {
  if (opts$diff.type==1) {
    
    # Mesoderm-specific
    diff.rna.mes <- fread(sprintf("%s/E7.5Mesoderm_vs_E7.5EndodermEctoderm.txt.gz",io$diff.rna)) %>%
      .[,lineage:="Mesoderm"] %>% .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff]
    
    # Ectoderm-specific
    diff.rna.ect <- fread(sprintf("%s/E7.5Ectoderm_vs_E7.5MesodermEndoderm.txt.gz",io$diff.rna)) %>%
      .[,lineage:="Ectoderm"] %>% .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff]
    
    # Endoderm-specific
    diff.rna.end <- fread(sprintf("%s/E7.5Endoderm_vs_E7.5MesodermEctoderm.txt.gz",io$diff.rna)) %>%
      .[,lineage:="Endoderm"] %>% .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff]
    
  } else if (opts$diff.type==2) {
    
    # Mesoderm-specific
    diff.rna.mes <- rbind(
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Endoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Ectoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Ectoderm"]
      ) %>% .[,lineage1:="Mesoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff] %>% 
      data.table::dcast(id+symbol+lineage1~lineage2, value.var=c("logFC","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Ectoderm==T & sign(logFC_Ectoderm)==sign(logFC_Endoderm)] %>%
      .[,logFC:=(logFC_Endoderm+logFC_Ectoderm)/2] %>%
      .[,c("id","symbol","lineage1","logFC","sig")] %>% setnames("lineage1","lineage")
    
    
    # Ectoderm-specific
    diff.rna.ect <- rbind(
      fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Endoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Endoderm"],
      fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Mesoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Mesoderm"]
    ) %>% .[,lineage1:="Ectoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff] %>% 
      data.table::dcast(id+symbol+lineage1~lineage2, value.var=c("logFC","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Endoderm==T & sign(logFC_Endoderm)==sign(logFC_Mesoderm)] %>%
      .[,logFC:=(logFC_Endoderm+logFC_Mesoderm)/2] %>%
      .[,c("id","symbol","lineage1","logFC","sig")] %>% setnames("lineage1","lineage")
    
    # Endoderm-specific
    diff.rna.end <- rbind(
      fread(sprintf("%s/E7.5Endoderm_vs_E7.5Ectoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Ectoderm"],
      fread(sprintf("%s/E7.5Endoderm_vs_E7.5Mesoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Mesoderm"]
    ) %>% .[,lineage1:="Endoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff] %>% 
      data.table::dcast(id+symbol+lineage1~lineage2, value.var=c("logFC","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Ectoderm==T & sign(logFC_Ectoderm)==sign(logFC_Mesoderm)] %>%
      .[,logFC:=(logFC_Ectoderm+logFC_Mesoderm)/2] %>%
      .[,c("id","symbol","lineage1","logFC","sig")] %>% setnames("lineage1","lineage")
  }
  
  diff.rna <- do.call("rbind", list(diff.rna.mes,diff.rna.end,diff.rna.ect))
  rm(diff.rna.mes,diff.rna.ect,diff.rna.end)
}

###############################################
## Load differential DNA methylation results ##
###############################################

if (!is.null(io$diff.met)) {
  
  # Mesoderm-specific
  if (opts$diff.type==1) {
    diff.met.mes <- lapply(names(opts$met.annos), function(x)
    # diff.met.mes <- lapply(list("H3K27ac_distal_E7.5_Mes_intersect12"), function(x)
      fread(sprintf("%s/E7.5Mesoderm_vs_E7.5EndodermEctoderm_%s.txt.gz",io$diff.met,x))
    ) %>% rbindlist() %>% .[,lineage:="Mesoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff]
  
    # Ectoderm-specific
    diff.met.ect <- lapply(names(opts$met.annos), function(x)
    # diff.met.ect <- lapply(list("H3K27ac_distal_E7.5_Ect_intersect12"), function(x)
      fread(sprintf("%s/E7.5Ectoderm_vs_E7.5MesodermEndoderm_%s.txt.gz",io$diff.met,x))
    ) %>% rbindlist() %>% .[,lineage:="Ectoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff]
  
    # Endoderm-specific
    diff.met.end <- lapply(names(opts$met.annos), function(x)
    # diff.met.end <- lapply(list("H3K27ac_distal_E7.5_End_intersect12"), function(x)
      fread(sprintf("%s/E7.5Endoderm_vs_E7.5MesodermEctoderm_%s.txt.gz",io$diff.met,x))
    ) %>% rbindlist() %>% .[,lineage:="Endoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff]
  }
  
  if (opts$diff.type==2) {
    
    # Mesoderm-specific
    diff.met.mes <- lapply(list("H3K27ac_distal_E7.5_Mes_intersect12"), function(x)
    # diff.met.mes <- lapply(names(opts$met.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist %>% .[,lineage1:="Mesoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Ectoderm==T & sign(diff_Ectoderm)==sign(diff_Endoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Ectoderm-specific
    diff.met.ect <- lapply(list("H3K27ac_distal_E7.5_Ect_intersect12"), function(x)
    # diff.met.ect <- lapply(names(opts$met.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Mesoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Ectoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Mesoderm==T & sign(diff_Endoderm)==sign(diff_Mesoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Mesoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Endoderm-specific
    diff.met.end <- lapply(list("H3K27ac_distal_E7.5_End_intersect12"), function(x)
    # diff.met.end <- lapply(names(opts$met.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Mesoderm"],
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Endoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Ectoderm==T & sign(diff_Ectoderm)==sign(diff_Mesoderm)] %>%
      .[,diff:=(diff_Mesoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
    }
  
  diff.met <- do.call("rbind", list(diff.met.mes,diff.met.end,diff.met.ect))
  rm(diff.met.mes,diff.met.ect,diff.met.end)
}

######################################################
## Load differential chromatin accessiblity results ##
######################################################



if (!is.null(io$diff.acc)) {
  
  if (opts$diff.type==1) {
    
    # Mesoderm-specific
    diff.acc.mes <- lapply(names(opts$acc.annos), function(x)
    # diff.acc.mes <- lapply(list("H3K27ac_distal_E7.5_Mes_intersect12"), function(x)
      fread(sprintf("%s/E7.5Mesoderm_vs_E7.5EndodermEctoderm_%s.txt.gz",io$diff.acc,x))
    ) %>% rbindlist() %>% .[,lineage:="Mesoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff]
  
    # Ectoderm-specific
    diff.acc.ect <- lapply(names(opts$acc.annos), function(x)
    # diff.acc.ect <- lapply(list("H3K27ac_distal_E7.5_Ect_intersect12"), function(x)
      fread(sprintf("%s/E7.5Ectoderm_vs_E7.5MesodermEndoderm_%s.txt.gz",io$diff.acc,x))
    ) %>% rbindlist() %>% .[,lineage:="Ectoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff]
  
    # Endoderm-specific
    diff.acc.end <- lapply(names(opts$acc.annos), function(x)
    # diff.acc.end <- lapply(list("H3K27ac_distal_E7.5_End_intersect12"), function(x)
      fread(sprintf("%s/E7.5Endoderm_vs_E7.5MesodermEctoderm_%s.txt.gz",io$diff.acc,x))
    ) %>% rbindlist() %>% .[,lineage:="Endoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff]
  
  }
  
  if (opts$diff.type==2) {
    # Mesoderm-specific
    # diff.acc.mes <- lapply(names(opts$acc.annos), function(x)
    diff.acc.mes <- lapply(list("H3K27ac_distal_E7.5_Mes_intersect12"), function(x)
      rbind(
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Mesoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Ectoderm==T & sign(diff_Endoderm)==sign(diff_Ectoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Ectoderm-specific
    diff.acc.ect <- lapply(list("H3K27ac_distal_E7.5_Ect_intersect12"), function(x)
    # diff.acc.ect <- lapply(names(opts$acc.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Mesoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Ectoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Mesoderm==T & sign(diff_Endoderm)==sign(diff_Mesoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Mesoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Endoderm-specific
    diff.acc.end <- lapply(list("H3K27ac_distal_E7.5_End_intersect12"), function(x)
    # diff.acc.end <- lapply(names(opts$acc.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Mesoderm"],
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Endoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Ectoderm==T & sign(diff_Mesoderm)==sign(diff_Ectoderm)] %>%
      .[,diff:=(diff_Mesoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  }
  
  diff.acc <- do.call("rbind", list(diff.acc.mes,diff.acc.end,diff.acc.ect))
  rm(diff.acc.mes,diff.acc.ect,diff.acc.end)
}
