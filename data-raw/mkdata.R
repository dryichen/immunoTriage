suppressMessages({library(matrixStats)})
V <- "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
P <- "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"

## --- built-in signature sets ---
gastro_signatures <- list(
  IFN_Ayers            = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG"),
  Tcell_inflamed_Ayers = c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
                           "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
                           "PSMB10","STAT1","TIGIT"),
  CYT_Rooney           = c("GZMA","PRF1"),
  Chemokine_Messina    = c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21",
                           "CXCL9","CXCL10","CXCL11","CXCL13"),
  APM                  = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","TAPBP","PSMB8","PSMB9","NLRC5"),
  CD8_Tcell            = c("CD8A","CD8B","CD3D","CD3E","GZMB","NKG7","PRF1","GNLY"),
  TLS_Cabrita          = c("CCL19","CCL21","CXCL13","CCL17","CCL22","CXCL9","CXCL10","CXCL11",
                           "CD79B","CD1D","CCR6","LAT","SKAP1","CETP","RBP5","PTGDS"),
  CAF                  = c("FAP","PDGFRA","PDGFRB","THY1","COL1A1","COL1A2","COL3A1",
                           "ACTA2","TAGLN","POSTN","LUM","DCN"),
  M2_TAM               = c("CD163","MRC1","MSR1","TREM2","APOE","SPP1","LGALS3","IL10","CCL18"),
  STAT1_axis           = c("STAT1","IRF1","GBP1","GBP4","PSMB9","TAP1","CXCL9","CXCL10"),
  Cytotoxic            = c("GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7","KLRD1"),
  Checkpoints          = c("PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","IDO1"),
  TGFB_Mariathasan     = c("ACTA2","ACTG2","ADAM12","ADAM19","CNN1","COL4A1","CTPS1","FSTL3",
                           "HSPB1","IGFBP3","JUNB","MYL9","MYLK","PALLD","PDLIM7","POSTN",
                           "SGK1","TAGLN","TGFB1I1","TGFBR2","TNS1","TPM1"))

# controls a genuinely immune-specific result must outperform
gastro_controls <- list(
  CTRL_proliferation = c("MKI67","TOP2A","BIRC5","CCNB1","BUB1","PLK1","AURKA","TYMS",
                         "RRM2","UBE2C","CDK1","CENPF","TPX2","ZWINT"),
  CTRL_hypoxia       = c("VEGFA","SLC2A1","PGAM1","ENO1","LDHA","CA9","ADM","P4HA1","ALDOA","MIF"),
  CTRL_ribosome      = c("RPL5","RPL11","RPS6","RPS3","RPL13A","RPS19","RPL23","RPS4X","RPL10","RPS8"))

imm <- c("IFN_Ayers","Tcell_inflamed_Ayers","CYT_Rooney","Chemokine_Messina","APM",
         "CD8_Tcell","TLS_Cabrita","STAT1_axis","Cytotoxic","Checkpoints")
str_ <- c("CAF","M2_TAM","TGFB_Mariathasan")
gastro_families <- c(setNames(rep("immune",length(imm)), imm),
                     setNames(rep("stroma",length(str_)), str_),
                     setNames(rep("confounder",length(gastro_controls)), names(gastro_controls)))


## --- gastric ICI cohort (Kim 2018, pembrolizumab) ---
e <- as.matrix(read.delim(file.path(V,"mGC_Pembrolizumab_Kim2018.TPM"), row.names=1, check.names=FALSE))
b <- read.delim(file.path(V,"benchmark_scores.tsv"), row.names=1, check.names=FALSE)
e <- log2(as.matrix(read.delim(file.path(V,"mGC_Pembrolizumab_Kim2018.TPM"),row.names=1,check.names=FALSE))+1)
e <- e[, rownames(b), drop=FALSE]
e <- e[rowSums(is.finite(e))==ncol(e) & rowSds(e)>0, , drop=FALSE]
sig_all <- unique(unlist(c(gastro_signatures, gastro_controls)))
keepg <- union(order(rowSds(e),decreasing=TRUE)[1:6000], which(rownames(e) %in% sig_all))
gc_ici <- list(expr = round(e[sort(keepg), , drop=FALSE], 3),
               response = as.integer(b$Responder),
               msi = b$MSI_H %in% c("True","TRUE",TRUE),
               ebv = b$EBV_pos %in% c("True","TRUE",TRUE))
cat("gc_ici:", dim(gc_ici$expr), " responders:", sum(gc_ici$response), "\n")

## --- gastric TCGA prognostic cohort ---
x <- as.matrix(read.delim(gzfile(file.path(P,"expr/STAD.gz")), row.names=1, check.names=FALSE))
colnames(x) <- gsub("\\.","-",colnames(x)); x <- x[, substr(colnames(x),14,15)=="01", drop=FALSE]
s <- read.delim(file.path(P,"Survival_SupplementalTable_S1_20171025_xena_sp"), check.names=FALSE)
m <- s[match(colnames(x), s$sample), ]
ok <- !is.na(m$OS) & !is.na(m$OS.time) & m$OS.time > 0
x <- x[, ok, drop=FALSE]; m <- m[ok, ]
x <- x[rowSums(is.finite(x))==ncol(x) & rowSds(x)>0, , drop=FALSE]
keepg <- union(order(rowSds(x),decreasing=TRUE)[1:6000], which(rownames(x) %in% sig_all))
gc_tcga <- list(expr = round(x[sort(keepg), , drop=FALSE], 3),
                time = as.numeric(m$OS.time), event = as.integer(m$OS),
                stage = as.character(m$ajcc_pathologic_tumor_stage),
                age = suppressWarnings(as.numeric(m$age_at_initial_pathologic_diagnosis)))
cat("gc_tcga:", dim(gc_tcga$expr), " events:", sum(gc_tcga$event), "\n")

save(gc_ici,  file="data/gc_ici.rda",  compress="xz")
save(gc_tcga, file="data/gc_tcga.rda", compress="xz")
save(gastro_signatures, file="data/gastro_signatures.rda", compress="xz")
save(gastro_controls,   file="data/gastro_controls.rda",   compress="xz")
save(gastro_families,   file="data/gastro_families.rda",   compress="xz")
cat("\nsaved\n")
