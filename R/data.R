#' Gastric cancer checkpoint-inhibitor cohort (Kim et al. 2018)
#'
#' Pembrolizumab-treated metastatic gastric cancer, PRJEB25780. Expression is
#' log2(TPM+1) restricted to the 6,000 most variable genes to keep the package small.
#'
#' @format A list with `expr` (genes x samples), `response` (1 = CR/PR),
#'   `msi` and `ebv` logical vectors.
#' @source European Nucleotide Archive PRJEB25780
"gc_ici"

#' Gastric cancer prognostic cohort (TCGA-STAD)
#'
#' Primary tumours with overall survival. Expression is log2 RSEM restricted to the
#' 6,000 most variable genes.
#'
#' @format A list with `expr`, `time`, `event`, `stage`, `age`.
#' @source TCGA via UCSC Xena; endpoints from TCGA-CDR
"gc_tcga"

#' Published immune and stromal signatures
#' @format Named list of character vectors.
"gastro_signatures"

#' Confounder control signatures
#'
#' Proliferation, hypoxia and ribosomal programmes. An immune-specific result must
#' outperform these; if it does not, the result reflects generic transcriptional
#' activity rather than immune biology.
#' @format Named list of character vectors.
"gastro_controls"

#' Family labels for the built-in signatures
#' @format Named character vector.
"gastro_families"
