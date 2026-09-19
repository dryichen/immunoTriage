#' Predict how a signature will behave, before touching a patient cohort
#'
#' Validating a signature normally means finding a cohort, running it, and seeing
#' what happens. That is slow, and — as the benchmark in this package shows — usually
#' run against the wrong endpoint. Which cells produce a signature's signal turns out
#' to determine how it behaves, and that is measurable in a single-cell atlas with no
#' outcome data at all.
#'
#' `attribute_signature()` computes where a signature's signal comes from.
#' `triage_signature()` turns that into a predicted effect for checkpoint response and
#' a recommendation.
#'
#' The supplied model was fitted on 50 published signatures whose effects were measured
#' across 20 cancer types and 12 checkpoint cohorts, and validated by nested
#' leave-one-signature-out: the held-out signature contributed nothing to either
#' fitting or model selection. Expect a rank correlation near 0.6 between predicted
#' and measured response effect — useful for ordering candidates, not for replacing
#' a trial.
#'
#' @param expr Genes x cells matrix from a single-cell reference, log scale.
#' @param cell_type Character vector, one label per cell (column of `expr`).
#' @param signatures Named list of gene sets.
#' @param min_genes Minimum genes that must be present.
#'
#' @examples
#' # a toy atlas: NK cells carry signature A, malignant epithelium carries signature B
#' set.seed(1)
#' types <- rep(c("Tcell", "NK", "Macro", "Epi", "Fibro"), each = 40)
#' expr  <- matrix(rnorm(300 * length(types)), 300,
#'                 dimnames = list(paste0("G", 1:300), NULL))
#' expr[1:20,  types == "NK"]  <- expr[1:20,  types == "NK"]  + 2
#' expr[21:40, types == "Epi"] <- expr[21:40, types == "Epi"] + 2
#' sigs <- list(A = paste0("G", 1:20), B = paste0("G", 21:40))
#'
#' att <- attribute_signature(expr, types, sigs)
#' triage_signature(att)   # A ranks above B: the NK-derived signature is the one to test
#'
#' # your own labels can be remapped onto the model's vocabulary
#' names(triage_model$coef)
#'
#' @return A data frame with one row per signature: the z-scored enrichment in each
#'   cell type, the cell type contributing most, and that z.
#' @export
attribute_signature <- function(expr, cell_type, signatures, min_genes = 3L) {
  stopifnot(is.matrix(expr), length(cell_type) == ncol(expr))
  Z <- .zrow(expr)
  types <- sort(unique(as.character(cell_type)))
  out <- lapply(names(signatures), function(nm) {
    g <- intersect(.as_signature(signatures[[nm]])$genes, rownames(expr))
    if (length(g) < min_genes) return(NULL)
    per_cell <- colMeans(Z[g, , drop = FALSE], na.rm = TRUE)
    mu <- vapply(types, function(t) mean(per_cell[cell_type == t], na.rm = TRUE), 1)
    z  <- (mu - mean(mu)) / (stats::sd(mu) + 1e-9)
    data.frame(signature = nm, top_cell_type = types[which.max(z)],
               top_z = max(z), t(z), check.names = FALSE)
  })
  out <- out[!vapply(out, is.null, TRUE)]
  if (!length(out)) stop("no signature had enough genes present")
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

#' @rdname attribute_signature
#'
#' @param attribution Output of [attribute_signature()].
#' @param cell_map Named character vector remapping your cell-type labels onto the
#'   model's vocabulary, e.g. `c("CD8T" = "Tcell", "Malignant" = "Epi")`. The model
#'   was fitted on the ten types in `names(triage_model$coef)`; any type you
#'   cannot map is set to zero, which is the neutral value for a z-score.
#'
#' @return A data frame with the predicted checkpoint-response effect and a
#'   recommendation, ordered best first.
#' @export
triage_signature <- function(attribution, cell_map = NULL) {
  stopifnot(is.data.frame(attribution), "signature" %in% names(attribution))
  m <- triage_model
  z <- attribution[, setdiff(names(attribution), c("signature", "top_cell_type", "top_z")),
                   drop = FALSE]
  if (!is.null(cell_map)) {
    hit <- names(z) %in% names(cell_map)
    names(z)[hit] <- unname(cell_map[names(z)[hit]])
  }
  want <- names(m$coef)
  missing <- setdiff(want, names(z))
  if (length(missing) == length(want))
    stop("none of the model's cell types were found; use `cell_map` to relabel. Model expects: ",
         paste(want, collapse = ", "))
  if (length(missing))
    warning("cell types not found, treated as neutral: ", paste(missing, collapse = ", "),
            call. = FALSE)

  # z-scores are already on a common scale by construction in attribute_signature(),
  # so the fitted coefficients apply directly; an absent type contributes nothing.
  Z <- matrix(0, nrow(z), length(want), dimnames = list(NULL, want))
  present <- intersect(want, names(z))
  Z[, present] <- as.matrix(z[, present, drop = FALSE])
  pred <- as.numeric(m$intercept + Z %*% m$coef)

  call <- ifelse(pred >= m$cut_use,   "test for treatment response",
          ifelse(pred >= m$cut_maybe, "weak; test only with a strong prior",
                                      "not worth a response cohort"))
  out <- data.frame(signature = attribution$signature,
                    top_cell_type = attribution$top_cell_type,
                    predicted_response_effect = pred,
                    recommendation = call,
                    stringsAsFactors = FALSE)
  out[order(-pred), ]
}

#' Fitted triage model
#'
#' Ridge coefficients relating a signature's cell-type attribution to the checkpoint
#' response effect it will show, fitted on the 50 signatures in [sig_behavior].
#' Held-out performance (leave-one-signature-out) is Spearman rho = 0.68 against the
#' measured effect, AUC 0.92 for identifying a signature that works. At `cut_use` the
#' precision is 0.85 against a base rate of 0.42, at 0.81 recall.
#'
#' @format A list: `intercept`, a named `coef` vector over ten cell types, the two
#'   decision thresholds, and the validation statistics.
"triage_model"

#' Signature-behaviour reference table
#'
#' Measured prognostic and checkpoint-response effects for 50 published signatures,
#' with the cell type producing each signature's signal and the leave-one-signature-out
#' prediction. The training data behind [triage_signature()], provided so the
#' calibration can be inspected or refitted.
#'
#' @format A data frame with 50 rows.
"sig_behavior"

#' Repair a signature by dropping the genes that report on the wrong cells
#'
#' A signature is a set of genes assumed to report on one population. Attributing each
#' gene separately in a reference atlas shows that some of them do not: an immune
#' signature can contain genes whose expression is dominated by epithelium or stroma,
#' and those genes dilute the quantity the set is meant to measure.
#'
#' Removing them improves the signature. Across 318 signature-by-cohort comparisons in
#' 12 checkpoint cohorts, attribution-guided removal raised the response AUC by +0.008 on
#' average (paired P = 1 x 10^-4), and 23 of 27 signatures improved. Because shorter gene
#' sets score better in general, the comparison that matters is against dropping the same
#' number of genes at random from the same signature: random removal made signatures
#' slightly worse (-0.002), and guided removal beat it at P = 4 x 10^-7.
#'
#' @param expr Genes x cells matrix from a single-cell reference, log scale.
#' @param cell_type Character vector, one label per cell.
#' @param signatures Named list of gene sets.
#' @param keep_types Cell types the signature is meant to report on. Genes attributed
#'   elsewhere are dropped. Defaults to the immune compartment.
#' @param min_keep Refuse to return a set shorter than this.
#'
#' @return A list with the same names as `signatures`, each the retained genes, plus a
#'   `report` attribute giving what was dropped and why.
#' @export
refine_signature <- function(expr, cell_type, signatures,
                             keep_types = c("Tcell", "NK", "Bcell", "Plasma",
                                            "Macro", "DC"),
                             min_keep = 3L) {
  stopifnot(is.matrix(expr), length(cell_type) == ncol(expr))
  Z <- .zrow(expr)
  types <- sort(unique(as.character(cell_type)))
  if (!any(keep_types %in% types))
    stop("none of `keep_types` are present in `cell_type`")
  # each gene's own cell-type profile, standardised across types
  prof <- vapply(types, function(t) rowMeans(Z[, cell_type == t, drop = FALSE], na.rm = TRUE),
                 numeric(nrow(Z)))
  gz <- (prof - rowMeans(prof)) / (apply(prof, 1, stats::sd) + 1e-9)
  top <- types[max.col(gz, ties.method = "first")]
  names(top) <- rownames(expr)

  rep_rows <- list()
  out <- lapply(names(signatures), function(nm) {
    g <- intersect(.as_signature(signatures[[nm]])$genes, rownames(expr))
    if (!length(g)) return(character(0))
    off <- !(top[g] %in% keep_types)
    keep <- g[!off]
    if (length(keep) < min_keep) keep <- g          # too little left; leave it alone
    rep_rows[[nm]] <<- data.frame(signature = nm, n_input = length(g),
                                  n_kept = length(keep),
                                  dropped = paste(g[off], collapse = ","),
                                  stringsAsFactors = FALSE)
    keep
  })
  names(out) <- names(signatures)
  attr(out, "report") <- do.call(rbind, rep_rows)
  out
}
