#' Run the whole assessment of a signature in one call
#'
#' The three questions a signature needs answering are usually asked in separate
#' analyses and often in the wrong order: what will it do, can it be improved, and how
#' much of the available signal does it capture. This runs them as one pipeline.
#'
#' The first two stages need only a reference atlas and no patient data at all, which is
#' the point: they can be run before committing a cohort. The third stage needs a cohort
#' and is skipped if none is supplied.
#'
#' 1. **Attribute** — which cells produce the signature's signal
#'    ([attribute_signature()]).
#' 2. **Triage and repair** — the predicted checkpoint-response effect and a
#'    recommendation ([triage_signature()]), and the set with off-target genes removed
#'    ([refine_signature()]).
#' 3. **Measure, if a cohort is given** — the original and refined sets against a matched
#'    random null ([benchmark_signatures()]), and the gap to what a full-transcriptome
#'    model reaches on the same data ([estimate_ceiling()]).
#'
#' @param signatures Named list of gene sets.
#' @param atlas_expr Genes x cells matrix from a single-cell reference, log scale.
#' @param atlas_cell_type One cell-type label per cell.
#' @param expr Optional genes x samples bulk matrix, to measure rather than predict.
#' @param outcome Optional outcome for `expr`: a binary vector or a [survival::Surv].
#' @param cell_map Passed to [triage_signature()] if your atlas labels differ from the
#'   model's vocabulary.
#' @param keep_types Passed to [refine_signature()].
#' @param ceiling If `TRUE` and a cohort is supplied, also fit the full-transcriptome
#'   model. This is the slow step; set `FALSE` to skip it.
#'
#' @return A data frame with one row per signature, and the intermediate objects in
#'   attributes `attribution`, `refined` and `ceiling`.
#' @export
signature_workflow <- function(signatures, atlas_expr, atlas_cell_type,
                               expr = NULL, outcome = NULL, cell_map = NULL,
                               keep_types = c("Tcell", "NK", "Bcell", "Plasma",
                                              "Macro", "DC"),
                               ceiling = TRUE) {
  stopifnot(is.list(signatures), length(signatures) > 0)

  att <- attribute_signature(atlas_expr, atlas_cell_type, signatures)
  tri <- triage_signature(att, cell_map = cell_map)
  ref <- refine_signature(atlas_expr, atlas_cell_type, signatures,
                          keep_types = keep_types)
  rep <- attr(ref, "report")

  out <- merge(tri, rep[, c("signature", "n_input", "n_kept")],
               by = "signature", all.x = TRUE)
  out$genes_dropped <- out$n_input - out$n_kept

  cl <- NULL
  if (!is.null(expr) && !is.null(outcome)) {
    b0 <- benchmark_signatures(expr, outcome, signatures)
    b1 <- benchmark_signatures(expr, outcome, ref[lengths(ref) >= 3])
    m <- merge(b0[, c("signature", "performance", "random_median", "beats_null")],
               b1[, c("signature", "performance")],
               by = "signature", suffixes = c("", "_refined"), all.x = TRUE)
    out <- merge(out, m, by = "signature", all.x = TRUE)
    out$refinement_gain <- out$performance_refined - out$performance
    if (isTRUE(ceiling)) {
      cl <- estimate_ceiling(t(.zrow(as.matrix(expr))), outcome)
      out$ceiling <- cl$ceiling
      # what fraction of the reachable signal, above the random floor, was captured
      out$fraction_of_reachable <-
        (out$performance - out$random_median) /
        pmax(cl$ceiling - out$random_median, 1e-9)
    }
  }
  out <- out[order(-out$predicted_response_effect), ]
  rownames(out) <- NULL
  attr(out, "attribution") <- att
  attr(out, "refined") <- ref
  attr(out, "ceiling") <- cl
  out
}
