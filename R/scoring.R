#' Score gene signatures with the method used in their original publication
#'
#' Signatures published with GSVA, ssGSEA or weighted-mean scoring are scored that
#' way rather than being forced into a single convention, because scoring method is
#' a common reviewer objection and the choice is cheap to respect.
#'
#' @param expr Numeric matrix, genes x samples, on a log scale.
#' @param signatures Named list. Each element is either a character vector of genes,
#'   or a list with `genes` and optional `weights`.
#' @param methods Named character vector giving the method per signature, one of
#'   `"gsva"`, `"ssgsea"`, `"weighted"`, `"zmean"`. Missing entries default to
#'   `default_method`.
#' @param default_method Method for signatures not named in `methods`.
#' @param min_genes Minimum genes that must be present for a signature to be scored.
#'
#' @return Numeric matrix, samples x signatures. Columns are `NA`-free; signatures
#'   with too few matched genes are dropped with a warning.
#' @export
score_signatures <- function(expr, signatures,
                             methods = character(0),
                             default_method = c("zmean", "gsva", "ssgsea"),
                             min_genes = 3L) {
  default_method <- match.arg(default_method)
  stopifnot(is.matrix(expr), !is.null(rownames(expr)))
  Z <- .zrow(expr)
  out <- list(); gsva_sets <- list(); ssgsea_sets <- list(); dropped <- character(0)

  for (nm in names(signatures)) {
    s <- .as_signature(signatures[[nm]])
    g <- intersect(s$genes, rownames(expr))
    if (length(g) < min_genes) { dropped <- c(dropped, nm); next }
    m <- if (nm %in% names(methods)) methods[[nm]] else default_method
    m <- tolower(m)
    if (m == "weighted") {
      w <- s$weights[match(g, s$genes)]
      if (all(is.na(w))) w <- rep(1, length(g))
      w[is.na(w)] <- 0
      if (sum(abs(w)) == 0) { dropped <- c(dropped, nm); next }
      out[[nm]] <- as.numeric(crossprod(Z[g, , drop = FALSE], w)) / sum(abs(w))
    } else if (m == "gsva")   gsva_sets[[nm]]   <- g
    else if (m == "ssgsea")   ssgsea_sets[[nm]] <- g
    else                      out[[nm]] <- colMeans(Z[g, , drop = FALSE], na.rm = TRUE)
  }

  if (length(gsva_sets) || length(ssgsea_sets)) {
    if (!requireNamespace("GSVA", quietly = TRUE))
      stop("GSVA is required for gsva/ssgsea scoring; install it or use method 'zmean'.")
    if (length(gsva_sets)) {
      r <- GSVA::gsva(GSVA::gsvaParam(expr, gsva_sets, minSize = min_genes), verbose = FALSE)
      for (n in rownames(r)) out[[n]] <- as.numeric(r[n, ])
    }
    if (length(ssgsea_sets)) {
      r <- GSVA::gsva(GSVA::ssgseaParam(expr, ssgsea_sets, minSize = min_genes), verbose = FALSE)
      for (n in rownames(r)) out[[n]] <- as.numeric(r[n, ])
    }
  }
  if (length(dropped))
    warning(sprintf("%d signature(s) dropped for having <%d matched genes: %s",
                    length(dropped), min_genes, paste(utils::head(dropped, 5), collapse = ", ")))
  if (!length(out)) stop("no signature could be scored")
  M <- do.call(cbind, out)
  rownames(M) <- colnames(expr)
  M
}

.as_signature <- function(x) {
  if (is.character(x)) return(list(genes = x, weights = rep(NA_real_, length(x))))
  g <- x$genes
  w <- if (!is.null(x$weights)) x$weights else rep(NA_real_, length(g))
  list(genes = g, weights = w)
}

.zrow <- function(m) {
  s <- matrixStats::rowSds(m, na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- 1
  (m - rowMeans(m, na.rm = TRUE)) / s
}

.zs <- function(v) {
  s <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- 1
  (v - mean(v, na.rm = TRUE)) / s
}
