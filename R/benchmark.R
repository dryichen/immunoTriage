#' Benchmark signatures against a matched random-gene-set null
#'
#' The null draws random gene sets of the same size as the signature and allows the
#' direction of association to be chosen post hoc — the same latitude a published
#' signature receives in practice. This makes the null conservative rather than
#' lenient, which is the point: a signature that cannot clear it has not been shown
#' to carry information specific to the biology it names.
#'
#' @param expr Genes x samples matrix, log scale.
#' @param outcome Either a binary vector (response; AUC is used) or a
#'   [survival::Surv] object (C-index is used).
#' @param signatures Named list of gene sets (see [score_signatures()]).
#' @param methods,default_method Passed to [score_signatures()].
#' @param n_random Number of random gene sets.
#' @param set_size Size of each random set; defaults to the median signature size.
#' @param seed RNG seed.
#'
#' @return A data frame with one row per signature: the performance metric, the
#'   random-null median and 95th percentile, and whether the signature exceeds it.
#' @export
benchmark_signatures <- function(expr, outcome, signatures,
                                 methods = character(0), default_method = "zmean",
                                 n_random = 300L, set_size = NULL, seed = 1L) {
  set.seed(seed)
  S <- score_signatures(expr, signatures, methods, default_method)
  metric <- .metric_for(outcome)
  if (is.null(set_size))
    set_size <- max(3L, stats::median(vapply(signatures, function(x) length(.as_signature(x)$genes), 1L)))

  Z <- .zrow(expr)
  obs <- vapply(colnames(S), function(s) metric(outcome, S[, s]), numeric(1))
  rnd <- vapply(seq_len(n_random), function(i)
    metric(outcome, colMeans(Z[sample(rownames(Z), set_size), , drop = FALSE])), numeric(1))
  rnd <- rnd[is.finite(rnd)]

  data.frame(
    signature   = colnames(S),
    n_genes     = vapply(colnames(S), function(s)
                    length(intersect(.as_signature(signatures[[s]])$genes, rownames(expr))), 1L),
    performance = unname(obs),
    random_median = stats::median(rnd),
    random_p95    = unname(stats::quantile(rnd, 0.95)),
    delta         = unname(obs) - stats::median(rnd),
    beats_null    = unname(obs) > stats::quantile(rnd, 0.95),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

#' Test whether a signature family behaves differently in two settings
#'
#' The prognostic-versus-predictive dissociation. Supplying a `confounder` family is
#' strongly recommended: proliferation, hypoxia and chromosomal-instability scores
#' should NOT shift between settings, and if they do, the shift observed for the
#' immune family is not specific and should not be interpreted.
#'
#' @param bench_a,bench_b Outputs of [benchmark_signatures()] from the two settings
#'   (for example prognosis and treatment response), aggregated across cohorts.
#' @param family Named character vector mapping signature to family.
#' @param target,confounder Family labels for the family of interest and the control.
#'
#' @return A list with the per-family shift table and the interaction test.
#' @export
test_dissociation <- function(bench_a, bench_b, family,
                              target = "immune", confounder = "confounder") {
  m <- merge(
    stats::aggregate(delta ~ signature, bench_a, stats::median),
    stats::aggregate(delta ~ signature, bench_b, stats::median),
    by = "signature", suffixes = c("_a", "_b"))
  m$shift  <- m$delta_b - m$delta_a
  m$family <- unname(family[m$signature])
  m <- m[!is.na(m$family), ]

  tab <- do.call(rbind, lapply(split(m, m$family), function(g) data.frame(
    family = g$family[1], n = nrow(g),
    delta_a = stats::median(g$delta_a), delta_b = stats::median(g$delta_b),
    shift = stats::median(g$shift),
    p_shift = if (nrow(g) > 5) stats::wilcox.test(g$shift)$p.value else NA_real_,
    row.names = NULL)))

  a <- m$shift[m$family == target]; b <- m$shift[m$family == confounder]
  inter <- if (length(a) >= 3 && length(b) >= 3)
    stats::wilcox.test(a, b, alternative = "greater")$p.value else NA_real_
  list(per_signature = m, per_family = tab, interaction_p = inter,
       note = "A significant target shift is interpretable only if the confounder family does not shift.")
}

#' Treatment-interaction test in a randomised trial
#'
#' The formal definition of a predictive biomarker. The same random gene sets are
#' pushed through the identical model so that the false-positive rate is measured on
#' these data rather than assumed to be 5%.
#'
#' @param expr Genes x samples matrix.
#' @param time,event Survival time and event indicator.
#' @param arm Binary treatment indicator (1 = experimental arm).
#' @param signatures Named list of gene sets.
#' @param n_random,set_size,seed As in [benchmark_signatures()].
#' @param methods,default_method Passed to [score_signatures()].
#'
#' @return A list with the per-signature interaction table, the empirical random-set
#'   false-positive rate, and a binomial test of the signature hit rate against it.
#' @export
test_interaction <- function(expr, time, event, arm, signatures,
                             methods = character(0), default_method = "zmean",
                             n_random = 300L, set_size = NULL, seed = 1L) {
  set.seed(seed)
  S <- score_signatures(expr, signatures, methods, default_method)
  if (is.null(set_size)) set_size <- max(3L, ncol(S) * 0 + 48L)
  y <- survival::Surv(time, event)

  fit1 <- function(v) {
    d <- data.frame(v = .zs(v), arm = arm)
    d$vx <- d$v * d$arm
    f <- try(survival::coxph(y ~ v + arm + vx, data = d), silent = TRUE)
    if (inherits(f, "try-error")) return(c(NA_real_, NA_real_))
    s <- summary(f)$coefficients
    if (!("vx" %in% rownames(s))) return(c(NA_real_, NA_real_))
    c(s["vx", "coef"], s["vx", ncol(s)])
  }
  res <- t(vapply(colnames(S), function(s) fit1(S[, s]), numeric(2)))
  tab <- data.frame(signature = rownames(res), interaction_coef = res[, 1],
                    p = res[, 2], row.names = NULL)
  tab$q <- stats::p.adjust(tab$p, "BH")

  Z <- .zrow(expr)
  rp <- vapply(seq_len(n_random), function(i)
    fit1(colMeans(Z[sample(rownames(Z), set_size), , drop = FALSE]))[2], numeric(1))
  rp <- rp[is.finite(rp)]
  rate <- mean(rp < 0.05)
  hits <- sum(tab$p < 0.05, na.rm = TRUE)
  bt <- stats::binom.test(hits, sum(is.finite(tab$p)),
                          p = max(rate, 1e-3), alternative = "greater")$p.value
  list(per_signature = tab[order(tab$p), ], random_fp_rate = rate,
       n_hits = hits, n_tested = sum(is.finite(tab$p)), binomial_p = bt)
}

.metric_for <- function(outcome) {
  if (inherits(outcome, "Surv")) {
    function(o, v) { np <- -as.numeric(v)
      c <- as.numeric(survival::concordance(o ~ np)$concordance); max(c, 1 - c) }
  } else {
    function(o, v) { ok <- is.finite(v) & !is.na(o); o <- o[ok]; v <- v[ok]
      if (length(unique(o)) < 2) return(NA_real_)
      r <- rank(v); n1 <- sum(o == 1); n0 <- sum(o == 0)
      a <- (sum(r[o == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0); max(a, 1 - a) }
  }
}
