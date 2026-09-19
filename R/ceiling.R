#' Estimate the information ceiling of a dataset
#'
#' A signature's performance is only interpretable against what *any* model could
#' achieve on the same data. This fits several machine-learning models to the full
#' feature matrix under honest validation and returns the best out-of-fold
#' performance, which upper-bounds what a hand-curated score could reach.
#'
#' The models are deliberately heterogeneous — a sparse linear model, a random
#' forest, and gradient boosting — because a ceiling claimed from one model class is
#' only a statement about that class. Validation is **grouped** when `group` is
#' supplied: whole cohorts are held out, which is the only honest way to estimate
#' transfer across datasets and typically much lower than within-cohort estimates.
#'
#' Interpretation matters more than the number. If the ceiling sits close to the best
#' signature, the task is information-limited and no amount of modelling will help;
#' reporting that is more useful than another marginally different score.
#'
#' @param x Samples x features matrix (genes or signature scores). Standardise per
#'   cohort before calling if cohorts are pooled.
#' @param y Binary outcome vector, or a [survival::Surv] object.
#' @param group Optional cohort labels for leave-one-group-out validation. If `NULL`,
#'   k-fold cross-validation is used and the result is labelled as within-cohort.
#' @param models Which model families to fit. Any of `"enet"`, `"rf"`, `"xgb"`.
#' @param nfolds Folds when `group` is `NULL`.
#' @param seed RNG seed.
#'
#' @return A list with per-model and per-fold performance, the ceiling estimate, and
#'   the validation scheme actually used.
#' @export
estimate_ceiling <- function(x, y, group = NULL,
                             models = c("enet", "rf", "xgb"),
                             nfolds = 5L, seed = 1L) {
  set.seed(seed)
  models <- match.arg(models, several.ok = TRUE)
  surv <- inherits(y, "Surv")
  x <- as.matrix(x)
  keep <- apply(is.finite(x), 2, all) & (matrixStats::colSds(x) > 0)
  x <- x[, keep, drop = FALSE]

  folds <- if (!is.null(group)) split(seq_along(group), group)
           else split(seq_len(nrow(x)), sample(rep(seq_len(nfolds), length.out = nrow(x))))
  scheme <- if (!is.null(group)) "leave-one-group-out" else sprintf("%d-fold (within-cohort)", nfolds)

  metric <- if (surv) function(o, p) { np <- -as.numeric(p)
      as.numeric(survival::concordance(o ~ np)$concordance)
    } else function(o, p) { r <- rank(p); n1 <- sum(o == 1); n0 <- sum(o == 0)
      if (n1 == 0 || n0 == 0) return(NA_real_)
      (sum(r[o == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

  res <- list()
  for (m in models) {
    perf <- vapply(names(folds), function(f) {
      te <- folds[[f]]; tr <- setdiff(seq_len(nrow(x)), te)
      p <- try(.fit_predict(m, x[tr, , drop = FALSE], y[tr], x[te, , drop = FALSE], surv), silent = TRUE)
      if (inherits(p, "try-error") || all(!is.finite(p))) return(NA_real_)
      v <- metric(y[te], p)
      if (!surv && is.finite(v)) v <- max(v, 1 - v)     # direction is not the question here
      v
    }, numeric(1))
    res[[m]] <- perf
  }
  R <- do.call(rbind, res)
  best <- apply(R, 1, stats::median, na.rm = TRUE)
  list(per_fold = R,
       per_model = data.frame(model = names(best), performance = unname(best), row.names = NULL),
       ceiling = max(best, na.rm = TRUE),
       best_model = names(which.max(best)),
       scheme = scheme,
       note = paste("A ceiling close to the best signature means the task is",
                    "information-limited; model complexity will not close the gap."))
}

.fit_predict <- function(model, xtr, ytr, xte, surv) {
  if (model == "enet") {
    fam <- if (surv) "cox" else "binomial"
    cv <- glmnet::cv.glmnet(xtr, ytr, family = fam, alpha = 0.5, nfolds = 5, standardize = FALSE)
    as.numeric(stats::predict(cv, xte, s = "lambda.min"))
  } else if (model == "rf") {
    if (!requireNamespace("ranger", quietly = TRUE)) stop("install 'ranger' for model 'rf'")
    if (surv) {
      d <- data.frame(time = ytr[, 1], status = ytr[, 2], xtr, check.names = FALSE)
      f <- ranger::ranger(survival::Surv(time, status) ~ ., data = d, num.trees = 500,
                          verbose = FALSE)
      pr <- stats::predict(f, data.frame(xte, check.names = FALSE))
      rowSums(pr$chf)                                    # cumulative hazard = risk
    } else {
      d <- data.frame(y = factor(ytr), xtr, check.names = FALSE)
      f <- ranger::ranger(y ~ ., data = d, num.trees = 500, probability = TRUE, verbose = FALSE)
      stats::predict(f, data.frame(xte, check.names = FALSE))$predictions[, 2]
    }
  } else if (model == "xgb") {
    if (!requireNamespace("xgboost", quietly = TRUE)) stop("install 'xgboost' for model 'xgb'")
    # xgboost >= 2.1 renamed data/label/eta and dropped verbose from this entry point;
    # the old names still run but warn that they will become errors. Call the training
    # interface directly, which is stable across both.
    lab <- if (surv) ifelse(ytr[, 2] == 1, ytr[, 1], -ytr[, 1]) else as.numeric(ytr)
    obj <- if (surv) "survival:cox" else "binary:logistic"
    b <- xgboost::xgb.train(
      params = list(objective = obj, max_depth = 3, eta = 0.05, subsample = 0.8,
                    colsample_bytree = 0.5, nthread = 2),
      data = xgboost::xgb.DMatrix(data = xtr, label = lab),
      nrounds = 200, verbose = 0)
    as.numeric(stats::predict(b, xte))
  } else stop("unknown model: ", model)
}

#' Compare signatures against the ceiling
#'
#' @param bench Output of [benchmark_signatures()].
#' @param ceiling Output of [estimate_ceiling()].
#' @return A data frame with the gap between each signature and the ceiling, and the
#'   gap between the ceiling and the random null.
#' @export
compare_to_ceiling <- function(bench, ceiling) {
  data.frame(
    signature = bench$signature,
    performance = bench$performance,
    ceiling = ceiling$ceiling,
    gap_to_ceiling = ceiling$ceiling - bench$performance,
    random_median = bench$random_median,
    headroom_above_random = ceiling$ceiling - bench$random_median,
    row.names = NULL
  )[order(-bench$performance), ]
}
