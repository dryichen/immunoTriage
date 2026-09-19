#!/usr/bin/env Rscript
# Rebuild the two fitted data objects the package ships, from the v3 analysis outputs.
#
#   triage_model   the ridge on ten cell-type attribution z-scores, its two
#                         decision cut-offs, and the validation figures quoted in the docs
#   sig_behavior         the measured prognostic and response effects the model was fitted on
#
# Both were previously created inline and had no build script, so they could not be
# regenerated when the underlying tables changed. They now come from icb/triage_v3.py:
# 50 signatures (up from 47, after three signatures lost to a CSV parsing error were
# restored) with both sides of the benchmark against a size- and engine-matched null.
#
#   Rscript data-raw/mkmodel.R      (run from the package root)
suppressMessages(library(jsonlite))
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"

M <- fromJSON(file.path(D, "v3_triage_model.json"))
S <- read.csv(file.path(D, "v3_sigtriage_final.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(S) == M$n_train)

works <- S$d_resp > 0.05
sel   <- S$pred_loso >= M$cut_use
triage_model <- list(
  intercept = M$intercept,
  coef      = setNames(as.numeric(M$coef), M$features),
  cut_use   = M$cut_use,
  cut_maybe = M$cut_maybe,
  validation = list(
    loso_rho             = M$loso_rho,
    loso_auc             = M$loso_auc,
    precision_at_cut_use = mean(works[sel]),
    recall_at_cut_use    = sum(works & sel) / sum(works),
    base_rate            = mean(works),
    n_signatures         = as.integer(M$n_train)))

sig_behavior <- data.frame(
  signature                 = S$signature,
  family                    = S$family,
  top_cell_type             = S$top_cell_type,
  prognostic_effect         = S$d_prog,
  response_effect           = S$d_resp,
  shift                     = S$shift,
  predicted_response_effect = S$pred_loso,
  stringsAsFactors = FALSE)
sig_behavior <- sig_behavior[order(-sig_behavior$response_effect), ]
rownames(sig_behavior) <- NULL

save(triage_model, file = "data/triage_model.rda", compress = "xz")
save(sig_behavior,       file = "data/sig_behavior.rda",       compress = "xz")
cat(sprintf("triage_model: %d features, cut_use %.5f, LOSO rho %.3f, AUC %.3f\n",
            length(triage_model$coef), triage_model$cut_use,
            triage_model$validation$loso_rho, triage_model$validation$loso_auc))
cat(sprintf("  precision at cut_use %.2f, recall %.2f, base rate %.2f, n = %d\n",
            triage_model$validation$precision_at_cut_use,
            triage_model$validation$recall_at_cut_use,
            triage_model$validation$base_rate,
            triage_model$validation$n_signatures))
cat(sprintf("sig_behavior: %d signatures\n", nrow(sig_behavior)))
