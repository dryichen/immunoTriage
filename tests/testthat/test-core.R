test_that("scoring returns one column per scorable signature", {
  data(gc_ici, gastro_signatures, package = "immunoTriage")
  S <- score_signatures(gc_ici$expr, gastro_signatures[c(1,2,4,5)], default_method = "zmean")
  expect_equal(nrow(S), ncol(gc_ici$expr))
  expect_true(all(is.finite(S)))
})

test_that("benchmark reports a null and a delta", {
  data(gc_ici, gastro_signatures, package = "immunoTriage")
  b <- benchmark_signatures(gc_ici$expr, gc_ici$response,
                            gastro_signatures[c(1,2,4,5)], n_random = 30, seed = 1)
  expect_true(all(c("performance","random_median","delta","beats_null") %in% names(b)))
  expect_true(all(b$performance >= 0.5))          # direction-agnostic metric
})

test_that("random signatures do not systematically beat their own null", {
  data(gc_ici, package = "immunoTriage")
  set.seed(7)
  fake <- lapply(1:6, function(i) sample(rownames(gc_ici$expr), 20))
  names(fake) <- paste0("random", 1:6)
  b <- benchmark_signatures(gc_ici$expr, gc_ici$response, fake, n_random = 100, seed = 2)
  expect_lt(mean(b$beats_null), 0.5)
})

test_that("ceiling estimation runs and is bounded", {
  data(gc_ici, package = "immunoTriage")
  x <- t(gc_ici$expr[1:300, ])
  cl <- suppressWarnings(estimate_ceiling(x, gc_ici$response, models = "enet", nfolds = 4, seed = 1))
  expect_true(is.finite(cl$ceiling))
  expect_lte(cl$ceiling, 1)
})
