library(testthat)

source("../../R/gkrr.R")
source("../../R/gkrr_boot.R")

# ── Shared data helpers ───────────────────────────────────────────────────────

make_clean <- function(n = 80, seed = 1) {
  set.seed(seed)
  x <- runif(n, 0, 10)
  y <- 2 + 3 * x + rnorm(n)
  data.frame(y = y, x = x)
}

make_outlier_y <- function(n = 80, pct = 0.15, seed = 1) {
  set.seed(seed)
  x <- runif(n, 0, 10)
  y <- 2 + 3 * x + rnorm(n)
  k <- floor(n * pct); y[seq_len(k)] <- y[seq_len(k)] + 30
  data.frame(y = y, x = x)
}

# =============================================================================
# BLOCK 1 -- gkrr()
# =============================================================================

test_that("gkrr: clean data converges and coefficients close to OLS (S3)", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  ols <- lm(y ~ x, data = df)
  expect_true(fit$converged)
  expect_equal(coef(fit)[["x"]], coef(ols)[["x"]], tolerance = 0.02)
})

test_that("gkrr: objective function S is non-increasing", {
  df   <- make_clean()
  fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
  expect_true(all(diff(fit$criterion) <= 1e-9))
})

test_that("gkrr: GKRR more robust than OLS under Y-space outliers (S3)", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  ols <- lm(y ~ x, data = df)
  expect_lt(abs(coef(fit)[["x"]] - 3), abs(coef(ols)[["x"]] - 3))
})

test_that("gkrr: kernel weights are in (0, 1]", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  expect_true(all(fit$weights > 0 & fit$weights <= 1))
})

test_that("gkrr: residuals = y - fitted.values", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s2")
  expect_equal(fit$residuals, df$y - fit$fitted.values, tolerance = 1e-10)
})

test_that("gkrr: predict without newdata returns fitted.values", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df)
  expect_equal(predict(fit), fit$fitted.values)
})

test_that("gkrr: predict with newdata uses correct coefficients", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df)
  nd  <- data.frame(x = c(1, 5, 9))
  expect_equal(predict(fit, newdata = nd),
               coef(fit)[1] + coef(fit)[2] * nd$x, tolerance = 1e-10)
})

test_that("gkrr: two predictors -- correct coefficient count", {
  set.seed(7)
  n  <- 100
  df <- data.frame(x1 = runif(n), x2 = runif(n),
                   y  = 1 + 2 * runif(n) - runif(n) + rnorm(n, 0, 0.3))
  fit <- gkrr(y ~ x1 + x2, data = df, sigma_method = "s3")
  expect_length(coef(fit), 3L)
  expect_true(fit$converged)
})

test_that("gkrr: fixed numeric gamma^2 is accepted", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df, sigma_method = 1.5)
  expect_equal(fit$gamma2, 1.5)
})

test_that("gkrr: error when n < number of parameters", {
  df <- data.frame(y = 1:2, x = 1:2)
  expect_error(gkrr(y ~ x, data = df), "n < number")
})

test_that("gkrr: error on NAs", {
  df      <- make_clean()
  df$y[1] <- NA
  expect_error(gkrr(y ~ x, data = df))
})

test_that("sigma2_s1: positive and finite", {
  y <- rnorm(50, 10); yhat <- y + rnorm(50, 0, .5)
  g <- sigma2_s1(y, yhat)
  expect_gt(g, 0); expect_true(is.finite(g))
})

test_that("sigma2_s2: positive and finite", {
  y <- rnorm(50, 10); yhat <- y + rnorm(50, 0, .5)
  g <- sigma2_s2(y, yhat)
  expect_gt(g, 0); expect_true(is.finite(g))
})

test_that("sigma2_s3: positive and finite", {
  y <- rnorm(50, 10); yhat <- y + rnorm(50, 0, .5)
  g <- sigma2_s3(y, yhat, p = 1L)
  expect_gt(g, 0); expect_true(is.finite(g))
})

# =============================================================================
# BLOCK 2 -- gkrr() boot argument
# =============================================================================

test_that("gkrr: boot = FALSE stores NULL", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3", boot = FALSE)
  expect_null(fit$boot)
})

test_that("gkrr: boot = TRUE stores a gkrr_boot object", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3",
              boot = TRUE, boot_args = list(B = 49, seed = 1))
  expect_s3_class(fit$boot, "gkrr_boot")
})

test_that("gkrr: boot = TRUE with boot_args passes parameters correctly", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3",
              boot = TRUE,
              boot_args = list(B = 49, type = "percentile", seed = 99))
  expect_equal(fit$boot$type, "percentile")
  expect_lte(fit$boot$B, 49L)
})

test_that("gkrr: boot accepts a pre-computed gkrr_boot object", {
  df   <- make_outlier_y()
  fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
  boot <- gkrr_boot(fit, B = 49, seed = 2)
  fit2 <- gkrr(y ~ x, data = df, sigma_method = "s3", boot = boot)
  expect_identical(fit2$boot, boot)
})

test_that("gkrr: invalid boot argument raises error", {
  df <- make_clean()
  expect_error(gkrr(y ~ x, data = df, boot = "yes"), "'boot' must be")
})

# =============================================================================
# BLOCK 3 -- summary.gkrr
# =============================================================================

test_that("summary.gkrr: runs without boot and returns invisibly", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df)
  out <- capture.output(s <- summary(fit))
  expect_type(s, "list")
  expect_named(s, c("coefficients", "r_squared"))
})

test_that("summary.gkrr: coefficient table has sandwich columns without boot", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df)
  s   <- suppressMessages(capture.output(sm <- summary(fit)))
  # Sandwich is always computed: Estimate, Std. Error, CI lower, CI upper, p-value
  expect_true(all(c("Estimate", "Std. Error", "p-value") %in%
                    colnames(sm$coefficients)))
})

test_that("summary.gkrr: with boot adds SE, CI and p-value columns", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3",
              boot = TRUE, boot_args = list(B = 49, seed = 1))
  s   <- suppressMessages(capture.output(sm <- summary(fit)))
  # Table must have Estimate, Boot.SE, two CI cols, p-value
  expect_gte(ncol(sm$coefficients), 5L)
})

test_that("summary.gkrr: boot argument overrides stored boot", {
  df   <- make_outlier_y()
  fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
  boot <- gkrr_boot(fit, B = 49, type = "normal", seed = 5)
  s    <- suppressMessages(capture.output(sm <- summary(fit, boot = boot)))
  expect_gte(ncol(sm$coefficients), 5L)
})

test_that("summary.gkrr: r_squared values are in [0, 1]", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df)
  s   <- suppressMessages(capture.output(sm <- summary(fit)))
  expect_gte(sm$r_squared["r2"], 0)
  expect_lte(sm$r_squared["r2"], 1)
})

# =============================================================================
# BLOCK 4 -- gkrr_boot()
# =============================================================================

test_that("gkrr_boot: returns a 'gkrr_boot' object with required fields", {
  df   <- make_outlier_y()
  fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
  boot <- gkrr_boot(fit, B = 49, type = "percentile", seed = 1)
  expect_s3_class(boot, "gkrr_boot")
  expect_named(boot, c("t0","t","se","bias","ci","pval",
                        "conf","type","B","B_failed","call"))
})

test_that("gkrr_boot: t matrix dimensions are B_ok x p_coef", {
  df   <- make_outlier_y()
  fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
  boot <- gkrr_boot(fit, B = 49, type = "percentile", seed = 2)
  expect_equal(ncol(boot$t), length(boot$t0))
  expect_lte(nrow(boot$t), 49L)
})

test_that("gkrr_boot: SE is positive and finite", {
  df   <- make_outlier_y()
  fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
  boot <- gkrr_boot(fit, B = 49, seed = 3)
  expect_true(all(boot$se > 0 & is.finite(boot$se)))
})

test_that("gkrr_boot: CI lower < upper", {
  df   <- make_outlier_y()
  fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
  boot <- gkrr_boot(fit, B = 49, seed = 4)
  expect_true(all(boot$ci["lower", ] < boot$ci["upper", ]))
})

test_that("gkrr_boot: p-values in (0, 1]", {
  df   <- make_outlier_y()
  fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
  boot <- gkrr_boot(fit, B = 49, seed = 4)
  expect_true(all(boot$pval > 0 & boot$pval <= 1))
})

test_that("gkrr_boot: all three CI types produce lower < upper", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  for (tp in c("normal", "percentile", "bca")) {
    boot <- gkrr_boot(fit, B = 49, type = tp, seed = 5)
    expect_true(all(boot$ci["lower",] < boot$ci["upper",]), info = tp)
  }
})

test_that("gkrr_boot: BCa CI valid for small sample", {
  df   <- make_outlier_y(n = 40, pct = 0.10, seed = 10)
  fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
  boot <- gkrr_boot(fit, B = 49, type = "bca", seed = 6)
  expect_true(all(boot$ci["lower",] < boot$ci["upper",]))
})


test_that("gkrr_boot: seed guarantees reproducibility", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  b1  <- gkrr_boot(fit, B = 49, type = "percentile", seed = 99)
  b2  <- gkrr_boot(fit, B = 49, type = "percentile", seed = 99)
  expect_equal(b1$t, b2$t)
  expect_equal(b1$ci, b2$ci)
})

test_that("gkrr_boot: error when object is not 'gkrr'", {
  expect_error(gkrr_boot(lm(y ~ x, data = make_clean())), "classe|class")
})

test_that("gkrr_boot: coverage reasonable in small Monte Carlo", {
  # 200 MC runs with n=60, 10% outliers; expect empirical coverage in [85%, 100%]
  set.seed(42)
  B_mc  <- 200L
  cover <- logical(B_mc)
  for (mc in seq_len(B_mc)) {
    x  <- runif(60, 0, 10)
    y  <- 2 + 3 * x + rnorm(60)
    y[1:6] <- y[1:6] + 25
    df   <- data.frame(y = y, x = x)
    fit  <- gkrr(y ~ x, data = df, sigma_method = "s3")
    # suppressWarnings: occasional non-convergence in resampled data is
    # expected and tested separately; suppress here to keep output clean.
    boot <- tryCatch(
      suppressWarnings(gkrr_boot(fit, B = 99L, type = "bca", seed = mc)),
      error = function(e) NULL
    )
    if (!is.null(boot))
      cover[mc] <- boot$ci["lower", "x"] <= 3 & 3 <= boot$ci["upper", "x"]
  }
  coverage <- mean(cover)
  expect_gte(coverage, 0.85)
  expect_lte(coverage, 1.00)
})

# =============================================================================
# BLOCK 5 -- plot.gkrr() and plot.gkrr_boot()
# =============================================================================

.plot_fixture <- function() {
  set.seed(42)
  n  <- 80
  x  <- runif(n, 0, 10)
  y  <- 2 + 3 * x + rnorm(n)
  y[1:12] <- y[1:12] + 25
  df  <- data.frame(y = y, x = x)
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  list(fit = fit, df = df)
}

test_that("plot.gkrr: all 6 panels render without error", {
  obj <- .plot_fixture()
  pdf(NULL); on.exit(dev.off())
  for (w in 1:6)
    expect_no_error(plot(obj$fit, which = w, ask = FALSE))
})

test_that("plot.gkrr: returns the object invisibly", {
  obj    <- .plot_fixture()
  pdf(NULL); on.exit(dev.off())
  result <- plot(obj$fit, which = 1, ask = FALSE)
  expect_identical(result, obj$fit)
})

test_that("plot.gkrr: n_id = 0 renders without error", {
  obj <- .plot_fixture()
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot(obj$fit, which = 1, n_id = 0L, ask = FALSE))
})

test_that("plot.gkrr_boot panel 1: histograms render without error", {
  obj  <- .plot_fixture()
  boot <- gkrr_boot(obj$fit, B = 49, type = "bca", seed = 1)
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot(boot, which = 1, ask = FALSE))
})

test_that("plot.gkrr_boot panel 2: scatter matrix with 2 predictors", {
  set.seed(7)
  df2 <- data.frame(x1 = runif(80), x2 = runif(80),
                    y  = 1 + 2*runif(80) + rnorm(80, 0, .3))
  df2$y[1:12] <- df2$y[1:12] + 10
  fit2  <- gkrr(y ~ x1 + x2, data = df2, sigma_method = "s3")
  boot2 <- gkrr_boot(fit2, B = 99, type = "percentile", seed = 2)
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot(boot2, which = 2, ask = FALSE))
})

test_that("plot.gkrr_boot panel 2: emits message with p_coef < 2", {
  obj   <- .plot_fixture()
  fit1  <- gkrr(y ~ x - 1, data = obj$df, sigma_method = "s3")
  boot1 <- gkrr_boot(fit1, B = 49, seed = 3)
  expect_equal(ncol(boot1$t), 1L)
  pdf(NULL); on.exit(dev.off())
  expect_message(plot(boot1, which = 2, ask = FALSE), "requires")
})

test_that("plot.gkrr_boot: returns object invisibly", {
  obj  <- .plot_fixture()
  boot <- gkrr_boot(obj$fit, B = 49, seed = 4)
  pdf(NULL); on.exit(dev.off())
  expect_identical(plot(boot, which = 1, ask = FALSE), boot)
})

# =============================================================================
# BLOCK 6 -- Sandwich variance estimator
# =============================================================================

test_that("gkrr: sandwich object is always present in the fit", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df)
  expect_true(!is.null(fit$sandwich))
  expect_named(fit$sandwich, c("vcov","se","tval","pval","ci_lo","ci_hi","conf"))
})

test_that("gkrr: sandwich SE is positive and finite", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  expect_true(all(fit$sandwich$se > 0))
  expect_true(all(is.finite(fit$sandwich$se)))
})

test_that("gkrr: sandwich p-values in (0, 1]", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  expect_true(all(fit$sandwich$pval > 0 & fit$sandwich$pval <= 1))
})

test_that("gkrr: sandwich CI lower < upper", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df)
  expect_true(all(fit$sandwich$ci_lo < fit$sandwich$ci_hi))
})

test_that("gkrr: sandwich CI contains true parameter in clean data", {
  # With clean data and n=100, sandwich 95% CI should contain beta1=3
  df  <- make_clean(n = 100)
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  expect_gte(fit$sandwich$ci_hi["x"], 3)
  expect_lte(fit$sandwich$ci_lo["x"], 3)
})

test_that("gkrr: vcov() returns symmetric positive-definite matrix", {
  df  <- make_clean()
  fit <- gkrr(y ~ x, data = df)
  V   <- vcov(fit)
  expect_equal(V, t(V), tolerance = 1e-10)         # symmetric
  expect_true(all(eigen(V)$values > 0))             # positive definite
})

test_that("gkrr: sandwich vcov dimensions match number of coefficients", {
  set.seed(5)
  df  <- data.frame(x1 = runif(80), x2 = runif(80),
                    y  = 1 + 2*runif(80) + rnorm(80, 0, .3))
  fit <- gkrr(y ~ x1 + x2, data = df, sigma_method = "s3")
  p   <- length(coef(fit))
  expect_equal(dim(vcov(fit)), c(p, p))
  expect_equal(length(fit$sandwich$se), p)
})

test_that("gkrr: summary uses sandwich when boot = FALSE", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  s   <- suppressMessages(capture.output(sm <- summary(fit)))
  # Column names should include Std. Error (sandwich), not Boot.SE
  expect_true("Std. Error" %in% colnames(sm$coefficients))
  expect_false("Boot.SE"   %in% colnames(sm$coefficients))
})

test_that("gkrr: summary uses bootstrap when boot object present", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3",
              boot = TRUE, boot_args = list(B = 49, seed = 1))
  s   <- suppressMessages(capture.output(sm <- summary(fit)))
  # With bootstrap, Std. Error column comes from boot$se
  expect_true("Std. Error" %in% colnames(sm$coefficients))
  # Bootstrap SE and sandwich SE should differ (bootstrap accounts for gamma^2)
  expect_false(isTRUE(all.equal(
    sm$coefficients[, "Std. Error"],
    fit$sandwich$se, tolerance = 1e-6
  )))
})

test_that("gkrr: conf argument changes sandwich CI width", {
  df   <- make_clean()
  fit1 <- gkrr(y ~ x, data = df, conf = 0.90)
  fit2 <- gkrr(y ~ x, data = df, conf = 0.99)
  # Wider conf -> wider CI
  width1 <- fit1$sandwich$ci_hi["x"] - fit1$sandwich$ci_lo["x"]
  width2 <- fit2$sandwich$ci_hi["x"] - fit2$sandwich$ci_lo["x"]
  expect_lt(width1, width2)
})

test_that("gkrr: sandwich coverage reasonable in Monte Carlo", {
  # 200 MC runs, n=80, clean data; expect empirical coverage near 95%
  set.seed(99)
  cover <- logical(200L)
  for (mc in seq_len(200L)) {
    x   <- runif(80, 0, 10)
    y   <- 2 + 3*x + rnorm(80)
    df  <- data.frame(y = y, x = x)
    fit <- suppressWarnings(
      gkrr(y ~ x, data = df, sigma_method = "s1")
    )
    cover[mc] <- fit$sandwich$ci_lo["x"] <= 3 & 3 <= fit$sandwich$ci_hi["x"]
  }
  coverage <- mean(cover)
  expect_gte(coverage, 0.85)
  expect_lte(coverage, 1.00)
})

test_that("summary.gkrr: no warning when only sandwich is available", {
  df  <- make_outlier_y()
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3")
  expect_no_warning(
    capture.output(summary(fit))
  )
})

test_that("summary.gkrr: no divergence warning when SEs agree", {
  # n=200, mild outliers: sandwich and bootstrap should agree
  set.seed(1)
  df <- data.frame(x = runif(200, 0, 10))
  df$y <- 2 + 3*df$x + rnorm(200); df$y[1:10] <- df$y[1:10] + 5
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3",
              boot = TRUE, boot_args = list(B = 199, seed = 1))
  expect_no_warning(capture.output(summary(fit)))
})

test_that("summary.gkrr: divergence warning fired when SEs differ a lot", {
  # n=30, 30% severe outliers: sandwich should underestimate SE
  set.seed(7)
  df <- data.frame(x = runif(30, 0, 10))
  df$y <- 2 + 3*df$x + rnorm(30); df$y[1:9] <- df$y[1:9] + 40
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3",
              boot = TRUE, boot_args = list(B = 199, seed = 7))
  expect_warning(
    capture.output(summary(fit)),
    regexp = "diverge"
  )
})

test_that("summary.gkrr: se_tol = Inf suppresses divergence warning", {
  set.seed(7)
  df <- data.frame(x = runif(30, 0, 10))
  df$y <- 2 + 3*df$x + rnorm(30); df$y[1:9] <- df$y[1:9] + 40
  fit <- gkrr(y ~ x, data = df, sigma_method = "s3",
              boot = TRUE, boot_args = list(B = 199, seed = 7))
  expect_no_warning(
    capture.output(summary(fit, se_tol = Inf))
  )
})

# =============================================================================
# BLOCK 7 -- sigma2_s4() and sigma_method = "auto"
# =============================================================================

test_that("sigma2_s4: returns a positive finite value", {
  set.seed(1)
  y    <- rnorm(100, 10); yhat <- y + rnorm(100, 0, 0.5)
  g    <- sigma2_s4(y, yhat)
  expect_gt(g, 0)
  expect_true(is.finite(g))
})

test_that("gkrr: sigma_method = 's4' produces a valid fit", {
  df  <- make_outlier_y(n = 100)
  fit <- gkrr(y ~ x, data = df, sigma_method = "s4")
  expect_s3_class(fit, "gkrr")
  expect_true(fit$converged)
  expect_true(all(is.finite(coef(fit))))
  expect_match(fit$sigma_method, "^s4$")
})

test_that("gkrr: sigma_method = 's4' gamma2 is positive", {
  df  <- make_outlier_y(n = 100)
  fit <- gkrr(y ~ x, data = df, sigma_method = "s4")
  expect_gt(fit$gamma2, 0)
})

test_that("gkrr: sigma_method = 'auto' selects one of s1/s2/s3/s4", {
  df  <- make_outlier_y(n = 80, seed = 5)
  fit <- suppressMessages(
    gkrr(y ~ x, data = df, sigma_method = "auto",
         auto_args = list(B = 19, seed = 1))
  )
  expect_s3_class(fit, "gkrr")
  expect_true(fit$converged)
  expect_match(fit$sigma_method, "^auto\\((s1|s2|s3|s4)\\)$")
})

test_that("gkrr: 'auto' emits message about selection", {
  df <- make_outlier_y(n = 80, seed = 5)
  expect_message(
    gkrr(y ~ x, data = df, sigma_method = "auto",
         auto_args = list(B = 19, seed = 1)),
    "auto"
  )
})

test_that("gkrr: 'auto' with auto_args$seed is reproducible", {
  df   <- make_outlier_y(n = 80, seed = 5)
  fit1 <- suppressMessages(
    gkrr(y ~ x, data = df, sigma_method = "auto",
         auto_args = list(B = 19, seed = 42))
  )
  fit2 <- suppressMessages(
    gkrr(y ~ x, data = df, sigma_method = "auto",
         auto_args = list(B = 19, seed = 42))
  )
  expect_equal(fit1$sigma_method, fit2$sigma_method)
  expect_equal(coef(fit1), coef(fit2))
})

test_that("gkrr: auto + boot = TRUE works without error", {
  df  <- make_outlier_y(n = 80, seed = 5)
  # This was failing because "auto(sX)" was not recognized by .gkrr_fit
  expect_no_error(suppressMessages(
    gkrr(y ~ x, data = df,
         sigma_method = "auto",
         boot = TRUE,
         boot_args = list(B = 19, seed = 1),
         auto_args  = list(B = 19, seed = 1))
  ))
})

test_that("gkrr_boot: works when sigma_method is 'auto(sX)'", {
  df   <- make_outlier_y(n = 80, seed = 5)
  fit  <- suppressMessages(
    gkrr(y ~ x, data = df, sigma_method = "auto",
         auto_args = list(B = 19, seed = 1))
  )
  # sigma_method stored as "auto(s?)" -- gkrr_boot must resolve it
  expect_match(fit$sigma_method, "^auto\\(s[1-4]\\)$")
  boot <- expect_no_error(gkrr_boot(fit, B = 49, seed = 2))
  expect_gte(boot$B, 10L)
})
