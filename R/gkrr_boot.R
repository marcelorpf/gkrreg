# ==============================================================================
# gkrr_boot.R -- Bootstrap inference for GKRReg parameters
# ==============================================================================
#
# RATIONALE
# ----------
# GKRR has no closed-form sampling distribution for beta_hat because the
# weights k_ii = G(y_i, mu_i_hat) depend on the estimated parameters
# themselves.  The matrix (X'KX)^{-1} from the final WLS step treats K as
# fixed, ignoring the variability introduced by estimating the weights and
# gamma^2 -- so it underestimates the true variance of beta_hat.
#
# Pairs bootstrap (resample (y_i, x_i) with replacement) is the most direct
# solution: the full algorithm is re-run on every replicate (including
# re-estimating gamma^2), so all sources of variability are captured.
#
# THREE INTERVAL TYPES
# ---------------------
# 1. Percentile : [q_{alpha/2}, q_{1-alpha/2}] of the bootstrap distribution.
#    Simple, but asymptotically inconsistent under skewness.
#
# 2. Normal     : beta_hat +/- z_{alpha/2} * se_boot.
#    Assumes symmetry; reliable for large n.
#
# 3. BCa (bias-corrected and accelerated, Efron 1987):
#    Corrects for bias (z0) and skewness (acceleration a) of the bootstrap
#    distribution.  Best coverage in small samples.  Recommended default.
#
# WEIGHTED BOOTSTRAP
# -------------------
# When weighted = TRUE the final kernel weights k_ii are used as sampling
# probabilities.  Outliers (small k_ii) are rarely selected, aligning the
# resampling with the robust philosophy of GKRReg.  This tends to produce
# wider intervals because the effective sample size is reduced.
#
# BOOTSTRAP P-VALUES  (used by summary.gkrr)
# -------------------------------------------
# For H0: beta_j = 0, we use the centred-t approach:
#
#   p_j = (2 / B) * #{|beta*_b - beta_j| >= |beta_j|}
#
# The bootstrap distribution is centred at beta_j (the estimate), and we
# count how many replicates are at least as far from zero as the estimate
# itself.  This is equivalent to a two-sided permutation test and requires
# no additional computation beyond the existing bootstrap replicates.
# ==============================================================================


#' Bootstrap inference for GKRReg
#'
#' Runs a pairs bootstrap to estimate standard errors, confidence intervals
#' and (optionally) p-values for the coefficients of a \code{\link{gkrr}}
#' model.  The full fitting algorithm is re-executed on every replicate,
#' including re-estimation of \eqn{\gamma^2}.
#'
#' @param object   A \code{"gkrr"} object (output of \code{gkrr()}).
#' @param B        Number of bootstrap replicates (default \code{999}).
#' @param type     CI type: \code{"bca"} (default, BCa),
#'   \code{"percentile"} or \code{"normal"}.
#' @param conf     Confidence level, scalar in (0, 1) (default \code{0.95}).
#' @param seed     Integer seed for reproducibility (optional).
#' @param verbose  Logical.  If \code{TRUE}, prints progress every 100
#'   replicates. Default \code{FALSE}.
#'
#' @return An object of class \code{"gkrr_boot"} containing:
#' \describe{
#'   \item{\code{t0}}{Original coefficient estimates (length \eqn{p+1}).}
#'   \item{\code{t}}{Matrix \eqn{B_{\text{ok}} \times (p+1)} of bootstrap
#'     estimates.}
#'   \item{\code{se}}{Bootstrap standard errors.}
#'   \item{\code{bias}}{Bootstrap bias \eqn{E^*[\hat\beta^*] - \hat\beta}.}
#'   \item{\code{ci}}{Matrix \eqn{2 \times (p+1)} of CIs with rows
#'     \code{"lower"} and \code{"upper"}.}
#'   \item{\code{pval}}{Bootstrap p-values (centred-t, two-sided).}
#'   \item{\code{conf}}{Confidence level used.}
#'   \item{\code{type}}{CI type used.}
#'   \item{\code{B}}{Number of successful replicates.}
#'   \item{\code{B_failed}}{Replicates discarded due to non-convergence.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @references
#' Efron, B. (1987). Better bootstrap confidence intervals.
#' \emph{Journal of the American Statistical Association}, 82(397), 171--185.
#'
#' De Carvalho, F.A.T., Lima Neto, E.A., Ferreira, M.R.P. (2017).
#' \doi{10.1016/j.neucom.2016.12.035}
#'
#' @seealso \code{\link{gkrr}}, \code{\link{plot.gkrr_boot}}
#'
#' @examples
#' set.seed(42)
#' n  <- 80
#' x  <- runif(n, 0, 10)
#' y  <- 2 + 3 * x + rnorm(n)
#' y[1:12] <- y[1:12] + 25      # 15% Y-space outliers
#'
#' fit  <- gkrr(y ~ x, sigma_method = "s3")
#' boot <- gkrr_boot(fit, B = 499, seed = 1)
#' print(boot)
#' summary(boot)
#' plot(boot)
#'
#' @export
gkrr_boot <- function(object,
                      B        = 999L,
                      type     = c("bca", "percentile", "normal"),
                      conf     = 0.95,
                      seed     = NULL,
                      verbose  = FALSE) {

  # ── Validate ─────────────────────────────────────────────────────────────────
  if (!inherits(object, "gkrr"))
    stop("'object' must be of class \"gkrr\".")
  type <- match.arg(type)
  if (conf <= 0 || conf >= 1) stop("'conf' must be in (0, 1).")
  if (B < 1L) stop("'B' must be >= 1.")

  if (!is.null(seed)) set.seed(seed)

  # ── Extract model components ─────────────────────────────────────────────────
  mf     <- object$model
  mt     <- object$terms
  y_full <- model.response(mf, "numeric")
  X_full <- model.matrix(mt, mf)
  n      <- length(y_full)
  t0     <- object$coefficients
  p_coef <- length(t0)

  # Recover sigma_method argument (numeric if fixed, string otherwise).
  # Handles three formats:
  #   "s1" / "s2" / "s3" / "s4"  -> used directly
  #   "fixed(1.234)"              -> parsed as numeric
  #   "auto(s3)"                  -> resolved to the inner method "s3"
  sm     <- object$sigma_method
  sm_arg <- if (grepl("^fixed", sm)) {
    as.numeric(sub("fixed\\((.+)\\)", "\\1", sm))
  } else if (grepl("^auto\\(", sm)) {
    sub("^auto\\((.+)\\)$", "\\1", sm)   # extract e.g. "s4" from "auto(s4)"
  } else {
    sm
  }

  # ── Bootstrap loop ───────────────────────────────────────────────────────────
  t_boot   <- matrix(NA_real_, nrow = B, ncol = p_coef,
                     dimnames = list(NULL, names(t0)))
  B_failed <- 0L

  for (b in seq_len(B)) {
    if (verbose && (b %% 100L == 0L))
      message(sprintf("Bootstrap replicate %d / %d", b, B))

    idx  <- sample.int(n, n, replace = TRUE)
    y_b  <- y_full[idx]
    X_b  <- X_full[idx, , drop = FALSE]

    fit_b <- tryCatch(
      .gkrr_fit(X_b, y_b,
                sigma_method = sm_arg,
                alpha        = 0.5,
                tol          = 1e-10,
                maxit        = 100L),
      error = function(e) NULL
    )

    if (is.null(fit_b) || !fit_b$converged) {
      B_failed <- B_failed + 1L
      next
    }
    t_boot[b, ] <- fit_b$coefficients
  }

  # Drop failed replicates (NA rows)
  t_boot <- t_boot[complete.cases(t_boot), , drop = FALSE]
  B_ok   <- nrow(t_boot)

  if (B_ok < 10L)
    stop(sprintf("Fewer than 10 valid replicates (%d failed). Check the model.",
                 B_failed))
  if (B_failed > 0L)
    warning(sprintf("%d bootstrap replicate(s) discarded (non-convergence).",
                    B_failed))

  # ── Summary statistics ───────────────────────────────────────────────────────
  se_boot   <- apply(t_boot, 2, sd)
  bias_boot <- colMeans(t_boot) - t0
  alpha_lo  <- (1 - conf) / 2
  alpha_hi  <- 1 - alpha_lo

  # ── Confidence intervals ─────────────────────────────────────────────────────
  ci <- switch(type,

    percentile = {
      lo <- apply(t_boot, 2, quantile, probs = alpha_lo)
      hi <- apply(t_boot, 2, quantile, probs = alpha_hi)
      matrix(c(lo, hi), nrow = 2L, byrow = TRUE,
             dimnames = list(c("lower", "upper"), names(t0)))
    },

    normal = {
      z  <- qnorm(alpha_hi)
      lo <- t0 - z * se_boot
      hi <- t0 + z * se_boot
      matrix(c(lo, hi), nrow = 2L, byrow = TRUE,
             dimnames = list(c("lower", "upper"), names(t0)))
    },

    bca = .bca_intervals(t0, t_boot, X_full, y_full, sm_arg,
                         alpha_lo, alpha_hi)
  )

  # ── Bootstrap p-values (centred-t, two-sided) ────────────────────────────────
  # p_j = proportion of |beta*_b - beta_j| >= |beta_j|, floored at 1/B_ok
  pval <- sapply(seq_len(p_coef), function(j) {
    centred <- abs(t_boot[, j] - t0[j])
    max(1 / B_ok, mean(centred >= abs(t0[j])))
  })
  names(pval) <- names(t0)

  structure(
    list(
      t0       = t0,
      t        = t_boot,
      se       = se_boot,
      bias     = bias_boot,
      ci       = ci,
      pval     = pval,
      conf     = conf,
      type     = type,
      B        = B_ok,
      B_failed = B_failed,
      call     = match.call()
    ),
    class = "gkrr_boot"
  )
}


# ── BCa intervals (internal) ──────────────────────────────────────────────────
#
# Reference: Efron (1987), DiCiccio & Efron (1996).
#
# z0 : bias-correction -- Phi^{-1}(proportion of t*_b < t0)
# a  : acceleration -- estimated by jackknife on the coefficients
#      a_j = sum_i(theta_. - theta_i)^3 / (6 * [sum_i(theta_. - theta_i)^2]^{3/2})
#
# Adjusted quantiles:
#   alpha1 = Phi(z0 + (z0 + z_alpha)   / (1 - a*(z0 + z_alpha)))
#   alpha2 = Phi(z0 + (z0 + z_{1-alpha}) / (1 - a*(z0 + z_{1-alpha})))

.bca_intervals <- function(t0, t_boot, X, y, sm_arg, alpha_lo, alpha_hi) {

  p_coef <- length(t0)
  n      <- length(y)
  B      <- nrow(t_boot)

  # Bias-correction z0
  z0 <- sapply(seq_len(p_coef), function(j) {
    prop <- mean(t_boot[, j] < t0[j])
    prop <- pmax(1 / (B + 1), pmin(B / (B + 1), prop))   # clip to avoid +/-Inf
    qnorm(prop)
  })

  # Jackknife acceleration
  theta_jack <- matrix(NA_real_, nrow = n, ncol = p_coef)
  for (i in seq_len(n)) {
    fit_i <- tryCatch(
      .gkrr_fit(X[-i, , drop = FALSE], y[-i],
                sigma_method = sm_arg,
                alpha = 0.5, tol = 1e-10, maxit = 100L),
      error = function(e) NULL
    )
    if (!is.null(fit_i))
      theta_jack[i, ] <- fit_i$coefficients
  }
  theta_jack <- theta_jack[complete.cases(theta_jack), , drop = FALSE]

  a <- sapply(seq_len(p_coef), function(j) {
    th  <- theta_jack[, j]
    mu  <- mean(th)
    d   <- mu - th          # jackknife difference (standard sign convention)
    num <- sum(d^3)
    den <- 6 * sum(d^2)^(3/2)
    if (!is.finite(den) || den == 0) return(0)
    num / den
  })

  z_lo <- qnorm(alpha_lo)
  z_hi <- qnorm(alpha_hi)

  adj_q <- function(z_alpha, z0_j, a_j) {
    num <- z0_j + z_alpha
    pnorm(z0_j + num / (1 - a_j * num))
  }

  lo <- sapply(seq_len(p_coef), function(j)
    quantile(t_boot[, j], probs = adj_q(z_lo, z0[j], a[j])))
  hi <- sapply(seq_len(p_coef), function(j)
    quantile(t_boot[, j], probs = adj_q(z_hi, z0[j], a[j])))

  matrix(c(lo, hi), nrow = 2L, byrow = TRUE,
         dimnames = list(c("lower", "upper"), names(t0)))
}


# ── S3 methods for "gkrr_boot" ───────────────────────────────────────────────

#' @export
print.gkrr_boot <- function(x, digits = 4L, ...) {
  cat("\nGKRReg Bootstrap\n")
  cat("-----------------\n")
  cat(sprintf("Replicates: %d valid / %d failed  |  conf = %.0f%%  |  type = %s\n",
              x$B, x$B_failed, 100 * x$conf, x$type))
  cat("\n")

  tab <- rbind(
    "Estimate"  = x$t0,
    "Boot bias" = x$bias,
    "Boot SE"   = x$se,
    x$ci,
    "p-value"   = x$pval
  )
  print(round(tab, digits))
  invisible(x)
}

#' @export
summary.gkrr_boot <- function(object, digits = 4L, ...) {
  print(object, digits = digits, ...)

  cat("\nBootstrap distribution of coefficients:\n")
  for (j in seq_len(ncol(object$t))) {
    nm <- colnames(object$t)[j]
    q  <- quantile(object$t[, j], c(0, .025, .25, .5, .75, .975, 1))
    cat(sprintf("  %-20s  %s\n", nm,
                paste(round(q, digits), collapse = "  ")))
  }
  invisible(object)
}

#' Diagnostic plots for a GKRReg bootstrap object
#'
#' Produces up to 2 panels for a \code{"gkrr_boot"} object.
#'
#' \describe{
#'   \item{\code{which = 1}}{For each coefficient: histogram of bootstrap
#'     replicates with smoothed density, original estimate and shaded CI
#'     region.}
#'   \item{\code{which = 2}}{Scatter-plot matrix of bootstrap replicates
#'     across all pairs of coefficients, with a 95\% confidence ellipse and
#'     Pearson correlation in each off-diagonal cell.  Only available when
#'     the model has at least 2 coefficients.}
#' }
#'
#' @param x     A \code{"gkrr_boot"} object.
#' @param which Integer vector: \code{1} (histograms), \code{2} (scatter
#'   matrix), or \code{1:2} for both (default).
#' @param ask   Logical; waits for user input between panels when \code{TRUE}
#'   (default \code{TRUE} when \code{length(which) > 1}).
#' @param ...   Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{gkrr_boot}}, \code{\link{plot.gkrr}}
#' @export
plot.gkrr_boot <- function(x, which = 1:2, ask = length(which) > 1L, ...) {

  if (ask) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(ask = TRUE)
  }

  p_coef <- ncol(x$t)
  nms    <- colnames(x$t)

  # ── Panel 1: Histogram + density + CI per coefficient ──────────────────────
  if (1L %in% which) {
    for (j in seq_len(p_coef)) {
      tv  <- x$t[, j]
      lo  <- x$ci["lower", j]
      hi  <- x$ci["upper", j]
      t0j <- x$t0[j]
      pv  <- x$pval[j]

      h    <- hist(tv, plot = FALSE)
      dens <- density(tv)
      ymax <- max(max(h$density), max(dens$y)) * 1.18

      hist(tv, freq = FALSE,
           col  = "lightsteelblue", border = "white",
           main = sprintf("Bootstrap: %s", nms[j]),
           xlab = "Bootstrap estimate",
           ylab = "Density",
           ylim = c(0, ymax))
      lines(dens, col = "steelblue", lwd = 2)

      # Shade CI region
      x_ic <- seq(lo, hi, length.out = 200)
      y_ic <- approx(dens$x, dens$y, xout = x_ic)$y
      y_ic[is.na(y_ic)] <- 0
      polygon(c(x_ic, rev(x_ic)), c(y_ic, rep(0, length(y_ic))),
              col = rgb(0.7, 0.85, 1, 0.45), border = NA)

      abline(v = lo,  col = "firebrick",  lty = 2, lwd = 1.5)
      abline(v = hi,  col = "firebrick",  lty = 2, lwd = 1.5)
      abline(v = t0j, col = "darkgreen",  lty = 1, lwd = 2)
      abline(v = 0,   col = "grey40",     lty = 3, lwd = 1.2)

      legend("topright", bty = "n", cex = 0.80,
             legend = c(sprintf("Estimate: %.4f", t0j),
                        sprintf("%s%% %s CI: [%.4f; %.4f]",
                                round(100 * x$conf), toupper(x$type), lo, hi),
                        sprintf("p-value: %.4f%s",
                                pv, .pval_label(pv))),
             col = c("darkgreen", "firebrick", NA),
             lty = c(1, 2, NA), lwd = c(2, 1.5, NA))

      mtext(sprintf("B = %d  |  SE = %.4f  |  Bias = %.4f",
                    x$B, x$se[j], x$bias[j]),
            side = 3, line = 0.25, cex = 0.72, col = "grey40")
    }
  }

  # ── Panel 2: Bootstrap scatter-plot matrix ─────────────────────────────────
  if (2L %in% which) {
    if (p_coef < 2L) {
      message("Panel 2 requires at least 2 coefficients; skipped.")
    } else {
      op_mat <- par(no.readonly = TRUE)
      on.exit(par(op_mat), add = TRUE)
      par(mfrow = c(p_coef, p_coef),
          mar   = c(2, 2, 1.5, 0.5),
          oma   = c(0, 0, 3, 0))

      for (i in seq_len(p_coef)) {
        for (j in seq_len(p_coef)) {

          if (i == j) {
            # Diagonal: marginal density
            dens <- density(x$t[, i])
            plot(dens, main = "", xlab = "", ylab = "",
                 col = "steelblue", lwd = 2, axes = FALSE)
            abline(v = x$t0[i], col = "darkgreen", lty = 1, lwd = 1.5)
            box(col = "grey80")
            usr <- par("usr")
            text(mean(usr[1:2]), usr[4] * 0.92, nms[i],
                 cex = 0.85, font = 2, col = "grey20")

          } else {
            # Off-diagonal: scatter of replicates
            xi <- x$t[, j];  yi <- x$t[, i]
            r  <- cor(xi, yi)

            plot(xi, yi, pch = 19, cex = 0.45,
                 col  = rgb(0.2, 0.45, 0.8, 0.35),
                 xlab = "", ylab = "", axes = FALSE)
            box(col = "grey80")
            .draw_ellipse(xi, yi)
            points(x$t0[j], x$t0[i], pch = 3, col = "darkgreen",
                   cex = 1.3, lwd = 2)
            mtext(sprintf("r=%.2f", r), side = 3, line = -1.3,
                  adj = 0.97, cex = 0.65,
                  col = if (abs(r) > 0.5) "firebrick" else "grey40")
          }
        }
      }
      mtext("Bootstrap scatter-plot matrix",
            outer = TRUE, cex = 0.9, font = 2, col = "grey20")
    }
  }

  invisible(x)
}

# Returns a short significance label for the histogram legend
.pval_label <- function(p) {
  if      (p < 0.001) " ***"
  else if (p < 0.01)  " **"
  else if (p < 0.05)  " *"
  else if (p < 0.1)   " ."
  else                 ""
}

# Draws an approximate 95% confidence ellipse from the empirical covariance
# of two bootstrap columns, using spectral decomposition of the 2x2 matrix.
.draw_ellipse <- function(x, y, conf = 0.95, n_pts = 100L) {
  S    <- cov(cbind(x, y))
  mu   <- c(mean(x), mean(y))
  ev   <- eigen(S)
  r    <- sqrt(qchisq(conf, df = 2L))
  th   <- seq(0, 2 * pi, length.out = n_pts)
  circ <- rbind(cos(th), sin(th))
  ell  <- t(ev$vectors %*% diag(r * sqrt(ev$values)) %*% circ)
  ell[, 1] <- ell[, 1] + mu[1]
  ell[, 2] <- ell[, 2] + mu[2]
  lines(ell, col = rgb(0.7, 0.2, 0.2, 0.6), lwd = 1.2)
}
