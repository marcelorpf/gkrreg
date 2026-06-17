# ==============================================================================
# gkrreg: Gaussian Kernel Robust Regression (GKRReg / GKRR)
# De Carvalho, Lima Neto & Ferreira (2017), Neurocomputing 234:58-74
# ==============================================================================
#
# Why only the Gaussian kernel?
# ------------------------------
# GKRR solves the normal equations
#
#   dS/dbeta_r = 2 * sum_i [ G(y_i, mu_i) * x_ir * (y_i - mu_i) / gamma^2 ] = 0
#
# where S = 2 * sum_i [1 - G(y_i, mu_i)]  and  G(a,b) = exp(-(a-b)^2 / gamma^2).
#
# The key property of the Gaussian kernel is:
#
#   dG/d(mu_i) = (2/gamma^2)(y_i - mu_i) * G(y_i, mu_i)          ... (*)
#
# which allows G to cancel in the normal equations, reducing the problem to an
# iterative WLS with weights k_ii = G(y_i, mu_i^{(t-1)}).  The Gaussian kernel
# is the ONLY one among the three kernels in Table 1 of the paper that is
# differentiable everywhere (including at residual = 0), which formally
# underpins Propositions 4.1 and 4.2 (guaranteed convergence).  For this
# reason we adopt exclusively the Gaussian kernel and rename the method GKRReg.
# ==============================================================================


# ── Gaussian kernel (internal) ────────────────────────────────────────────────

.gkern <- function(a, b, gamma2) exp(-(a - b)^2 / gamma2)


# ── Width hyper-parameter estimators ─────────────────────────────────────────

#' Width hyper-parameter estimators for GKRReg
#'
#' Four estimators for the Gaussian kernel width parameter \eqn{\gamma^2}.
#' All functions receive \code{y} (observed response) and \code{yhat}
#' (OLS fitted values) and return a positive scalar \eqn{\hat\gamma^2}.
#'
#' \describe{
#'   \item{\strong{S1 — Caputo}}{Mean of the 0.1 and 0.9 quantiles of the
#'     squared residuals on a sub-sample of size
#'     \eqn{n^* = \lfloor\alpha n\rfloor}.
#'     Recommended for clean data.}
#'   \item{\strong{S2 — Pairwise median}}{Median of \eqn{(y_i -
#'     \hat\mu_j^{\text{OLS}})^2},\ \eqn{i \neq j}.
#'     Recommended for Y-space outliers up to 10\% and X-space outliers
#'     up to 15\%.}
#'   \item{\strong{S3 — Residual variance}}{
#'     \eqn{\hat\gamma^2 = \sum(y_i - \hat\mu_i)^2 / (n - p - 1)}.
#'     Recommended for Y-space outliers \eqn{\geq 15\%} and leverage points.}
#'   \item{\strong{S4 — AICc bandwidth}}{Squares the bandwidth \eqn{h}
#'     selected by minimising the corrected Akaike information criterion in
#'     a nonparametric regression of \eqn{y} on \eqn{\hat\mu^{\text{OLS}}},
#'     via \code{sm::h.select(..., method = "aicc")}.
#'     Recommended for large samples (\eqn{n \geq 200}).
#'     Requires the \pkg{sm} package.}
#' }
#'
#' @param y     Numeric vector of observed responses.
#' @param yhat  Numeric vector of OLS fitted values.
#' @param p     Number of predictors excluding the intercept (used in S3 only).
#' @param alpha Fraction of the sample used in S1 (default \code{0.5}).
#'
#' @return A positive scalar \eqn{\hat\gamma^2}.
#' @name sigma2_estimators
#' @references De Carvalho et al. (2017) \doi{10.1016/j.neucom.2016.12.035}
NULL

#' @rdname sigma2_estimators
#' @export
sigma2_s1 <- function(y, yhat, alpha = 0.5) {
  n   <- length(y)
  m   <- max(2L, floor(n * alpha))
  idx <- sample.int(n, m, replace = FALSE)
  r2  <- (y[idx] - yhat[idx])^2
  r2  <- r2[r2 > 0]
  if (length(r2) < 2L) return(var(y - yhat))
  mean(quantile(r2, probs = c(0.1, 0.9)))
}

#' @rdname sigma2_estimators
#' @export
sigma2_s2 <- function(y, yhat) {
  # r_ij = (y_i - mu_j^OLS)^2, i != j — vectorised via outer() (O(n^2) memory,
  # no explicit R loop; acceptable for n <= ~2000)
  mat       <- outer(y, yhat, function(yi, muj) (yi - muj)^2)
  diag(mat) <- NA
  median(mat, na.rm = TRUE)
}

#' @rdname sigma2_estimators
#' @export
sigma2_s3 <- function(y, yhat, p) {
  n <- length(y)
  sum((y - yhat)^2) / (n - p - 1L)
}

#' @rdname sigma2_estimators
#' @param y     Numeric vector of observed responses (S4 only).
#' @param yhat  Numeric vector of OLS fitted values (S4 only).
#' @export
sigma2_s4 <- function(y, yhat) {
  # S4: bandwidth selected by AICc in the nonparametric regression of y on
  # yhat, via sm::h.select().  The bandwidth h is squared to obtain a
  # variance-scale estimate of gamma^2.
  h <- sm::h.select(y, yhat, method = "aicc", display = "none")
  h^2
}


# ── Automatic sigma_method selection (not exported) ──────────────────────────
#
# Selects among S1, S2, S3 and S4 by out-of-bag (OOB) bootstrap MSE.
# For each candidate method, B bootstrap replicates are drawn.  In each
# replicate, the model is fitted on the in-bag sample and the OOB MSE
# (mean squared prediction error on the left-out observations) is computed.
# The method with the lowest average OOB MSE is returned.
#
# Arguments:
#   X, y        -- design matrix and response (already parsed)
#   p           -- number of predictors (excl. intercept)
#   alpha       -- passed to sigma2_s1
#   tol, maxit  -- passed to .gkrr_fit
#   auto_args   -- list with optional element B (default 99)

.gkrr_auto_select <- function(X, y, p, alpha, tol, maxit, auto_args) {

  B_auto    <- if (!is.null(auto_args$B)) as.integer(auto_args$B) else 99L
  seed_auto <- auto_args$seed
  if (!is.null(seed_auto)) set.seed(seed_auto)

  n         <- nrow(X)
  candidates <- c("s1", "s2", "s3", "s4")
  oob_mse   <- setNames(rep(0, length(candidates)), candidates)
  oob_count <- setNames(rep(0L, length(candidates)), candidates)

  for (b in seq_len(B_auto)) {
    idx_in  <- sample.int(n, n, replace = TRUE)
    idx_out <- setdiff(seq_len(n), unique(idx_in))
    if (length(idx_out) == 0L) next

    X_in <- X[idx_in, , drop = FALSE]
    y_in <- y[idx_in]
    X_out <- X[idx_out, , drop = FALSE]
    y_out <- y[idx_out]

    for (sm in candidates) {
      fit_b <- tryCatch(
        .gkrr_fit(X_in, y_in, sigma_method = sm,
                  alpha = alpha, tol = tol, maxit = maxit),
        error = function(e) NULL
      )
      if (is.null(fit_b) || !fit_b$converged) next

      pred <- as.vector(X_out %*% fit_b$coefficients)
      oob_mse[sm]   <- oob_mse[sm]   + mean((y_out - pred)^2)
      oob_count[sm] <- oob_count[sm] + 1L
    }
  }

  # Average OOB MSE; if a method never converged, set to Inf
  avg_mse <- ifelse(oob_count > 0L, oob_mse / oob_count, Inf)
  best    <- names(which.min(avg_mse))

  message(sprintf(
    "sigma_method = \"auto\": \"%s\" selected (avg. OOB MSE: %s)",
    best,
    paste(sprintf("%s=%.4f", names(avg_mse), avg_mse), collapse = ", ")
  ))

  best
}


# ── Algorithm core (not exported) ─────────────────────────────────────────────
#
# Receives pre-built matrices/vectors so it can be reused by the bootstrap
# without re-parsing the formula on every replicate.

.gkrr_fit <- function(X, y, sigma_method, alpha, tol, maxit) {

  p <- ncol(X) - 1L

  # Step 0: OLS initialisation
  ols  <- lm.fit(X, y)
  beta <- ols$coefficients
  yhat <- as.vector(X %*% beta)

  # Estimate gamma^2
  if (is.numeric(sigma_method)) {
    gamma2 <- sigma_method
  } else {
    gamma2 <- switch(sigma_method,
      s1 = sigma2_s1(y, yhat, alpha = alpha),
      s2 = sigma2_s2(y, yhat),
      s3 = sigma2_s3(y, yhat, p = p),
      s4 = sigma2_s4(y, yhat)
    )
    if (!is.finite(gamma2) || gamma2 <= 0)
      gamma2 <- sigma2_s3(y, yhat, p = p)
  }

  # Iterative WLS
  weights  <- .gkern(y, yhat, gamma2)
  S        <- sum(2 - 2 * weights)
  crit_seq <- S
  it       <- 1L

  repeat {
    it       <- it + 1L
    wls      <- lm.wfit(X, y, w = weights)
    beta     <- wls$coefficients
    yhat     <- as.vector(X %*% beta)
    weights  <- .gkern(y, yhat, gamma2)
    S_new    <- sum(2 - 2 * weights)
    crit_seq <- c(crit_seq, S_new)
    if (abs(S_new - crit_seq[it - 1L]) <= tol || it >= maxit) break
  }

  list(
    coefficients  = beta,
    fitted.values = yhat,
    residuals     = unname(y - yhat),
    weights       = as.vector(weights),
    gamma2        = gamma2,
    criterion     = crit_seq,
    iterations    = it,
    converged     = abs(crit_seq[it] - crit_seq[it - 1L]) <= tol
  )
}


# ── Sandwich variance estimator (not exported) ────────────────────────────────
#
# At convergence, GKRR satisfies the estimating equations
#
#   Psi(beta) = sum_i k_ii(beta) * x_i * (y_i - x_i' beta) = 0
#
# Treating this as a generalised M-estimator, the asymptotic sandwich
# covariance is V = A^{-1} B A^{-1}, where:
#
#   A = (1/n) * sum_i k_ii * x_i x_i' * (1 - 2*e_i^2 / gamma^2)
#   B = (1/n) * sum_i k_ii^2 * e_i^2 * x_i x_i'
#
# The derivative dk_ii/d(beta) = (2/gamma^2) * e_i * k_ii * x_i (from the
# Gaussian kernel property (*) in the file header) is what gives the
# correction term (1 - 2*e_i^2/gamma^2) in A.
#
# Standard errors: se_j = sqrt(V_jj / n)
# Wald p-values:   p_j  = 2 * Phi(-|beta_j / se_j|)
# CIs:             beta_j +/- z_{alpha/2} * se_j

.gkrr_sandwich <- function(X, y, beta, weights, gamma2, conf = 0.95) {

  n      <- nrow(X)
  e      <- as.vector(y - X %*% beta)   # residuals
  k      <- weights                      # kernel weights k_ii
  z_half <- qnorm(1 - (1 - conf) / 2)

  # ── Bread: A = (1/n) X' diag(k * (1 - 2e^2/gamma^2)) X ────────────────────
  d_A <- k * (1 - 2 * e^2 / gamma2)    # length-n diagonal multiplier
  A   <- crossprod(X * d_A, X) / n     # p x p

  # ── Meat: B = (1/n) X' diag(k^2 * e^2) X ──────────────────────────────────
  d_B <- k^2 * e^2
  B   <- crossprod(X * d_B, X) / n

  # ── Sandwich: V = A^{-1} B A^{-1} ──────────────────────────────────────────
  # Use solve() with a fallback to pseudoinverse if A is singular
  A_inv <- tryCatch(solve(A), error = function(e) MASS::ginv(A))
  V     <- A_inv %*% B %*% A_inv

  se   <- sqrt(diag(V) / n)
  tval <- beta / se
  pval <- pmax(2e-16, 2 * pnorm(-abs(tval)))

  ci_lo <- beta - z_half * se
  ci_hi <- beta + z_half * se

  list(
    vcov  = V / n,     # full covariance matrix (scale: 1/n already in A, B)
    se    = se,
    tval  = tval,
    pval  = pval,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    conf  = conf
  )
}



#' Gaussian Kernel Robust Regression (GKRReg)
#'
#' Fits a robust linear regression model using the Gaussian kernel to
#' down-weight poorly fitted observations (outliers and leverage points).
#' Weights \eqn{k_{ii} = G(y_i, \hat\mu_i) \in (0, 1]} are updated at each
#' iteration of an IRWLS until the objective function
#' \eqn{S(\beta) = 2\sum_{i=1}^n [1 - G(y_i, \hat\mu_i)]} converges.
#'
#' Convergence is guaranteed by Propositions 4.1 and 4.2 of
#' De Carvalho et al. (2017).
#'
#' When \code{boot = TRUE} a \code{\link{gkrr_boot}} object is computed and
#' stored inside the fitted model, making it available automatically to
#' \code{\link{summary.gkrr}} for inference (confidence intervals and
#' bootstrap p-values).
#'
#' @param formula      Model formula, e.g. \code{y ~ x1 + x2}.
#' @param data         Optional \code{data.frame} containing the model variables.
#' @param sigma_method Estimator for \eqn{\gamma^2}: \code{"s1"} (Caputo
#'   sub-sample quantiles), \code{"s2"} (pairwise median),
#'   \code{"s3"} (residual variance), \code{"s4"} (AICc bandwidth via
#'   \code{sm::h.select()}, recommended for \eqn{n \geq 200}), or
#'   \code{"auto"} (selects among S1--S4 by out-of-bag bootstrap MSE;
#'   incurs additional computational cost — see \code{auto_args}).
#'   A single positive numeric value fixes \eqn{\gamma^2} directly.
#' @param auto_args Named list of arguments controlling the \code{"auto"}
#'   selection bootstrap.  Recognised elements: \code{B} (number of
#'   replicates, default \code{99}) and \code{seed} (integer seed for
#'   reproducibility).  Ignored when \code{sigma_method != "auto"}.
#' @param alpha        Fraction of the sample used in estimator S1
#'   (default \code{0.5}).
#' @param tol          Convergence tolerance (default \code{1e-10}).
#' @param maxit        Maximum number of iterations (default \code{100}).
#' @param boot         Logical or a \code{"gkrr_boot"} object.
#'   \itemize{
#'     \item \code{FALSE} (default): sandwich standard errors, confidence
#'       intervals and Wald p-values are computed analytically and stored in
#'       \code{fit$sandwich}.  \code{summary()} uses them automatically.
#'     \item \code{TRUE}: runs \code{gkrr_boot()} with default arguments
#'       (\eqn{B = 999}, \code{type = "bca"}) and stores the result inside
#'       the fitted object.  Bootstrap inference takes precedence over the
#'       sandwich in \code{summary()}.
#'     \item A pre-computed \code{"gkrr_boot"} object: stored directly,
#'       avoiding redundant computation.
#'   }
#' @param boot_args    Named list of extra arguments passed to
#'   \code{\link{gkrr_boot}} when \code{boot = TRUE} (e.g.
#'   \code{list(B = 499, type = "percentile", seed = 1)}).
#' @param conf         Confidence level for sandwich intervals (default
#'   \code{0.95}; ignored when \code{boot != FALSE}).
#'
#' @return An object of class \code{"gkrr"} containing:
#' \describe{
#'   \item{\code{coefficients}}{Named vector of estimated coefficients.}
#'   \item{\code{fitted.values}}{Fitted values \eqn{\hat\mu}.}
#'   \item{\code{residuals}}{Residuals \eqn{y - \hat\mu}.}
#'   \item{\code{weights}}{Final kernel weights \eqn{k_{ii}}.}
#'   \item{\code{gamma2}}{Value of \eqn{\hat\gamma^2} used.}
#'   \item{\code{criterion}}{Sequence of \eqn{S} values across iterations.}
#'   \item{\code{iterations}}{Number of iterations performed.}
#'   \item{\code{converged}}{Logical; was convergence achieved?}
#'   \item{\code{sigma_method}}{Label of the \eqn{\gamma^2} estimator.}
#'   \item{\code{sandwich}}{A list with sandwich inference components:
#'     \code{vcov}, \code{se}, \code{tval}, \code{pval}, \code{ci_lo},
#'     \code{ci_hi}, \code{conf}.  Always computed.}
#'   \item{\code{boot}}{A \code{"gkrr_boot"} object if \code{boot != FALSE},
#'     otherwise \code{NULL}.}
#'   \item{\code{call}, \code{terms}, \code{model}}{Model metadata.}
#' }
#'
#' @references
#' De Carvalho, F.A.T., Lima Neto, E.A., Ferreira, M.R.P. (2017).
#' A robust regression method based on exponential-type kernel functions.
#' \emph{Neurocomputing}, 234, 58--74. \doi{10.1016/j.neucom.2016.12.035}
#'
#' @seealso \code{\link{gkrr_boot}} for bootstrap inference,
#'   \code{\link{summary.gkrr}} for the inference table.
#'
#' @examples
#' set.seed(1)
#' n  <- 80
#' x  <- runif(n, 0, 10)
#' y  <- 2 + 3 * x + rnorm(n)
#' y[1:12] <- y[1:12] + 25      # 15% Y-space outliers
#'
#' # Basic fit with sandwich inference (default)
#' fit <- gkrr(y ~ x, sigma_method = "s3")
#' summary(fit)
#'
#' # Using S4 estimator (AICc bandwidth)
#' \donttest{
#' fit_s4 <- gkrr(y ~ x, sigma_method = "s4")
#' summary(fit_s4)
#' }
#'
#' # Automatic estimator selection
#' \donttest{
#' fit_auto <- gkrr(y ~ x, sigma_method = "auto",
#'                  auto_args = list(B = 49, seed = 1))
#' summary(fit_auto)
#' }
#'
#' # Fit with bootstrap inference (BCa, B = 999 by default)
#' \donttest{
#' fit_b <- gkrr(y ~ x, sigma_method = "s3", boot = TRUE,
#'               boot_args = list(B = 499, seed = 1))
#' summary(fit_b)
#' plot(fit_b)
#' }
#'
#' @export
gkrr <- function(formula,
                 data         = NULL,
                 sigma_method = c("s1", "s2", "s3", "s4", "auto"),
                 alpha        = 0.5,
                 tol          = 1e-10,
                 maxit        = 100L,
                 boot         = FALSE,
                 boot_args    = list(),
                 auto_args    = list(),
                 conf         = 0.95) {

  cl <- match.call()

  # ── Parse model ─────────────────────────────────────────────────────────────
  mf                    <- match.call(expand.dots = FALSE)
  m                     <- match(c("formula", "data"), names(mf), 0L)
  mf                    <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action          <- quote(stats::na.fail)
  mf[[1L]]              <- quote(stats::model.frame)
  mf                    <- eval(mf, parent.frame())
  mt                    <- attr(mf, "terms")
  y                     <- model.response(mf, "numeric")
  X                     <- model.matrix(mt, mf)

  n <- length(y);  p <- ncol(X) - 1L

  # ── Validate inputs ──────────────────────────────────────────────────────────
  if (n < ncol(X) + 1L) stop("n < number of parameters.")
  if (anyNA(y) || anyNA(X)) stop("NAs are not allowed.")
  if (tol   <= 0) stop("'tol' must be positive.")
  if (maxit  < 1) stop("'maxit' must be >= 1.")
  if (conf  <= 0 || conf >= 1) stop("'conf' must be in (0, 1).")

  # ── Normalise sigma_method ───────────────────────────────────────────────────
  if (is.numeric(sigma_method)) {
    if (sigma_method <= 0) stop("Numeric 'sigma_method' must be > 0.")
    sm_label <- paste0("fixed(", round(sigma_method, 4), ")")
    sm_arg   <- sigma_method
  } else {
    sm_arg <- match.arg(sigma_method, c("s1", "s2", "s3", "s4", "auto"))

    if (sm_arg == "auto") {
      # Resolve "auto" before fitting: select best estimator by OOB MSE
      message("sigma_method = \"auto\": running bootstrap selection ",
              "(B = ", if (!is.null(auto_args$B)) auto_args$B else 99L,
              " replicates). Use auto_args = list(B = ...) to change.")
      sm_arg <- .gkrr_auto_select(X, y, p = p,
                                   alpha    = alpha,
                                   tol      = tol,
                                   maxit    = maxit,
                                   auto_args = auto_args)
      sm_label <- paste0("auto(", sm_arg, ")")
    } else {
      sm_label <- sm_arg
    }
  }

  # ── Fit ──────────────────────────────────────────────────────────────────────
  fit <- .gkrr_fit(X, y, sigma_method = sm_arg,
                   alpha = alpha, tol = tol, maxit = maxit)

  if (!fit$converged)
    warning(sprintf("Algorithm did not converge in %d iterations.", maxit))

  names(fit$coefficients) <- colnames(X)

  # ── Sandwich variance estimator (always computed) ────────────────────────────
  sw <- .gkrr_sandwich(X, y,
                       beta    = fit$coefficients,
                       weights = fit$weights,
                       gamma2  = fit$gamma2,
                       conf    = conf)
  names(sw$se)    <- colnames(X)
  names(sw$tval)  <- colnames(X)
  names(sw$pval)  <- colnames(X)
  names(sw$ci_lo) <- colnames(X)
  names(sw$ci_hi) <- colnames(X)
  rownames(sw$vcov) <- colnames(sw$vcov) <- colnames(X)

  # ── Bootstrap ────────────────────────────────────────────────────────────────
  boot_obj <- NULL

  if (isTRUE(boot)) {
    # The temporary object passed to gkrr_boot must store sm_arg (the resolved
    # method, e.g. "s4") in sigma_method -- not sm_label (e.g. "auto(s4)") --
    # because gkrr_boot recovers the gamma^2 estimator from that field to
    # re-estimate on each bootstrap replicate.
    tmp <- structure(
      c(fit, list(sigma_method = sm_arg, call = cl, terms = mt, model = mf,
                  boot = NULL)),
      class = "gkrr"
    )
    boot_obj <- do.call(gkrr_boot, c(list(object = tmp), boot_args))

  } else if (inherits(boot, "gkrr_boot")) {
    boot_obj <- boot

  } else if (!isFALSE(boot)) {
    stop("'boot' must be FALSE, TRUE, or a \"gkrr_boot\" object.")
  }

  structure(
    c(fit,
      list(sigma_method = sm_label,
           sandwich     = sw,
           boot         = boot_obj,
           call         = cl,
           terms        = mt,
           model        = mf)),
    class = "gkrr"
  )
}


# ── S3 methods for "gkrr" ─────────────────────────────────────────────────────

#' @export
coef.gkrr <- function(object, ...) object$coefficients

#' @export
fitted.gkrr <- function(object, ...) object$fitted.values

#' @export
residuals.gkrr <- function(object, ...) object$residuals

#' @export
vcov.gkrr <- function(object, ...) object$sandwich$vcov

#' @export
print.gkrr <- function(x, digits = 4L, ...) {
  cat("\nGKRReg -- Gaussian Kernel Robust Regression\n")
  cat("--------------------------------------------\n")
  cat(sprintf("gamma^2 : %.4g  (method: %s)\n", x$gamma2, x$sigma_method))
  cat(sprintf("Iterations: %d  [%s]\n\n", x$iterations,
              if (x$converged) "converged" else "DID NOT converge"))
  cat("Coefficients:\n")
  print(round(x$coefficients, digits))
  if (!is.null(x$boot)) {
    cat("\n(Bootstrap inference available -- use summary() for the full table)\n")
  } else {
    cat("\n(Sandwich inference available -- use summary() for the full table)\n")
  }
  invisible(x)
}

#' Summary of a GKRReg fit
#'
#' Prints a coefficient table modelled after \code{summary.lm}.  By default
#' inference is based on the sandwich variance estimator (always computed at
#' fit time).  When a \code{\link{gkrr_boot}} object is available it takes
#' precedence, replacing the sandwich standard errors and p-values with their
#' bootstrap counterparts.
#'
#' \strong{Sandwich inference} (default, \code{boot = FALSE} in
#' \code{\link{gkrr}}):
#' Standard errors are derived from the asymptotic sandwich covariance matrix
#' \eqn{\hat V = A^{-1} B A^{-1} / n}, where
#' \eqn{A = n^{-1} X^\top K \,\mathrm{diag}(1 - 2e_i^2/\hat\gamma^2) X}
#' and \eqn{B = n^{-1} X^\top K^2 \,\mathrm{diag}(e_i^2) X}.
#' P-values use the two-sided Wald z-test
#' \eqn{p_j = 2\,\Phi(-|\hat\beta_j / \mathrm{se}_j|)}.
#' Confidence intervals are \eqn{\hat\beta_j \pm z_{\alpha/2}\,\mathrm{se}_j}.
#'
#' \strong{Note:} the sandwich estimator is asymptotically valid but does not
#' account for the variability introduced by estimating \eqn{\hat\gamma^2}.
#' For small samples or when \eqn{\gamma^2} estimation has high variance,
#' bootstrap inference (\code{boot = TRUE}) is more reliable.
#'
#' \strong{Bootstrap inference} (when \code{boot != FALSE}):
#' Bootstrap p-values use the centred-t approach:
#' \deqn{p_j = \frac{2}{B}\#\{|\hat\beta^*_b - \hat\beta_j| \geq |\hat\beta_j|\}}
#' which counts how many bootstrap replicates are at least as extreme as the
#' observed estimate relative to zero.
#'
#' @param object       A \code{"gkrr"} object.
#' @param boot         Optional \code{"gkrr_boot"} object.  If \code{NULL}
#'   (default), the method checks whether \code{object$boot} already contains
#'   one (set via \code{boot = TRUE} in \code{\link{gkrr}}).  If neither is
#'   available, the sandwich estimator is used.
#' @param digits       Number of significant digits (default \code{4}).
#' @param signif_stars Logical; print significance stars (default \code{TRUE}).
#' @param se_tol       Relative divergence threshold between sandwich and
#'   bootstrap standard errors.  When both are available and
#'   \eqn{|\text{SE}_{\text{sw}} - \text{SE}_{\text{boot}}| /
#'   \text{SE}_{\text{boot}} > \text{se\_tol}} for any coefficient, a warning
#'   is issued listing which coefficients diverge and in which direction.
#'   Set to \code{Inf} to suppress the check.  Default \code{0.25} (25\%).
#' @param ...          Currently ignored.
#'
#' @return Invisibly returns a list with components \code{coefficients}
#'   (the printed table as a matrix) and \code{r_squared}.
#' @seealso \code{\link{gkrr}}, \code{\link{gkrr_boot}}
#' @export
summary.gkrr <- function(object, boot = NULL, digits = 4L,
                         signif_stars = TRUE, se_tol = 0.25, ...) {

  # Resolve inference source: bootstrap argument > stored > sandwich
  if (is.null(boot)) boot <- object$boot
  has_boot <- inherits(boot, "gkrr_boot")
  sw       <- object$sandwich     # always present

  res  <- object$residuals
  wt   <- object$weights
  yobs <- object$fitted.values + res
  n    <- length(res)

  # ── R-squared ───────────────────────────────────────────────────────────────
  ss_res <- sum(res^2)
  ss_tot <- sum((yobs - mean(yobs))^2)
  r2     <- 1 - ss_res / ss_tot
  yw     <- sum(wt * yobs) / sum(wt)
  r2_w   <- 1 - sum(wt * res^2) / sum(wt * (yobs - yw)^2)

  # ── Coefficient table ────────────────────────────────────────────────────────
  coefs  <- object$coefficients
  p_coef <- length(coefs)
  conf   <- sw$conf

  if (has_boot) {
    # Bootstrap inference
    se_use  <- boot$se
    ci_lo   <- boot$ci["lower", ]
    ci_hi   <- boot$ci["upper", ]
    pval    <- sapply(seq_len(p_coef), function(j) {
      centred <- abs(boot$t[, j] - coefs[j])
      max(1 / boot$B, mean(centred >= abs(coefs[j])))
    })
    conf_pct  <- sprintf("%.0f%%", 100 * boot$conf)
    ci_label  <- sprintf("CI %s [%s]", conf_pct, toupper(boot$type))
    inf_label <- sprintf("Bootstrap p-values: centred-t, B = %d, %s CI",
                         boot$B, toupper(boot$type))
  } else {
    # Sandwich inference (default)
    se_use   <- sw$se
    ci_lo    <- sw$ci_lo
    ci_hi    <- sw$ci_hi
    pval     <- sw$pval
    conf_pct <- sprintf("%.0f%%", 100 * conf)
    ci_label <- sprintf("CI %s [Sandwich]", conf_pct)
    inf_label <- "Sandwich SE (asymptotic Wald z-test)"
  }

  tab <- cbind(
    Estimate = coefs,
    "Std. Error" = se_use,
    setNames(data.frame(ci_lo, ci_hi),
             c(paste("CI", conf_pct, "lower"),
               paste("CI", conf_pct, "upper"))),
    "p-value" = pval
  )

  # ── SE divergence check (only when both sandwich and bootstrap are present) ──
  if (has_boot && is.finite(se_tol)) {
    rel_diff <- abs(sw$se - boot$se) / boot$se
    diverged <- which(rel_diff > se_tol)
    if (length(diverged) > 0L) {
      nms_div <- names(coefs)[diverged]
      dir_str <- ifelse(sw$se[diverged] < boot$se[diverged],
                        "sandwich < bootstrap (sandwich may underestimate SE)",
                        "sandwich > bootstrap (sandwich may overestimate SE)")
      msg <- paste0(
        "Sandwich and bootstrap SEs diverge by more than ",
        round(se_tol * 100), "% for: ",
        paste(sprintf("%s [%.0f%%, %s]",
                      nms_div,
                      rel_diff[diverged] * 100,
                      dir_str),
              collapse = "; "),
        ".\n",
        "  This may indicate small n, heavy outlier contamination, or high\n",
        "  variability in the gamma^2 estimator. Bootstrap inference is\n",
        "  recommended in this case."
      )
      warning(msg, call. = FALSE)
    }
  }

  # ── Print ────────────────────────────────────────────────────────────────────
  cat("\nCall:\n"); print(object$call)

  cat("\nResiduals:\n")
  res_q <- quantile(res, c(0, .25, .5, .75, 1))
  names(res_q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(round(res_q, digits))

  cat("\nCoefficients:\n")

  if (signif_stars) {
    # Round all columns except p-value (rounding kills small p-values)
    tab_num <- tab
    tab_fmt <- apply(round(tab_num[, colnames(tab_num) != "p-value",
                                   drop = FALSE], digits), 2, format)
    # Format p-values: show "< 2e-16" for very small values
    pval_num <- tab_num[, "p-value"]
    pval_fmt <- ifelse(pval_num < 0.001,
                       formatC(pval_num, format = "e", digits = 3),
                       format(round(pval_num, digits)))
    pval_fmt <- ifelse(pval_num < 2e-16, "< 2e-16", pval_fmt)
    tab_fmt  <- cbind(tab_fmt, "p-value" = pval_fmt)
    # as.character() converts the factor levels to their labels
    stars   <- as.character(.pval_stars(pval_num))
    tab_fmt <- cbind(tab_fmt, " " = stars)
    print(tab_fmt, quote = FALSE, right = TRUE)
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    cat(sprintf("(%s)\n", inf_label))
  } else {
    print(round(tab, digits))
  }

  cat(sprintf(
    "\ngamma^2: %.4g  |  method: %s  |  iterations: %d  |  converged: %s\n",
    object$gamma2, object$sigma_method, object$iterations, object$converged))
  cat(sprintf("R-squared: %.4f  |  Weighted R-squared: %.4f\n", r2, r2_w))
  cat(sprintf("Kernel weights -- min: %.4f  mean: %.4f  max: %.4f\n",
              min(wt), mean(wt), max(wt)))

  # ── Proactive note when only sandwich is available ───────────────────────────
  # Suggest bootstrap if n is small or if the GKRR assigned near-zero weight
  # to more than 10% of observations (wt < 0.01 corresponds to a residual of
  # at least 2.15 * sqrt(gamma^2), which is a genuinely large deviation).
  # Using wt < 0.01 rather than wt < 0.1 avoids false positives when gamma^2
  # is small and ordinary residuals happen to produce moderate low weights.
  # Emitted as message(), not warning() -- it is a suggestion, not a problem.
  if (!has_boot) {
    small_n      <- n < 50L
    heavy_contam <- mean(wt < 0.01) > 0.10
    if (small_n || heavy_contam) {
      reasons <- c(
        if (small_n)      sprintf("n = %d (small sample)", n),
        if (heavy_contam) sprintf(
          "%.0f%% of observations received near-zero kernel weight (< 0.01)",
          mean(wt < 0.01) * 100)
      )
      message(
        "Note: sandwich inference may be less reliable here (",
        paste(reasons, collapse = "; "), ").\n",
        "  Consider bootstrap inference via boot = TRUE in gkrr() or\n",
        "  summary(fit, boot = gkrr_boot(fit))."
      )
    }
  }

  invisible(list(
    coefficients = tab,
    r_squared    = c(r2 = r2, r2_weighted = r2_w)
  ))
}

# Significance stars (same thresholds as summary.lm)
.pval_stars <- function(p) {
  cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
      labels = c("***", "**", "*", ".", " "))
}

#' @export
predict.gkrr <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  tt <- delete.response(object$terms)
  mf <- model.frame(tt, newdata, xlev = .getXlevels(tt, object$model))
  X  <- model.matrix(tt, mf)
  as.vector(X %*% object$coefficients)
}

#' Diagnostic plots for a GKRReg fit
#'
#' Produces up to 6 diagnostic panels for a \code{"gkrr"} object.
#' Point size is inversely proportional to the kernel weight \eqn{k_{ii}},
#' so outliers (small weights) appear large and red while well-fitted
#' observations appear small and blue.
#'
#' \describe{
#'   \item{\code{which = 1}}{Residuals vs. fitted values.}
#'   \item{\code{which = 2}}{Observed vs. fitted (\eqn{y} vs. \eqn{\hat\mu}).}
#'   \item{\code{which = 3}}{Kernel weight vs. residual, with the theoretical
#'     curve \eqn{G(e) = \exp(-e^2/\hat\gamma^2)} overlaid.}
#'   \item{\code{which = 4}}{Kernel weight vs. observation index.}
#'   \item{\code{which = 5}}{Normal QQ-plot of residuals, coloured by weight.}
#'   \item{\code{which = 6}}{Objective function \eqn{S(\beta)} by iteration.}
#' }
#'
#' @param x     A \code{"gkrr"} object.
#' @param which Integer vector selecting panels to draw (default \code{1:5}).
#' @param n_id  Number of extreme observations to label in panels 1--5
#'   (default \code{3}).
#' @param ask   Logical; if \code{TRUE} waits for user input between panels
#'   (default \code{TRUE} when \code{length(which) > 1}).
#' @param ...   Additional arguments (ignored; included for S3 compatibility).
#'
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{gkrr}}, \code{\link{plot.gkrr_boot}}
#' @export
plot.gkrr <- function(x, which = 1:5, n_id = 3L,
                      ask = length(which) > 1L, ...) {

  if (ask) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(ask = TRUE)
  }

  res  <- x$residuals
  fit  <- x$fitted.values
  wt   <- x$weights
  yobs <- fit + res
  n    <- length(res)

  # Point size: inversely proportional to weight.
  # weight = 1  ->  cex_min (small point, well-fitted observation)
  # weight -> 0 ->  cex_max (large point, highlighted outlier)
  cex_min <- 0.4
  cex_max <- 2.8
  .cex <- function(w) cex_max - (cex_max - cex_min) * w

  # Colour ramp: blue (well-fitted) -> red (outlier).
  # Size and colour are redundant but complementary channels.
  .col <- function(w) {
    r <- 1 - w
    rgb(r, 0.15 * (1 - r), 1 - r, alpha = 0.75)
  }

  # Indices of the n_id observations with the lowest weights
  n_id  <- min(max(0L, n_id), n)
  low_i <- if (n_id > 0L) order(wt)[seq_len(n_id)] else integer(0)

  # Helper: label the n_id extreme points (no-op when n_id = 0)
  .label <- function(ix, px, py, ...) {
    if (length(ix) == 0L) return(invisible(NULL))
    text(px[ix], py[ix], labels = ix, pos = 3, cex = 0.65,
         col = "firebrick", font = 2)
  }

  # ── Panel 1: Residuals vs. Fitted ──────────────────────────────────────────
  if (1L %in% which) {
    plot(fit, res,
         xlab = expression(hat(mu)),
         ylab = expression(y - hat(mu)),
         main = "Residuals vs. Fitted",
         pch  = 19, cex = .cex(wt), col = .col(wt))
    abline(h = 0, lty = 2, col = "grey50")
    .label(low_i, fit, res)
    legend("topright", bty = "n", cex = 0.75,
           legend = c("high weight (well fitted)", "low weight (outlier)"),
           pch = 19, col = c(.col(1), .col(0)),
           pt.cex = c(.cex(1), .cex(0)))
  }

  # ── Panel 2: Observed vs. Fitted ───────────────────────────────────────────
  if (2L %in% which) {
    lim <- range(c(yobs, fit))
    plot(fit, yobs,
         xlab = expression(hat(mu)),
         ylab = "y  (observed)",
         main = "Observed vs. Fitted",
         xlim = lim, ylim = lim,
         pch  = 19, cex = .cex(wt), col = .col(wt))
    abline(0, 1, lty = 2, col = "grey50")
    .label(low_i, fit, yobs)
    legend("topleft", bty = "n", cex = 0.75,
           legend = c("high weight", "low weight"),
           pch = 19, col = c(.col(1), .col(0)),
           pt.cex = c(.cex(1), .cex(0)))
  }

  # ── Panel 3: Kernel Weight vs. Residual + theoretical curve ────────────────
  if (3L %in% which) {
    e_seq <- seq(min(res) * 1.15, max(res) * 1.15, length.out = 300)
    g_seq <- exp(-e_seq^2 / x$gamma2)

    plot(res, wt,
         xlab = expression(y - hat(mu)),
         ylab = expression(k[ii] == G(y[i], hat(mu)[i])),
         main = "Kernel Weight vs. Residual",
         ylim = c(0, 1),
         pch  = 19, cex = .cex(wt), col = .col(wt))
    lines(e_seq, g_seq, col = "steelblue", lwd = 2)
    .label(low_i, res, wt)
    legend("top", bty = "n", cex = 0.75,
           legend = sprintf("G(e) = exp(-e^2 / %.3g)", x$gamma2),
           col = "steelblue", lwd = 2)
  }

  # ── Panel 4: Kernel Weight vs. Index ───────────────────────────────────────
  if (4L %in% which) {
    idx <- seq_len(n)
    plot(idx, wt,
         xlab = "Observation index",
         ylab = expression(k[ii]),
         main = "Kernel Weight by Observation",
         ylim = c(0, 1),
         pch  = 19, cex = .cex(wt), col = .col(wt))
    abline(h = 1, lty = 2, col = "grey50")
    abline(h = mean(wt), lty = 3, col = "grey40")
    mtext(sprintf("mean weight = %.3f", mean(wt)),
          side = 1, line = -1.2, adj = 0.98, cex = 0.7, col = "grey40")
    segments(idx, 1, idx, wt, col = rgb(0.7, 0.7, 0.7, 0.4), lwd = 0.8)
    .label(low_i, idx, wt)
  }

  # ── Panel 5: Normal QQ-plot of residuals ───────────────────────────────────
  if (5L %in% which) {
    q_teo  <- qnorm(ppoints(n))
    q_emp  <- sort(res)
    wt_ord <- wt[order(res)]

    plot(q_teo, q_emp,
         xlab = "Theoretical N(0,1) quantiles",
         ylab = expression(y - hat(mu)),
         main = "Normal Q-Q Plot of Residuals",
         pch  = 19, cex = .cex(wt_ord), col = .col(wt_ord))
    # Reference line through quartiles (more robust than OLS line)
    q25   <- quantile(res, 0.25);  q75 <- quantile(res, 0.75)
    slope <- (q75 - q25) / (qnorm(0.75) - qnorm(0.25))
    abline(a = q25 - slope * qnorm(0.25), b = slope, lty = 2, col = "grey50")
    low_qq <- order(wt_ord)[seq_len(n_id)]
    .label(low_qq, q_teo, q_emp)
  }

  # ── Panel 6: Convergence of S(beta) ────────────────────────────────────────
  if (6L %in% which) {
    its <- seq_along(x$criterion)
    plot(its, x$criterion,
         type = "b", pch = 19, col = "steelblue", lwd = 1.5,
         xlab = "Iteration",
         ylab = expression(S(beta)),
         main = "Convergence of the Objective Function")
    points(its[length(its)], x$criterion[length(x$criterion)],
           pch = 21, bg = "steelblue", col = "white", cex = 1.8)
    mtext(sprintf("Final S = %.6g  |  %d iterations  |  converged: %s",
                  x$criterion[length(x$criterion)],
                  x$iterations, x$converged),
          side = 3, line = 0.2, cex = 0.72, col = "grey40")
  }

  invisible(x)
}
