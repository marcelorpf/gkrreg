# gkrreg 0.4.0

## New features

* `sigma_method = "s4"`: new width hyper-parameter estimator based on the
  AICc-selected bandwidth from `sm::h.select()` (squared to obtain a
  variance-scale estimate of `gamma^2`).  Recommended for large samples
  (`n >= 200`).  Requires the **sm** package (now in `Imports`).

* `sigma_method = "auto"`: automatic selection among S1, S2, S3 and S4 by
  out-of-bag bootstrap MSE.  The number of replicates and seed are
  configurable via the new `auto_args` argument (e.g.,
  `auto_args = list(B = 99, seed = 1)`).  A `message()` is emitted
  informing the user of the selected method and comparative OOB MSEs.

* New `auto_args` argument in `gkrr()` for controlling the `"auto"`
  selection bootstrap (default `B = 99`).

## Removed features

* The `weighted` argument has been removed from `gkrr_boot()`.  The
  weighted bootstrap was found to produce wider confidence intervals than
  the standard pairs bootstrap in all tested scenarios, because the
  robustness of GKRReg already resides in the kernel weights — the
  bootstrap itself does not need to replicate this.  The standard pairs
  bootstrap is the recommended and only available option.

## Other changes

* **sm** added to `Imports` (previously not a dependency).

---

# gkrreg 0.3.0

## New features

* Sandwich variance estimator (HC0) is now computed automatically at fit
  time and used by `summary()` for inference (standard errors, confidence
  intervals and Wald z-tests) when no bootstrap object is available.
* `vcov()` method added, returning the sandwich covariance matrix.
* `summary()` gains a `se_tol` argument controlling the threshold for
  divergence warnings between sandwich and bootstrap standard errors.
* `summary()` emits a proactive note suggesting bootstrap inference when
  small sample size (`n < 50`) or heavy contamination is detected.

## Bug fixes

* `par()` settings are now properly restored in all plot methods and
  vignette chunks using `oldpar <- par(no.readonly = TRUE)` /
  `on.exit(par(oldpar))`.

---

# gkrreg 0.2.0

First public release.

## New features

* `gkrr()` fits a Gaussian Kernel Robust Regression model via IRWLS.
  Three estimators for the kernel width hyper-parameter: `"s1"` (Caputo),
  `"s2"` (pairwise median) and `"s3"` (residual variance).
* `gkrr_boot()` runs a pairs bootstrap to produce standard errors,
  confidence intervals (percentile, normal, BCa) and p-values.
* `boot` argument in `gkrr()`: set `boot = TRUE` to compute bootstrap
  inference at fit time.
* `summary.gkrr()` prints a coefficient table modelled after `summary.lm()`.
* `plot.gkrr()` provides six diagnostic panels where point size is
  inversely proportional to the kernel weight.
* `plot.gkrr_boot()` provides histogram and scatter-plot matrix panels.
* Six real datasets bundled: `belgium_calls`, `cloud_point`, `kootenay`,
  `delivery`, `mammals`, `stars_cyg`.
