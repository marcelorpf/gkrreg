## Resubmission (0.4.0)

This is a resubmission. The following changes have been made since
version 0.3.0:

### New features

* `sigma_method = "s4"`: AICc-based bandwidth estimator via `sm::h.select()`.
* `sigma_method = "auto"`: automatic selection among S1--S4 by out-of-bag
  bootstrap MSE, controlled by the new `auto_args` argument.
* The **sm** package has been added to `Imports`.

### Removed

* The `weighted` argument has been removed from `gkrr_boot()` (the weighted
  bootstrap was found to be inferior to the standard pairs bootstrap for
  GKRReg in all tested scenarios).

## R CMD check results

0 errors | 0 warnings | 0 notes
