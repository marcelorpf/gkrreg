# ==============================================================================
# data.R -- Documentation for datasets bundled with gkrreg
# ==============================================================================
# All six datasets are used in the real-data applications of
# De Carvalho, Lima Neto & Ferreira (2017), Sections 6.1--6.6, except
# 'mammals' and 'stars_cyg' which are added to illustrate the multivariate
# case and the leverage-point scenario, respectively.
# ==============================================================================


#' International Telephone Calls from Belgium (1950--1973)
#'
#' Number of international telephone calls made from Belgium between 1950 and
#' 1973, taken from the Belgian Statistical Survey published by the Ministry of
#' Economy.  The data contain a block of Y-space outliers (observations 15--20)
#' caused by a recording error in the original source: calls were recorded in
#' total minutes instead of number of calls for those years.
#'
#' This dataset is used in De Carvalho et al. (2017), Section 6.2, to
#' illustrate the robustness of GKRR to Y-space outliers.  The recommended
#' width estimator for this scenario is \code{sigma_method = "s3"}.
#'
#' @format A data frame with 24 rows and 2 columns:
#' \describe{
#'   \item{year}{Last two digits of the calendar year (50 = 1950, \ldots,
#'     73 = 1973).}
#'   \item{calls}{Number of international calls originating from Belgium,
#'     in tens of millions.}
#' }
#'
#' @section Outlier structure:
#' Observations 15--20 (years 1964--1969) are Y-space outliers.  Their
#' recorded values are anomalously large relative to the linear trend of the
#' remaining years because calls were measured in total minutes rather than
#' number of calls.
#'
#' @source
#' P.J. Rousseeuw and A.M. Leroy (1987).
#' \emph{Robust Regression and Outlier Detection}.
#' John Wiley & Sons, Table 6, p. 26.
#'
#' Available in \pkg{robustbase} as \code{telef}.
#'
#' @references
#' De Carvalho, F.A.T., Lima Neto, E.A., Ferreira, M.R.P. (2017).
#' A robust regression method based on exponential-type kernel functions.
#' \emph{Neurocomputing}, 234, 58--74.
#' \doi{10.1016/j.neucom.2016.12.035}
#'
#' @examples
#' data(belgium_calls)
#'
#' # Fit competing methods and compare
#' fit_ols  <- lm(calls ~ year, data = belgium_calls)
#' fit_gkrr <- gkrr(calls ~ year, data = belgium_calls, sigma_method = "s3")
#'
#' plot(belgium_calls$year + 1900, belgium_calls$calls,
#'      xlab = "Year", ylab = "Calls (tens of millions)",
#'      main = "Belgium International Calls",
#'      pch = 19, col = "grey60")
#' abline(fit_ols,  col = "black",     lwd = 2, lty = 2)
#' abline(fit_gkrr, col = "firebrick", lwd = 2)
#' legend("topleft", c("OLS", "GKRR"), col = c("black","firebrick"),
#'        lty = c(2,1), lwd = 2, bty = "n")
"belgium_calls"


#' Cloud Point of a Liquid
#'
#' Measurements on the cloud point of a liquid, which is a measure of the
#' degree of crystallization in a stock.  The dataset contains three leverage
#' points (baseline measurements at \code{percentage_i8 = 0}) that pull
#' poorly-specified models away from the underlying trend.
#'
#' This dataset is used in De Carvalho et al. (2017), Section 6.3, to
#' illustrate robustness to leverage points.  The recommended width estimator
#' is \code{sigma_method = "s3"}.
#'
#' @format A data frame with 19 rows and 2 columns:
#' \describe{
#'   \item{percentage_i8}{Percentage of I-8 in the liquid mixture (0--10).}
#'   \item{cloud_point}{Cloud point temperature of the liquid (degrees Celsius).}
#' }
#'
#' @section Outlier structure:
#' Observations 1, 10 and 16 have \code{percentage_i8 = 0} and act as leverage
#' points because they correspond to a baseline condition that is far from the
#' centroid of the predictor space.
#'
#' @source
#' N.R. Draper and R. Craig Smith (1998).
#' \emph{Applied Regression Analysis}, 3rd edition.
#' John Wiley & Sons, Exercise 4, p. 629.
#'
#' Available in \pkg{robustbase} as \code{cloud}.
#'
#' @references
#' De Carvalho, F.A.T., Lima Neto, E.A., Ferreira, M.R.P. (2017).
#' \doi{10.1016/j.neucom.2016.12.035}
#'
#' @examples
#' data(cloud_point)
#'
#' fit_ols  <- lm(cloud_point ~ percentage_i8, data = cloud_point)
#' fit_gkrr <- gkrr(cloud_point ~ percentage_i8, data = cloud_point,
#'                  sigma_method = "s3")
#'
#' plot(cloud_point$percentage_i8, cloud_point$cloud_point,
#'      xlab = "Percentage I-8", ylab = "Cloud point (C)",
#'      main = "Cloud Point Data",
#'      pch = 19, col = "grey60")
#' abline(fit_ols,  col = "black",     lwd = 2, lty = 2)
#' abline(fit_gkrr, col = "firebrick", lwd = 2)
#' legend("topleft", c("OLS", "GKRR"), col = c("black","firebrick"),
#'        lty = c(2,1), lwd = 2, bty = "n")
"cloud_point"


#' Kootenay River Water-Flow Measurements
#'
#' Annual water-flow measurements (in units of \eqn{10^6} cubic feet) at two
#' gauging stations on the Kootenay River: Libby (Montana, USA) and Newgate
#' (British Columbia, Canada).  The dataset contains a single X-space outlier
#' (year 1934) where the Libby measurement is anomalously large, likely due to
#' an upstream dam operation.
#'
#' This dataset is used in De Carvalho et al. (2017), Section 6.5, to
#' illustrate robustness to X-space outliers.  The recommended width estimator
#' is \code{sigma_method = "s1"}.
#'
#' @format A data frame with 13 rows and 2 columns, with row names giving the
#'   year of measurement (1931--1943):
#' \describe{
#'   \item{libby}{Water flow at Libby, Montana (\eqn{10^6} cubic feet).}
#'   \item{newgate}{Water flow at Newgate, British Columbia
#'     (\eqn{10^6} cubic feet).}
#' }
#'
#' @section Outlier structure:
#' The observation for year 1934 has an anomalously large \code{libby} value
#' (77.6 vs. a typical range of 17--39), making it an X-space outlier that
#' strongly influences naive regression fits.
#'
#' @source
#' J. Neter, M.H. Kutner, C.J. Nachtsheim, and W. Wasserman (1996).
#' \emph{Applied Linear Statistical Models}, 4th edition.
#' Irwin/McGraw-Hill, p. 414.
#'
#' Available in \pkg{robustbase} as \code{kootenay}.
#'
#' @references
#' De Carvalho, F.A.T., Lima Neto, E.A., Ferreira, M.R.P. (2017).
#' \doi{10.1016/j.neucom.2016.12.035}
#'
#' @examples
#' data(kootenay)
#'
#' fit_ols  <- lm(newgate ~ libby, data = kootenay)
#' fit_gkrr <- gkrr(newgate ~ libby, data = kootenay, sigma_method = "s1")
#'
#' plot(kootenay$libby, kootenay$newgate,
#'      xlab = "Libby flow (10^6 ft^3)", ylab = "Newgate flow (10^6 ft^3)",
#'      main = "Kootenay River Water Flow",
#'      pch = 19, col = "grey60")
#' # Highlight the X-space outlier (1934)
#' points(kootenay["1934", "libby"], kootenay["1934", "newgate"],
#'        col = "red", pch = 17, cex = 1.5)
#' abline(fit_ols,  col = "black",     lwd = 2, lty = 2)
#' abline(fit_gkrr, col = "firebrick", lwd = 2)
#' legend("topleft", c("OLS", "GKRR", "1934 (outlier)"),
#'        col = c("black","firebrick","red"),
#'        lty = c(2,1,NA), pch = c(NA,NA,17), lwd = 2, bty = "n")
"kootenay"


#' Soft-Drink Delivery Time Data
#'
#' Data on the time required by a route driver to service vending machines at
#'25 stops.  Two predictors (number of products delivered and distance
#' travelled) are used to model delivery time.  Observation 9 is a bad leverage
#' point: it has both a very large number of products (30) and a very large
#' distance (1460 feet), pulling regression fits away from the main trend.
#'
#' This dataset is used in De Carvalho et al. (2017), Section 6.4, to
#' illustrate robustness to leverage points in multiple regression.  The
#' recommended width estimator is \code{sigma_method = "s3"}.
#'
#' @format A data frame with 25 rows and 3 columns:
#' \describe{
#'   \item{n_products}{Number of products stocked (cases).}
#'   \item{distance}{Distance walked by the driver (feet).}
#'   \item{delivery_time}{Total service time at the stop (minutes).}
#' }
#'
#' @section Outlier structure:
#' Observation 9 (\code{n_products = 30}, \code{distance = 1460}) is a
#' leverage point that simultaneously deviates in the predictor space and
#' has an unusually large response, making it a bad leverage point.
#'
#' @source
#' D.C. Montgomery and E.A. Peck (1992).
#' \emph{Introduction to Linear Regression Analysis}, 2nd edition.
#' John Wiley & Sons, Table B.7.
#'
#' Available in \pkg{robustbase} as \code{delivery}.
#'
#' @references
#' De Carvalho, F.A.T., Lima Neto, E.A., Ferreira, M.R.P. (2017).
#' \doi{10.1016/j.neucom.2016.12.035}
#'
#' @examples
#' data(delivery)
#'
#' fit_ols  <- lm(delivery_time ~ n_products + distance, data = delivery)
#' fit_gkrr <- gkrr(delivery_time ~ n_products + distance, data = delivery,
#'                  sigma_method = "s3")
#'
#' # Compare coefficients
#' rbind(OLS  = coef(fit_ols),
#'       GKRR = coef(fit_gkrr))
#'
#' # Diagnostic plot: kernel weights highlight observation 9
#' plot(fit_gkrr, which = 4, ask = FALSE)
"delivery"


#' Brain and Body Mass of 62 Mammal Species
#'
#' Body mass (kg) and brain mass (g) for 62 species of land mammals.
#' On the natural logarithm scale the relationship is approximately linear, but
#' the dataset contains several leverage points (very large mammals such as
#' African and Asian elephants) that strongly influence OLS estimates.
#'
#' @format A data frame with 62 rows and 5 columns:
#' \describe{
#'   \item{species}{Common name of the mammal species (character).}
#'   \item{body_mass}{Body mass in kilograms.}
#'   \item{brain_mass}{Brain mass in grams.}
#'   \item{log_body}{Natural logarithm of body mass.}
#'   \item{log_brain}{Natural logarithm of brain mass.}
#' }
#'
#' @section Outlier structure:
#' On the log scale, African elephant (\code{log_body} \eqn{\approx} 8.8) and
#' Asian elephant (\code{log_body} \eqn{\approx} 7.8) are high-leverage
#' observations.  Additionally, some primates (e.g., humans) deviate from the
#' overall trend (Y-space outliers relative to the bulk regression line).
#'
#' @source
#' P.J. Rousseeuw and A.M. Leroy (1987).
#' \emph{Robust Regression and Outlier Detection}.
#' John Wiley & Sons.
#'
#' Originally from: T. Allison and D.V. Cicchetti (1976).
#' Sleep in mammals: ecological and constitutional correlates.
#' \emph{Science}, 194, 732--734.
#'
#' Available in \pkg{MASS} as \code{mammals}.
#'
#' @examples
#' data(mammals)
#'
#' fit_ols  <- lm(log_brain ~ log_body, data = mammals)
#' fit_gkrr <- gkrr(log_brain ~ log_body, data = mammals,
#'                  sigma_method = "s3")
#'
#' plot(mammals$log_body, mammals$log_brain,
#'      xlab = "log(body mass, kg)", ylab = "log(brain mass, g)",
#'      main = "Brain vs. Body Mass (62 Mammal Species)",
#'      pch = 19, col = "grey60", cex = 0.8)
#' # Label the leverage points
#' elephants <- mammals$species %in% c("African elephant","Asian elephant")
#' points(mammals$log_body[elephants], mammals$log_brain[elephants],
#'        col = "red", pch = 17, cex = 1.4)
#' abline(fit_ols,  col = "black",     lwd = 2, lty = 2)
#' abline(fit_gkrr, col = "firebrick", lwd = 2)
#' legend("topleft", c("OLS","GKRR","Elephants"),
#'        col = c("black","firebrick","red"),
#'        lty = c(2,1,NA), pch = c(NA,NA,17), lwd = 2, bty = "n")
"mammals"


#' Hertzsprung-Russell Diagram of Star Cluster CYG OB1
#'
#' The Hertzsprung-Russell diagram plots the log-luminosity against the
#' log-effective-temperature of 47 stars in the star cluster CYG OB1.
#' The dataset contains 4 giant stars (observations 11, 20, 30 and 34) that
#' are well off the main sequence and act as leverage points, severely
#' distorting OLS estimates of the linear trend in the main-sequence stars.
#'
#' @format A data frame with 47 rows and 2 columns:
#' \describe{
#'   \item{log_temp}{Logarithm (base 10) of the effective surface temperature
#'     of the star (Kelvin).  Higher values correspond to hotter, bluer stars.}
#'   \item{log_light}{Logarithm (base 10) of the light intensity of the star
#'     relative to that of the Sun.}
#' }
#'
#' @section Outlier structure:
#' Observations 11, 20, 30 and 34 are giant stars that lie far from the main
#' sequence (low temperature, high luminosity).  They are leverage points
#' because their \code{log_temp} values are much lower than the bulk of the
#' data, and they also deviate from the linear relationship of the main sequence
#' (bad leverage points).
#'
#' @source
#' Humphreys, R.M. (1978).
#' Studies of luminous stars in nearby galaxies.
#' \emph{Astrophysical Journal Supplement}, 38, 309--350.
#'
#' Cited by: P.J. Rousseeuw and A.M. Leroy (1987).
#' \emph{Robust Regression and Outlier Detection}.
#' John Wiley & Sons, p. 27.
#'
#' Available in \pkg{robustbase} as \code{starsCYG}.
#'
#' @references
#' Rousseeuw, P.J. and Leroy, A.M. (1987).
#' \emph{Robust Regression and Outlier Detection}.
#' John Wiley & Sons.
#'
#' @examples
#' data(stars_cyg)
#'
#' fit_ols  <- lm(log_light ~ log_temp, data = stars_cyg)
#' fit_gkrr <- gkrr(log_light ~ log_temp, data = stars_cyg,
#'                  sigma_method = "s3")
#'
#' plot(stars_cyg$log_temp, stars_cyg$log_light,
#'      xlab = "log10(Temperature, K)", ylab = "log10(Luminosity, Sun = 1)",
#'      main = "Hertzsprung-Russell Diagram (CYG OB1)",
#'      pch = 19, col = "grey60")
#' # Highlight the 4 giant stars
#' giants <- c(11, 20, 30, 34)
#' points(stars_cyg$log_temp[giants], stars_cyg$log_light[giants],
#'        col = "red", pch = 17, cex = 1.5)
#' abline(fit_ols,  col = "black",     lwd = 2, lty = 2)
#' abline(fit_gkrr, col = "firebrick", lwd = 2)
#' legend("topleft", c("OLS","GKRR","Giant stars"),
#'        col = c("black","firebrick","red"),
#'        lty = c(2,1,NA), pch = c(NA,NA,17), lwd = 2, bty = "n")
"stars_cyg"
