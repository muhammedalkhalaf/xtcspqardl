#' Cross-Sectionally Augmented Panel Quantile ARDL
#'
#' Estimates the Cross-Sectionally Augmented Panel Quantile ARDL (CS-PQARDL)
#' model or the Quantile Common Correlated Effects Mean Group (QCCEMG/QCCEPMG)
#' estimator for panel data with cross-sectional dependence.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 | z1 + z2} where
#'   variables before \code{|} are short-run regressors and variables after
#'   \code{|} are long-run regressors (for CS-PQARDL). For QCCEMG, use
#'
#'   \code{y ~ x1 + x2}.
#' @param data A data frame containing panel data.
#' @param id Character string naming the cross-sectional unit identifier.
#' @param time Character string naming the time variable.
#' @param tau Numeric vector of quantiles to estimate (between 0 and 1).
#' @param estimator Character string specifying the estimator:
#'   \code{"qccemg"} for Quantile CCE Mean Group (default),
#'   \code{"qccepmg"} for Quantile CCE Pooled Mean Group,
#'   \code{"cspqardl"} for CS Panel Quantile ARDL.
#' @param p Integer specifying the number of lags for the dependent variable
#'   (default 1, for CS-PQARDL).
#' @param q Integer or vector specifying the number of lags for each regressor
#'   (default 1, for CS-PQARDL).
#' @param cr_lags Integer specifying the number of lags for cross-sectional
#'   averages. Default is \code{floor(T^(1/3))} following Chudik and Pesaran
#'   (2015).
#' @param constant Logical; if \code{TRUE} (default), include a constant term.
#' @param model Character string for CS-PQARDL pooling:
#'   \code{"pmg"} for Pooled Mean Group (default),
#'   \code{"mg"} for Mean Group,
#'   \code{"dfe"} for Dynamic Fixed Effects.
#'
#' @return An object of class \code{"xtcspqardl"} containing:
#' \item{coefficients}{Mean group coefficients across panels.}
#' \item{se}{Standard errors using mean group variance.}
#' \item{vcov}{Variance-covariance matrix.}
#' \item{individual}{List of unit-specific estimates.}
#' \item{long_run}{Long-run coefficient estimates.}
#' \item{speed_adj}{Speed of adjustment coefficients.}
#' \item{half_life}{Half-life of adjustment.}
#' \item{tau}{Quantiles estimated.}
#' \item{call}{The matched call.}
#' \item{formula}{The formula used.}
#' \item{n_panels}{Number of panels.}
#' \item{n_obs}{Total observations.}
#' \item{avg_T}{Average time periods per panel.}
#'
#' @details
#' The package implements two main estimators for panel quantile regression
#' with cross-sectional dependence:
#'
#' \strong{QCCEMG (Quantile CCE Mean Group):}
#' Estimates unit-by-unit quantile regressions augmented with cross-sectional
#' averages, then aggregates using mean group estimator. The model is:
#' \deqn{y_{it} = \lambda_i y_{i,t-1} + \beta_i' x_{it} + \delta_i' \bar{z}_t + u_{it}}
#' where \eqn{\bar{z}_t} contains cross-sectional averages of \eqn{y} and
#' \eqn{x}.
#'
#' \strong{CS-PQARDL (CS Panel Quantile ARDL):}
#' Extends the ARDL approach to cointegration (Pesaran, Shin & Smith, 2001)
#' to quantile regression with CCE augmentation. Estimates error-correction
#' form with long-run relationships.
#'
#' Cross-sectional dependence is handled through the CCE approach (Pesaran,
#' 2006), which augments regressions with cross-sectional averages of all
#' variables. Lagged CSA follow Chudik & Pesaran (2015) with default
#' \code{floor(T^{1/3})} lags.
#'
#' @references
#' Chudik, A. and Pesaran, M.H. (2015). Common Correlated Effects Estimation
#' of Heterogeneous Dynamic Panel Data Models with Weakly Exogenous Regressors.
#' \emph{Journal of Econometrics}, 188(2), 393-420.
#' \doi{10.1016/j.jeconom.2015.03.007}
#'
#' Harding, M., Lamarche, C., and Pesaran, M.H. (2018). Common Correlated
#' Effects Estimation of Heterogeneous Dynamic Panel Quantile Regression
#' Models. \emph{Journal of Applied Econometrics}, 35(3), 294-314.
#' \doi{10.1016/j.jeconom.2018.07.010}
#'
#' Pesaran, M.H. (2006). Estimation and Inference in Large Heterogeneous
#' Panels with a Multifactor Error Structure. \emph{Econometrica}, 74(4),
#' 967-1012. \doi{10.1002/jae.890}
#'
#' Pesaran, M.H., Shin, Y., and Smith, R.J. (2001). Bounds Testing Approaches
#' to the Analysis of Level Relationships. \emph{Journal of Applied
#' Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' @examples
#' \donttest{
#' # Generate example panel data
#' set.seed(123)
#' N <- 20  # panels
#' T <- 50  # time periods
#' data <- data.frame(
#'   id = rep(1:N, each = T),
#'   time = rep(1:T, N),
#'   x = rnorm(N * T),
#'   y = rnorm(N * T)
#' )
#' # Add dynamics
#' for (i in 1:N) {
#'   idx <- ((i-1)*T + 2):(i*T)
#'   data$y[idx] <- 0.5 * data$y[idx-1] + 0.3 * data$x[idx] + rnorm(T-1, sd=0.5)
#' }
#'
#' # QCCEMG estimation
#' fit <- xtcspqardl(y ~ x, data = data, id = "id", time = "time",
#'                   tau = c(0.25, 0.50, 0.75), estimator = "qccemg")
#' summary(fit)
#' }
#'
#' @export
xtcspqardl <- function(formula, data, id, time,
                       tau = 0.5,
                       estimator = c("qccemg", "qccepmg", "cspqardl"),
                       p = 1L, q = 1L,
                       cr_lags = NULL,
                       constant = TRUE,
                       model = c("pmg", "mg", "dfe")) {

  # Match arguments
  cl <- match.call()
  estimator <- match.arg(estimator)
  model <- match.arg(model)

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!id %in% names(data)) {
    stop(sprintf("id variable '%s' not found in data", id))
  }
  if (!time %in% names(data)) {
    stop(sprintf("time variable '%s' not found in data", time))
  }
  if (any(tau <= 0 | tau >= 1)) {
    stop("'tau' must contain values strictly between 0 and 1")
  }

  # Parse formula
  formula_parts <- parse_formula(formula)
  depvar <- formula_parts$depvar
  sr_vars <- formula_parts$sr_vars
  lr_vars <- formula_parts$lr_vars

  # Validate variables exist
  all_vars <- c(depvar, sr_vars, lr_vars)
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variables not found in data: %s",
                 paste(missing_vars, collapse = ", ")))
  }

  # Ensure proper ordering
  data <- data[order(data[[id]], data[[time]]), ]

  # Get panel info
  ids <- unique(data[[id]])
  n_panels <- length(ids)
  n_obs <- nrow(data)
  avg_T <- n_obs / n_panels

  # Default cr_lags following Chudik & Pesaran (2015)
  if (is.null(cr_lags)) {
    cr_lags <- floor(avg_T^(1/3))
  }
  cr_lags <- as.integer(max(0, cr_lags))

  # Compute cross-sectional averages
  csa_result <- compute_csa(data, id, time, depvar, sr_vars, cr_lags)
  data_aug <- csa_result$data
  csa_vars <- csa_result$csa_vars

  # Dispatch to appropriate estimator
  if (estimator %in% c("qccemg", "qccepmg")) {
    result <- estimate_qccemg(
      data = data_aug,
      id = id,
      time = time,
      depvar = depvar,
      indepvars = sr_vars,
      csa_vars = csa_vars,
      tau = tau,
      constant = constant,
      pooled = (estimator == "qccepmg")
    )
  } else {
    # CS-PQARDL
    if (length(lr_vars) == 0) {
      stop("CS-PQARDL requires long-run variables specified after '|' in formula")
    }
    result <- estimate_cspqardl(
      data = data_aug,
      id = id,
      time = time,
      depvar = depvar,
      sr_vars = sr_vars,
      lr_vars = lr_vars,
      csa_vars = csa_vars,
      tau = tau,
      p = p,
      q = q,
      constant = constant,
      model = model
    )
  }

  # Build output object
  out <- list(
    coefficients = result$coefficients,
    se = result$se,
    vcov = result$vcov,
    individual = result$individual,
    long_run = result$long_run,
    speed_adj = result$speed_adj,
    half_life = result$half_life,
    csa_coefficients = result$csa_coefficients,
    tau = tau,
    call = cl,
    formula = formula,
    estimator = estimator,
    n_panels = n_panels,
    valid_panels = result$valid_panels,
    n_obs = n_obs,
    avg_T = avg_T,
    cr_lags = cr_lags,
    depvar = depvar,
    indepvars = sr_vars,
    lr_vars = lr_vars,
    id = id,
    time = time,
    model = model
  )

  class(out) <- "xtcspqardl"
  out
}


#' Parse formula for xtcspqardl
#'
#' @param formula Formula object
#' @return List with depvar, sr_vars, lr_vars
#' @keywords internal
parse_formula <- function(formula) {
  # Get terms
  tf <- terms(formula, specials = NULL)
  vars <- all.vars(formula)
  depvar <- vars[1]

  # Check for | separator (long-run variables)
  formula_str <- deparse(formula)
  if (grepl("\\|", formula_str)) {
    # Split at |
    parts <- strsplit(formula_str, "\\|")[[1]]
    lhs_rhs <- strsplit(parts[1], "~")[[1]]
    depvar <- trimws(lhs_rhs[1])
    sr_vars <- trimws(unlist(strsplit(lhs_rhs[2], "\\+")))
    sr_vars <- sr_vars[sr_vars != ""]
    lr_vars <- trimws(unlist(strsplit(parts[2], "\\+")))
    lr_vars <- lr_vars[lr_vars != ""]
  } else {
    sr_vars <- vars[-1]
    lr_vars <- character(0)
  }

  list(depvar = depvar, sr_vars = sr_vars, lr_vars = lr_vars)
}
