#' Compute Cross-Sectional Averages
#'
#' Computes cross-sectional averages (CSA) of all variables at each time
#' period, following the CCE approach of Pesaran (2006).
#'
#' @param data Data frame with panel data.
#' @param id Character string naming the cross-sectional identifier.
#' @param time Character string naming the time variable.
#' @param depvar Character string naming the dependent variable.
#' @param indepvars Character vector of independent variable names.
#' @param cr_lags Integer number of lags for CSA (Chudik & Pesaran, 2015).
#'
#' @return A list containing:
#' \item{data}{Data frame with CSA columns added.}
#' \item{csa_vars}{Character vector of CSA variable names.}
#'
#' @details
#' Cross-sectional averages are computed as:
#' \deqn{\bar{z}_t = \frac{1}{N} \sum_{i=1}^{N} z_{it}}
#' for each variable \eqn{z} in \eqn{\{y, x_1, x_2, ...\}}.
#'
#' Following Chudik and Pesaran (2015), lagged CSA are included with
#' default lag order \code{floor(T^{1/3})}.
#'
#' @references
#' Chudik, A. and Pesaran, M.H. (2015). Common Correlated Effects Estimation
#' of Heterogeneous Dynamic Panel Data Models with Weakly Exogenous Regressors.
#' \emph{Journal of Econometrics}, 188(2), 393-420.
#' \doi{10.1016/j.jeconom.2015.03.007}
#'
#' Pesaran, M.H. (2006). Estimation and Inference in Large Heterogeneous
#' Panels with a Multifactor Error Structure. \emph{Econometrica}, 74(4),
#' 967-1012. \doi{10.1002/jae.890}
#'
#' @export
compute_csa <- function(data, id, time, depvar, indepvars, cr_lags = 0) {

  # All variables for CSA
  all_vars <- c(depvar, indepvars)
  csa_vars <- character(0)

  # Compute CSA for each time period
  time_vals <- unique(data[[time]])

  for (v in all_vars) {
    csa_name <- paste0("csa_", v)
    data[[csa_name]] <- NA_real_

    for (t in time_vals) {
      idx <- data[[time]] == t
      data[[csa_name]][idx] <- mean(data[[v]][idx], na.rm = TRUE)
    }
    csa_vars <- c(csa_vars, csa_name)
  }

  # Compute lagged CSA within each panel
  if (cr_lags > 0) {
    ids <- unique(data[[id]])
    base_csa <- csa_vars  # CSA at lag 0

    for (lag in 1:cr_lags) {
      for (csa_v in base_csa) {
        lag_name <- paste0("L", lag, "_", csa_v)
        data[[lag_name]] <- NA_real_

        for (i in ids) {
          idx_i <- which(data[[id]] == i)
          n_i <- length(idx_i)
          if (n_i > lag) {
            data[[lag_name]][idx_i[(lag + 1):n_i]] <-
              data[[csa_v]][idx_i[1:(n_i - lag)]]
          }
        }
        csa_vars <- c(csa_vars, lag_name)
      }
    }
  }

  list(data = data, csa_vars = csa_vars)
}


#' Compute lagged variable within panel
#'
#' @param x Numeric vector.
#' @param id Panel identifier vector.
#' @param lag Number of lags.
#'
#' @return Numeric vector with lagged values.
#' @keywords internal
panel_lag <- function(x, id, lag = 1) {
  if (lag == 0) return(x)

  ids <- unique(id)
  result <- rep(NA_real_, length(x))

  for (i in ids) {
    idx <- which(id == i)
    n_i <- length(idx)
    if (n_i > lag) {
      result[idx[(lag + 1):n_i]] <- x[idx[1:(n_i - lag)]]
    }
  }
  result
}


#' Compute first difference within panel
#'
#' @param x Numeric vector.
#' @param id Panel identifier vector.
#'
#' @return Numeric vector with first differences.
#' @keywords internal
panel_diff <- function(x, id) {
  ids <- unique(id)
  result <- rep(NA_real_, length(x))

  for (i in ids) {
    idx <- which(id == i)
    n_i <- length(idx)
    if (n_i > 1) {
      result[idx[2:n_i]] <- diff(x[idx])
    }
  }
  result
}
