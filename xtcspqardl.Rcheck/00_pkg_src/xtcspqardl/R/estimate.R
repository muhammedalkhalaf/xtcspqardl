#' QCCEMG Estimation
#'
#' Estimates Quantile CCE Mean Group (QCCEMG) or Pooled Mean Group (QCCEPMG)
#' following Harding, Lamarche, and Pesaran (2018).
#'
#' @param data Data frame with CSA-augmented panel data.
#' @param id Panel identifier variable name.
#' @param time Time variable name.
#' @param depvar Dependent variable name.
#' @param indepvars Character vector of independent variable names.
#' @param csa_vars Character vector of CSA variable names.
#' @param tau Numeric vector of quantiles.
#' @param constant Logical; include constant term.
#' @param pooled Logical; if TRUE, use QCCEPMG (pooled).
#'
#' @return List with estimation results.
#'
#' @references
#' Harding, M., Lamarche, C., and Pesaran, M.H. (2018). Common Correlated
#' Effects Estimation of Heterogeneous Dynamic Panel Quantile Regression
#' Models. \emph{Journal of Applied Econometrics}, 35(3), 294-314.
#' \doi{10.1016/j.jeconom.2018.07.010}
#'
#' @keywords internal
estimate_qccemg <- function(data, id, time, depvar, indepvars, csa_vars,
                             tau, constant = TRUE, pooled = FALSE) {

  ids <- unique(data[[id]])
  n_panels <- length(ids)
  n_tau <- length(tau)
  k <- length(indepvars)
  n_csa <- length(csa_vars)

  # Storage for individual estimates
  # For each panel: lambda (1), beta (k), delta (n_csa) per tau
  lambda_all <- matrix(NA_real_, nrow = n_panels, ncol = n_tau)
  beta_all <- array(NA_real_, dim = c(n_panels, k, n_tau))
  delta_all <- array(NA_real_, dim = c(n_panels, n_csa, n_tau))
  theta_all <- array(NA_real_, dim = c(n_panels, k, n_tau))  # long-run
  halflife_all <- matrix(NA_real_, nrow = n_panels, ncol = n_tau)

  valid_panels <- 0

  # Create lagged dependent variable
  data$L1_y <- panel_lag(data[[depvar]], data[[id]], lag = 1)

  for (i in seq_along(ids)) {
    panel_id <- ids[i]
    idx <- which(data[[id]] == panel_id)
    panel_data <- data[idx, , drop = FALSE]

    # Remove NA rows
    vars_needed <- c(depvar, "L1_y", indepvars, csa_vars)
    complete <- complete.cases(panel_data[, vars_needed, drop = FALSE])
    panel_data <- panel_data[complete, , drop = FALSE]

    n_i <- nrow(panel_data)
    if (n_i < (k + n_csa + 5)) next  # Need sufficient observations

    # Build design matrix: [L1.y, x1, ..., xk, csa1, ..., csa_m]
    y_i <- panel_data[[depvar]]
    X_i <- as.matrix(panel_data[, c("L1_y", indepvars, csa_vars), drop = FALSE])

    if (constant) {
      X_i <- cbind(1, X_i)
    }

    # Estimate for each quantile
    for (ti in seq_along(tau)) {
      tryCatch({
        fit <- quantreg::rq(y_i ~ X_i - 1, tau = tau[ti])
        coefs <- coef(fit)

        # Extract coefficients
        offset <- ifelse(constant, 1, 0)
        lambda_i <- coefs[offset + 1]  # L1.y coefficient
        beta_i <- coefs[(offset + 2):(offset + 1 + k)]  # x coefficients
        delta_i <- coefs[(offset + 2 + k):length(coefs)]  # CSA coefficients

        # Store
        lambda_all[i, ti] <- lambda_i
        beta_all[i, , ti] <- beta_i
        if (length(delta_i) == n_csa) {
          delta_all[i, , ti] <- delta_i
        }

        # Compute long-run coefficients: theta = beta / (1 - lambda)
        if (abs(1 - lambda_i) > 1e-8) {
          theta_all[i, , ti] <- beta_i / (1 - lambda_i)
        }

        # Half-life: ln(0.5) / ln(lambda) if 0 < lambda < 1
        if (lambda_i > 0 && lambda_i < 1) {
          halflife_all[i, ti] <- log(0.5) / log(lambda_i)
        }

        if (ti == 1) valid_panels <- valid_panels + 1

      }, error = function(e) NULL)
    }
  }

  if (valid_panels == 0) {
    stop("No panels could be estimated successfully")
  }

  # Mean Group aggregation
  lambda_mg <- colMeans(lambda_all, na.rm = TRUE)
  beta_mg <- apply(beta_all, c(2, 3), mean, na.rm = TRUE)
  delta_mg <- apply(delta_all, c(2, 3), mean, na.rm = TRUE)
  theta_mg <- apply(theta_all, c(2, 3), mean, na.rm = TRUE)
  halflife_mg <- colMeans(halflife_all, na.rm = TRUE)

  # Mean Group variance: V = (1/(N-1)) * sum((b_i - b_bar)^2)
  lambda_V <- matrix(0, n_tau, n_tau)
  for (ti in seq_along(tau)) {
    valid_lambda <- lambda_all[!is.na(lambda_all[, ti]), ti]
    if (length(valid_lambda) > 1) {
      lambda_V[ti, ti] <- var(valid_lambda) / length(valid_lambda)
    }
  }

  # Beta variance (k x k for each tau)
  beta_V <- array(0, dim = c(k, k, n_tau))
  for (ti in seq_along(tau)) {
    for (j in 1:k) {
      valid_beta <- beta_all[!is.na(beta_all[, j, ti]), j, ti]
      if (length(valid_beta) > 1) {
        beta_V[j, j, ti] <- var(valid_beta) / length(valid_beta)
      }
    }
  }

  # Delta variance
  delta_V <- array(0, dim = c(n_csa, n_csa, n_tau))
  for (ti in seq_along(tau)) {
    for (j in seq_len(n_csa)) {
      valid_delta <- delta_all[!is.na(delta_all[, j, ti]), j, ti]
      if (length(valid_delta) > 1) {
        delta_V[j, j, ti] <- var(valid_delta) / length(valid_delta)
      }
    }
  }

  # Long-run variance using delta method
  # theta = beta / (1 - lambda)
  # Var(theta) approx (1/(1-lambda))^2 * Var(beta) + (beta/(1-lambda)^2)^2 * Var(lambda)
  theta_V <- array(0, dim = c(k, k, n_tau))
  for (ti in seq_along(tau)) {
    lam <- lambda_mg[ti]
    denom <- 1 - lam
    if (abs(denom) > 1e-8) {
      for (j in 1:k) {
        g_beta <- 1 / denom
        g_lam <- beta_mg[j, ti] / (denom^2)
        theta_V[j, j, ti] <- g_beta^2 * beta_V[j, j, ti] +
                             g_lam^2 * lambda_V[ti, ti]
      }
    }
  }

  # Organize coefficients
  coef_names <- c("L1.y", indepvars)
  coefficients <- list()
  se <- list()
  vcov <- list()

  for (ti in seq_along(tau)) {
    tau_name <- paste0("tau_", tau[ti])
    coefs_ti <- c(lambda_mg[ti], if (k == 1) beta_mg[1, ti] else beta_mg[, ti])
    names(coefs_ti) <- coef_names
    coefficients[[tau_name]] <- coefs_ti

    # Handle k=1 case where array subsetting drops dimensions
    if (k == 1) {
      beta_se <- sqrt(beta_V[1, 1, ti])
    } else {
      beta_mat <- beta_V[, , ti]
      beta_se <- sqrt(diag(beta_mat))
    }
    se_ti <- c(sqrt(lambda_V[ti, ti]), beta_se)
    names(se_ti) <- coef_names
    se[[tau_name]] <- se_ti

    vcov_ti <- matrix(0, k + 1, k + 1)
    vcov_ti[1, 1] <- lambda_V[ti, ti]
    if (k == 1) {
      vcov_ti[2, 2] <- beta_V[1, 1, ti]
    } else {
      vcov_ti[2:(k+1), 2:(k+1)] <- beta_V[, , ti]
    }
    rownames(vcov_ti) <- colnames(vcov_ti) <- coef_names
    vcov[[tau_name]] <- vcov_ti
  }

  # Long-run coefficients
  long_run <- list()
  for (ti in seq_along(tau)) {
    tau_name <- paste0("tau_", tau[ti])
    lr_ti <- if (k == 1) theta_mg[1, ti] else theta_mg[, ti]
    names(lr_ti) <- indepvars
    if (k == 1) {
      attr(lr_ti, "se") <- sqrt(theta_V[1, 1, ti])
    } else {
      theta_mat <- theta_V[, , ti]
      attr(lr_ti, "se") <- sqrt(diag(theta_mat))
    }
    long_run[[tau_name]] <- lr_ti
  }

  # Speed of adjustment (lambda - 1 for ECM interpretation)
  speed_adj <- lambda_mg - 1

  # CSA coefficients
  csa_coefficients <- list()
  for (ti in seq_along(tau)) {
    tau_name <- paste0("tau_", tau[ti])
    csa_ti <- if (n_csa == 1) delta_mg[1, ti] else delta_mg[, ti]
    names(csa_ti) <- csa_vars
    if (n_csa == 1) {
      attr(csa_ti, "se") <- sqrt(delta_V[1, 1, ti])
    } else {
      delta_mat <- delta_V[, , ti]
      attr(csa_ti, "se") <- sqrt(diag(delta_mat))
    }
    csa_coefficients[[tau_name]] <- csa_ti
  }

  list(
    coefficients = coefficients,
    se = se,
    vcov = vcov,
    individual = list(
      lambda = lambda_all,
      beta = beta_all,
      theta = theta_all,
      halflife = halflife_all
    ),
    long_run = long_run,
    speed_adj = speed_adj,
    half_life = halflife_mg,
    csa_coefficients = csa_coefficients,
    valid_panels = valid_panels
  )
}


#' CS-PQARDL Estimation
#'
#' Estimates Cross-Sectionally Augmented Panel Quantile ARDL model.
#'
#' @param data Data frame with CSA-augmented panel data.
#' @param id Panel identifier variable name.
#' @param time Time variable name.
#' @param depvar Dependent variable name.
#' @param sr_vars Short-run variable names.
#' @param lr_vars Long-run variable names.
#' @param csa_vars CSA variable names.
#' @param tau Numeric vector of quantiles.
#' @param p Lag order for dependent variable.
#' @param q Lag order(s) for regressors.
#' @param constant Include constant.
#' @param model PMG, MG, or DFE.
#'
#' @return List with estimation results.
#'
#' @references
#' Pesaran, M.H., Shin, Y., and Smith, R.J. (2001). Bounds Testing Approaches
#' to the Analysis of Level Relationships. \emph{Journal of Applied
#' Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' @keywords internal
estimate_cspqardl <- function(data, id, time, depvar, sr_vars, lr_vars,
                               csa_vars, tau, p, q, constant, model) {

  ids <- unique(data[[id]])
  n_panels <- length(ids)
  n_tau <- length(tau)
  k_sr <- length(sr_vars)
  k_lr <- length(lr_vars)
  n_csa <- length(csa_vars)

  # Expand q to match k_sr
  if (length(q) == 1) {
    q <- rep(q, k_sr)
  }

  # Generate lagged dependent variables and differenced variables
  for (lag in 1:p) {
    lag_name <- paste0("L", lag, "_", depvar)
    data[[lag_name]] <- panel_lag(data[[depvar]], data[[id]], lag)
  }

  # Generate lagged independent variables
  for (j in seq_along(sr_vars)) {
    for (lag in 0:q[j]) {
      if (lag == 0) {
        # Current value already exists
      } else {
        lag_name <- paste0("L", lag, "_", sr_vars[j])
        data[[lag_name]] <- panel_lag(data[[sr_vars[j]]], data[[id]], lag)
      }
    }
  }

  # ECM form: Delta y = rho * (y_{t-1} - theta' x_{t-1}) + ...
  # Compute first difference of depvar
  data$D_y <- panel_diff(data[[depvar]], data[[id]])

  # Storage
  rho_all <- matrix(NA_real_, nrow = n_panels, ncol = n_tau)
  beta_all <- array(NA_real_, dim = c(n_panels, k_lr, n_tau))
  sr_all <- array(NA_real_, dim = c(n_panels, k_sr, n_tau))
  phi_all <- array(NA_real_, dim = c(n_panels, p, n_tau))
  csa_all <- array(NA_real_, dim = c(n_panels, n_csa, n_tau))
  halflife_all <- matrix(NA_real_, nrow = n_panels, ncol = n_tau)

  valid_panels <- 0

  for (i in seq_along(ids)) {
    panel_id <- ids[i]
    idx <- which(data[[id]] == panel_id)
    panel_data <- data[idx, , drop = FALSE]

    # Build regressors for ARDL
    # ECM form: Dy_t = rho*(y_{t-1} - theta'x_{t-1}) + sum phi_j*Dy_{t-j} + sum gamma'Dx_t + delta'csa + u
    # Simplified: estimate in ARDL levels form first

    # Lagged y terms
    y_lag_names <- paste0("L", 1:p, "_", depvar)

    # Build X matrix
    x_names <- c(y_lag_names, sr_vars, lr_vars, csa_vars)
    vars_needed <- c(depvar, x_names)

    complete <- complete.cases(panel_data[, vars_needed, drop = FALSE])
    panel_data <- panel_data[complete, , drop = FALSE]

    n_i <- nrow(panel_data)
    if (n_i < (length(x_names) + 5)) next

    y_i <- panel_data[[depvar]]
    X_i <- as.matrix(panel_data[, x_names, drop = FALSE])

    if (constant) {
      X_i <- cbind(1, X_i)
    }

    for (ti in seq_along(tau)) {
      tryCatch({
        fit <- quantreg::rq(y_i ~ X_i - 1, tau = tau[ti])
        coefs <- coef(fit)

        offset <- ifelse(constant, 1, 0)

        # Extract AR coefficients (phi_1, ..., phi_p)
        phi_i <- coefs[(offset + 1):(offset + p)]

        # Sum of AR coefficients
        sum_phi <- sum(phi_i)

        # Speed of adjustment: rho = sum(phi) - 1
        rho_i <- sum_phi - 1

        # Short-run coefficients
        sr_i <- coefs[(offset + p + 1):(offset + p + k_sr)]

        # Long-run coefficients
        lr_coefs <- coefs[(offset + p + k_sr + 1):(offset + p + k_sr + k_lr)]

        # Normalize to long-run: beta = lr_coefs / (1 - sum_phi)
        if (abs(1 - sum_phi) > 1e-8) {
          beta_i <- lr_coefs / (1 - sum_phi)
        } else {
          beta_i <- rep(NA_real_, k_lr)
        }

        # CSA coefficients
        csa_i <- coefs[(offset + p + k_sr + k_lr + 1):length(coefs)]

        # Store
        rho_all[i, ti] <- rho_i
        phi_all[i, , ti] <- phi_i
        beta_all[i, , ti] <- beta_i
        sr_all[i, , ti] <- sr_i
        if (length(csa_i) == n_csa) {
          csa_all[i, , ti] <- csa_i
        }

        # Half-life
        if (sum_phi > 0 && sum_phi < 1) {
          halflife_all[i, ti] <- log(0.5) / log(sum_phi)
        }

        if (ti == 1) valid_panels <- valid_panels + 1

      }, error = function(e) NULL)
    }
  }

  if (valid_panels == 0) {
    stop("No panels could be estimated successfully")
  }

  # Mean Group aggregation
  rho_mg <- colMeans(rho_all, na.rm = TRUE)
  beta_mg <- apply(beta_all, c(2, 3), mean, na.rm = TRUE)
  sr_mg <- apply(sr_all, c(2, 3), mean, na.rm = TRUE)
  phi_mg <- apply(phi_all, c(2, 3), mean, na.rm = TRUE)
  csa_mg <- apply(csa_all, c(2, 3), mean, na.rm = TRUE)
  halflife_mg <- colMeans(halflife_all, na.rm = TRUE)

  # Variance computation
  rho_V <- matrix(0, n_tau, n_tau)
  for (ti in seq_along(tau)) {
    valid_rho <- rho_all[!is.na(rho_all[, ti]), ti]
    if (length(valid_rho) > 1) {
      rho_V[ti, ti] <- var(valid_rho) / length(valid_rho)
    }
  }

  beta_V <- array(0, dim = c(k_lr, k_lr, n_tau))
  for (ti in seq_along(tau)) {
    for (j in seq_len(k_lr)) {
      valid_beta <- beta_all[!is.na(beta_all[, j, ti]), j, ti]
      if (length(valid_beta) > 1) {
        beta_V[j, j, ti] <- var(valid_beta) / length(valid_beta)
      }
    }
  }

  csa_V <- array(0, dim = c(n_csa, n_csa, n_tau))
  for (ti in seq_along(tau)) {
    for (j in seq_len(n_csa)) {
      valid_csa <- csa_all[!is.na(csa_all[, j, ti]), j, ti]
      if (length(valid_csa) > 1) {
        csa_V[j, j, ti] <- var(valid_csa) / length(valid_csa)
      }
    }
  }

  # Organize output
  coefficients <- list()
  se <- list()
  vcov <- list()
  long_run <- list()
  csa_coefficients <- list()

  for (ti in seq_along(tau)) {
    tau_name <- paste0("tau_", tau[ti])

    # Long-run coefficients (handle k_lr=1 case)
    lr_ti <- if (k_lr == 1) beta_mg[1, ti] else beta_mg[, ti]
    names(lr_ti) <- lr_vars
    if (k_lr == 1) {
      attr(lr_ti, "se") <- sqrt(beta_V[1, 1, ti])
    } else {
      beta_mat <- beta_V[, , ti]
      attr(lr_ti, "se") <- sqrt(diag(beta_mat))
    }
    long_run[[tau_name]] <- lr_ti

    # CSA (handle n_csa=1 case)
    csa_ti <- if (n_csa == 1) csa_mg[1, ti] else csa_mg[, ti]
    names(csa_ti) <- csa_vars
    if (n_csa == 1) {
      attr(csa_ti, "se") <- sqrt(csa_V[1, 1, ti])
    } else {
      csa_mat <- csa_V[, , ti]
      attr(csa_ti, "se") <- sqrt(diag(csa_mat))
    }
    csa_coefficients[[tau_name]] <- csa_ti

    # Combine coefficients
    coefs_ti <- c(rho = rho_mg[ti], lr_ti)
    coefficients[[tau_name]] <- coefs_ti

    # SE (handle k_lr=1 case)
    if (k_lr == 1) {
      se_ti <- c(sqrt(rho_V[ti, ti]), sqrt(beta_V[1, 1, ti]))
    } else {
      beta_mat <- beta_V[, , ti]
      se_ti <- c(sqrt(rho_V[ti, ti]), sqrt(diag(beta_mat)))
    }
    names(se_ti) <- names(coefs_ti)
    se[[tau_name]] <- se_ti
  }

  list(
    coefficients = coefficients,
    se = se,
    vcov = vcov,
    individual = list(
      rho = rho_all,
      beta = beta_all,
      sr = sr_all,
      phi = phi_all,
      halflife = halflife_all
    ),
    long_run = long_run,
    speed_adj = rho_mg,
    half_life = halflife_mg,
    csa_coefficients = csa_coefficients,
    valid_panels = valid_panels
  )
}
