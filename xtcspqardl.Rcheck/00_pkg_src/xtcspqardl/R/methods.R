#' Print method for xtcspqardl objects
#'
#' @param x An object of class \code{"xtcspqardl"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.xtcspqardl <- function(x, ...) {
  cat("\n")
  cat("==========================================================================\n")

  if (x$estimator %in% c("qccemg", "qccepmg")) {
    est_label <- ifelse(x$estimator == "qccemg",
                        "QCCEMG - Quantile CCE Mean Group",
                        "QCCEPMG - Quantile CCE Pooled Mean Group")
  } else {
    est_label <- "CS-PQARDL - Cross-Sectionally Augmented Panel Quantile ARDL"
  }

  cat(sprintf("  %s\n", est_label))
  cat("  Harding, Lamarche & Pesaran (2018)\n")
  cat("==========================================================================\n\n")

  cat(sprintf("  Dependent variable:  %s\n", x$depvar))
  cat(sprintf("  Regressors:          %s\n", paste(x$indepvars, collapse = ", ")))
  if (length(x$lr_vars) > 0) {
    cat(sprintf("  Long-run vars:       %s\n", paste(x$lr_vars, collapse = ", ")))
  }
  cat(sprintf("  Panels (N):          %d\n", x$n_panels))
  cat(sprintf("  Valid panels:        %d\n", x$valid_panels))
  cat(sprintf("  Avg. T:              %.1f\n", x$avg_T))
  cat(sprintf("  Observations:        %d\n", x$n_obs))
  cat(sprintf("  CSA lags:            %d\n", x$cr_lags))
  cat(sprintf("  Quantiles:           %s\n",
              paste(sprintf("%.2f", x$tau), collapse = ", ")))
  cat("\n")

  invisible(x)
}


#' Summary method for xtcspqardl objects
#'
#' @param object An object of class \code{"xtcspqardl"}.
#' @param ... Additional arguments (ignored).
#'
#' @return An object of class \code{"summary.xtcspqardl"}.
#'
#' @export
summary.xtcspqardl <- function(object, ...) {
  # Build summary tables
  tables <- list()

  for (ti in seq_along(object$tau)) {
    tau_val <- object$tau[ti]
    tau_name <- paste0("tau_", tau_val)

    coefs <- object$coefficients[[tau_name]]
    ses <- object$se[[tau_name]]

    if (is.null(coefs) || length(coefs) == 0) next

    z_stat <- coefs / ses
    p_val <- 2 * (1 - stats::pnorm(abs(z_stat)))

    tbl <- data.frame(
      Estimate = coefs,
      `Std.Error` = ses,
      `z.value` = z_stat,
      `Pr(>|z|)` = p_val,
      check.names = FALSE
    )
    rownames(tbl) <- names(coefs)

    tables[[tau_name]] <- tbl
  }

  # Long-run table
  lr_tables <- list()
  for (ti in seq_along(object$tau)) {
    tau_val <- object$tau[ti]
    tau_name <- paste0("tau_", tau_val)

    lr <- object$long_run[[tau_name]]
    if (is.null(lr)) next

    lr_se <- attr(lr, "se")
    if (is.null(lr_se)) lr_se <- rep(NA, length(lr))

    z_stat <- lr / lr_se
    p_val <- 2 * (1 - stats::pnorm(abs(z_stat)))

    tbl <- data.frame(
      `LR.Coef` = lr,
      `Std.Error` = lr_se,
      `z.value` = z_stat,
      `Pr(>|z|)` = p_val,
      check.names = FALSE
    )
    rownames(tbl) <- names(lr)

    lr_tables[[tau_name]] <- tbl
  }

  # Speed of adjustment table
  ecm_table <- data.frame(
    Quantile = object$tau,
    `Speed.of.Adj` = object$speed_adj,
    `Half.Life` = object$half_life,
    check.names = FALSE
  )

  out <- list(
    call = object$call,
    estimator = object$estimator,
    n_panels = object$n_panels,
    valid_panels = object$valid_panels,
    n_obs = object$n_obs,
    avg_T = object$avg_T,
    cr_lags = object$cr_lags,
    tau = object$tau,
    depvar = object$depvar,
    indepvars = object$indepvars,
    tables = tables,
    lr_tables = lr_tables,
    ecm_table = ecm_table
  )

  class(out) <- "summary.xtcspqardl"
  out
}


#' Print method for summary.xtcspqardl
#'
#' @param x An object of class \code{"summary.xtcspqardl"}.
#' @param digits Number of significant digits.
#' @param signif.stars Logical; print significance stars.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.summary.xtcspqardl <- function(x, digits = 4, signif.stars = TRUE, ...) {
  cat("\n")
  cat("==========================================================================\n")

  if (x$estimator %in% c("qccemg", "qccepmg")) {
    est_label <- ifelse(x$estimator == "qccemg",
                        "QCCEMG - Quantile CCE Mean Group",
                        "QCCEPMG - Quantile CCE Pooled Mean Group")
  } else {
    est_label <- "CS-PQARDL - Cross-Sectionally Augmented Panel Quantile ARDL"
  }

  cat(sprintf("  %s\n", est_label))
  cat("==========================================================================\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat(sprintf("Panels: %d (%d valid)  |  Avg T: %.1f  |  Obs: %d  |  CSA lags: %d\n\n",
              x$n_panels, x$valid_panels, x$avg_T, x$n_obs, x$cr_lags))

  # Short-run coefficients
  cat("--------------------------------------------------------------------------\n")
  cat("                        SHORT-RUN COEFFICIENTS\n")
  cat("--------------------------------------------------------------------------\n\n")

  for (ti in seq_along(x$tau)) {
    tau_val <- x$tau[ti]
    tau_name <- paste0("tau_", tau_val)

    if (is.null(x$tables[[tau_name]])) next

    cat(sprintf("--- Quantile tau = %.2f ---\n\n", tau_val))
    print_coef_table(x$tables[[tau_name]], digits = digits,
                     signif.stars = signif.stars)
    cat("\n")
  }

  # Long-run coefficients
  if (length(x$lr_tables) > 0) {
    cat("--------------------------------------------------------------------------\n")
    cat("                        LONG-RUN COEFFICIENTS\n")
    cat("--------------------------------------------------------------------------\n\n")

    for (ti in seq_along(x$tau)) {
      tau_val <- x$tau[ti]
      tau_name <- paste0("tau_", tau_val)

      if (is.null(x$lr_tables[[tau_name]])) next

      cat(sprintf("--- Quantile tau = %.2f ---\n\n", tau_val))
      print_coef_table(x$lr_tables[[tau_name]], digits = digits,
                       signif.stars = signif.stars)
      cat("\n")
    }
  }

  # ECM speed of adjustment
  cat("--------------------------------------------------------------------------\n")
  cat("                     SPEED OF ADJUSTMENT & HALF-LIFE\n")
  cat("--------------------------------------------------------------------------\n\n")

  ecm_fmt <- x$ecm_table
  ecm_fmt[, 2] <- round(ecm_fmt[, 2], digits)
  ecm_fmt[, 3] <- round(ecm_fmt[, 3], 2)

  # Add status
  ecm_fmt$Status <- ifelse(ecm_fmt[, 2] < -0.5, "Strong",
                           ifelse(ecm_fmt[, 2] < -0.1, "Moderate",
                                  ifelse(ecm_fmt[, 2] < 0, "Weak", "No conv.")))

  print(ecm_fmt, row.names = FALSE)
  cat("\n")

  cat("--------------------------------------------------------------------------\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat("==========================================================================\n\n")

  invisible(x)
}


#' Print coefficient table with significance stars
#'
#' @param tbl Data frame with coefficient table.
#' @param digits Number of digits.
#' @param signif.stars Show significance stars.
#'
#' @keywords internal
print_coef_table <- function(tbl, digits = 4, signif.stars = TRUE) {
  # Format numbers
  tbl_fmt <- tbl
  tbl_fmt[, 1] <- sprintf(paste0("%.", digits, "f"), tbl[, 1])
  tbl_fmt[, 2] <- sprintf(paste0("%.", digits, "f"), tbl[, 2])
  tbl_fmt[, 3] <- sprintf("%.3f", tbl[, 3])
  tbl_fmt[, 4] <- sprintf("%.4f", tbl[, 4])

  if (signif.stars) {
    stars <- ifelse(tbl[, 4] < 0.001, "***",
                    ifelse(tbl[, 4] < 0.01, "**",
                           ifelse(tbl[, 4] < 0.05, "*",
                                  ifelse(tbl[, 4] < 0.1, ".", ""))))
    tbl_fmt$` ` <- stars
  }

  print(tbl_fmt, right = TRUE)
}


#' Extract coefficients from xtcspqardl object
#'
#' @param object An object of class \code{"xtcspqardl"}.
#' @param tau Optional quantile(s) to extract. If NULL, returns all.
#' @param type Character; \code{"short_run"} or \code{"long_run"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Named numeric vector or list of coefficients.
#'
#' @export
coef.xtcspqardl <- function(object, tau = NULL, type = "short_run", ...) {
  if (type == "long_run") {
    coefs <- object$long_run
  } else {
    coefs <- object$coefficients
  }

  if (is.null(tau)) {
    return(coefs)
  }

  tau_names <- paste0("tau_", tau)
  coefs[tau_names]
}


#' Extract variance-covariance matrix from xtcspqardl object
#'
#' @param object An object of class \code{"xtcspqardl"}.
#' @param tau Optional quantile to extract. If NULL, returns all.
#' @param ... Additional arguments (ignored).
#'
#' @return Variance-covariance matrix or list of matrices.
#'
#' @export
vcov.xtcspqardl <- function(object, tau = NULL, ...) {
  if (is.null(tau)) {
    return(object$vcov)
  }

  tau_name <- paste0("tau_", tau)
  object$vcov[[tau_name]]
}


#' Fitted values from xtcspqardl model
#'
#' @param object An object of class \code{"xtcspqardl"}.
#' @param ... Additional arguments (ignored).
#'
#' @return This function is not yet implemented for xtcspqardl objects.
#'
#' @export
fitted.xtcspqardl <- function(object, ...) {
  warning("fitted() is not yet implemented for xtcspqardl objects")
  NULL
}


#' Residuals from xtcspqardl model
#'
#' @param object An object of class \code{"xtcspqardl"}.
#' @param ... Additional arguments (ignored).
#'
#' @return This function is not yet implemented for xtcspqardl objects.
#'
#' @export
residuals.xtcspqardl <- function(object, ...) {
  warning("residuals() is not yet implemented for xtcspqardl objects")
  NULL
}
