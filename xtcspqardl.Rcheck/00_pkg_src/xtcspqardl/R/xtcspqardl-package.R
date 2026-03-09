#' xtcspqardl: Cross-Sectionally Augmented Panel Quantile ARDL
#'
#' Implements the Cross-Sectionally Augmented Panel Quantile Autoregressive
#' Distributed Lag (CS-PQARDL) model and the Quantile Common Correlated
#' Effects Mean Group (QCCEMG) estimator for panel data with cross-sectional
#' dependence.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{xtcspqardl}}: Main estimation function for CS-PQARDL
#'     and QCCEMG/QCCEPMG models.
#'   \item \code{\link{compute_csa}}: Compute cross-sectional averages for
#'     CCE augmentation.
#' }
#'
#' @section Estimators:
#' The package provides three main estimators:
#'
#' \strong{QCCEMG (Quantile CCE Mean Group):}
#' Estimates unit-by-unit quantile regressions augmented with cross-sectional
#' averages, then aggregates using the mean group estimator. This follows
#' Harding, Lamarche, and Pesaran (2018).
#'
#' \strong{QCCEPMG (Quantile CCE Pooled Mean Group):}
#' Similar to QCCEMG but pools the long-run coefficients across units while
#' allowing heterogeneous short-run dynamics.
#'
#' \strong{CS-PQARDL (CS Panel Quantile ARDL):}
#' Extends the ARDL approach to cointegration (Pesaran, Shin & Smith, 2001)
#' to quantile regression with CCE augmentation for handling cross-sectional
#' dependence.
#'
#' @section Cross-Sectional Dependence:
#' The package handles cross-sectional dependence through the Common
#' Correlated Effects (CCE) approach of Pesaran (2006). Cross-sectional
#' averages of all variables are computed and included as proxies for
#' unobserved common factors.
#'
#' Following Chudik and Pesaran (2015), lagged cross-sectional averages
#' are included with default lag order \code{floor(T^{1/3})}.
#'
#' @section References:
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
#' @docType package
#' @name xtcspqardl-package
#' @aliases xtcspqardl-package
#' @keywords internal
#' @importFrom quantreg rq
#' @importFrom stats coef var complete.cases terms pnorm setNames
"_PACKAGE"
