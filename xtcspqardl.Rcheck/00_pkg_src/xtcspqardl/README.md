# xtcspqardl: Cross-Sectionally Augmented Panel Quantile ARDL

## Overview

The `xtcspqardl` package implements the **Cross-Sectionally Augmented Panel Quantile ARDL** (CS-PQARDL) model and the **Quantile Common Correlated Effects Mean Group** (QCCEMG/QCCEPMG) estimator for panel data with cross-sectional dependence.

### Key Features

- **QCCEMG**: Quantile CCE Mean Group estimator (Harding, Lamarche & Pesaran, 2018)
- **QCCEPMG**: Quantile CCE Pooled Mean Group estimator
- **CS-PQARDL**: Cross-sectionally augmented Panel Quantile ARDL
- Handles **cross-sectional dependence** through CCE augmentation (Pesaran, 2006)
- **Lagged CSA** following Chudik & Pesaran (2015)
- **Long-run coefficients** with delta-method standard errors
- **Speed of adjustment** and half-life calculations

## Installation

### From CRAN (when available)

```r
install.packages("xtcspqardl")
```

### From GitHub

```r
# install.packages("devtools")
devtools::install_github("muhammedalkhalaf/xtcspqardl")
```

## Usage

### QCCEMG Estimation

```r
library(xtcspqardl)

# Generate panel data
set.seed(123)
N <- 20  # panels
T <- 50  # time periods

data <- data.frame(
  id = rep(1:N, each = T),
  time = rep(1:T, N),
  x = rnorm(N * T)
)

# Add dynamics
for (i in 1:N) {
  idx <- ((i-1)*T + 2):(i*T)
  data$y[idx] <- 0.5 * data$y[idx-1] + 0.3 * data$x[idx] + rnorm(T-1, sd=0.5)
}

# Estimate QCCEMG
fit <- xtcspqardl(
  formula = y ~ x,
  data = data,
  id = "id",
  time = "time",
  tau = c(0.25, 0.50, 0.75),
  estimator = "qccemg"
)

# View results
summary(fit)
```

### CS-PQARDL Estimation

```r
# With long-run variables
fit_ardl <- xtcspqardl(
  formula = y ~ dx | x,  # dx is short-run, x is long-run
  data = data,
  id = "id",
  time = "time",
  tau = c(0.25, 0.50, 0.75),
  estimator = "cspqardl",
  p = 1,
  q = 1
)

summary(fit_ardl)
```

## Methodology

### Cross-Sectional Averages (CSA)

The CCE approach augments individual regressions with cross-sectional averages:

$$\bar{z}_t = \frac{1}{N} \sum_{i=1}^{N} z_{it}$$

where $z$ includes both the dependent and independent variables.

### QCCEMG Model

$$y_{it} = \lambda_i(\tau) y_{i,t-1} + \beta_i(\tau)' x_{it} + \delta_i(\tau)' \bar{z}_t + u_{it}(\tau)$$

### Long-Run Coefficients

Long-run effects are computed as:

$$\theta(\tau) = \frac{\beta(\tau)}{1 - \lambda(\tau)}$$

### CSA Lag Order

Following Chudik & Pesaran (2015), the default lag order for CSA is:

$$p_T = \lfloor T^{1/3} \rfloor$$

## References

- Chudik, A. and Pesaran, M.H. (2015). Common Correlated Effects Estimation of Heterogeneous Dynamic Panel Data Models with Weakly Exogenous Regressors. *Journal of Econometrics*, 188(2), 393-420. [DOI: 10.1016/j.jeconom.2015.03.007](https://doi.org/10.1016/j.jeconom.2015.03.007)

- Harding, M., Lamarche, C., and Pesaran, M.H. (2018). Common Correlated Effects Estimation of Heterogeneous Dynamic Panel Quantile Regression Models. *Journal of Applied Econometrics*, 35(3), 294-314. [DOI: 10.1016/j.jeconom.2018.07.010](https://doi.org/10.1016/j.jeconom.2018.07.010)

- Pesaran, M.H. (2006). Estimation and Inference in Large Heterogeneous Panels with a Multifactor Error Structure. *Econometrica*, 74(4), 967-1012. [DOI: 10.1002/jae.890](https://doi.org/10.1002/jae.890)

- Pesaran, M.H., Shin, Y., and Smith, R.J. (2001). Bounds Testing Approaches to the Analysis of Level Relationships. *Journal of Applied Econometrics*, 16(3), 289-326. [DOI: 10.1002/jae.616](https://doi.org/10.1002/jae.616)

## Author

**Merwan Roudane** - merwanroudane920@gmail.com

## License

GPL (>= 3)
