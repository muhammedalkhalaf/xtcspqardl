# Tests for xtcspqardl package

test_that("xtcspqardl QCCEMG estimation works", {
  skip_if_not_installed("quantreg")

  # Generate simple panel data
  set.seed(123)
  N <- 10
  T <- 30

  data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    x = rnorm(N * T)
  )

  # Generate y with AR(1) dynamics and dependence on x
  data$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T + 1):(i * T)
    data$y[idx[1]] <- rnorm(1)
    for (t in 2:T) {
      data$y[idx[t]] <- 0.5 * data$y[idx[t-1]] + 0.3 * data$x[idx[t]] + rnorm(1, sd = 0.5)
    }
  }

  # Test QCCEMG
  fit <- xtcspqardl(
    formula = y ~ x,
    data = data,
    id = "id",
    time = "time",
    tau = c(0.25, 0.5, 0.75),
    estimator = "qccemg"
  )

  expect_s3_class(fit, "xtcspqardl")
  expect_equal(fit$estimator, "qccemg")
  expect_equal(length(fit$tau), 3)
  expect_true(fit$valid_panels > 0)
  expect_true(!is.null(fit$coefficients))
  expect_true(!is.null(fit$long_run))
})


test_that("xtcspqardl handles missing data correctly", {
  skip_if_not_installed("quantreg")

  set.seed(456)
  N <- 8
  T <- 25

  data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    x = rnorm(N * T)
  )

  data$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T + 1):(i * T)
    data$y[idx[1]] <- rnorm(1)
    for (t in 2:T) {
      data$y[idx[t]] <- 0.4 * data$y[idx[t-1]] + 0.2 * data$x[idx[t]] + rnorm(1, sd = 0.3)
    }
  }

  # Introduce some missing values
  data$y[sample(1:nrow(data), 10)] <- NA

  fit <- xtcspqardl(
    formula = y ~ x,
    data = data,
    id = "id",
    time = "time",
    tau = 0.5,
    estimator = "qccemg"
  )

  expect_s3_class(fit, "xtcspqardl")
  expect_true(fit$valid_panels > 0)
})


test_that("compute_csa generates correct variables", {
  set.seed(789)
  N <- 5
  T <- 20

  data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = rnorm(N * T),
    x = rnorm(N * T)
  )

  result <- compute_csa(data, "id", "time", "y", "x", cr_lags = 2)

  expect_true("csa_y" %in% names(result$data))
  expect_true("csa_x" %in% names(result$data))
  expect_true("L1_csa_y" %in% names(result$data))
  expect_true("L2_csa_y" %in% names(result$data))

  # Check CSA is the mean across panels at each time
  for (t in unique(data$time)) {
    idx <- data$time == t
    expect_equal(
      unique(result$data$csa_y[idx]),
      mean(data$y[idx]),
      tolerance = 1e-10
    )
  }
})


test_that("panel_lag works correctly", {
  data <- data.frame(
    id = rep(1:2, each = 5),
    x = c(1, 2, 3, 4, 5, 10, 20, 30, 40, 50)
  )

  lagged <- panel_lag(data$x, data$id, lag = 1)

  # First obs of each panel should be NA
  expect_true(is.na(lagged[1]))
  expect_true(is.na(lagged[6]))

  # Check correct lag values
  expect_equal(lagged[2:5], c(1, 2, 3, 4))
  expect_equal(lagged[7:10], c(10, 20, 30, 40))
})


test_that("print and summary methods work", {
  skip_if_not_installed("quantreg")

  set.seed(321)
  N <- 8
  T <- 25

  data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    x = rnorm(N * T)
  )

  data$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T + 1):(i * T)
    data$y[idx[1]] <- rnorm(1)
    for (t in 2:T) {
      data$y[idx[t]] <- 0.5 * data$y[idx[t-1]] + 0.3 * data$x[idx[t]] + rnorm(1, sd = 0.5)
    }
  }

  fit <- xtcspqardl(
    formula = y ~ x,
    data = data,
    id = "id",
    time = "time",
    tau = c(0.5),
    estimator = "qccemg"
  )

  # Test print
  expect_output(print(fit), "QCCEMG")

  # Test summary
  smry <- summary(fit)
  expect_s3_class(smry, "summary.xtcspqardl")
  expect_output(print(smry), "SHORT-RUN")
})


test_that("coef and vcov methods work", {
  skip_if_not_installed("quantreg")

  set.seed(654)
  N <- 8
  T <- 25

  data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    x = rnorm(N * T)
  )

  data$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T + 1):(i * T)
    data$y[idx[1]] <- rnorm(1)
    for (t in 2:T) {
      data$y[idx[t]] <- 0.5 * data$y[idx[t-1]] + 0.3 * data$x[idx[t]] + rnorm(1, sd = 0.5)
    }
  }

  fit <- xtcspqardl(
    formula = y ~ x,
    data = data,
    id = "id",
    time = "time",
    tau = c(0.25, 0.5, 0.75),
    estimator = "qccemg"
  )

  # Test coef
  all_coefs <- coef(fit)
  expect_type(all_coefs, "list")
  expect_equal(length(all_coefs), 3)

  # Test coef with specific tau
  coefs_50 <- coef(fit, tau = 0.5)
  expect_type(coefs_50, "list")
  expect_true("tau_0.5" %in% names(coefs_50))

  # Test long-run coefs
  lr_coefs <- coef(fit, type = "long_run")
  expect_type(lr_coefs, "list")

  # Test vcov
  vcov_all <- vcov(fit)
  expect_type(vcov_all, "list")
})


test_that("input validation works", {
  data <- data.frame(id = 1:10, time = 1:10, y = rnorm(10), x = rnorm(10))

  # Wrong id column

  expect_error(
    xtcspqardl(y ~ x, data, id = "wrong", time = "time", tau = 0.5),
    "not found"
  )

  # Invalid tau
  expect_error(
    xtcspqardl(y ~ x, data, id = "id", time = "time", tau = 1.5),
    "tau"
  )

  expect_error(
    xtcspqardl(y ~ x, data, id = "id", time = "time", tau = 0),
    "tau"
  )
})


test_that("formula parsing works", {
  # Simple formula
  result <- xtcspqardl:::parse_formula(y ~ x1 + x2)
  expect_equal(result$depvar, "y")
  expect_equal(result$sr_vars, c("x1", "x2"))
  expect_equal(length(result$lr_vars), 0)

  # Formula with long-run variables
  result2 <- xtcspqardl:::parse_formula(y ~ x1 + x2 | z1 + z2)
  expect_equal(result2$depvar, "y")
  expect_equal(result2$sr_vars, c("x1", "x2"))
  expect_equal(result2$lr_vars, c("z1", "z2"))
})
