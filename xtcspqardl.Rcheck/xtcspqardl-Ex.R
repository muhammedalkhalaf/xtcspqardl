pkgname <- "xtcspqardl"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "xtcspqardl-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('xtcspqardl')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("xtcspqardl")
### * xtcspqardl

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: xtcspqardl
### Title: Cross-Sectionally Augmented Panel Quantile ARDL
### Aliases: xtcspqardl

### ** Examples

## No test: 
# Generate example panel data
set.seed(123)
N <- 20  # panels
T <- 50  # time periods
data <- data.frame(
  id = rep(1:N, each = T),
  time = rep(1:T, N),
  x = rnorm(N * T),
  y = rnorm(N * T)
)
# Add dynamics
for (i in 1:N) {
  idx <- ((i-1)*T + 2):(i*T)
  data$y[idx] <- 0.5 * data$y[idx-1] + 0.3 * data$x[idx] + rnorm(T-1, sd=0.5)
}

# QCCEMG estimation
fit <- xtcspqardl(y ~ x, data = data, id = "id", time = "time",
                  tau = c(0.25, 0.50, 0.75), estimator = "qccemg")
summary(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("xtcspqardl", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
