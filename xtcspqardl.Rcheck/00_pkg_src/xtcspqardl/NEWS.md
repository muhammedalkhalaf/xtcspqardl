# xtcspqardl 1.0.0

## Initial Release

* Implements Cross-Sectionally Augmented Panel Quantile ARDL (CS-PQARDL) estimation
* QCCEMG (Quantile CCE Mean Group) estimator following Harding, Lamarche & Pesaran (2018)
* QCCEPMG (Quantile CCE Pooled Mean Group) estimator
* Cross-sectional dependence handling via CCE approach (Pesaran, 2006)
* Automatic CSA lag selection following Chudik & Pesaran (2015): floor(T^{1/3})
* Long-run coefficient estimation with delta-method standard errors
* Speed of adjustment and half-life calculations
* Mean group variance estimation
* S3 methods: print, summary, coef, vcov
* Comprehensive documentation with DOI references
* Test suite using testthat
