# Regression Analysis with Average Hazard

The `ahreg` function performs a regression analysis for the average
hazard (AH).

## Usage

``` r
ahreg(formula, tau, data, link="log", conf.int=0.95, 
             cens_strata=NULL, cens_covs=NULL)
```

## Arguments

- formula:

  A formula object, with the response on the left of a ~ operator, and
  the terms on the right. The response must be a survival object as
  returned by the Surv function. For a multi-state model the formula may
  be a list of formulas.

- tau:

  A scalar value to specify a time point for calculating the average
  hazard.

- data:

  A data.frame in which to interpret the variables named in the formula.

- link:

  A link function to be used, either "log" (default) or "identity".

- conf.int:

  A confidence coefficient for calculating confidence intervals. The
  default is `conf.int=0.95`.

- cens_strata:

  A variable name for specifying group-specific censoring. Only one of
  `cens_strata` or `cens_cov` can be specified. The default is `NULL`.

- cens_covs:

  A set of variable names used for modeling censoring time distribution.
  Only one of `cens_strata` or `cens_covs` can be specified. The default
  is `NULL`.

## Value

an object of class ahreg.

- result:

  A table containing the coefficient estimates, standard errors,
  confidence intervals, z-values, and two-sided p-values for each
  predicotor.

## Details

The function implements the average hazard regression.

## References

\#' Uno H, Tian L, Horiguchi M, Hattori S, Kehl KL. Regression models
for average hazard. Biometrics. 2024; 80(2):ujae037. \<doi:
10.1093/biomtc/ujae037\>

## Author

Hajime Uno, Miki Horiguchi

## Examples

``` r
#================================================================================
# ahreg.sample.data: Sample data from the pbc data in the survival package
#================================================================================
D = ahreg.sample.data()


#-- Independent censoring
a1 = ahreg(Surv(time,status) ~ arm + edema + bili, tau=7, data=D)
print(a1)
#> Call: ahreg() 
#> 
#> [1] "Link: log"
#> Surv(time, status) ~ arm + edema + bili
#> <environment: 0x10f458cd8>
#> 
#>                Est      SE low_0.95 upp_0.95         Z       p
#> Intercept -3.41306 0.20310 -3.81113 -3.01499 -16.80468 0.00000
#> arm        0.29691 0.21736 -0.12911  0.72293   1.36598 0.17194
#> edema      1.38853 0.35653  0.68974  2.08731   3.89457 0.00010
#> bili       0.11498 0.01630  0.08304  0.14692   7.05514 0.00000

#-- Group specific censoring
a2 = ahreg(Surv(time,status) ~ arm + edema + bili, tau=7, data=D, cens_strata="arm")
print(a2)
#> Call: ahreg() 
#> 
#> [1] "Link: log"
#> Surv(time, status) ~ arm + edema + bili
#> <environment: 0x10f458cd8>
#> 
#>                Est      SE low_0.95 upp_0.95         Z       p
#> Intercept -3.40649 0.20717 -3.81254 -3.00043 -16.44264 0.00000
#> arm        0.27832 0.22935 -0.17119  0.72783   1.21354 0.22492
#> edema      1.39166 0.35774  0.69050  2.09283   3.89011 0.00010
#> bili       0.11507 0.01634  0.08305  0.14709   7.04357 0.00000

#-- Covariate dependent censoring via Cox 
a3 = ahreg(Surv(time,status) ~ arm + edema + bili, tau=7, data=D, cens_covs=c("arm","edema"))
print(a3)
#> Call: ahreg() 
#> 
#> [1] "Link: log"
#> Surv(time, status) ~ arm + edema + bili
#> <environment: 0x10f458cd8>
#> 
#>                Est      SE low_0.95 upp_0.95         Z       p
#> Intercept -3.40686 0.21436 -3.82699 -2.98673 -15.89349 0.00000
#> arm        0.29992 0.23368 -0.15809  0.75793   1.28347 0.19933
#> edema      1.46953 0.39803  0.68941  2.24965   3.69203 0.00022
#> bili       0.11277 0.01738  0.07870  0.14684   6.48750 0.00000

#-- Covariate dependent censoring via Cox (identity link)
a4 = ahreg(Surv(time,status) ~ arm + edema + bili, tau=7, data=D, cens_covs=c("arm","edema"), 
           link="identity")
print(a4)
#> Call: ahreg() 
#> 
#> [1] "Link: identity"
#> Surv(time, status) ~ arm + edema + bili
#> <environment: 0x10f458cd8>
#> 
#>                Est      SE low_0.95 upp_0.95        Z       p
#> Intercept -0.00139 0.01212 -0.02514  0.02236 -0.11466 0.90872
#> arm        0.00442 0.01609 -0.02712  0.03597  0.27485 0.78343
#> edema      0.23775 0.11805  0.00638  0.46912  2.01402 0.04401
#> bili       0.02609 0.00475  0.01678  0.03541  5.49178 0.00000

```
