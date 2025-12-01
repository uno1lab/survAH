# Two-Sample Comparison of Average Hazard

The `ah2` function performs two-sample comparisons using the average
hazard (AH) as a summary measure of the survival time distribution. Two
kinds of between-group contrast metrics, the ratio of AH (RAH) and the
difference in AH (DAH), are calculated.

## Usage

``` r
ah2(time, status, arm, tau=NULL, conf.int=0.95, eta=0, strata=NULL)
```

## Arguments

- time:

  A numeric vector of follow-up times for right-censored data.

- status:

  A numeric vector indicating the event status; 1 = event occurred, 0 =
  right censored.

- arm:

  A binary vector indicating group assignment; elements should be either
  0 or 1. Typically, 0 = control group and 1 = treatment group.

- tau:

  A scalar specifying the end time point (`tau`) for calculating the
  average hazard. If `tau = NULL`, the default is the maximum time point
  at which the risk set size in both groups is at least 10.

- conf.int:

  A numeric value specifying the confidence level for confidence
  intervals. The default is `0.95`.

- eta:

  A scalar specifying the start time point (`eta`) for calculating the
  average hazard over the time window \[eta, tau\]. The default is `0`.

- strata:

  An optional numeric vector specifying a stratification factor for
  stratified analysis. If `strata = NULL` (default), an unstratified
  analysis is performed.

## Value

An object of class `ah2`, which contains the following components:

- note:

  A note indicating the time window used in the analysis, specified as
  \[eta, tau\].

- n.obs:

  A summary of the number of observations, including the total number,
  number of events by tau, number of censored observations by tau, and
  the size of the risk set at tau.

- ah:

  Estimated average hazard in each arm.

- rah:

  Ratio of average hazards (RAH), calculated as treatment over control.

- dah:

  Difference of average hazards (DAH), calculated as treatment minus
  control.

- conventional_rah:

  Ratio of average hazards based on the conventional stratified analysis
  method.

- conventional_dah:

  Difference of average hazards based on the conventional stratified
  analysis method.

- stratified_ah:

  Estimated average hazard by arm based on the proposed stratified
  analysis method.

- stratified_rah:

  Ratio of average hazards based on the proposed stratified analysis
  method.

- stratified_dah:

  Difference of average hazards based on the proposed stratified
  analysis method.

## Details

The function provides the AH for each of the two groups, the absolute
difference and the absolute ratio of AH (DAH and RAH) between the two
groups, and the corresponding confidence intervals. It also calculates
p-values for the two-sided tests based on the RAH and DAH.

## References

Uno H and Horiguchi M. Ratio and difference of average hazard with
survival weight: new measures to quantify survival benefit of new
therapy. Statistics in Medicine. 2023; 42(7):936-952.
\<doi:10.1002/sim.9651\> Horiguchi M, Tian L, Kehl KL, Uno H. Assessing
delayed treatment benefits of immunotherapy using long-term average
hazard: a novel test/estimation approach. Lifetime Data Anal. 2025;
31(4):784-809. \<doi:10.1007/s10985-025-09671-0\> Qian Z, Tian L,
Horiguchi M, Uno H. A novel stratified analysis method for testing and
estimating overall treatment effects on time-to-event outcomes using
average hazard with survival weight. Stat Med. 2025; 44(7):e70056.
\<doi:10.1002/sim.70056\>

## Author

Hajime Uno, Miki Horiguchi, Zihan Qian

## Examples

``` r
#====================================================================
# cm214_pfs: The sample reconstructed data of the CheckMate214 study.
#====================================================================
# The code below reproduces the results reported by
# Uno H and Horiguchi M (StatMed; 2023) in Table 6.
D      = cm214_pfs
time   = D$time
status = D$status
arm    = D$arm
tau    = 21

a = ah2(time=time, status=status, arm=arm, tau=tau, conf.int=0.95, eta=0)
print(a, digits=3)
#> 
#> The time window: [eta, tau] = [0, 21] was specified.
#> 
#> Number of observations: 
#>      Total N Event by tau Censor by tau At risk at tau
#> arm0     422          225           163             34
#> arm1     425          219           160             46
#> 
#> 
#> Average Hazard (AH) by arm: 
#>            Est. Lower 0.95 Upper 0.95
#> AH (arm0) 0.066      0.057      0.076
#> AH (arm1) 0.049      0.042      0.057
#> 
#> 
#> Between-group contrast: 
#>                                Est. Lower 0.95 Upper 0.95 P-value
#> Ratio of AH (arm1/arm0)       0.747      0.608      0.917   0.005
#> Difference of AH (arm1-arm0) -0.017     -0.029     -0.005   0.006


# The code below reproduces the results reported by
# Horiguchi M, Tian L, Kehl K.L, Uno H (arXiv; 2024) in Table 3.
b = ah2(time=time, status=status, arm=arm, tau=21, conf.int=0.95, eta=7)
print(b, digits=3)
#> 
#> The time window: [eta, tau] = [7, 21] was specified.
#> 
#> Number of observations: 
#>      Total N Event by tau Censor by tau At risk at tau
#> arm0     171           67            70             34
#> arm1     213           61           106             46
#> 
#> 
#> Average Hazard (AH) by arm: 
#>            Est. Lower 0.95 Upper 0.95
#> AH (arm0) 0.051      0.040      0.065
#> AH (arm1) 0.028      0.022      0.037
#> 
#> 
#> Between-group contrast: 
#>                                Est. Lower 0.95 Upper 0.95 P-value
#> Ratio of AH (arm1/arm0)       0.553      0.387      0.791   0.001
#> Difference of AH (arm1-arm0) -0.023     -0.037     -0.008   0.002


#====================================================================
# Stratified analysis example
#====================================================================
D      = myeloid
time   = D$futime/365.25
status = D$death
arm    = as.numeric(D$trt=="A")
tau    = 3
strata = as.numeric(D$flt3)

b = ah2(time=time, status=status, arm=arm, strata=strata, tau=tau)
print(b, digits=3)
#> 
#> The time window: [eta, tau] = [0, 3] was specified.
#> 
#> Number of observations: 
#>         total arm0 arm1
#> strata1   149   75   74
#> strata2   319  165  154
#> strata3   178   89   89
#> total     646  329  317
#> 
#> 
#>      Total N Event by tau Censor by tau At risk at tau
#> arm0     329          142            18            169
#> arm1     317          160            28            129
#> 
#> 
#> <Unstratified analysis> Average Hazard (AH) by arm: 
#>            Est. Lower 0.95 Upper 0.95
#> AH (arm0) 0.207      0.175      0.246
#> AH (arm1) 0.290      0.245      0.343
#> 
#> 
#> <Unstratified analysis> Between-group contrast: 
#>                               Est. Lower 0.95 Upper 0.95 P-value
#> Ratio of AH (arm1/arm0)      1.398      1.099      1.777   0.006
#> Difference of AH (arm1-arm0) 0.082      0.022      0.143   0.007
#> 
#> 
#> <Stratified analysis> Average Hazard (AH) by arm: 
#>            Est. Lower 0.95 (orginal scale) Upper 0.95 (orginal scale)
#> AH (arm0) 0.207                      0.170                      0.243
#> AH (arm1) 0.286                      0.235                      0.337
#>           Lower 0.95 (based on log scale) Upper 0.95 (based on log scale)
#> AH (arm0)                           0.173                           0.247
#> AH (arm1)                           0.239                           0.342
#> 
#> 
#> <Stratified analysis> Between-group contrast: 
#>                               Est. Lower 0.95 Upper 0.95 P-value
#> Ratio of AH (arm1/arm0)      1.383      1.076      1.778   0.011
#> Difference of AH (arm1-arm0) 0.079      0.016      0.142   0.013
 
```
