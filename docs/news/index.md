# Changelog

## survAH version 1.2.0 (2025-12-01)

#### New features

- The *ahreg()* function has been introduced to enable regression
  analysis for the average hazard.
- *print.ahreg()* method has been added to provide formatted output from
  *ahreg()*.

------------------------------------------------------------------------

## survAH version 1.1.2 (2025-06-15)

#### Changes

- The *ah1()* and *ah2()* functions check the input data to make sure
  time does not involve 0 or negative numbers.
- An example for stratified analysis is added for *ah2()*.
- *print.ah2()* is updated to produce a nicer output.

------------------------------------------------------------------------

## survAH version 1.1.1 (2025-04-30)

#### Changes

- The *ah1()* function is now provide the variance for RMST and
  cumulative incidence probability as well as averaage hazard.

------------------------------------------------------------------------

## survAH version 1.1.0 (2025-04-04)

#### New features

- Added the *eta* argument to *ah2()* to allow calculation of the
  average hazard within a specified time window \[eta, tau\].
- Added the *strata* argument to *ah2()* to enable stratified analysis.

#### Changes

- The *ah1()* function is now exported and available for estimating the
  average hazard in a single sample (previously a hidden function).

------------------------------------------------------------------------

## survAH version 1.0.0 (2023-01-16)

CRAN release: 2023-01-17

- Initial release.
