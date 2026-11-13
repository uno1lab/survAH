# survAH version 1.2.1 (2026-11-07)

### Changes
* *ah1()* now sets *need_stop = TRUE*, and *ah2()* stops, whenever no subject remains at risk at *tau* in a group or stratum, whether the last observation is an event or a censoring.
* The error messages from *ah2()* name the arm and the stratum, state the reason, and report the largest observed time.
* *ahreg()* now defines how an observation recorded exactly at *tau* is handled. The method assumes continuous time, so that no event or censoring coincides with *tau*; such ties do occur, however, when times are recorded in days or months. An event exactly at *tau* is now counted as an event in [0, tau], following the convention that events precede censorings at tied times, with person-time *tau* and weight 1/G(tau); a censoring exactly at *tau* is event-free through *tau* with the same weight. This is the convention used by *ah1()* and *ah2()*, and it applies to all three censoring-adjustment approaches. Previously the three approaches did not handle this case consistently.
* *ahreg()* gains *maxiter* and *tol* arguments, reports non-convergence, uses a finite starting value for the identity link when no events occur before *tau*, and checks that *tau* does not exceed the largest observed time (overall and within censoring strata), that the risk set at *tau* is not too small, that the censoring weights are finite, and that the covariates are not collinear. It returns the new attributes *ipcw_max* and *n_risk_tau*, and the element *itration* of *convergence_information* is renamed *iterations* (the old name is kept for one release).
* The YAML header of the vignette was put on single lines so that it parses with Pandoc 3.13 (the CRAN check ERROR of version 1.2.0 on the Linux flavors).

---

# survAH version 1.2.0 (2025-12-01)

### New features
* The *ahreg()* function has been introduced to enable regression analysis for the average hazard.
* *print.ahreg()* method has been added to provide formatted output from *ahreg()*.


---

# survAH version 1.1.2 (2025-06-15)

### Changes
* The *ah1()* and *ah2()* functions check the input data to make sure time does not involve 0 or negative numbers. 
* An example for stratified analysis is added for *ah2()*.
* *print.ah2()* is updated to produce a nicer output.

---

# survAH version 1.1.1 (2025-04-30)

### Changes
* The *ah1()* function is now provide the variance for RMST and cumulative incidence probability as well as averaage hazard. 


---


# survAH version 1.1.0 (2025-04-04)

### New features
* Added the *eta* argument to *ah2()* to allow calculation of the average hazard within a specified time window [eta, tau].
* Added the *strata* argument to *ah2()* to enable stratified analysis.

### Changes
* The *ah1()* function is now exported and available for estimating the average hazard in a single sample (previously a hidden function). 


---

# survAH version 1.0.0 (2023-01-16)
* Initial release.
