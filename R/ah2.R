#' @name ah2
#' @aliases ah2
#' @title Two-Sample Comparison of Average Hazard
#' @description The \code{ah2} function performs two-sample comparisons using the average hazard (AH) as a summary measure of the survival time distribution.
#' Two kinds of between-group contrast metrics, the ratio of AH (RAH) and the difference in AH (DAH), are calculated.
#' @details The function provides the AH for each of the two groups, the absolute difference and the absolute ratio of AH (DAH and RAH) between the two groups, and the corresponding confidence intervals.
#' It also calculates p-values for the two-sided tests based on the RAH and DAH.
#' @author Hajime Uno, Miki Horiguchi, Zihan Qian
#' @references
#' Uno H and Horiguchi M. Ratio and difference of average hazard with survival weight: new measures to quantify survival benefit of new therapy. Statistics in Medicine. 2023;1-17. <doi:10.1002/sim.9651>

#' Horiguchi M, Tian L, Kehl K.L., and Uno H. Assessing delayed treatment benefits of immunotherapy using long-term average hazard: a novel test/estimation approach. <arXiv:2403.10742>

#' Qian Z, Tian L, Horiguchi M, and Uno H. A novel stratified analysis method for testing and estimating overall treatment effects on time-to-event outcomes using average hazard with survival weight. Statistics in Medicine. 2025. <arXiv:2404.00788>
#' @usage ah2(time, status, arm, tau=NULL, conf.int=0.95, eta=0, strata=NULL)
#' @param time A numeric vector of follow-up times for right-censored data.
#' @param status A numeric vector indicating the event status; 1 = event occurred, 0 = right censored.
#' @param arm A binary vector indicating group assignment; elements should be either 0 or 1. Typically, 0 = control group and 1 = treatment group.
#' @param tau A scalar specifying the end time point (\code{tau}) for calculating the average hazard. If \code{tau = NULL}, the default is the maximum time point at which the risk set size in both groups is at least 10.
#' @param conf.int A numeric value specifying the confidence level for confidence intervals. The default is \code{0.95}.
#' @param eta A scalar specifying the start time point (\code{eta}) for calculating the average hazard over the time window [eta, tau]. The default is \code{0}.
#' @param strata An optional numeric vector specifying a stratification factor for stratified analysis. If \code{strata = NULL} (default), an unstratified analysis is performed.
#' @return An object of class \code{ah2}, which contains the following components:
#' @return \item{note}{A note indicating the time window used in the analysis, specified as [eta, tau].}
#' @return \item{n.obs}{A summary of the number of observations, including the total number, number of events by tau, number of censored observations by tau, and the size of the risk set at tau.}
#' @return \item{ah}{Estimated average hazard in each arm.}
#' @return \item{rah}{Ratio of average hazards (RAH), calculated as treatment over control.}
#' @return \item{dah}{Difference of average hazards (DAH), calculated as treatment minus control.}
#' @return \item{conventional_rah}{Ratio of average hazards based on the conventional stratified analysis method.}
#' @return \item{conventional_dah}{Difference of average hazards based on the conventional stratified analysis method.}
#' @return \item{stratified_ah}{Estimated average hazard by arm based on the proposed stratified analysis method.}
#' @return \item{stratified_rah}{Ratio of average hazards based on the proposed stratified analysis method.}
#' @return \item{stratified_dah}{Difference of average hazards based on the proposed stratified analysis method.}


#' @examples
#' #====================================================================
#' # cm214_pfs: The sample reconstructed data of the CheckMate214 study.
#' #====================================================================
#' # The code below reproduces the results reported by
#' # Uno H and Horiguchi M (StatMed; 2023) in Table 6.
#' D      = cm214_pfs
#' time   = D$time
#' status = D$status
#' arm    = D$arm
#' tau    = 21
#'
#' a = ah2(time=time, status=status, arm=arm, tau=tau, conf.int=0.95, eta=0)
#' print(a, digits=3)
#' 
#' 
#' # The code below reproduces the results reported by
#' # Horiguchi M, Tian L, Kehl K.L, Uno H (arXiv; 2024) in Table 3.
#' b = ah2(time=time, status=status, arm=arm, tau=21, conf.int=0.95, eta=7)
#' print(b, digits=3)

#'@export
ah2 <- function(time, status, arm, tau=NULL, conf.int=0.95, eta=0, strata=NULL){
  #---------
  #-- tau --
  #---------
  #--default_tau: time point where each arm has at least 10 at risk
  wk_rs       = find_t_nrisk(time=time, status=status, arm=arm, n.risk=10)
  default_tau = wk_rs$time_at_nrisk
  
  #---
  NOTE = NULL
  if(is.null(tau)){
    tau = default_tau
    NOTE = paste0("The time window: [eta, tau] = [", eta, ", ", tau, "] was specified.")
  }
  if(!is.null(tau)){
    if(tau > default_tau){
      tau  = tau
      NOTE = paste0("The time window: [eta, tau] = [", eta, ", ", tau, "] was specified. Warning: The normal approximation may be questionable with the specified tau. A smaller value of tau would be recommended for this data.")
    }else{
      tau  = tau
      NOTE = paste0("The time window: [eta, tau] = [", eta, ", ", tau, "] was specified.")
    }
  }
  
  #---tau based on the original data
  tau_org = tau
  
  #---
  if(eta==0){
    #---------------------------
    #--data for AH [0, tau]
    #---------------------------
    if(is.null(strata)){
      indat = data.frame(time=time, status=status, arm=arm)
    }else{
      indat = data.frame(time=time, status=status, arm=arm, strata=strata)
    }
    
    tau   = tau_org
  }else{
    #------------------------------------
    #---data for long-term AH [eta, tau]
    #------------------------------------
    if(is.null(strata)){
      tmp   = data.frame(time=time, status=status, arm=arm)
      tmp2  = tmp[tmp$time-eta>0,]
      indat = data.frame(time=tmp2$time-eta, status=tmp2$status, arm=tmp2$arm)
    }else{
      tmp   = data.frame(time=time, status=status, arm=arm, strata=strata)
      tmp2  = tmp[tmp$time-eta>0,]
      indat = data.frame(time=tmp2$time-eta, status=tmp2$status, arm=tmp2$arm, strata=tmp2$strata)
    }
    
    tau   = tau_org - eta
  }
  
  if(is.null(strata)){
    Z2 = ah2_core1(indat=indat, tau=tau, conf.int=conf.int)
  }else{
    #---If !is.null(strata), perform a stratified analysis
    Z2 = ah2_strata(time=indat$time, status=indat$status, arm=indat$arm, strata=indat$strata, tau=tau, conf.int=conf.int) 
  }
  
  #----
  class(Z2) = "ah2"
  
  #----
  Z2$note = NOTE
  return(Z2)
}

