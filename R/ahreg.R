#' @name ahreg
#' @aliases ahreg
#' @title Regression Analysis with Average Hazard
#' @description The \code{ahreg} function performs a regression analysis for the average hazard (AH).
#' @author Hajime Uno, Miki Horiguchi
#' @references #' Uno H, Tian L, Horiguchi M, Hattori S, Kehl KL. Regression models for average hazard. Biometrics. 2024; 80(2):ujae037. <doi: 10.1093/biomtc/ujae037>
#' @usage ahreg(formula, tau, data, link="log", conf.int=0.95, 
#'              cens_strata=NULL, cens_covs=NULL)
#' @param formula A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the Surv function. For a multi-state model the formula may be a list of formulas.
#' @param tau A scalar value to specify a time point for calculating the average hazard.
#' @param data A data.frame in which to interpret the variables named in the formula.
#' @param link A link function to be used, either "log" (default) or "identity". 
#' @param conf.int A confidence coefficient for calculating confidence intervals. The default is \code{conf.int=0.95}.
#' @param cens_strata A variable name for specifying group-specific censoring. Only one of \code{cens_strata} or \code{cens_cov} can be specified. The default is \code{NULL}.
#' @param cens_covs A set of variable names used for modeling censoring time distribution. Only one of \code{cens_strata} or \code{cens_covs} can be specified. The default is \code{NULL}. 
#' @details The function implements the average hazard regression. 
#' @return an object of class ahreg.
#' @return \item{result}{A table containing the coefficient estimates, standard errors, confidence intervals, z-values, and two-sided p-values for each predicotor.}


#' @examples
#' #================================================================================
#' # ahreg.sample.data: Sample data from the pbc data in the survival package
#' #================================================================================
#' D = ahreg.sample.data()
#'
#'
#' #-- Independent censoring
#' a1 = ahreg(Surv(time,status) ~ arm + edema + bili, tau=7, data=D)
#' print(a1)
#'
#' #-- Group specific censoring
#' a2 = ahreg(Surv(time,status) ~ arm + edema + bili, tau=7, data=D, cens_strata="arm")
#' print(a2)
#' 
#' #-- Covariate dependent censoring via Cox 
#' a3 = ahreg(Surv(time,status) ~ arm + edema + bili, tau=7, data=D, cens_covs=c("arm","edema"))
#' print(a3)
#' 
#' #-- Covariate dependent censoring via Cox (identity link)
#' a4 = ahreg(Surv(time,status) ~ arm + edema + bili, tau=7, data=D, cens_covs=c("arm","edema"), 
#'            link="identity")
#' print(a4)
#' 
#' 
#'@export
ahreg <- function(formula, tau, data, link="log", conf.int=0.95, cens_strata=NULL, cens_covs=NULL){

  #---
  wt_ipcw = NULL
  
  #--- Initial check ---
  if(missing(formula)){
    stop("a formula argument is required")
  }
  
  if(is.null(tau)){
    stop("Please specify tau")
  }
  
  if(!is.null(cens_strata) & !is.null(cens_covs)){
    stop("It is not permitted to specify both cens_strata and cens_covs simultaneously. Only one of these parameters can be specified.")
  }

  #--- Formula  ---
  fmla     = as.formula(formula)
  terms    = terms(formula)
  response = model.response(model.frame(terms, data))
  
  if(!inherits(response, "Surv")){
    stop("The response must be a Surv object.")
  }

  #--- Data extraction (outcome model) ---
  time   =  response[, 1]
  status = response[, 2]
  
  if(length(all.vars(terms))==2){
    covariates = NULL
  }
  
  if(length(all.vars(terms))>2){
    covariates = model.matrix(~ . - 1, data = data[, all.vars(terms)[c(-1,-2)], drop=FALSE])
  }
  
  #--- Strata for censoring ---
  if(is.null(cens_strata)){
    strata        = rep(1, length(time))
    unique_strata = 1
    nstrata       = 1
  }
  
  if(!is.null(cens_strata)){
    strata        = data[,cens_strata]
    unique_strata = sort(unique(strata))
    nstrata       = length(unique_strata)
  }
  
  #--- Covariates for censoring ---
  if(is.null(cens_covs)){
    covariates4cens = NULL
  }
  
  if(!is.null(cens_covs)){
    covariates4cens = as.matrix(data[,cens_covs])
  }

  #--- Initial value for beta -- 
  ah_rs = sum(status*as.numeric(time<tau))/sum(pmin(time, tau))
  
  if(is.null(covariates)){
    b_ini = log(ah_rs)
  }else{
    b_ini = c(log(ah_rs), rep(0, ncol(covariates)))
  }

  design_mat = cbind(rep(1,length(time)), covariates) 
  
  #--- When the number of covariates is 1, this line is required ---
  if(is.null(colnames(design_mat))){
    colnames(design_mat) = paste0("X", 1:ncol(design_mat)-1)
  }
  
  if(!is.null(colnames(design_mat))){
    colnames(design_mat)[1] = "Intercept"
  }

  
  #=====================================================
  # IPCW (ref. t-year survivor)
  #=====================================================
  if(!is.null(wt_ipcw)){
    wt = wt_ipcw
  }else{

    #---------------------------------------------------
    #------ Independent or group-specific censoring ----
    #---------------------------------------------------
    if(is.null(covariates4cens)){
      wt        = rep(0,length(time))
      kmc_list  = list()
      psii_list = list()
      
      for(i in 1:nstrata){
        idx = strata==unique_strata[i]
        kmc = kmcens2(time[idx], status[idx], tau=tau)
        
        wt[idx]        = kmc$ipcw
        kmc_list[[i]]  = kmc
        psii_list[[i]] = kmc$psii
      }
    }

    #----------------------------------------------------
    #------ Cox censoring -------------------------------
    #----------------------------------------------------
    if(!is.null(covariates4cens)){
      censoring = 1-status
      
    	#---Ignore censoring at and after tau
    	censoring[time>=tau] = 0 
      
    	fc_cox = cox1(time, censoring, covariates4cens, tau=tau)
      wt     = fc_cox$ipcw
    }
    
  }
  
  
  #=====================================================
  # Estimate beta using the Newton-Raphson method
  #=====================================================
  hx            = c()
  ii            = 0
  qq            = 99999
  beta          = b_ini
  max_iteration = 20
  criterion     = 1E-15
  convergence   = 9
  
  
  if(link=="log"){
    while(ii < max_iteration & qq > criterion){
      ii = ii+1
      
      gbeta      = exp(design_mat %*% beta)
      yy         = (as.numeric(time < tau) - gbeta*pmin(time, tau))*wt
      yy_mat     = t(matrix(yy, ncol=length(yy), nrow=ncol(design_mat), byrow=T))
      S_beta     = apply(design_mat * yy_mat, 2, mean)
      zz         = gbeta*pmin(time, tau)*wt
      zz_mat     = t(matrix(zz, ncol=length(zz), nrow=ncol(design_mat), byrow=T))
      A          = -t(design_mat * zz_mat) %*% design_mat / nrow(design_mat)
      A_inv      = solve(A)
      difference = A_inv %*% S_beta
      qq         = sum(difference^2)
      hx         = rbind(hx, c(ii, qq, beta))
      beta       = beta - difference
      
      if(qq<=criterion){
        convergence = 0
      }
    }
  }
  
  if(link=="identity"){
    while(ii < max_iteration & qq > criterion){
      ii = ii+1
      
      gbeta      = design_mat %*% beta
      yy         = (as.numeric(time < tau) - gbeta*pmin(time, tau))*wt
      yy_mat     = t(matrix(yy, ncol=length(yy), nrow=ncol(design_mat), byrow=T))
      S_beta     = apply(design_mat * yy_mat,2,mean)
      zz         = pmin(time, tau)*wt
      zz_mat     = t(matrix(zz, ncol=length(zz), nrow=ncol(design_mat), byrow=T))
      A          = -t(design_mat * zz_mat) %*% design_mat / nrow(design_mat)
      A_inv      = solve(A)
      difference = A_inv %*% S_beta
      qq         = sum(difference^2)
      hx         = rbind(hx, c(ii, qq, beta))
      beta       = beta - difference
      
      if(qq<=criterion){
        convergence = 0
      }
    }
  }
  
  #--- Save convergence info ---
  convergence_information             = list()
  convergence_information$convergence = 0
  convergence_information$itration    = ii
  convergence_information$qq          = hx[,2]
  
  history                             = data.frame(hx)
  colnames(history)                   = c("Iteration", "Q", paste0("beta",1:(ncol(hx)-2)))
  convergence_information$history     = history

  
  #=====================================================
  # Asymptotic variance 
  #=====================================================
  if(is.null(wt_ipcw)){
      
    if(link=="log"){
      gbeta   = exp(design_mat %*% beta)
      yy      = (as.numeric(time < tau) - gbeta*pmin(time, tau))*wt
      yy_mat  = t(vtm(yy, ncol(design_mat)))
      S_beta  = apply(design_mat * yy_mat, 2, mean)
      zz      = gbeta*pmin(time, tau)*wt
    }
    
    if(link=="identity"){
      gbeta   = design_mat %*% beta
      yy      = (as.numeric(time < tau) - gbeta*pmin(time, tau))*wt
      yy_mat  = t(vtm(yy, ncol(design_mat)))
      S_beta  = apply(design_mat * yy_mat, 2, mean)
      zz      = pmin(time, tau)*wt
    }
  
    zz_mat  = t(vtm(zz, ncol(design_mat)))
    A       = -t(design_mat * zz_mat) %*% design_mat / nrow(design_mat)


  #------------------------------------------------------
  #------ Independent or group-specific censoring -------
  #------------------------------------------------------
  if(is.null(covariates4cens)){

    #--- K(beta,S) ---
    psii_ks = matrix(NA, nrow=length(time), ncol=ncol(design_mat))
    
    for(i in 1:nstrata){
      idx  = strata==unique_strata[i]
      psii = psii_list[[i]]

      #--- Distinct time points ---
      distinct = unique(sort(c(time[idx],tau)))
      
      wk1      = t(vtm(time[idx], length(distinct)))
      wk2      = vtm(distinct, sum(idx))
      wk3      = (wk1 <= wk2)
      
      K_mat    = matrix(NA, nrow=ncol(design_mat), ncol=length(distinct)) 
      dK_mat   = matrix(NA, nrow=ncol(design_mat), ncol=length(distinct)) 
    
      for(jj in 1:ncol(design_mat)){
        yy2     = (as.numeric(time[idx] < tau) - gbeta[idx]*pmin(time[idx], tau))*wt[idx]*design_mat[idx,jj]
        yy2_mat = t(vtm(yy2, ncol(wk3)))
        
        wk4        = wk3 * yy2_mat
        wk5        = apply(wk4, 2, mean)
        K_mat[jj,] = apply(wk3 * yy2_mat, 2, mean)
        
        wk6                 = K_mat[jj,]
        wk6[distinct > tau] = 0
        wk7                 = diff(c(0,wk6))
        dK_mat[jj,]         = wk7
      }
      
      psii_ks[idx,] = psii%*%t(dK_mat)
    }
  
    Ui = design_mat*(yy%*%t(rep(1,length(beta)))) - psii_ks
  }

  #------------------------------------------------------
  #------ Cox censoring ---------------------------------
  #------------------------------------------------------
  if(!is.null(covariates4cens)){
    
    #--- 1st term ---
    wk1a  = (as.numeric(time < tau) - gbeta*pmin(time, tau))*wt
    term1 = t(vtm(wk1a,ncol(design_mat)))*design_mat
    
    #--- 2nd term ---
    wk2a    = (as.numeric(time < tau) - gbeta*pmin(time, tau))*wt*fc_cox$Lam0.ebz
    wk2b    = t(vtm(wk2a,ncol(design_mat)))*design_mat
    K_gamma = t(wk2b) %*% covariates4cens / length(time)
    term2   = t(K_gamma %*% t(fc_cox$eta_beta_i))
    
    #--- 3rd term ---
      #--- K_lambda(S) ---
	    distinct = fc_cox$distinct
    	wk1      = t(vtm(time, length(distinct)))
    	wk2      = vtm(distinct, length(time))
    	wk3      = (wk1 >= wk2)*1
    	K_mat    = matrix(NA, nrow=ncol(design_mat), ncol=length(distinct)) 
    	dK_mat   = matrix(NA, nrow=ncol(design_mat), ncol=length(distinct)) 

      for(jj in 1:ncol(design_mat)){
        yy2      = (as.numeric(time < tau) - gbeta*pmin(time, tau))*wt*fc_cox$ebz*design_mat[,jj]
        yy2_mat  = t(vtm(yy2, ncol(wk3))) ; 
        wk4 = wk3 * yy2_mat
        dim(yy2_mat)
        dim(wk3)
        dim(wk4)
        wk5=apply(wk4,2,mean)
        K_mat[jj,] = wk5
        wk6 = K_mat[jj,]
        wk6[distinct > tau]=0
        wk7 = diff(c(0,wk6))
        dK_mat[jj,] = wk7
      }
    	
    term3 = fc_cox$eta_lam_i %*%t(dK_mat)
    	
    #--- Sum of i.i.d.  ---
    Ui = term1 + term2 + term3
  } 
    
    
  #=====================================================
  # Results
  #=====================================================
  B       = (t(Ui) %*% Ui)/nrow(design_mat)
  V       = solve(A) %*% B %*% solve(A)
  se      = sqrt(diag(V)/nrow(design_mat))
  low     = beta - se*abs(qnorm((1-conf.int)/2))
  upp     = beta + se*abs(qnorm((1-conf.int)/2))
  zstat   = beta/se
  pval    = (1-pnorm(abs(zstat)))*2

  result           = data.frame(cbind(beta, se, low, upp, zstat, pval))
  colnames(result) = c("Est", "SE", paste0("low_",conf.int), paste0("upp_",conf.int), "Z", "p")
  
  #--- Save variance info ---
  variance_information   = list()
  variance_information$A = A
  variance_information$B = B
  variance_information$V = V
  variance_information$n = nrow(design_mat)
  
}
  
  #=====================================================
  # Output 
  #=====================================================
  Z = list()
  
  if(is.null(wt_ipcw)){
    Z$result = result
  }else{
    Z$result = beta
  }

  class(Z)          = "ahreg"
  attr(Z,"formula") = fmla
  attr(Z,"beta")    = beta
  
  if(is.null(wt_ipcw)){
    attr(Z,"beta.var")             = V/nrow(design_mat)
    attr(Z,"Ainv")                 = solve(A)
    attr(Z,"B")                    = B
    attr(Z,"predicted")            = gbeta
    attr(Z,"variance_information") = variance_information
    
    if(is.null(covariates4cens)){
      attr(Z,"ipcw_information") = kmc_list   
      attr(Z,"nstrata")          = nstrata
    }
    
    if(!is.null(covariates4cens)){
      attr(Z,"ipcw_information") = fc_cox
      attr(Z,"nstrata")          = NA
    }
  }
  
  attr(Z,"convergence_information") = convergence_information
  attr(Z,"time")                    = time
  attr(Z,"status")                  = status
  attr(Z,"design_mat")              = design_mat
  attr(Z,"link")                    = link
  attr(Z,"tau")                     = tau
  attr(Z,"ipcw")                    = wt
  attr(Z,"cens_strata")             = cens_strata
  attr(Z,"cens_covs")               = cens_covs
  attr(Z,"data")                    = data
  
  return(Z)
}

