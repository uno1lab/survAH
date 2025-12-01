#============================================================
# cox1: Cox Regression Models (To get ipcw via Cox)
#============================================================

cox1 = function(time, status, covariates, tau=NULL){

  #--- Number of covaraites ---
  ncov = ncol(as.matrix(covariates))

  #--- Fit Cox and get beta ---
  tmpD = data.frame(cbind(time, status, covariates))
  vars = colnames(tmpD)[c(-1,-2)]
  fmla = as.formula(paste("Surv(time, status)~", paste(vars, collapse="+")))
  ft   = coxph(fmla, data=tmpD)
  beta = summary(ft)$coef[,1]
  
  if(ncov==1){ 
    bz = as.vector(unlist(covariates)) * rep(beta,length(time))
  }else{
    bz = covariates%*%beta
  }    
  
  ebz  = exp(bz)
  
  #--- Initial ---
  if(is.null(tau)){
    tau = max(time)
  }
  
  distinct = unique(sort(c(0, time, tau)))
  t        = length(distinct)
  n        = length(time)

  #--- n x t matrix ---
  wk1   = vtm(distinct,n)
  wk2   = t(vtm(time, t))
  dNi   = as.matrix(wk2==wk1 & t(vtm(status, t))==1)
  Ni    = t(sapply(data.frame(t(dNi)), cumsum))
  Yi    = as.matrix(wk2>=wk1)*1
  Ybar  = apply(Yi, 2, sum)
  dNbar = apply(dNi, 2, sum)

  #--- S0bar, S1bar, S2bar ---
  Yiebz = Yi*t(vtm(ebz,t))
  S0bar = apply(Yiebz, 2, mean)
  
 
  #---------------------------
  if(ncov!=1){
    wk1   = data.frame(covariates)
    tmp1  = function(x){apply(Yiebz*t(vtm(x,t)), 2, mean)}
    S1bar = t(sapply(wk1, tmp1))

    wk2   = data.frame(t(covariates))
    tmp2  = function(x){as.vector(x%*%t(x))}
  
    wk3             = data.frame(t(sapply(wk2, tmp2)))
    S2bar           = t(sapply(wk3, tmp1))
    S2bar_per_S0bar = S2bar/vtm(S0bar, ncov*ncov)
    S1bar_per_S0bar = S1bar/vtm(S0bar, ncov)
  
    wk4                = data.frame(S1bar_per_S0bar)
    S1bar_per_S0bar_2  = sapply(wk4, tmp2)
    v_beta_t           = S2bar_per_S0bar - S1bar_per_S0bar_2
  
    wk5  = data.frame(t(v_beta_t))
    tmp5 = function(x){apply(vtm(x,n) * dNi, 1, sum)}
  
    wk6       = sapply(wk5, tmp5)
    Ibeta1    = matrix(apply(wk6, 2, mean), nrow=ncov, ncol=ncov)
    InvIbeta1 = solve(Ibeta1)
  }

  if(ncov==1){
    S1bar = apply(Yiebz * t(vtm(as.matrix(covariates),t)), 2, mean)
    S2bar = apply(Yiebz * t(vtm(as.matrix(covariates) * as.matrix(covariates),t)), 2, mean)
  
    S2bar_per_S0bar    = S2bar/S0bar
    S1bar_per_S0bar    = S1bar/S0bar
    S1bar_per_S0bar_2  = S1bar_per_S0bar * S1bar_per_S0bar
    v_beta_t           = S2bar_per_S0bar - S1bar_per_S0bar_2
    Ibeta1             = mean(dNi %*% v_beta_t) 
    InvIbeta1          = 1/Ibeta1
  }
  #---------------------------
  
     
  #--- Bleslow Estimator (t x 1) ---
  breslow = cumsum(apply(dNi/vtm(S0bar,n),2,mean))
  
  #--- Martingale (n x t) ---
  dAi = Yiebz * vtm(diff(c(0,breslow)), n)
  Ai  = cumsum(apply(dAi, 2, mean))
  Mi  = Ni - Ai
  dMi = dNi - dAi     

  #--- Lambda hat (sum of iid terms) ---
  wk1       = data.frame(t(vtm(1/S0bar, n)*dMi))
  wk2       = sapply(wk1, cumsum)
  eta_lam_i = t(wk2)

  
  #---------------------------
  if(ncov!=1){
    #--- Information matrix ---
    wk7      = v_beta_t * vtm(S0bar, ncov*ncov) * vtm(diff(c(0,breslow)), ncov*ncov)
    Ibeta    = matrix(apply(wk7, 1, sum), nrow=ncov, ncol=ncov)
    InvIbeta = solve(Ibeta)
    
    tmp2 = c()
    for(i in 1:ncov){
      tmp1 = apply((t(vtm(covariates[,i], t)) - vtm(S1bar_per_S0bar[i,], n)) * dMi, 1, sum)
      tmp2 = rbind(tmp2,tmp1)
    }
    
    #--- Beta hat (sum of i.i.d. terms) ---
    eta_beta_i = t(-solve(Ibeta)%*%tmp2)
  }
  
  if(ncov==1){
    #--- Information matrix ---
    wk7      = v_beta_t * S0bar * diff(c(0,breslow))
    Ibeta    = sum(wk7)
    InvIbeta = 1/Ibeta
    
    tmp1 = apply((t(vtm(as.matrix(covariates), t)) -  vtm(S1bar_per_S0bar, n)) * dMi, 1, sum)
    
    #--- Beta hat (sum of i.i.d. terms) ---
    eta_beta_i = -1/(Ibeta)*tmp1
  }
  #---------------------------

  
  #-------------------------------------
  #--- IPCW part (status=0 is event) ---
  #-------------------------------------
  event            = 1-status
  event[time>=tau] = 2
  
  wk1              = vtm(distinct, n)
  wk2              = t(vtm(pmin(time, tau), t))
  dNCi             = as.matrix(wk2==wk1)
  dNCi[dNCi==TRUE] = 1
  
  wk3            = vtm(breslow,n)
  Lam0           = apply(wk3*dNCi, 1, max)
  Ghat           = exp(-Lam0*ebz)
  ipcw           = rep(0, n)
  ipcw           = 1/Ghat
  ipcw[event==0] = 0
 
 #--------------
 #--- Output ---
 #--------------
 Z = list()
  
 Z$fit        = ft
 Z$distinct   = distinct
 Z$breslow    = breslow
 Z$ebz        = ebz 
 Z$InvIbeta   = InvIbeta1
 Z$ipcw       = ipcw
 Z$Lam0.ebz   = Lam0*ebz
 Z$ebz        = ebz
 Z$eta_lam_i  = eta_lam_i  
 Z$eta_beta_i = eta_beta_i
 
 class(Z) = "cox1"
 return(Z)
}
