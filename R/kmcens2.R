#======================================================================================
# kmcens2: Kaplan-Meier (To get ipcw and get i.i.d. form for KM for censoring)
#======================================================================================

kmcens2 = function(time, status, tau=NULL){
  
  #---Ininitial ---
  if(is.null(tau)){
    tau = max(time)
  }
  
  distinct = unique(sort(c(time,tau)))
  t        = length(distinct)
  n        = length(time)
  surv     = rep(0,t)
  nel.wk   = rep(0,t)
  nelson   = rep(0,t)
  
  #--- i=1 --- 
  yi        = sum(as.numeric(time>=distinct[1]))
  di        = sum(as.numeric(time==distinct[1] & status==0))
  surv[1]   = 1 * (1-di/yi)
  nel.wk[1] = di/yi
  
  #--- i=2 to t ---
  for(i in 2:t){
    yi        = sum(as.numeric(time>=distinct[i]))
    di        = sum(as.numeric(time==distinct[i] & status==0))
    surv[i]   = surv[i-1]*(1-di/yi)
    nel.wk[i] = di/yi
  }
  
  #--- Shift ---
  surv[2:t]   = surv[1:(t-1)]
  surv[1]     = 1
  nel.wk[2:t] = nel.wk[1:(t-1)]
  nel.wk[1]   = 0
  
  #---
  nelson = cumsum(nel.wk)
  
  #---
  pi_0 = rep(0,t)
  pi_X = rep(0,t)
  pi_T = rep(0,t)
  Fn   = rep(0,t)
  
  for(i in 1:t){
    yi      = as.numeric(time >= distinct[i])
    ni      = as.numeric(time <= distinct[i] & status==1)
    pi_0[i] = mean(yi)
    pi_X[i] = mean(yi)*surv[i]
    pi_T[i] = mean(yi)/surv[i]
    Fn[i]   = mean(ni)
  }
  
  #---------------------------------
  #--- Mi(t): (n x t) Martingale ---
  #---------------------------------
  Mi  = matrix(0, n, t)
  wk1 = matrix(0, n, t)
  wk2 = matrix(0, n, t)
  wk3 = matrix(0, n, t)
  
  for(j in 1:t){ 
    wk1[,j] = as.numeric(time <= distinct[j] & status==0)
    wk2[,j] = as.numeric(time >= distinct[j])
  }
  
  for(k in 1:n){
    Mi[k,] = wk1[k,] - wk2[k,]*nelson
  }
  
  #-------------------------------
  #--- psii(t): (n x t matrix) ---
  #-------------------------------
  psii = matrix(0, n, t) 
  wk1  = matrix(0, n, t)
  wk2  = matrix(0, n, t)
  
  for(j in 1:t){ 
    wk1[,j] = as.numeric(time == distinct[j] & status==0)
    wk2[,j] = as.numeric(time >= distinct[j])
  }
  
  for(i in 1:n){ 
    psii[i,] = (cumsum(wk1[i,]/pi_0) - cumsum(nel.wk*wk2[i,]/pi_0))*surv
  }
  
  #-----------------
  #--- IPCW part ---
  #-----------------
  
  #--- ghat(tau) ---
  ghat.tau = surv[distinct==tau]
  
  #--- ghat.X --- 
  ghat.x = rep(0,n)
  
  for(i in 1:n){
    ghat.x[i] = surv[distinct==time[i]]
  }
  
  ghat.x[is.na(ghat.x)] = 0
  
  #--- pi(t) ---
  pit                        = rep(0,n)
  pit[time>=tau]             = ghat.tau
  pit[time<=tau & status==1] = ghat.x[time<=tau & status==1]
  
  #--- Indicator ---
  vit = as.numeric(time<=tau)*status + as.numeric(time>=tau)
  
  #--- IPCW ---
  ipcw         = rep(0,n)
  ipcw[vit==1] = 1/pit[vit==1]

  
  #--------------
  #--- Output ---
  #--------------
  Z = list()
  
  Z$surv     = surv
  Z$nelson   = nelson
  Z$distinct = distinct
  Z$pi_0     = pi_0
  Z$pi_X     = pi_X
  Z$pi_T     = pi_T
  Z$Mi       = Mi
  Z$psii     = psii
  Z$Fn       = Fn
  Z$ipcw     = ipcw
  
  class(Z) = "kmcens2"
  return(Z)
}

