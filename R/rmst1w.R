#=====================
# funcKM2w -- hidden
#=====================
funcKM2w <- function(time, status, tau, weight=NULL){

  if(is.null(weight)){
    weight = rep(1, length(time))
  }

  wk = data.frame(time=time, status=status, weight=weight)
  wk = wk[order(wk$time),]

  t_idx = unique(sort(c(wk$time, 0, tau)))

  #----
  Y = rep(0,length(t_idx)) #n.risk
  N = rep(0,length(t_idx)) #n.event
  C = rep(0,length(t_idx)) #n.censor
  S = rep(0,length(t_idx)) #surv
  H = rep(0,length(t_idx)) #forSE
  D = rep(0,length(t_idx)) #forSE
  E = rep(0,length(t_idx)) #SE

  #i=1
  Y[1] = sum(wk$weight)
  N[1] = 0
  C[1] = 0
  S[1] = 1
  H[1] = 0
  D[1] = 0
  E[1] = 0

  #i>=2
  for(i in 2:length(t_idx)){
    Y[i] = Y[i-1] - N[i-1] - C[i-1]
    if(sum(wk$time==t_idx[i] & wk$status==1)>0){
      N[i] = sum(wk[which(wk$time==t_idx[i] & wk$status==1),"weight"])
    }else{
      N[i] = 0
    }

    if(sum(wk$time==t_idx[i] & wk$status==0)>0){
      C[i] = sum(wk[which(wk$time==t_idx[i] & wk$status==0),"weight"])
    }else{
      C[i] = 0
    }

    if(Y[i]<0){
      Y[i] = 0
    }
    if(Y[i]==0){
      S[i] = S[i-1]
    }else{
      S[i] = S[i-1]*(1-(N[i]/Y[i]))
    }

    if(Y[i]*(Y[i]-N[i])==0){
      H[i] = 0
    }else{
      H[i] = N[i]/(Y[i]*(Y[i]-N[i]))
    }

    if(S[i]<0){
      S[i] = 0
    }

    D[i] = sum(H[2:i])

    E[i] = sqrt((S[i]**2)*D[i])

    if(is.na(S[i])){
      S[i] = 0
    }
    if(is.na(E[i])){
      E[i] = 0
    }
  }

  #---
  out = data.frame(t_idx=t_idx, n_risk=Y, n_event=N, n_censor=C, surv=S, SE=E)
  return(out)
}


#=============================
# rmst1w (one-arm) -- hidden
#=============================
rmst1w <- function(obj_km, tau, alpha=0.05, var.method="Greenwood"){
  #-- obj_KM -- object (dataframe) from funcKM2w()
  #-- alpha -- gives (1-alpha) confidence interval

  ft  = data.frame(time=obj_km$t_idx, surv=obj_km$surv, n.risk=obj_km$n_risk, n.event=obj_km$n_event)

  idx = ft$time<=tau

  wk.time    = sort(c(ft$time[idx],tau))
  wk.surv    = ft$surv[idx]
  wk.n.risk  = ft$n.risk[idx]
  wk.n.event = ft$n.event[idx]

  #--rmst
  time.diff = diff(c(0, wk.time))
  areas     = time.diff * c(1, wk.surv)
  rmst      = sum(areas)
  #rmst

  if(var.method=="Greenwood"){
    wk.var = ifelse((wk.n.risk-wk.n.event)==0, 0, wk.n.event /(wk.n.risk * (wk.n.risk - wk.n.event)))
  }
  if(var.method=="Aalen"){
    wk.var = ifelse(wk.n.risk==0, 0, wk.n.event /(wk.n.risk * wk.n.risk))
  }

  wk.var   = c(wk.var,0)
  rmst.var = sum(cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se  = sqrt(rmst.var)

  #--- check ---
  # print(ft, rmean=tau)

  #--- output for rmst ---
  out = matrix(0,2,4)
  out[1,] = c(rmst, rmst.se, rmst-qnorm(1-alpha/2)*rmst.se, rmst+qnorm(1-alpha/2)*rmst.se)
  out[2,] = c(tau-out[1,1], rmst.se, tau-out[1,4], tau-out[1,3])
  rownames(out) = c("RMST","RMTL")
  colnames(out) = c("Est.", "se", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))

  #-- average hazard
  df = 1-ft[ft$time==tau,]$surv
  ah = df/rmst

  #--- output
  Z = list()
  Z$result          = out
  Z$rmst            = out[1,]
  Z$rmtl            = out[2,]
  Z$rmst.var        = rmst.var
  Z$fit             = ft
  Z$tau             = tau
  Z$dist_func       = df
  Z$average_hazard  = ah

  class(Z) = "rmst1w"

  return(Z)
}
