#=================================================================
# ah1_var4strata.R: Calculate Variance of Average hazard (hidden)
#=================================================================
# time: the follow-up time for right censored data
# status: the event indicator (1:event, 0:censor)
# tau: the truncation time point for calculating F, R, and AH
# F_tau_bar: Weighted average of Event rate, F(tau), across strata
# R_tau_bar: Weighted average of RMST, R(tau), across strata 
#=================================================================
#============================================
# This function calculates the variance of average hazard with the survival weight
# using the given F_tau_bar and R_tau_bar 
#============================================
ah1_var4strata = function(time, status, tau, F_tau_bar, R_tau_bar){

#--- initialize ---
n     = length(time)
indx  = order(time)
x     = time[indx]
d     = status[indx]

#-- distinct failure time (up to tau)
timegrd = unique(sort(c(0,(x[d==1 & x<=tau]), tau)))
k       = length(timegrd)

#-- working matrix to create counting process, at-risk process, etc.
wkN = t(matrix(x, ncol=length(x), nrow=k, byrow=T))
wkD = t(matrix(d, ncol=length(d), nrow=k, byrow=T))
wkT = matrix(timegrd, ncol=length(timegrd), nrow=n, byrow=T)

#-- Ni(t) = I(Xi<=t, delta_i=1) -- counting process for each
Nt  = (wkN <= wkT & wkD == 1)*1
#-- dNi(t) = I(Xi==t, delta_i=1)
dNt = (wkN == wkT & wkD == 1)*1

#-- Yi(t) = I(Xi>=t)  -- at risk process for each
Yt  = (wkN >= wkT)*1

#--- D(t) = cumulative number of events at time t
#---dD(t) = number of events at time t
Dt = apply(Nt, 2, sum)
#dDt = diff(c(0,Dt))
dDt = apply(dNt, 2, sum)

#--- H(t): cumulative hazard function (NA estimator)
Yt_bar = apply(Yt, 2, sum)
tmp = dDt/Yt_bar; tmp[Yt_bar==0]=0 ; 
Ht = cumsum(tmp)
dHt    = diff(c(0,Ht))

#--- S(t): survival function (KM estimator)
ft = funcKM2w(x, d, tau=tau)
St = ft$surv[ft$t_idx %in% timegrd]

#--- RMST(t) = RMST function
RMSTt = rep(0,k)
for(j in 2:k){
  	wk_time  = timegrd[1:j]
  	wk_surv  = St[1:(j-1)]
  	RMSTt[j] = diff(wk_time)%*%wk_surv
}

#--- G(t)=Pr(X>=t) --
Gt = apply(Yt, 2, mean)

#--- Mi(t) = Ni(t) - \int Yi(u)dHt(u) -- Martingle for each
#--- dMi(t) = dNi(t) - Yi(t)dHt(t) -- derivative
Mt = dMt = matrix(0, nrow=n, ncol=length(timegrd))
for(i in 1:n){
  	 Mt[i,] =  Nt[i,] - cumsum(Yt[i,] * dHt)
   	dMt[i,] = diff(c(0, Mt[i,]))
}

#-- make sure the sum is zero
chk = apply(Mt,2,sum);chk

#-----
F_tau_bar_vec = rep(F_tau_bar, length(RMSTt))
R_tau_bar_vec = rep(R_tau_bar, length(RMSTt))

#--- variance estimate of AH 
tmp   = dHt*((1/R_tau_bar_vec - F_tau_bar_vec*RMSTt/R_tau_bar_vec^2)^2/Gt)
v_Ujk = sum(tmp[dHt!=0])

tmp    = dHt*((1/F_tau_bar_vec - RMSTt/R_tau_bar_vec)^2/Gt)
v_Xijk = sum(tmp[dHt!=0])

#----------------
#--- output ----
#----------------
Z = data.frame(var_ah=v_Ujk, var_log_ah=v_Xijk)

return(Z)
}