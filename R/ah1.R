#====================================================
# ah1.R: Calculate Average hazard (one-arm) --hidden
#====================================================
# time: the follow-up time for right censored data
# status: the event indicator (1:event, 0:censor)
# tau: the truncation time point for calculating F, R, and AH
# conf.int: a confidence coefficient for calculating the confidence interval
#============================================
#============================================
# This function calculates average hazard with the survival weight and its variance.
# It also provides t-year event rate and RMST and their variances based on the
# *LOG-TRANSFORMATION* as intermediate products. Note that the CI for t-year event rate and
# that for RMST do not exactly match the ones from the surv2sampleComp or survRM2 pacakges.
#============================================
ah1 = function(time, status, tau, conf.int){

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

#--- D(t) = cummulative number of events at time t
#---dD(t) = number of events at time t
Dt = apply(Nt, 2, sum)
#dDt = diff(c(0,Dt))
dDt = apply(dNt, 2, sum)

#--- H(t): cumulative hazard function (NA estimator)
Yt_bar = apply(Yt, 2, sum)
Ht     = cumsum(dDt/Yt_bar)
dHt    = diff(c(0,Ht))

#--- S(t): survival function (KM estimator)
ft = funcKM2w(x, d, tau=tau)
St = ft$surv[ft$t_idx %in% timegrd]

#====================================
# S(tau) needs to be away from 0 (>0)
#====================================
if(sum(St%in%0)>0){
  #stop("Stop")
  need_stop = TRUE
}else{
  need_stop = FALSE
}


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

#--- integrand for W^F　---
F_tau     = 1-St[k]
int1_W_Ft = (1-F_tau)/F_tau/Gt

#--- integrand for W^R　---
RMST_tau  = RMSTt[k]
int2_W_Rt = (RMST_tau - RMSTt)/RMST_tau/Gt

#--- integration ---
wk_F = wk_R = wk_Q = rep(0, n)
for(i in 1:n){
  	wk_F[i] = dMt[i,] %*% (int1_W_Ft)
  	wk_R[i] = dMt[i,] %*% (int2_W_Rt)
  	wk_Q[i] = dMt[i,] %*% (int1_W_Ft + int2_W_Rt)
}

#--- variance estimates for W^F and W^R (log-transformed)
v_F = wk_F%*%wk_F/n
v_R = wk_R%*%wk_R/n

#--- variance estimates for W^Q (log-transformed)
v_Q = dHt %*% ((1/F_tau - RMSTt/RMST_tau)*(1/F_tau - RMSTt/RMST_tau)/Gt)

#--- variance estimates for W^U (NOT log-transformed)
v_U = dHt %*% ((1/RMST_tau - F_tau*RMSTt/RMST_tau^2)^2/Gt)
#===============
#--- output  ---
#===============
z_alpha = qnorm(1-(1-conf.int)/2)

out = c()

#--F(tau)--
m = F_tau; se = sqrt(v_F/n)
tmp = c(log(m), se, log(m)-z_alpha*se, log(m)+z_alpha*se, m, exp(log(m)-z_alpha*se), exp(log(m)+z_alpha*se)); out=rbind(out,tmp)

#--RMST(tau)--
m = RMST_tau; se = sqrt(v_R/n)
tmp = c(log(m), se, log(m)-z_alpha*se, log(m)+z_alpha*se, m, exp(log(m)-z_alpha*se), exp(log(m)+z_alpha*se)); out=rbind(out,tmp)

#--AH(tau)--
m   = F_tau / RMST_tau; se = sqrt(v_Q/n)
tmp = c(log(m), se, log(m)-z_alpha*se, log(m)+z_alpha*se, m, exp(log(m)-z_alpha*se), exp(log(m)+z_alpha*se)); out=rbind(out,tmp)

rownames(out)=c("F(tau)","RMST(tau)","AH(tau)")
colnames(out)=c("log(Est)", "SE(log(Est))", paste0(c("log(low_","log(upp_"), conf.int*100, ")"), "Est", paste0(c("low_","upp_"), conf.int*100))


#-----------
Z=list()
Z$need_stop = need_stop
#Z$time_grid   = timegrd
#Z$cum_haz_fun = Ht
#Z$surv_fun    = St
#Z$rmst_fun    = RMSTt
#Z$ah_fun      = (1-St)/RMSTt

#Z$v_F         = v_F
#Z$v_R         = v_R
Z$v_Q         = v_Q
Z$v_U         = v_U
Z$tau         = tau
Z$result1     = out

class(Z) = "ah1"
return(Z)
}
