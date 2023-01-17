#==========================
# core_process.R  --hidden
#==========================
core_process <- function(time_idx, status_idx, time_length, tau_start=tau_start, tau=tau, method="KM"){

#-- distinct failure time (up to tau)
timegrd = unique(sort(c(0,(time_idx[status_idx==1 & time_idx<=tau]), tau_start, tau)))
k = length(timegrd)
n = time_length

#-- working matrix to create counting process, at-risk process, etc.
wkN = t(matrix(time_idx, ncol=length(time_idx), nrow=k, byrow=T))
wkD = t(matrix(status_idx, ncol=length(status_idx), nrow=k, byrow=T))
wkT = matrix(timegrd, ncol=length(timegrd), nrow=n, byrow=T)

#-- Ni(t) = I(Xi<=t, delta_i=1) -- counting process for each
Nt = (wkN <= wkT & wkD == 1)*1
#-- dNi(t) = I(Xi==t, delta_i=1)
dNt = (wkN == wkT & wkD == 1)*1

#-- Yi(t) = I(Xi>=t)  -- at risk process for each
Yt = (wkN >= wkT)*1

#--- D(t) = cummulative number of events at time t
#---dD(t) = number of events at time t
Dt = apply(Nt, 2, sum)
#dDt = diff(c(0,Dt))
dDt = apply(dNt, 2, sum)

#--- H(t)= cumulative hazard function (NA estimator)
Yt_bar = apply(Yt, 2, sum)
Ht = cumsum(dDt/Yt_bar)
dHt = diff(c(0,Ht))

#--- S(t) = survival function (NA estimator or KM estimator)
if(method=="NA"){
  St = exp(-Ht)
}
if(method=="KM"){
  ft = funcKM2w(time_idx, status_idx, tau=c(tau_start, tau))
  St = ft$surv[ft$t_idx %in% timegrd]
}

#--- RMST(t) = RMST function (depend on the method to estimate St)
RMSTt=rep(0,k)
for(j in 2:k){
  wk_time=timegrd[1:j]
  wk_surv=St[1:(j-1)]
  RMSTt[j]=diff(wk_time)%*%wk_surv
}

#--- G(t)=Pr(X>=t) --
Gt = apply(Yt, 2, mean)

#--- Mi(t) = Ni(t) - \int Yi(u)dHt(u) -- Martingle for each
#--- dMi(t) = dNi(t) - Yi(t)dHt(t) -- derivative
Mt=dMt=matrix(0, nrow=n, ncol=length(timegrd))
for(i in 1:n){
  Mt[i,] =  Nt[i,] - cumsum(Yt[i,] * dHt)
  dMt[i,] = diff(c(0, Mt[i,]))
}
#-- make sure the sum is zero
chk=apply(Mt,2,sum)

#--output
zout = list()
zout$timegrd = timegrd
zout$k       = length(timegrd)
zout$Nt      = Nt
zout$dNt     = dNt
zout$Yt      = Yt
zout$Dt      = Dt
zout$dDt     = dDt
zout$Ht      = Ht
zout$dHt     = dHt
zout$St      = St
zout$RMSTt   = RMSTt
zout$Gt      = Gt
zout$Mt      = Mt
zout$dMt     = dMt
zout$check   = round(chk)

class(zout) = "core_process"
return(zout)
}
