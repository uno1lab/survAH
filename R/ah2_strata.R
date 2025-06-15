#========================================================
# ah2_strata (AH stratified analysis main code) (hidden)
#========================================================
ah2_strata = function(time, status, arm, strata, tau, conf.int, weights=NULL){ 

  #--- number of strata ---
  ustrata = sort(unique(strata))
  nstrata = length(ustrata)

  if(is.null(strata) | nstrata<=1){
       stop("There is only one stratum.")
  }

  #--- strata and weight --
  if(!is.null(strata) & nstrata!=1 & !is.null(weights)){
    if(nstrata!=length(weights)){
       stop("Give a weight for each stattum.")
    }else{
     	 print(data.frame(strata=ustrata, wt=weights/sum(weights)))}
  }
    
  #--- initialization --- 	
  F0tau = F1tau = rep(NA,nstrata) #-- event rate, F(tau) for each stratum and each group
  R0tau = R1tau = rep(NA,nstrata) #-- RMTS, R(tau) for each stratum and each group
  eta0tau = eta1tau = rep(NA,nstrata) #-- Average Hazard, eta(tau) for each stratum and each group
  n0 = n1 = nn = rep(NA,nstrata)  #-- Sample size for each stratum and each group, and total
  log_rah = log_rah_var = rep(NA,nstrata) #--- log(rah) and its variance
  dah = dah_var = rep(NA,nstrata) #-- dah and its variance

  for (i in 1:nstrata){
	  idx       = strata==ustrata[i]
	  st_time   = time[idx] 
	  st_status = status[idx] 
    st_arm    = arm[idx]
    n1[i]     = sum(st_arm==1)
    n0[i]     = sum(st_arm==0) 
    nn[i]     = n1[i] + n0[i]
    
    #--- run AH2 ---
    ft = ah2_core1(indat=data.frame(time=st_time, status=st_status, arm=st_arm), tau=tau, conf.int=conf.int)
    
    F0tau[i]   = attributes(ft)$F0tau
    F1tau[i]   = attributes(ft)$F1tau
    R0tau[i]   = attributes(ft)$R0tau
    R1tau[i]   = attributes(ft)$R1tau
    eta1tau[i] = attributes(ft)$eta1tau
    eta0tau[i] = attributes(ft)$eta0tau

    log_rah[i]     = attributes(ft)$log_RAH_obs
    log_rah_var[i] = attributes(ft)$log_RAH_var
    dah[i]         = attributes(ft)$DAH_obs
    dah_var[i]     = attributes(ft)$DAH_var
  }
    
  #-- RAH/DAH for each stratum
  out1 = data.frame(ustrata, nn, log_rah, log_rah_var, exp(log_rah), dah, dah_var) 

#-----------------------------------------------------------------
#---- Direct standardization 
#-----------------------------------------------------------------
    #--- the default weights -- use marginal sample size as the weight
    if(is.null(weights)){
      wt = nn/sum(nn)
    }else{
      wt = weights/sum(weights)
    }
    		
    F0_bar = F0tau %*% wt ; R0_bar = R0tau %*% wt ; eta0 = F0_bar/R0_bar
    F1_bar = F1tau %*% wt ; R1_bar = R1tau %*% wt ; eta1 = F1_bar/R1_bar
    rah_standerdized = eta1 / eta0
    dah_standerdized = eta1 - eta0

    #--------------------------------------------------------
    #--- variance calculations for eat0, eat1, dah, and rah
    #--------------------------------------------------------
    eta0tau_var     = eta1tau_var     = rep(NA,nstrata) #-- variance of eta0 and eta1
    log_eta0tau_var = log_eta1tau_var = rep(NA,nstrata) #-- variance of log(eta0) and log(eta1)
    
    for (i in 1:nstrata){
 	    idx0 = strata==ustrata[i] & arm==0 ; getvar0=ah1_var4strata(time[idx0], status[idx0], tau, F0_bar, R0_bar)
 	    idx1 = strata==ustrata[i] & arm==1 ; getvar1=ah1_var4strata(time[idx1], status[idx1], tau, F1_bar, R1_bar)
	    eta0tau_var[i]     = getvar0$var_ah
	    eta1tau_var[i]     = getvar1$var_ah
	    log_eta0tau_var[i] = getvar0$var_log_ah
	    log_eta1tau_var[i] = getvar1$var_log_ah
    }

    #--- F, R, eta by strata
    out2 = data.frame(ustrata, nn, n0, n1, F0tau, F1tau, R0tau, R1tau, eta1tau, eta0tau, eta0tau_var, eta1tau_var)
    
    eta1_var = eta1tau_var %*% ((wt^2)/(n1/sum(n1))) / sum(n1)
    eta0_var = eta0tau_var %*% ((wt^2)/(n0/sum(n0))) / sum(n0)

    log_eta1_var = log_eta1tau_var %*% ((wt^2)/(n1/sum(n1))) / sum(n1)
    log_eta0_var = log_eta0tau_var %*% ((wt^2)/(n0/sum(n0))) / sum(n0)

    #--- standardized summaries
    out3 = data.frame(F0_bar, R0_bar, F1_bar, R1_bar, 
                      eta0, eta1, eta0_var, eta1_var, 
                      rah_standerdized, dah_standerdized)
       
   #-----------------------------------------------------
   # Average Hazard for each arm
   #-----------------------------------------------------
    z_alpha = qnorm(1-(1-conf.int)/2)

    #--- based on the original scale
    low1 = eta1 - sqrt(eta1_var)* z_alpha
    upp1 = eta1 + sqrt(eta1_var)* z_alpha
    low0 = eta0 - sqrt(eta0_var)* z_alpha
    upp0 = eta0 + sqrt(eta0_var)* z_alpha
   
    #--- based on log(AH)
    low1a = exp(log(eta1) - sqrt(log_eta1_var)* z_alpha)
    upp1a = exp(log(eta1) + sqrt(log_eta1_var)* z_alpha)
    low0a = exp(log(eta0) - sqrt(log_eta0_var)* z_alpha)
    upp0a = exp(log(eta0) + sqrt(log_eta0_var)* z_alpha)

    result1     = matrix(NA, 2, 5)
    result1[1,] = c(eta0, low0, upp0, low0a, upp0a)
    result1[2,] = c(eta1, low1, upp1, low1a, upp1a)
    colnames(result1) = c("Est.", 
                          paste0("Lower ", conf.int," (orginal scale)"), 
                          paste0("Upper ", conf.int," (orginal scale)"),
                          paste0("Lower ", conf.int," (based on log scale)"), 
                          paste0("Upper ", conf.int," (based on log scale)"))
    rownames(result1) = c(paste0("AH (arm0)"), paste0("AH (arm1)"))

    #-----------
    #--- RAH ---
    #-----------
    log_rah_standerdized_se = sqrt(log_eta1_var + log_eta0_var)
    
    low_rah = exp(log(rah_standerdized) - log_rah_standerdized_se*z_alpha)
    upp_rah = exp(log(rah_standerdized) + log_rah_standerdized_se*z_alpha)

    pval2 = pnorm(-abs(log(rah_standerdized))/log_rah_standerdized_se)*2

    result2 = matrix(NA, 1, 4)
    result2[1,c(1:3)] = c(rah_standerdized, low_rah, upp_rah)
    result2[1,4]      = pval2
    colnames(result2) = c("Est.", paste0("Lower ", conf.int), paste0("Upper ", conf.int), "P-value")
    rownames(result2) = paste0("Ratio of AH (arm1/arm0)")
    out2RAH = result2

    #-----------
    #--- DAH ---
    #-----------
    dah_standerdized_se = sqrt(eta1_var + eta0_var)
    
    low_dah = dah_standerdized - dah_standerdized_se*z_alpha
    upp_dah = dah_standerdized + dah_standerdized_se*z_alpha
    pval2   = pnorm(-abs(dah_standerdized)/dah_standerdized_se)*2

    result2 = matrix(NA, 1, 4)
    result2[1,c(1:3)] = c(dah_standerdized, low_dah, upp_dah)
    result2[1,4]      = pval2
    colnames(result2) = c("Est.", paste0("Lower ", conf.int), paste0("Upper ", conf.int), "P-value")
    rownames(result2) = paste0("Difference of AH (arm1-arm0)")
    out2DAH = result2


#---------------------------------------------------------------------------------
#---- Conventional Stratified Analysis (using 1/variance weight) -- for reference
#---------------------------------------------------------------------------------

    #-----------
    #--- RAH ---
    #-----------
    wt = 1/log_rah_var 
    log_rah_integrated     = log_rah %*% wt /sum(wt)
    log_rah_integrated_var = t(wt)%*% diag(log_rah_var) %*% wt/sum(wt)/sum(wt)
    log_rah_integrated_se  = sqrt(log_rah_integrated_var)

    low_rah = exp(log_rah_integrated - log_rah_integrated_se*z_alpha)
    upp_rah = exp(log_rah_integrated + log_rah_integrated_se*z_alpha)
    pval2 = pnorm(-abs(log_rah_integrated)/log_rah_integrated_se)*2
    result3 = matrix(NA, 1, 4)
    result3[1,c(1:3)] = c(exp(log_rah_integrated), low_rah, upp_rah)
    result3[1,4]      = pval2
    colnames(result3) = c("Est.", paste0("Lower ", conf.int), paste0("Upper ", conf.int), "P-value")
    rownames(result3) = paste0("Ratio of AH (arm1/arm0)")
    out3RAH = result3

    #-----------
    #--- DAH ---
    #-----------
    wt = 1/dah_var 
    dah_integrated     = dah %*% wt /sum(wt)
    dah_integrated_var = t(wt)%*% diag(dah_var) %*% wt/sum(wt)/sum(wt)
    dah_integrated_se  = sqrt(dah_integrated_var)
    
    low_dah = dah_integrated - dah_integrated_se*z_alpha
    upp_dah = dah_integrated + dah_integrated_se*z_alpha
    pval2 = pnorm(-abs(dah_integrated)/dah_integrated_se)*2

    result3 = matrix(NA, 1, 4)
    result3[1,c(1:3)] = c(dah_integrated, low_dah, upp_dah)
    result3[1,4]      = pval2
    colnames(result3) = c("Est.", paste0("Lower ", conf.int), paste0("Upper ", conf.int), "P-value")
    rownames(result3) = paste0("Difference of AH (arm1-arm0)")
    out3DAH = result3


#--------------------------------------------
#-- Results of the unstratified analysis 
#--------------------------------------------
out_unstrat = ah2_core1(indat=data.frame(time=time, status=status, arm=arm), tau=tau, conf.int=conf.int)
    
#---------------
#---- output ---
#---------------
  Z=list()
  
  Z$n.obs = out_unstrat$n.obs
  Z$ah    = out_unstrat$ah
  Z$rah   = out_unstrat$rah
  Z$dah   = out_unstrat$dah

    tmp=data.frame(total=nn, arm0=n0, arm1=n1)
    tmp1=rbind(tmp, apply(tmp,2,sum))
    rownames(tmp1)=c(paste0("strata",1:nrow(tmp)), "total")
  Z$nn = tmp1 
  Z$conventional_rah = out3RAH #-- results with the conventional stratified analysis (RAH)
  Z$conventional_dah = out3DAH #-- results with the conventional stratified analysis (DAH)
  Z$stratified_ah    = result1 #--- results of AH for each arm (proposed analysis)
  Z$stratified_rah   = out2RAH #--- results of RAH (proposed analysis)
  Z$stratified_dah   = out2DAH #--- results of DAH (proposed analysis)
  return(Z)
}
