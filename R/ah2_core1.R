#=========================
# ah2_core1.R (hidden)
#=========================

ah2_core1 <- function(indat, tau, conf.int){
  #-------------------------------------------------------
  # AH (average hazard), F (t-year event rate), R (rmst) - observed
  #-------------------------------------------------------
  #--Arm0:
  D0 = indat[indat$arm==0,]
  n0 = nrow(D0)
  
  wk0 = ah1(time=D0$time, status=D0$status, tau=tau, conf.int=conf.int)
  
  #--Check S_{0}(tau)>0
  if(wk0$need_stop==TRUE){
    stop("The survival probability at the specified tau for arm 0 needs to be >0. Please specify another tau.")
  }
  
  #--
  df0_obs    = wk0$result1["F(tau)", "Est"]
  mu0_obs    = wk0$result1["RMST(tau)", "Est"]
  ah0_obs    = wk0$result1["AH(tau)", "Est"]
  logdf0_obs = log(df0_obs)
  logmu0_obs = log(mu0_obs)
  logah0_obs = log(ah0_obs)
  
  
  #--Arm1:
  D1 = indat[indat$arm==1,]
  n1 = nrow(D1)
  
  wk1 = ah1(time=D1$time, status=D1$status, tau=tau, conf.int=conf.int)
  
  #--Check S_{1}(tau)>0
  if(wk1$need_stop==TRUE){
    stop("The survival probability at the specified tau for arm 1 needs to be >0. Please specify another tau.")
  }
  
  #--
  df1_obs    = wk1$result1["F(tau)", "Est"]
  mu1_obs    = wk1$result1["RMST(tau)", "Est"]
  ah1_obs    = wk1$result1["AH(tau)", "Est"]
  logdf1_obs = log(df1_obs)
  logmu1_obs = log(mu1_obs)
  logah1_obs = log(ah1_obs)
  
  
  #--Ratio (arm1/arm0)
  dfr_obs    = df1_obs/df0_obs
  logdfr_obs = log(dfr_obs)
  
  mur_obs    = mu1_obs/mu0_obs
  logmur_obs = log(mur_obs)
  
  #--RAH: Ratio of average hazards (arm1/arm0)
  rah_obs    = ah1_obs/ah0_obs
  lograh_obs = log(rah_obs)
  
  #--DAH: Difference of average hazards (arm1-arm0)
  dah_obs = ah1_obs-ah0_obs
  
  #-----------------------------------------------------
  # Average Hazard for each arm
  #-----------------------------------------------------
  arms    = c(0,1)
  label   = "AH"
  rowout  = "AH(tau)"
  
  #---
  out1     = matrix(NA, 2, 3)
  out1[1,] = wk0$result1[rowout, c("Est", paste0(c("low_","upp_"), conf.int*100))]
  out1[2,] = wk1$result1[rowout, c("Est", paste0(c("low_","upp_"), conf.int*100))]
  
  colnames(out1) = c("Est.", paste0("Lower ", conf.int), paste0("Upper ", conf.int))
  rownames(out1) = c(paste0(label, " (arm0)"), paste0(label, " (arm1)"))
  
  #-----------------------------------------------------
  # Ratio of AH (RAH)
  #     - delta method (using analytic variance inside)
  #-----------------------------------------------------
  V1      = attributes(wk1)$v_Q
  V0      = attributes(wk0)$v_Q
  log_obs = lograh_obs
  
  #--Ratio
  #--se_diff: se(log(XX1)-log(XX0))
  se_diff  = sqrt(V1/n1 + V0/n0)
  low_diff = exp(log_obs - se_diff*abs(qnorm((1-conf.int)/2)))
  upp_diff = exp(log_obs + se_diff*abs(qnorm((1-conf.int)/2)))
  
  pval2 = pnorm(-abs(log_obs)/se_diff)*2
  
  out2 = matrix(NA, 1, 4)
  out2[1,c(1:3)] = c(exp(log_obs), low_diff, upp_diff)
  out2[1,4]      = pval2
  colnames(out2) = c("Est.", paste0("Lower ", conf.int), paste0("Upper ", conf.int), "P-value")
  rownames(out2) = paste0("Ratio of ", label, " (arm1/arm0)")
  
  #---
  out1RAH = out1
  out2RAH = out2
  
  #-----------------------------------------------------
  # Difference of AH
  #     - delta method (using analytic variance inside)
  #-----------------------------------------------------
  V1  = attributes(wk1)$v_U
  V0  = attributes(wk0)$v_U
  obs = dah_obs
  
  #--Difference
  se  = sqrt(V1/n1 + V0/n0)
  low = obs - se*abs(qnorm((1-conf.int)/2))
  upp = obs + se*abs(qnorm((1-conf.int)/2))
  
  pval2 = pnorm(-abs(obs)/se)*2
  
  out2 = matrix(NA, 1, 4)
  out2[1,c(1:3)] = c(obs, low, upp)
  out2[1,4]      = pval2
  colnames(out2) = c("Est.", paste0("Lower ", conf.int), paste0("Upper ", conf.int), "P-value")
  rownames(out2) = paste0("Difference of ", label, " (arm1-arm0)")
  
  #---
  out2DAH = out2
  
  
  #-------------------------------------------
  # Table to summarize the number of patients
  #-------------------------------------------
  n_table = matrix(NA, nrow=2, ncol=4)
  
  #--number of total
  n_table[1,1] = nrow(D0)
  n_table[2,1] = nrow(D1)
  
  #--number of event by tau
  n_table[1,2] = sum(D0$status==1 & D0$time<tau)
  n_table[2,2] = sum(D1$status==1 & D1$time<tau)
  
  #--number of censoring by tau
  n_table[1,3] = sum(D0$status==0 & D0$time<tau)
  n_table[2,3] = sum(D1$status==0 & D1$time<tau)
  
  #--number at risk at tau
  n_table[1,4] = sum(D0$time>=tau)
  n_table[2,4] = sum(D1$time>=tau)
  
  #--labels
  colnames(n_table) = c("Total N", "Event by tau", "Censor by tau", "At risk at tau")
  rownames(n_table) = c("arm0", "arm1")
  
  
  #-------------------------------------------
  # Output
  #-------------------------------------------
  Z=list()
  Z$n.obs = n_table
  Z$ah    = out1RAH
  Z$rah   = out2RAH
  Z$dah   = out2DAH
  
  #----
  class(Z)="ah2_core1"
  
  attr(Z, "F0tau")       = df0_obs
  attr(Z, "F1tau")       = df1_obs
  attr(Z, "R0tau")       = mu0_obs
  attr(Z, "R1tau")       = mu1_obs
  attr(Z, "eta0tau")     = ah0_obs
  attr(Z, "eta1tau")     = ah1_obs
  attr(Z, "log_RAH_obs") = lograh_obs
  attr(Z, "log_RAH_var") = se_diff^2
  attr(Z, "DAH_obs")     = dah_obs
  attr(Z, "DAH_var")     = se^2
  
  return(Z)
}

