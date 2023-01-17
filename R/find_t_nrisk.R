#===============================================================================================================================================
# find_t_nrisk.R  --hidden
# The function calculates The maximum time point
# where the size of risk set for both groups remains at least \code{n.risk}.
# n.risk: number of at risk
#===============================================================================================================================================
find_t_nrisk <- function(time, status, arm, n.risk){
  D = data.frame(time=time, status=status, arm=arm)

  #-- Count number at risk
  timepoint = unique(sort(D$time))
  t.idx     = c(0, timepoint)
  rs1       = matrix(0, nrow=3, ncol=length(t.idx))

  for (r in 1:length(t.idx)){
    rs1[1,r] = t.idx[r]
    rs1[2,r] = sum(D$time[D$arm==0]>=t.idx[r]) # for arm0
    rs1[3,r] = sum(D$time[D$arm==1]>=t.idx[r]) # for arm1
  }

  rownames(rs1) = c("time", "arm0", "arm1")
  rs1 = data.frame(t(rs1))

  #--
  t_nr = min(max(rs1[rs1$arm0>=n.risk,]$time), max(rs1[rs1$arm1>=n.risk,]$time))

  #--
  zz=list()
  zz$time_at_nrisk = t_nr
  return(zz)
}
