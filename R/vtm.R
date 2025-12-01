#-------------------------------------
# Vector to Matrix (same as VTM())
#-------------------------------------
vtm<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}
