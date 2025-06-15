#' @name print.ah2
#' @aliases print.ah2
#' @title print.ah2
#' @description S3 method for class 'ah2'
#' @usage \method{print}{ah2}(x, digits=3, ...)
#' @param x Object to be printed.
#' @param digits Integer indicating the number of decimal places.
#' @param ... Further arguments ignored in this function.
#' @return returns summary output for class 'ah2'

#' @export
#' 
######################################
# print.ah2 (hidden)
######################################
print.ah2=function(x, digits=3,...){

  cat("\n")

  cat(x$note)

  cat("\n\n")

  #----------------------------------------------
  #-- when strata is not given (standard output)
  #----------------------------------------------
  if(is.null(x$stratified_ah)){
    
    cat ("Number of observations: \n")
    
    prmatrix(x$n.obs)
    
    cat("\n\n")
    
    cat ("Average Hazard (AH) by arm: \n")

    prmatrix(round(x$ah , digits=digits))

    cat("\n\n")

    cat ("Between-group contrast: \n")

    prmatrix(round(rbind(x$rah, x$dah), digits=digits))
  
  }
  
  #----------------------------------------------
  #-- strata is given
  #----------------------------------------------
  if(!is.null(x$stratified_ah)){
    cat ("Number of observations: \n")
    
    prmatrix(x$nn)
    
    cat("\n\n")

    prmatrix(x$n.obs)

    cat("\n\n")
    
    cat ("<Unstratified analysis> Average Hazard (AH) by arm: \n")
    
    prmatrix(round(x$ah , digits=digits))
    
    cat("\n\n")
    
    cat ("<Unstratified analysis> Between-group contrast: \n")
    
    prmatrix(round(rbind(x$rah, x$dah), digits=digits))
    #-----------------------------------------
    cat("\n\n")
    
    cat ("<Stratified analysis> Average Hazard (AH) by arm: \n")
    
    prmatrix(round(x$stratified_ah , digits=digits))
    
    cat("\n\n")
    
    cat ("<Stratified analysis> Between-group contrast: \n")
    
    prmatrix(round(rbind(x$stratified_rah, x$stratified_dah), digits=digits))
  }
  
  invisible(x)

}
