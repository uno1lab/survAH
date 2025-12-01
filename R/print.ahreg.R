#' @name print.ahreg
#' @aliases print.ahreg
#' @title print.ahreg
#' @description S3 method for class 'ahreg'
#' @usage \method{print}{ahreg}(x, digits=5, ...)
#' @param x Object to be printed.
#' @param digits Integer indicating the number of decimal places.
#' @param ... Further arguments ignored in this function.
#' @return returns summary output for class 'ahreg'
#' 
#' @export
#' 
#====================================================
# Printing ahreg object
#====================================================
print.ahreg=function(x, digits=5, ...){
  cat("Call: ahreg() \n")
  cat("\n")
  
  print(paste0("Link: ", attr(x,"link")))
  print(attr(x,"formula"))
  
  cat("\n")
  
  if(!is.null(digits)){
    prmatrix(round(x$result, digits=digits))
  }else{
    prmatrix(x$result)
  }
  
}