#' @name ahreg.sample.data
#' @aliases  ahreg.sample.data
#' @title Generate a sample data from the pbc data
#' @description This is a function to retrieve 312 randomized patients from the pbc data.
#' @usage ahreg.sample.data(t.unit="year")
#' @param t.unit Specify the time unit. It supports "year" (default), "month", and "day".
#' @details The function creates a sample dataset to illustrate the usage of the function \code{ahreg()} in this package.
#' The original pbc data in \code{survival} package consists of 418 patients data.
#' This function loads the pbc data, select the 312 patients who were randomized.
#' The status variable is edited, so that 1 indicates death and 0 indicates alive.
#' @return returns a data frame
#' @seealso \code{pbc} in survival package
#' @examples
#' D=ahreg.sample.data()
#' head(D)
#' @export

#=======================================
# ahreg sample data
#=======================================
ahreg.sample.data=function(t.unit="year"){
  
  data("pbc", envir = environment()) ;
  tmp <- get("pbc", envir  = environment()) ;
  
  D=tmp[1:312,c(2:5,10,11,13,19)]
  
  if(t.unit=="year"){
    D$time=D$time/365.25
  }
  if(t.unit=="month"){
    D$time=D$time/365.25*12
  }
  if(t.unit=="day"){
    D$time=D$time
  }
  
  D$status=as.numeric(D$status==2);
  D$arm=as.numeric(D$trt==1) ;
  DA=D[,c(1, 2, 9, 4:8)]
  DA
}
NULL