#' EVZINB tidy function
#'
#' @param x 
#' @param coef 
#' @param p_value 
#' @param confint 
#' @param component 
#' @param ... 
#'
#' @return An EVZINB glance function
#' @export
#'
#' @examples tidy(evzinb)
glance.evzinb <- function(x,...){
  

  res <- tibble::tibble(
  nobs = nrow(x$data$x.nb),
  npar = length(x$par.all),
  alpha = x$coef$Alpha.NB,
  parameter = x$coef$C,
  aic = x$AIC,
  bic = x$BIC,
  logLik = x$log.lik)
  

  return(res)
  
  
}

glance.evzinb <- function(x,...){
  

  res <- tibble::tibble(
  nobs = nrow(x$data$x.nb),
  npar = length(x$par.all),
  alpha = round(x$coef$Alpha.NB,2),
  parameter = x$coef$C,
  aic = x$AIC,
  bic = x$BIC,
  logLik = round(x$log.lik,2))
  

  return(res)
  
  
}
