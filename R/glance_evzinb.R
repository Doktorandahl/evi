#' EVZINB and EVINB glance functions
#'
#' @param x An EVZINB or EVINB object
#' @param ... Further arguments to be passed to glance()
#'
#' @return An EVZINB glance function
#' @aliases evzinb-glance glance, evzinb-method, 
#' @aliases evinb-glance glance, evinb-method
#' @seealso \code{\link[generics]{glance}}
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

#' @rdname glance.evzinb
glance.evinb <- function(x,...){
  

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
