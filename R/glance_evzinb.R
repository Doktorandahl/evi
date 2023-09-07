#' EVZINB and EVINB glance functions
#'
#' @param x An EVZINB or EVINB object
#' @param ... Further arguments to be passed to glance()
#'
#' @return An EVZINB glance function
#' @seealso \code{\link[generics]{glance}}
#' @export
#'
#' @examples 
#' \dontrun{
#' data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb, n_bootstraps = 10)
#' glance(model)
#' }
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

#' EVZINB and EVINB glance functions
#'
#' @param x An EVZINB or EVINB object
#' @param ... Further arguments to be passed to glance()
#'
#' @return An EVZINB glance function
#' @seealso \code{\link[generics]{glance}}
#' @export
#'
#' @examples 
#' \dontrun{
#' data(genevzinb)
#' model <- evinb(y~x1+x2+x3,data=genevzinb, n_bootstraps = 10)
#' glance(model)
#' }
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
