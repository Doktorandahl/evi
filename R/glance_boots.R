#' zinbboot and nboot glance functions
#'
#' @param x An nbboot or zinbboot object
#' @param ... Further arguments to be passed to glance()
#'
#' @return An nbboot glance function
#' @seealso \code{\link[generics]{glance}}
#' @export
#'
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' zinb_comp <- compare_models(model)
#' glance(zinb_comp$zinb)
glance.zinbboot <- function(x,...){
  
  
  res <- tibble::tibble(
    nobs = x$full_run$n,
    npar = nrow(x$full_run$vcov)+1,
    alpha = 1/x$full_run$theta,
    aic = AIC(x$full_run),
    bic = BIC(x$full_run),
    logLik = x$full_run$log.lik)
  
  
  return(res)
  
  
}

#' zinbboot and nboot glance functions
#'
#' @param x An nbboot or zinbboot object
#' @param ... Further arguments to be passed to glance()
#'
#' @return An nbboot glance function
#' @seealso \code{\link[generics]{glance}}
#' @export
#'
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' zinb_comp <- compare_models(model)
#' glance(zinb_comp$nb)
glance.nbboot <- function(x,...){
  
  
  res <- tibble::tibble(
    nobs = nrow(x$full_run$model),
    npar = x$full_run$rank,
    alpha = 1/x$full_run$theta,
    aic = AIC(x$full_run),
    bic = BIC(x$full_run),
    logLik = logLik(x$full_run))
  
  
  return(res)
  
  
}
