
#' Tidy function for zinbboot
#'
#' @param x A fitted bootstrapped zero-inflated model 
#' @param coef_type What type of coefficient should be reported, original, bootstrapped mean, or bootstrapped median
#' @param standard_error Should bootstrapped standard errors be reported?
#' @param p_value What type of p-value should be reported? Bootstrapped p_values, approximate p-values, or none?
#' @param confint What type of confidence intervals should be reported? Bootstrapped p_values, approximate p-values, or none?
#' @param conf_level Confidence level for confidence intervals
#' @param approx_t_value Should approximate t_values be reported
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#' @param ... Other arguments to be passed to tidy
#' @param component Which component should be shown?
#'
#' @return A tidy function for a bootstrapped zinb model
#' 
#' @export
#'
#' @examples 
#' \dontrun{
#' data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' zinb_comp <- compare_models(model)
#' tidy(zinb_comp$zinb)
#' }
tidy.zinbboot <- function(x, component = c('zi','count','all'),coef_type = c('original','bootstrap_mean','bootstrap_median'), standard_error=TRUE, p_value = c('bootstrapped','approx','none'), confint = c('none','bootstrapped','approx'),conf_level = 0.95,approx_t_value = TRUE,symmetric_bootstrap_p = TRUE,...){
  
  coef_type <- match.arg(coef_type, c('original','bootstrap_mean','bootstrap_median'))
  p_value <- match.arg(p_value, c('bootstrapped','approx','none'))
  confint <- match.arg(confint, c('none','bootstrapped','approx'))
  component <- match.arg(component, c('zi','count','all'))
  
  inv_leftjoin <- invisible(dplyr::left_join)
  
  nobs <- x$full_run$n
  npar <- nrow(x$full_run$vcov)+1
  if(!is.null(x$bootstraps)){
    
    n_bootstraps_org <- length(x$bootstraps)
    x$bootstraps <- x$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(x$bootstraps)
    
    nb_boot <- x$bootstraps %>% purrr::map('coefficients') %>% purrr::map('count') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'term')
    
    zi_boot <- x$bootstraps %>% purrr::map('coefficients') %>% purrr::map('zero') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'term')
    
    
  }
  
  
  if(coef_type == 'original'){
    nb <- dplyr::tibble(term = names(x$full_run$coefficients$count),estimate = x$full_run$coefficients$count)
    zi <- dplyr::tibble(term = names(x$full_run$coefficients$zero),estimate = x$full_run$coefficients$zero)
    
  }else if(coef_type == 'bootstrap_mean'){
    nb <- dplyr::tibble(term = names(x$full_run$coefficients$count)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                                            dplyr::summarize(estimate = mean(value)),by='term')
    zi <- dplyr::tibble(term = names(x$full_run$coefficients$zero)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(term) %>%
                                                                                           dplyr::summarize(estimate = mean(value)),by='term')
    
  }else if(coef_type == 'bootstrap_median'){
    nb <- dplyr::tibble(term = names(x$full_run$coefficients$count)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                                            dplyr::summarize(estimate = median(value)),by='term')
    zi <- dplyr::tibble(term = names(x$full_run$coefficients$zero)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(term) %>%
                                                                                           dplyr::summarize(estimate = median(value)),by='term')
    
  }
  if(standard_error){
    nb <- nb %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(std.error = sd(value)))
    zi <- zi %>% dplyr::left_join(zi_boot %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(std.error = sd(value)))
    
  }
  
  if(approx_t_value){
    nb <- nb %>% dplyr::mutate(statistic = estimate/std.error)
    zi <- zi %>% dplyr::mutate(statistic = estimate/std.error)
    
  }
  
  if(p_value == 'bootstrapped'){
    nb <- nb %>% dplyr::left_join(dplyr::left_join(nb_boot,nb) %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    zi <- zi %>% dplyr::left_join(dplyr::left_join(zi_boot,zi) %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    
  }else if(p_value == 'approx'){
    nb <- nb %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    zi <- zi %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    
  }
  
  if(component == 'zi'){
    return(zi)
  }else if(component == 'count'){
    return(nb)
  }else if(component == 'all'){
    return(dplyr::bind_rows(dplyr::mutate(zi,y.level='zi',.before=1),
                            dplyr::mutate(nb,y.level='count',.before=1)))
  }
  
  return(res)
  
  
}

#' Tidy function for nbboot
#'
#' @param x A fitted bootstrapped zero-inflated model 
#' @param coef_type What type of coefficient should be reported, original, bootstrapped mean, or bootstrapped median
#' @param standard_error Should bootstrapped standard errors be reported?
#' @param p_value What type of p-value should be reported? Bootstrapped p_values, approximate p-values, or none?
#' @param confint What type of confidence intervals should be reported? Bootstrapped p_values, approximate p-values, or none?
#' @param conf_level Confidence level for confidence intervals
#' @param approx_t_value Should approximate t_values be reported
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#' @param ... Other arguments to be passed to tidy
#' @param include_ylev Logical. Should y.lev be included in the tidy output? Makes for nicer tables when using modelsummary
#' 
#' @return A tidy function for a bootstrapped nb model
#' 
#' @export
#'
#' @examples 
#' \dontrun{
#' data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' zinb_comp <- compare_models(model)
#' tidy(zinb_comp$nb) 
#' }
tidy.nbboot <- function(x, coef_type = c('original','bootstrap_mean','bootstrap_median'), standard_error=TRUE, p_value = c('bootstrapped','approx','none'), confint = c('none','bootstrapped','approx'),conf_level = 0.95,approx_t_value = TRUE,symmetric_bootstrap_p = TRUE,include_ylev = FALSE,...){
  
  coef_type <- match.arg(coef_type, c('original','bootstrap_mean','bootstrap_median'))
  p_value <- match.arg(p_value, c('bootstrapped','approx','none'))
  confint <- match.arg(confint, c('none','bootstrapped','approx'))
  
  inv_leftjoin <- invisible(dplyr::left_join)
  
  nobs <- x$full_run$n
  npar <- x$full_run$rank
  if(!is.null(x$bootstraps)){
    
    n_bootstraps_org <- length(x$bootstraps)
    x$bootstraps <- x$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(x$bootstraps)
    
    nb_boot <- x$bootstraps %>% purrr::map('coefficients') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'term')
    
  }
  
  
  if(coef_type == 'original'){
    nb <- dplyr::tibble(term = names(x$full_run$coefficients),estimate = x$full_run$coefficients)
    
    
  }else if(coef_type == 'bootstrap_mean'){
    nb <- dplyr::tibble(term = names(x$full_run$coefficients)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                                      dplyr::summarize(estimate = mean(value)),by='term')
    
  }else if(coef_type == 'bootstrap_median'){
    nb <- dplyr::tibble(term = names(x$full_run$coefficients)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                                      dplyr::summarize(estimate = median(value)),by='term')
    
  }
  if(standard_error){
    nb <- nb %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(std.error = sd(value)))
    
    
  }
  
  if(approx_t_value){
    nb <- nb %>% dplyr::mutate(statistic = estimate/std.error)
    
    
  }
  
  if(p_value == 'bootstrapped'){
    nb <- nb %>% dplyr::left_join(dplyr::left_join(nb_boot,nb) %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    
    
  }else if(p_value == 'approx'){
    nb <- nb %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    
    
  }
  
  if(include_ylev){
    nb <- nb %>% dplyr::mutate(y.level = 'count')
  }
  
  return(nb)
}
