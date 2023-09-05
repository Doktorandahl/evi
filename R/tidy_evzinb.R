#' EVZINB tidy function
#'
#' @param x An evzinb object
#' @param coef_type Type of coefficients. Original are the coefficient estimates from the non-bootstrapped version of the model. 'bootstrapped_mean' are the mean coefficients across bootstraps, and 'bootstrapped_median' are the median coefficients across bootstraps
#' @param p_value What type of p_values should be computed? 'bootstrapped' are bootstrapped p_values through confidence interval inversion. 'approx' are p-values based on the t-value produced by dividing the coefficient with the standard error.
#' @param confint What type of confidence should be computed. Same options as p_value
#' @param component Which component should be shown?
#' @param ... Other arguments passsed to tidy function
#' @param standard_error Should standard errors be computed?
#' @param conf_level What confidence level should be used for the confidence interval
#' @param approx_t_value Should approximate t-values be returned
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#'
#' @return An EVZINB tidy function
#' @export
#'
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' tidy(model)
tidy.evzinb <- function(x,component = c('zi','evinf','count','pareto','all'), coef_type = c('original','bootstrap_mean','bootstrap_median'), standard_error=TRUE, p_value = c('bootstrapped','approx','none'), confint = c('none','bootstrapped','approx'),conf_level = 0.95,approx_t_value = TRUE,symmetric_bootstrap_p = TRUE,...){

  coef_type <- match.arg(coef_type, c('original','bootstrap_mean','bootstrap_median'))
  p_value <- match.arg(p_value, c('bootstrapped','approx','none'))
  confint <- match.arg(confint, c('none','bootstrapped','approx'))
  component <- match.arg(component, c('zi','evinf','count','pareto','all'))
  
  
  inv_leftjoin <- invisible(dplyr::left_join)
  
  nobs <- nrow(x$data$x.nb)
  npar <- length(x$par.all)
  if(!is.null(x$bootstraps)){
    n_bootstraps_org <- length(x$bootstraps)
    x$bootstraps <- x$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(x$bootstraps)
    
    nb_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term')
    
    zi_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.ZC') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term')
    
    evi_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term') 
    
    pareto_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term') 

  }
  
  
  if(coef_type == 'original'){
    nb <- dplyr::tibble(term = names(x$coef$Beta.NB),estimate = x$coef$Beta.NB)
    zi <- dplyr::tibble(term = names(x$coef$Beta.multinom.ZC),estimate = x$coef$Beta.multinom.ZC)
    evinf <- dplyr::tibble(term = names(x$coef$Beta.multinom.PL),estimate = x$coef$Beta.multinom.PL)
    pareto <- dplyr::tibble(term = names(x$coef$Beta.PL),estimate = x$coef$Beta.PL)
  }else if(coef_type == 'bootstrap_mean'){
    nb <- dplyr::tibble(term = names(x$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                                      dplyr::summarize(estimate = mean(value)),by='term')
    zi <- dplyr::tibble(term = names(x$coef$Beta.multinom.ZC)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(term) %>%
                                                                                               dplyr::summarize(estimate = mean(value)),by='term')
    evinf <- dplyr::tibble(term = names(x$coef$Beta.multinom.PL)) %>% 
      dplyr::left_join(evi_boot %>% dplyr::group_by(term) %>%
                         dplyr::summarize(estimate = mean(value)),by='term')
    pareto <- dplyr::tibble(term = names(x$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(term) %>%
                                                                                          dplyr::summarize(estimate = mean(value)),by='term')
  }else if(coef_type == 'bootstrap_median'){
    nb <- dplyr::tibble(term = names(x$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                                      dplyr::summarize(estimate = median(value)),by='term')
    zi <- dplyr::tibble(term = names(x$coef$Beta.multinom.ZC)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(term) %>%
                                                                                               dplyr::summarize(estimate = median(value)),by='term')
    evinf <- dplyr::tibble(term = names(x$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(term) %>%
                                                                                                dplyr::summarize(estimate = median(value)),by='term')
    pareto <- dplyr::tibble(term = names(x$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(term) %>%
                                                                                          dplyr::summarize(estimate = median(value)),by='term')
  }
  if(standard_error){
    nb <- nb %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(std.error = sd(value)))
    zi <- zi %>% dplyr::left_join(zi_boot %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(std.error = sd(value)))
    evinf <- evinf %>% dplyr::left_join(evi_boot %>% dplyr::group_by(term) %>%
                                      dplyr::summarize(std.error = sd(value)))
    pareto <- pareto %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(term) %>%
                                            dplyr::summarize(std.error = sd(value)))
  }
  
  if(approx_t_value){
    nb <- nb %>% dplyr::mutate(statistic = estimate/std.error)
    zi <- zi %>% dplyr::mutate(statistic = estimate/std.error)
    evinf <- evinf %>% dplyr::mutate(statistic = estimate/std.error)
    pareto <- pareto %>% dplyr::mutate(statistic = estimate/std.error)
  }
  
  if(p_value == 'bootstrapped'){
    nb <- nb %>% dplyr::left_join(dplyr::left_join(nb_boot,nb) %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    zi <- zi %>% dplyr::left_join(dplyr::left_join(zi_boot,zi) %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    evinf <- evinf %>% dplyr::left_join(dplyr::left_join(evi_boot,evinf) %>% dplyr::group_by(term) %>%
                                      dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    pareto <- pareto %>% dplyr::left_join(dplyr::left_join(pareto_boot,pareto) %>% dplyr::group_by(term) %>%
                                            dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    
  }else if(p_value == 'approx'){
    nb <- nb %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    zi <- zi %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    evinf <- evinf %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    pareto <- pareto %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
  }
  
  if(component == 'zi'){
    return(zi)
  }else if(component == 'evinf'){
    return(evinf)
  }else if(component == 'count'){
    return(nb)
  }else if(component == 'pareto'){
    return(pareto)
  }else if(component == 'all'){
    return(dplyr::bind_rows(dplyr::mutate(zi,y.level='zi',.before=1),
                            dplyr::mutate(evinf,y.level='evinf',.before=1),
                            dplyr::mutate(nb,y.level='count',.before=1),
                            dplyr::mutate(pareto,y.level='pareto',.before=1)))
  }
  
return(res)
  

}

#' EVINB tidy function
#'
#' @param x An evinb object
#' @param coef_type Type of coefficients. Original are the coefficient estimates from the non-bootstrapped version of the model. 'bootstrapped_mean' are the mean coefficients across bootstraps, and 'bootstrapped_median' are the median coefficients across bootstraps
#' @param p_value What type of p_values should be computed? 'bootstrapped' are bootstrapped p_values through confidence interval inversion. 'approx' are p-values based on the t-value produced by dividing the coefficient with the standard error.
#' @param confint What type of confidence should be computed. Same options as p_value
#' @param component Which component should be shown?
#' @param ... Other arguments passsed to tidy function
#' @param standard_error Should standard errors be computed?
#' @param conf_level What confidence level should be used for the confidence interval
#' @param approx_t_value Should approximate t-values be returned
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#'
#' @return An EVINB tidy function
#' @export
#'
#' @examples data(genevzinb)
#' model <- evinb(y~x1+x2+x3,data=genevzinb)
#' tidy(model)
tidy.evinb <- function(x,component = c('evinf','count','pareto','all'), coef_type = c('original','bootstrap_mean','bootstrap_median'), standard_error=TRUE, p_value = c('bootstrapped','approx','none'), confint = c('none','bootstrapped','approx'),conf_level = 0.95,approx_t_value = TRUE,symmetric_bootstrap_p = TRUE,...){
  
  coef_type <- match.arg(coef_type, c('original','bootstrap_mean','bootstrap_median'))
  p_value <- match.arg(p_value, c('bootstrapped','approx','none'))
  confint <- match.arg(confint, c('none','bootstrapped','approx'))
  component <- match.arg(component, c('evinf','count','pareto','all'))
  
  
  inv_leftjoin <- invisible(dplyr::left_join)
  
  nobs <- nrow(x$data$x.nb)
  npar <- length(x$par.all)
  if(!is.null(x$bootstraps)){
    n_bootstraps_org <- length(x$bootstraps)
    x$bootstraps <- x$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(x$bootstraps)
    
    nb_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term')
    
    evi_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term') 
    
    pareto_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term') 
    
  }
  
  
  if(coef_type == 'original'){
    nb <- dplyr::tibble(term = names(x$coef$Beta.NB),estimate = x$coef$Beta.NB)
      evinf <- dplyr::tibble(term = names(x$coef$Beta.multinom.PL),estimate = x$coef$Beta.multinom.PL)
    pareto <- dplyr::tibble(term = names(x$coef$Beta.PL),estimate = x$coef$Beta.PL)
  }else if(coef_type == 'bootstrap_mean'){
    nb <- dplyr::tibble(term = names(x$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                             dplyr::summarize(estimate = mean(value)),by='term')
 
    evinf <- dplyr::tibble(term = names(x$coef$Beta.multinom.PL)) %>% 
      dplyr::left_join(evi_boot %>% dplyr::group_by(term) %>%
                         dplyr::summarize(estimate = mean(value)),by='term')
    pareto <- dplyr::tibble(term = names(x$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(term) %>%
                                                                                 dplyr::summarize(estimate = mean(value)),by='term')
  }else if(coef_type == 'bootstrap_median'){
    nb <- dplyr::tibble(term = names(x$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                             dplyr::summarize(estimate = median(value)),by='term')
 
    evinf <- dplyr::tibble(term = names(x$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(term) %>%
                                                                                       dplyr::summarize(estimate = median(value)),by='term')
    pareto <- dplyr::tibble(term = names(x$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(term) %>%
                                                                                 dplyr::summarize(estimate = median(value)),by='term')
  }
  if(standard_error){
    nb <- nb %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(std.error = sd(value)))
  
    evinf <- evinf %>% dplyr::left_join(evi_boot %>% dplyr::group_by(term) %>%
                                      dplyr::summarize(std.error = sd(value)))
    pareto <- pareto %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(term) %>%
                                            dplyr::summarize(std.error = sd(value)))
  }
  
  if(approx_t_value){
    nb <- nb %>% dplyr::mutate(statistic = estimate/std.error)
  
    evinf <- evinf %>% dplyr::mutate(statistic = estimate/std.error)
    pareto <- pareto %>% dplyr::mutate(statistic = estimate/std.error)
  }
  
  if(p_value == 'bootstrapped'){
    nb <- nb %>% dplyr::left_join(dplyr::left_join(nb_boot,nb) %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
  
    evinf <- evinf %>% dplyr::left_join(dplyr::left_join(evi_boot,evinf) %>% dplyr::group_by(term) %>%
                                      dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    pareto <- pareto %>% dplyr::left_join(dplyr::left_join(pareto_boot,pareto) %>% dplyr::group_by(term) %>%
                                            dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    
  }else if(p_value == 'approx'){
    nb <- nb %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
     evinf <- evinf %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    pareto <- pareto %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
  }
  
if(component == 'evinf'){
    return(evinf)
  }else if(component == 'count'){
    return(nb)
  }else if(component == 'pareto'){
    return(pareto)
  }else if(component == 'all'){
    return(dplyr::bind_rows(
                            dplyr::mutate(evinf,y.level='evinf',.before=1),
                            dplyr::mutate(nb,y.level='count',.before=1),
                            dplyr::mutate(pareto,y.level='pareto',.before=1)))
  }
  
  return(res)
  
  
}


