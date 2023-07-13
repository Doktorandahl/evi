#' EVZINB tidy function
#'
#' @param x 
#' @param coef 
#' @param p_value 
#' @param confint 
#' @param component 
#' @param ... 
#'
#' @return An EVZINB tidy function
#' @export
#'
#' @examples tidy(evzinb)
tidy.evzinb <- function(x,component = c('zi','evi','count','pareto','all'), coef = c('original','bootstrapped_mean','bootstrapped_median'), standard_error=TRUE, p_value = c('bootstrapped','approx','none'), confint = c('none','bootstrapped','approx'),conf_level = 0.95,approx_t_value = TRUE,symmetric_bootstrap_p = TRUE,...){

  coef <- match.arg(coef, c('original','bootstrapped_mean','bootstrapped_median'))
  p_value <- match.arg(p_value, c('bootstrapped','approx','none'))
  confint <- match.arg(confint, c('none','bootstrapped','approx'))
  component <- match.arg(component, c('zi','evi','count','pareto','all'))
  
  
  inv_leftjoin <- invisible(dplyr::left_join)
  
  nobs <- nrow(x$data$x.nb)
  npar <- length(x$par.all)
  if(!is.null(x$bootstraps)){

    nb_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term')
    
    zi_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.ZC') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term')
    
    evi_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term') 
    
    pareto_boot <- x$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'term') 

  }
  
  
  if(coef == 'original'){
    nb <- dplyr::tibble(term = names(x$coef$Beta.NB),estimate = x$coef$Beta.NB)
    zi <- dplyr::tibble(term = names(x$coef$Beta.multinom.ZC),estimate = x$coef$Beta.multinom.ZC)
    evi <- dplyr::tibble(term = names(x$coef$Beta.multinom.PL),estimate = x$coef$Beta.multinom.PL)
    pareto <- dplyr::tibble(term = names(x$coef$Beta.PL),estimate = x$coef$Beta.PL)
  }else if(coef == 'bootstrapped_mean'){
    nb <- dplyr::tibble(term = names(x$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                                      dplyr::summarize(estimate = mean(value)),by='term')
    zi <- dplyr::tibble(term = names(x$coef$Beta.multinom.ZC)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(term) %>%
                                                                                               dplyr::summarize(estimate = mean(value)),by='term')
    evi <- dplyr::tibble(term = names(x$coef$Beta.multinom.PL)) %>% 
      dplyr::left_join(evi_boot %>% dplyr::group_by(term) %>%
                         dplyr::summarize(estimate = mean(value)),by='term')
    pareto <- dplyr::tibble(term = names(x$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(term) %>%
                                                                                          dplyr::summarize(estimate = mean(value)),by='term')
  }else if(coef == 'bootstrapped_median'){
    nb <- dplyr::tibble(term = names(x$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                                                                      dplyr::summarize(estimate = median(value)),by='term')
    zi <- dplyr::tibble(term = names(x$coef$Beta.multinom.ZC)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(term) %>%
                                                                                               dplyr::summarize(estimate = median(value)),by='term')
    evi <- dplyr::tibble(term = names(x$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(term) %>%
                                                                                                dplyr::summarize(estimate = median(value)),by='term')
    pareto <- dplyr::tibble(term = names(x$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(term) %>%
                                                                                          dplyr::summarize(estimate = median(value)),by='term')
  }
  if(standard_error){
    nb <- nb %>% dplyr::left_join(nb_boot %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(std.error = sd(value)))
    zi <- zi %>% dplyr::left_join(zi_boot %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(std.error = sd(value)))
    evi <- evi %>% dplyr::left_join(evi_boot %>% dplyr::group_by(term) %>%
                                      dplyr::summarize(std.error = sd(value)))
    pareto <- pareto %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(term) %>%
                                            dplyr::summarize(std.error = sd(value)))
  }
  
  if(approx_t_value){
    nb <- nb %>% dplyr::mutate(statistic = estimate/std.error)
    zi <- zi %>% dplyr::mutate(statistic = estimate/std.error)
    evi <- evi %>% dplyr::mutate(statistic = estimate/std.error)
    pareto <- pareto %>% dplyr::mutate(statistic = estimate/std.error)
  }
  
  if(p_value == 'bootstrapped'){
    nb <- nb %>% dplyr::left_join(dplyr::left_join(nb_boot,nb) %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    zi <- zi %>% dplyr::left_join(dplyr::left_join(zi_boot,zi) %>% dplyr::group_by(term) %>%
                                    dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    evi <- evi %>% dplyr::left_join(dplyr::left_join(evi_boot,evi) %>% dplyr::group_by(term) %>%
                                      dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    pareto <- pareto %>% dplyr::left_join(dplyr::left_join(pareto_boot,pareto) %>% dplyr::group_by(term) %>%
                                            dplyr::summarize(p.value = bootstrap_p_value_calculator(value,estimate[1],symmetric=symmetric_bootstrap_p)))
    
  }else if(p_value == 'approx'){
    nb <- nb %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    zi <- zi %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    evi <- evi %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    pareto <- pareto %>% dplyr::mutate(p.value = 2*(1-pt(abs(approx_t),df = nobs-npar)))
  }
  
  if(component == 'zi'){
    return(zi)
  }else if(component == 'evi'){
    return(evi)
  }else if(component == 'count'){
    return(nb)
  }else if(component == 'pareto'){
    return(pareto)
  }else if(component == 'all'){
    return(dplyr::bind_rows(dplyr::mutate(zi,y.level='zi',.before=1),
                            dplyr::mutate(evi,y.level='evi',.before=1),
                            dplyr::mutate(nb,y.level='count',.before=1),
                            dplyr::mutate(pareto,y.level='pareto',.before=1)))
  }
  
return(res)
  

}

