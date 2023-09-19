#' EVZINB summary function
#'
#' @param object an EVZINB object with bootstraps
#' @param coef Type of coefficients. Original are the coefficient estimates from the non-bootstrapped version of the model. 'bootstrapped_mean' are the mean coefficients across bootstraps, and 'bootstrapped_median' are the median coefficients across bootstraps
#' @param standard_error Should standard errors be computed?
#' @param p_value What type of p_values should be computed? 'bootstrapped' are bootstrapped p_values through confidence interval inversion. 'approx' are p-values based on the t-value produced by dividing the coefficient with the standard error.
#' @param bootstrapped_props Type of bootstrapped proportions of component proportions to be returned
#' @param approx_t_value Should approximate t-values be returned
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#' @param ... Additional arguments passed to the summary function
#'
#'
#' @return An EVZINB summary object
#' @export
#'
#' @examples 
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 10, multicore = TRUE, ncores = 2)
#' summary(model)
#' 
summary.evzinb <- function(object,coef = c('original','bootstrapped_mean','bootstrapped_median'),standard_error = TRUE, p_value = c('bootstrapped','approx','both','none'), bootstrapped_props = c('none','mean','median'),approx_t_value = TRUE,
                           symmetric_bootstrap_p = TRUE,...){

coef <- match.arg(coef,c('original','bootstrapped_mean','bootstrapped_median'))

p_value <- match.arg(p_value, c('bootstrapped','approx','both','none'))

bootstrapped_props <- match.arg(bootstrapped_props,c('none','mean','median'))


  nobs <- nrow(object$data$x.nb)
  npar <- length(object$par.all)
  
  props <- object$props %>% dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'state') %>%
    dplyr::group_by(.data$state) %>%
    dplyr::summarize(mean_prop = mean(.data$value)) %>%
    dplyr::slice(c(3,1,2))
  alpha_nb <- c(object$coef$Alpha.NB)
  names(alpha_nb) <- c('Alpha_NB')
  C_est <- c(object$coef$C)
  names(C_est) <- c('C')
  
  if(!is.null(object$bootstraps)){
    n_bootstraps_org <- length(object$bootstraps)
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(object$bootstraps)
    
  nb_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'Variable')

  zi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.ZC') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'Variable')

  evi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'Variable') 

  pareto_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'Variable') 

  props_boot <-  object$bootstraps %>% purrr::map('props') %>% purrr::map(colMeans) %>% purrr::reduce(rbind) %>%
    dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'state') %>%
    dplyr::group_by(.data$state) %>%
    dplyr::summarize(bootstrap_mean = mean(.data$value),bootstrap_median = median(.data$value), standard_error = sd(.data$value))
  
  if(bootstrapped_props == 'median'){
    props_boot <- props_boot %>% dplyr::select(.data$state,.data$bootstrap_median,.data$standard_error)
  }else if(bootstrapped_props == 'mean'){
    props_boot <- props_boot %>% dplyr::select(.data$state,.data$bootstrap_mean,.data$standard_error)
  }else{
    props_boot <- props_boot %>% dplyr::select(.data$state,.data$standard_error)
  }
  
    props <- dplyr::left_join(props,props_boot)

  alpha_nb_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Alpha.NB') %>% purrr::reduce(c)
  
  alpha_nb <- c(alpha_nb, mean(alpha_nb_boot),median(alpha_nb_boot),sd(alpha_nb_boot))
  names(alpha_nb) <- c('Alpha_NB','bootstrap_mean','bootstrap_median','standard_error')

  C_est_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('C')%>% purrr::reduce(c)
  C_est <- c(C_est, mean(C_est_boot),median(C_est_boot),sd(C_est_boot))
  names(C_est) <- c('C','bootstrap_mean','bootstrap_median','standard_error')
}


  if(coef == 'original'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB),Estimate = object$coef$Beta.NB)
    zi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.ZC),Estimate = object$coef$Beta.multinom.ZC)
    evinf <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL),Estimate = object$coef$Beta.multinom.PL)
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL),Estimate = object$coef$Beta.PL)
  }else if(coef == 'bootstrapped_mean'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(.data$Variable) %>%
      dplyr::summarize(Estimate = mean(.data$value)))
  zi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.ZC)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                             dplyr::summarize(Estimate = mean(.data$value)))
  evinf <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                             dplyr::summarize(Estimate = mean(.data$value)))
  pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                              dplyr::summarize(Estimate = mean(.data$value)))
  }else if(coef == 'bootstrapped_median'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                      dplyr::summarize(Estimate = median(.data$value)))
    zi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.ZC)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                               dplyr::summarize(Estimate = median(.data$value)))
    evinf <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                                dplyr::summarize(Estimate = median(.data$value)))
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                          dplyr::summarize(Estimate = median(.data$value)))
  }
  if(standard_error){
    nb <- nb %>% dplyr::left_join(nb_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                      dplyr::summarize(se = sd(.data$value)))
    zi <- zi %>% dplyr::left_join(zi_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                               dplyr::summarize(se = sd(.data$value)))
    evinf <- evinf %>% dplyr::left_join(evi_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                                dplyr::summarize(se = sd(.data$value)))
    pareto <- pareto %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                          dplyr::summarize(se = sd(.data$value)))
  }
  
  if(approx_t_value){
    nb <- nb %>% dplyr::mutate(approx_t = .data$Estimate/.data$se)
    zi <- zi %>% dplyr::mutate(approx_t = .data$Estimate/.data$se)
    evinf <- evinf %>% dplyr::mutate(approx_t = .data$Estimate/.data$se)
    pareto <- pareto %>% dplyr::mutate(approx_t = .data$Estimate/.data$se)
  }
  
  if(p_value %in% c('bootstrapped','both')){
    nb <- nb %>% dplyr::left_join(dplyr::left_join(nb_boot,nb) %>% dplyr::group_by(.data$Variable) %>%
                                    dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(.data$value,.data$Estimate[1],symmetric=symmetric_bootstrap_p)))
    zi <- zi %>% dplyr::left_join(dplyr::left_join(zi_boot,zi) %>% dplyr::group_by(.data$Variable) %>%
                                    dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(.data$value,.data$Estimate[1],symmetric=symmetric_bootstrap_p)))
    evinf <- evinf %>% dplyr::left_join(dplyr::left_join(evi_boot,evinf) %>% dplyr::group_by(.data$Variable) %>%
                                      dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(.data$value,.data$Estimate[1],symmetric=symmetric_bootstrap_p)))
    pareto <- pareto %>% dplyr::left_join(dplyr::left_join(pareto_boot,pareto) %>% dplyr::group_by(.data$Variable) %>%
                                            dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(.data$value,.data$Estimate[1],symmetric=symmetric_bootstrap_p)))
    
  }else if(p_value %in% c('approx','both')){
    nb <- nb %>% dplyr::mutate(approx_p = 2*(1-pt(abs(.data$approx_t),df = nobs-npar)))
    zi <- zi %>% dplyr::mutate(approx_p = 2*(1-pt(abs(.data$approx_t),df = nobs-npar)))
    evinf <- evinf %>% dplyr::mutate(approx_p = 2*(1-pt(abs(.data$approx_t),df = nobs-npar)))
    pareto <- pareto %>% dplyr::mutate(approx_p = 2*(1-pt(abs(.data$approx_t),df = nobs-npar)))
  }

  res <- list(coefficients = list(negative_binomial = nb, 
                                    zero_inflation = zi,
                                    extreme_value_inflation = evinf,
                                    pareto = pareto),
                model_statistics = list(Alpha_nb = alpha_nb,
                                        C = C_est,
                                        Obs = c(Obs = nobs,pars=npar,df=nobs-npar)),
                component_proportions = props,
                n_failed_bootstraps = n_failed_bootstraps)

  class(res) <- 'summary.evzinb'

  return(res)


}

#' EVINB summary function
#'
#' @param object an EVINB object with bootstraps
#' @param coef Type of coefficients. Original are the coefficient estimates from the non-bootstrapped version of the model. 'bootstrapped_mean' are the mean coefficients across bootstraps, and 'bootstrapped_median' are the median coefficients across bootstraps
#' @param standard_error Should standard errors be computed?
#' @param p_value What type of p_values should be computed? 'bootstrapped' are bootstrapped p_values through confidence interval inversion. 'approx' are p-values based on the t-value produced by dividing the coefficient with the standard error.
#' @param bootstrapped_props Type of bootstrapped proportions of component proportions to be returned
#' @param approx_t_value Should approximate t-values be returned
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#' @param ... Additional arguments passed to the summary function
#'
#'
#' @return An EVINB summary object
#' @export
#'
#' @examples 
#' 
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 10, multicore = TRUE, ncores = 2)
#' summary(model)
#' 
summary.evinb <- function(object,coef = c('original','bootstrapped_mean','bootstrapped_median'),standard_error = TRUE, p_value = c('bootstrapped','approx','both','none'), bootstrapped_props = c('none','mean','median'),approx_t_value = TRUE,
                          symmetric_bootstrap_p = TRUE,...){
  
  coef <- match.arg(coef,c('original','bootstrapped_mean','bootstrapped_median'))
  
  p_value <- match.arg(p_value, c('bootstrapped','approx','both','none'))
  
  bootstrapped_props <- match.arg(bootstrapped_props,c('none','mean','median'))
  
  
  nobs <- nrow(object$data$x.nb)
  npar <- length(object$par.all)
  
  props <- object$props %>% dplyr::as_tibble(.name_repair = ~c('negative_binomial','pareto')) %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'state') %>%
    dplyr::group_by(.data$state) %>%
    dplyr::summarize(mean_prop = mean(.data$value))
  alpha_nb <- c(object$coef$Alpha.NB)
  names(alpha_nb) <- c('Alpha_NB')
  C_est <- c(object$coef$C)
  names(C_est) <- c('C')
  
  if(!is.null(object$bootstraps)){
    n_bootstraps_org <- length(object$bootstraps)
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(object$bootstraps)
    
    nb_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'Variable')
    
    
    evi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'Variable') 
    
    pareto_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'Variable') 
    
    props_boot <-  object$bootstraps %>% purrr::map('props') %>% purrr::map(colMeans) %>% purrr::reduce(rbind) %>%
      dplyr::as_tibble(.name_repair = ~c('negative_binomial','pareto')) %>% tidyr::pivot_longer(dplyr::everything(),names_to = 'state') %>%
      dplyr::group_by(.data$state) %>%
      dplyr::summarize(bootstrap_mean = mean(.data$value),bootstrap_median = median(.data$value), standard_error = sd(.data$value))
    
    if(bootstrapped_props == 'median'){
      props_boot <- props_boot %>% dplyr::select(.data$state,.data$bootstrap_median,.data$standard_error)
    }else if(bootstrapped_props == 'mean'){
      props_boot <- props_boot %>% dplyr::select(.data$state,.data$bootstrap_mean,.data$standard_error)
    }else{
      props_boot <- props_boot %>% dplyr::select(.data$state,.data$standard_error)
    }
    
    props <- dplyr::left_join(props,props_boot)
    
    alpha_nb_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Alpha.NB') %>% purrr::reduce(c)
    
    alpha_nb <- c(alpha_nb, mean(alpha_nb_boot),median(alpha_nb_boot),sd(alpha_nb_boot))
    names(alpha_nb) <- c('Alpha_NB','bootstrap_mean','bootstrap_median','standard_error')
    
    C_est_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('C')%>% purrr::reduce(c)
    C_est <- c(C_est, mean(C_est_boot),median(C_est_boot),sd(C_est_boot))
    names(C_est) <- c('C','bootstrap_mean','bootstrap_median','standard_error')
  }
  
  
  if(coef == 'original'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB),Estimate = object$coef$Beta.NB)

    evinf <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL),Estimate = object$coef$Beta.multinom.PL)
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL),Estimate = object$coef$Beta.PL)
  }else if(coef == 'bootstrapped_mean'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                      dplyr::summarize(Estimate = mean(.data$value)))

    evinf <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                                dplyr::summarize(Estimate = mean(.data$value)))
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                          dplyr::summarize(Estimate = mean(.data$value)))
  }else if(coef == 'bootstrapped_median'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                      dplyr::summarize(Estimate = median(.data$value)))

    evinf <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                                dplyr::summarize(Estimate = median(.data$value)))
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(.data$Variable) %>%
                                                                                          dplyr::summarize(Estimate = median(.data$value)))
  }
  if(standard_error){
    nb <- nb %>% dplyr::left_join(nb_boot %>% dplyr::group_by(.data$Variable) %>%
                                    dplyr::summarize(se = sd(.data$value)))

    evinf <- evinf %>% dplyr::left_join(evi_boot %>% dplyr::group_by(.data$Variable) %>%
                                      dplyr::summarize(se = sd(.data$value)))
    pareto <- pareto %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(.data$Variable) %>%
                                            dplyr::summarize(se = sd(.data$value)))
  }
  
  if(approx_t_value){
    nb <- nb %>% dplyr::mutate(approx_t = .data$Estimate/.data$se)
    evinf <- evinf %>% dplyr::mutate(approx_t = .data$Estimate/.data$se)
    pareto <- pareto %>% dplyr::mutate(approx_t = .data$Estimate/.data$se)
  }
  
  if(p_value %in% c('bootstrapped','both')){
    nb <- nb %>% dplyr::left_join(dplyr::left_join(nb_boot,nb) %>% dplyr::group_by(.data$Variable) %>%
                                    dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(.data$value,.data$Estimate[1],symmetric=symmetric_bootstrap_p)))
    
    
    evinf <- evinf %>% dplyr::left_join(dplyr::left_join(evi_boot,evinf) %>% dplyr::group_by(.data$Variable) %>%
                                      dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(.data$value,.data$Estimate[1],symmetric=symmetric_bootstrap_p)))
    pareto <- pareto %>% dplyr::left_join(dplyr::left_join(pareto_boot,pareto) %>% dplyr::group_by(.data$Variable) %>%
                                            dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(.data$value,.data$Estimate[1],symmetric=symmetric_bootstrap_p)))
    
  }else if(p_value %in% c('approx','both')){
    nb <- nb %>% dplyr::mutate(approx_p = 2*(1-pt(abs(.data$approx_t),df = nobs-npar)))
    evinf <- evinf %>% dplyr::mutate(approx_p = 2*(1-pt(abs(.data$approx_t),df = nobs-npar)))
    pareto <- pareto %>% dplyr::mutate(approx_p = 2*(1-pt(abs(.data$approx_t),df = nobs-npar)))
  }
  
  res <- list(coefficients = list(negative_binomial = nb, 
                                  extreme_value_inflation = evinf,
                                  pareto = pareto),
              model_statistics = list(Alpha_nb = alpha_nb,
                                      C = C_est,
                                      Obs = c(Obs = nobs,pars=npar,df=nobs-npar)),
              component_proportions = props,
              n_failed_bootstraps = n_failed_bootstraps)
  
  class(res) <- 'summary.evinb'
  
  return(res)
  
  
}


error_remover <- function(object){
  if('try-error' %in% class(object)){
    return(NULL)
  }else{
    return(object)
  }
}

bootstrap_p_value_calculator <- function(x,estimate = NULL, estimate_fallback = c('median','mean'), symmetric = TRUE){
  if(is.null(estimate)){
    estimate_fallback <- match.arg(estimate_fallback,c('median','mean'))
    estimate <- do.call(estimate_fallback,list(x=x))
  }
  if(symmetric){
   if(estimate>=0){
     return(min(1,2*mean(x<0)))
   }else{
     return(min(1,2*mean(x>=0)))
   } 
  }else{
    return(mean(abs(x-estimate)>=abs(estimate)))
  }
}



