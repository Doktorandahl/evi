#' EVZINB summary function
#'
#' @param object an EVZINB object with bootstraps
#' @param boot_mean Logical: Should the bootstrapped mean of the estimates be computed?
#' @param boot_se Logical: Should the bootstrapped standard errors of the estimates be computed?
#' @param p_approx Logical: Should the normal approximated p-value be computed? Not recommended unless the number of bootstraps is too low for bootstrapped p-values
#'
#' @return An EVZINB summary object
#' @export
#'
#' @examples summary(evzinb)
summary.evzinb <- function(object,coef = c('original','bootstrapped_mean','bootstrapped_median'),standard_error = TRUE, p_value = c('bootstrapped','approx','both','none'), bootstrapped_props = c('none','mean','median'),approx_t_value = TRUE,
                           symmetric_bootstrap_p = TRUE,...){

coef <- match.arg(coef,c('original','bootstrapped_mean','bootstrapped_median'))

p_value <- match.arg(p_value, c('bootstrapped','approx','both','none'))

bootstrapped_props <- match.arg(bootstrapped_props,c('none','mean','median'))


  nobs <- nrow(object$data$x.nb)
  npar <- length(object$par.all)
  
  props <- object$props %>% dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(everything(),names_to = 'state') %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(mean_prop = mean(value)) %>%
    dplyr::slice(c(3,1,2))
  alpha_nb <- c(object$coef$Alpha.NB)
  names(alpha_nb) <- c('Alpha_NB')
  C_est <- c(object$coef$C)
  names(C_est) <- c('C')
  
  if(!is.null(object$bootstraps)){
    n_bootstraps_org <- length(object$bootstraps)
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(object$bootstraps)
    
  nb_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable')

  zi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.ZC') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable')

  evi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') 

  pareto_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') 

  props_boot <-  object$bootstraps %>% purrr::map('props') %>% purrr::map(colMeans) %>% purrr::reduce(rbind) %>%
    dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(everything(),names_to = 'state') %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(bootstrap_mean = mean(value),bootstrap_median = median(value), standard_error = sd(value))
  
  if(bootstrapped_props == 'median'){
    props_boot <- props_boot %>% dplyr::select(state,bootstrap_median,standard_error)
  }else if(bootstrapped_props == 'mean'){
    props_boot <- props_boot %>% dplyr::select(state,bootstrap_mean,standard_error)
  }else{
    props_boot <- props_boot %>% dplyr::select(state,standard_error)
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
    evi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL),Estimate = object$coef$Beta.multinom.PL)
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL),Estimate = object$coef$Beta.PL)
  }else if(coef == 'bootstrapped_mean'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(Variable) %>%
      dplyr::summarize(Estimate = mean(value)))
  zi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.ZC)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(Variable) %>%
                                                                                             dplyr::summarize(Estimate = mean(value)))
  evi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(Variable) %>%
                                                                                             dplyr::summarize(Estimate = mean(value)))
  pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(Variable) %>%
                                                                                              dplyr::summarize(Estimate = mean(value)))
  }else if(coef == 'bootstrapped_median'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(Variable) %>%
                                                                                      dplyr::summarize(Estimate = median(value)))
    zi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.ZC)) %>% dplyr::left_join(zi_boot %>% dplyr::group_by(Variable) %>%
                                                                                               dplyr::summarize(Estimate = median(value)))
    evi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(Variable) %>%
                                                                                                dplyr::summarize(Estimate = median(value)))
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(Variable) %>%
                                                                                          dplyr::summarize(Estimate = median(value)))
  }
  if(standard_error){
    nb <- nb %>% dplyr::left_join(nb_boot %>% dplyr::group_by(Variable) %>%
                                                                                      dplyr::summarize(se = sd(value)))
    zi <- zi %>% dplyr::left_join(zi_boot %>% dplyr::group_by(Variable) %>%
                                                                                               dplyr::summarize(se = sd(value)))
    evi <- evi %>% dplyr::left_join(evi_boot %>% dplyr::group_by(Variable) %>%
                                                                                                dplyr::summarize(se = sd(value)))
    pareto <- pareto %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(Variable) %>%
                                                                                          dplyr::summarize(se = sd(value)))
  }
  
  if(approx_t_value){
    nb <- nb %>% dplyr::mutate(approx_t = Estimate/se)
    zi <- zi %>% dplyr::mutate(approx_t = Estimate/se)
    evi <- evi %>% dplyr::mutate(approx_t = Estimate/se)
    pareto <- pareto %>% dplyr::mutate(approx_t = Estimate/se)
  }
  
  if(p_value %in% c('bootstrapped','both')){
    nb <- nb %>% dplyr::left_join(dplyr::left_join(nb_boot,nb) %>% dplyr::group_by(Variable) %>%
                                    dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(value,Estimate[1],symmetric=symmetric_bootstrap_p)))
    zi <- zi %>% dplyr::left_join(dplyr::left_join(zi_boot,zi) %>% dplyr::group_by(Variable) %>%
                                    dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(value,Estimate[1],symmetric=symmetric_bootstrap_p)))
    evi <- evi %>% dplyr::left_join(dplyr::left_join(evi_boot,evi) %>% dplyr::group_by(Variable) %>%
                                      dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(value,Estimate[1],symmetric=symmetric_bootstrap_p)))
    pareto <- pareto %>% dplyr::left_join(dplyr::left_join(pareto_boot,pareto) %>% dplyr::group_by(Variable) %>%
                                            dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(value,Estimate[1],symmetric=symmetric_bootstrap_p)))
    
  }else if(p_value %in% c('approx','both')){
    nb <- nb %>% dplyr::mutate(approx_p = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    zi <- zi %>% dplyr::mutate(approx_p = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    evi <- evi %>% dplyr::mutate(approx_p = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    pareto <- pareto %>% dplyr::mutate(approx_p = 2*(1-pt(abs(approx_t),df = nobs-npar)))
  }

  res <- list(coefficients = list(negative_binomial = nb, 
                                    zero_inflation = zi,
                                    extreme_value_inflation = evi,
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
#' @param object an EVZINB object with bootstraps
#' @param boot_mean Logical: Should the bootstrapped mean of the estimates be computed?
#' @param boot_se Logical: Should the bootstrapped standard errors of the estimates be computed?
#' @param p_approx Logical: Should the normal approximated p-value be computed? Not recommended unless the number of bootstraps is too low for bootstrapped p-values
#'
#' @return An EVZINB summary object
#' @export
#'
#' @examples summary(evinb)
summary.evinb <- function(object,coef = c('original','bootstrapped_mean','bootstrapped_median'),standard_error = TRUE, p_value = c('bootstrapped','approx','both','none'), bootstrapped_props = c('none','mean','median'),approx_t_value = TRUE,
                          symmetric_bootstrap_p = TRUE,...){
  
  coef <- match.arg(coef,c('original','bootstrapped_mean','bootstrapped_median'))
  
  p_value <- match.arg(p_value, c('bootstrapped','approx','both','none'))
  
  bootstrapped_props <- match.arg(bootstrapped_props,c('none','mean','median'))
  
  
  nobs <- nrow(object$data$x.nb)
  npar <- length(object$par.all)
  
  props <- object$props %>% dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(everything(),names_to = 'state') %>%
    dplyr::filter(state != 'zero') %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(mean_prop = mean(value)) %>%
  alpha_nb <- c(object$coef$Alpha.NB)
  names(alpha_nb) <- c('Alpha_NB')
  C_est <- c(object$coef$C)
  names(C_est) <- c('C')
  
  if(!is.null(object$bootstraps)){
    n_bootstraps_org <- length(object$bootstraps)
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(object$bootstraps)
    
    nb_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable')
    
    
    evi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') 
    
    pareto_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') 
    
    props_boot <-  object$bootstraps %>% purrr::map('props') %>% purrr::map(colMeans) %>% purrr::reduce(rbind) %>%
      dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(everything(),names_to = 'state') %>%
      dplyr::filter(state != 'zero') %>%
      dplyr::group_by(state) %>%
      dplyr::summarize(bootstrap_mean = mean(value),bootstrap_median = median(value), standard_error = sd(value))
    
    if(bootstrapped_props == 'median'){
      props_boot <- props_boot %>% dplyr::select(state,bootstrap_median,standard_error)
    }else if(bootstrapped_props == 'mean'){
      props_boot <- props_boot %>% dplyr::select(state,bootstrap_mean,standard_error)
    }else{
      props_boot <- props_boot %>% dplyr::select(state,standard_error)
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

    evi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL),Estimate = object$coef$Beta.multinom.PL)
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL),Estimate = object$coef$Beta.PL)
  }else if(coef == 'bootstrapped_mean'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(Variable) %>%
                                                                                      dplyr::summarize(Estimate = mean(value)))

    evi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(Variable) %>%
                                                                                                dplyr::summarize(Estimate = mean(value)))
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(Variable) %>%
                                                                                          dplyr::summarize(Estimate = mean(value)))
  }else if(coef == 'bootstrapped_median'){
    nb <- dplyr::tibble(Variable = names(object$coef$Beta.NB)) %>% dplyr::left_join(nb_boot %>% dplyr::group_by(Variable) %>%
                                                                                      dplyr::summarize(Estimate = median(value)))

    evi <- dplyr::tibble(Variable = names(object$coef$Beta.multinom.PL)) %>% dplyr::left_join(evi_boot %>% dplyr::group_by(Variable) %>%
                                                                                                dplyr::summarize(Estimate = median(value)))
    pareto <- dplyr::tibble(Variable = names(object$coef$Beta.PL)) %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(Variable) %>%
                                                                                          dplyr::summarize(Estimate = median(value)))
  }
  if(standard_error){
    nb <- nb %>% dplyr::left_join(nb_boot %>% dplyr::group_by(Variable) %>%
                                    dplyr::summarize(se = sd(value)))

    evi <- evi %>% dplyr::left_join(evi_boot %>% dplyr::group_by(Variable) %>%
                                      dplyr::summarize(se = sd(value)))
    pareto <- pareto %>% dplyr::left_join(pareto_boot %>% dplyr::group_by(Variable) %>%
                                            dplyr::summarize(se = sd(value)))
  }
  
  if(approx_t_value){
    nb <- nb %>% dplyr::mutate(approx_t = Estimate/se)
    evi <- evi %>% dplyr::mutate(approx_t = Estimate/se)
    pareto <- pareto %>% dplyr::mutate(approx_t = Estimate/se)
  }
  
  if(p_value %in% c('bootstrapped','both')){
    nb <- nb %>% dplyr::left_join(dplyr::left_join(nb_boot,nb) %>% dplyr::group_by(Variable) %>%
                                    dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(value,Estimate[1],symmetric=symmetric_bootstrap_p)))
    
    
    evi <- evi %>% dplyr::left_join(dplyr::left_join(evi_boot,evi) %>% dplyr::group_by(Variable) %>%
                                      dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(value,Estimate[1],symmetric=symmetric_bootstrap_p)))
    pareto <- pareto %>% dplyr::left_join(dplyr::left_join(pareto_boot,pareto) %>% dplyr::group_by(Variable) %>%
                                            dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(value,Estimate[1],symmetric=symmetric_bootstrap_p)))
    
  }else if(p_value %in% c('approx','both')){
    nb <- nb %>% dplyr::mutate(approx_p = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    evi <- evi %>% dplyr::mutate(approx_p = 2*(1-pt(abs(approx_t),df = nobs-npar)))
    pareto <- pareto %>% dplyr::mutate(approx_p = 2*(1-pt(abs(approx_t),df = nobs-npar)))
  }
  
  res <- list(coefficients = list(negative_binomial = nb, 
                                  extreme_value_inflation = evi,
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



