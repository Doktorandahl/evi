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
#' @examples summary(evinb)
summary.evzinb <- function(object,boot_mean=F,boot_se=F,p_approx=F,...){


  nb <- dplyr::tibble(Variable = names(object$par.mat$Beta.NB),Estimate = object$par.mat$Beta.NB)
  zi <- dplyr::tibble(Variable = names(object$par.mat$Beta.multinom.ZC),Estimate = object$par.mat$Beta.multinom.ZC)
  evi <- dplyr::tibble(Variable = names(object$par.mat$Beta.multinom.PL),Estimate = object$par.mat$Beta.multinom.PL)
  pareto <- dplyr::tibble(Variable = names(object$par.mat$Beta.PL),Estimate = object$par.mat$Beta.PL)

  nobs <- nrow(object$x.nb)
  npar <- length(object$par.all)
  
  mean_props <- object$par.mat$Props %>% dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(everything(),names_to = 'state') %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(mean_prop = mean(value)) %>%
    dplyr::slice(c(3,1,2))
  alpha_nb <- c(object$par.mat$Alpha.NB)
  names(alpha_nb) <- c('Alpha_NB')
  C_est <- c(object$par.mat$C)
  names(C_est) <- c('C')
  
  if(!is.null(object$bootstraps)){
    n_bootstraps_org <- length(object$bootstraps)
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(object$bootstraps)
    
  nb_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarize(boot_mean = mean(value),
              boot_median = median(value),
              boot_se = sd(value),
              p_positive = mean(value>0))

  zi_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Beta.multinom.ZC') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarize(boot_mean = mean(value),
              boot_median = median(value),
              boot_se = sd(value),
              p_positive = mean(value>0))

  evi_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarize(boot_mean = mean(value),
              boot_median = median(value),
              boot_se = sd(value),
              p_positive = mean(value>0))

  pareto_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarize(boot_mean = mean(value),
              boot_median = median(value),
              boot_se = sd(value),
              p_positive = mean(value>0))

  nb <- dplyr::left_join(nb,nb_boot) %>%
    dplyr::mutate(p_boot = dplyr::case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                              sign(Estimate)==-1 ~ p_positive/2),
           p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
           sig_boot = dplyr::case_when(p_boot < 0.001 ~ '***',
                                p_boot < 0.01 ~ '**',
                                p_boot < 0.05 ~ '*',
                                p_boot < 0.1 ~ '.',
                                T ~ ' '),
           sig_approx = dplyr::case_when(p_approx < 0.001 ~ '***',
                                  p_approx < 0.01 ~ '**',
                                  p_approx < 0.05 ~ '*',
                                  p_approx < 0.1 ~ '.',
                                  T ~ ' ')) %>%
    dplyr::select(-p_positive)


  zi <- dplyr::left_join(zi,zi_boot) %>%
    dplyr::mutate(p_boot = dplyr::case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                              sign(Estimate)==-1 ~ p_positive/2),
           p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
           sig_boot = dplyr::case_when(p_boot < 0.001 ~ '***',
                                p_boot < 0.01 ~ '**',
                                p_boot < 0.05 ~ '*',
                                p_boot < 0.1 ~ '.',
                                T ~ ' '),
           sig_approx = dplyr::case_when(p_approx < 0.001 ~ '***',
                                  p_approx < 0.01 ~ '**',
                                  p_approx < 0.05 ~ '*',
                                  p_approx < 0.1 ~ '.',
                                  T ~ ' ')) %>%
    dplyr::select(-p_positive)

  evi <- dplyr::left_join(evi,evi_boot) %>%
    dplyr::mutate(p_boot = dplyr::case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                              sign(Estimate)==-1 ~ p_positive/2),
           p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
           sig_boot = dplyr::case_when(p_boot < 0.001 ~ '***',
                                p_boot < 0.01 ~ '**',
                                p_boot < 0.05 ~ '*',
                                p_boot < 0.1 ~ '.',
                                T ~ ' '),
           sig_approx = dplyr::case_when(p_approx < 0.001 ~ '***',
                                  p_approx < 0.01 ~ '**',
                                  p_approx < 0.05 ~ '*',
                                  p_approx < 0.1 ~ '.',
                                  T ~ ' ')) %>%
    dplyr::select(-p_positive)

  pareto <- dplyr::left_join(pareto,pareto_boot) %>%
    dplyr::mutate(p_boot = dplyr::case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                              sign(Estimate)==-1 ~ p_positive/2),
           p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
           sig_boot = dplyr::case_when(p_boot < 0.001 ~ '***',
                                p_boot < 0.01 ~ '**',
                                p_boot < 0.05 ~ '*',
                                p_boot < 0.1 ~ '.',
                                T ~ ' '),
           sig_approx = dplyr::case_when(p_approx < 0.001 ~ '***',
                                  p_approx < 0.01 ~ '**',
                                  p_approx < 0.05 ~ '*',
                                  p_approx < 0.1 ~ '.',
                                  T ~ ' ')) %>%
    dplyr::select(-p_positive)
  
  
  mean_props_boot <-  object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Props') %>% purrr::map(colMeans) %>% purrr::reduce(rbind) %>%
    dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(everything(),names_to = 'state') %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(boot_mean_prop = mean(value),boot_median_prop = median(value), boot_se_prop = sd(value))
  
  mean_props <- dplyr::left_join(mean_props,mean_props_boot)

  alpha_nb_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Alpha.NB') %>% purrr::reduce(c)
  
  alpha_nb <- c(alpha_nb, mean(alpha_nb_boot),median(alpha_nb_boot),sd(alpha_nb_boot))
  names(alpha_nb) <- c('Alpha_NB','boot_mean_Alpha_NB','boot_median_Alpha_NB','boot_se_Alpha_NB')

  C_est_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('C')%>% purrr::reduce(c)
  C_est <- c(C_est, mean(C_est_boot),median(C_est_boot),sd(C_est_boot))
  names(C_est) <- c('C','boot_mean_C','boot_median_C','boot_se_C')
  
  if(!boot_mean){
    nb <- nb %>% dplyr::select(-boot_mean)
    zi <- zi %>% dplyr::select(-boot_mean)
    evi <- evi %>% dplyr::select(-boot_mean)
    pareto <- pareto %>% dplyr::select(-boot_mean)
    mean_props<- mean_props %>% dplyr::select(-boot_mean_prop)
  }
  
  if(!boot_se){
    nb <- nb %>% dplyr::select(-boot_se)
    zi <- zi %>% dplyr::select(-boot_se)
    evi <- evi %>% dplyr::select(-boot_se)
    pareto <- pareto %>% dplyr::select(-boot_se)
    mean_props<- mean_props %>% dplyr::select(-boot_se_prop)
  }
  if(!p_approx){
    nb <- nb %>% dplyr::select(-c(p_approx,sig_approx))
    zi <- zi %>% dplyr::select(-c(p_approx,sig_approx))
    evi <- evi %>% dplyr::select(-c(p_approx,sig_approx))
    pareto <- pareto %>% dplyr::select(-c(p_approx,sig_approx))
  }
  
}




  
  if(!is.null(object$bootstraps)){
    res <- list(negative_binomial = list(summary = nb, alpha_nb = alpha_nb),
                pareto = list(summary = pareto, C = C_est),
                zero_inflation = list(summary = zi, props = mean_props),
                extreme_value_inflation = list(summary = evi, props = mean_props),
                n_failed_bootstraps = n_failed_bootstraps)
  }else{
    res <- list(negative_binomial = list(summary = nb, alpha_nb = alpha_nb),
                pareto = list(summary = pareto, C = C_est),
                zero_inflation = list(summary = zi, props = mean_props),
                extreme_value_inflation = list(summary = evi, props = mean_props))
  }
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
summary.evinb <- function(object,boot_mean=F,boot_se=F,p_approx=F,...){
  
  
  nb <- dplyr::tibble(Variable = names(object$par.mat$Beta.NB),Estimate = object$par.mat$Beta.NB)
  zi <- dplyr::tibble(Variable = names(object$par.mat$Beta.multinom.ZC),Estimate = object$par.mat$Beta.multinom.ZC)
  evi <- dplyr::tibble(Variable = names(object$par.mat$Beta.multinom.PL),Estimate = object$par.mat$Beta.multinom.PL)
  pareto <- dplyr::tibble(Variable = names(object$par.mat$Beta.PL),Estimate = object$par.mat$Beta.PL)
  
  nobs <- nrow(object$x.nb)
  npar <- length(object$par.all)
  
  mean_props <- object$par.mat$Props %>% dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(everything(),names_to = 'state') %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(mean_prop = mean(value)) %>%
    dplyr::slice(c(3,1,2))
  alpha_nb <- c(object$par.mat$Alpha.NB)
  names(alpha_nb) <- c('Alpha_NB')
  C_est <- c(object$par.mat$C)
  names(C_est) <- c('C')
  n_failed_bootstraps <- NULL
  
  if(!is.null(object$bootstraps)){
    n_bootstraps_org <- length(object$bootstraps)
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org-length(object$bootstraps)
    
    nb_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarize(boot_mean = mean(value),
                       boot_median = median(value),
                       boot_se = sd(value),
                       p_positive = mean(value>0))
    
    zi_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Beta.multinom.ZC') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarize(boot_mean = mean(value),
                       boot_median = median(value),
                       boot_se = sd(value),
                       p_positive = mean(value>0))
    
    evi_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarize(boot_mean = mean(value),
                       boot_median = median(value),
                       boot_se = sd(value),
                       p_positive = mean(value>0))
    
    pareto_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() %>% tidyr::pivot_longer(everything(),names_to = 'Variable') %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarize(boot_mean = mean(value),
                       boot_median = median(value),
                       boot_se = sd(value),
                       p_positive = mean(value>0))
    
    nb <- dplyr::left_join(nb,nb_boot) %>%
      dplyr::mutate(p_boot = dplyr::case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                                              sign(Estimate)==-1 ~ p_positive/2),
                    p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
                    sig_boot = dplyr::case_when(p_boot < 0.001 ~ '***',
                                                p_boot < 0.01 ~ '**',
                                                p_boot < 0.05 ~ '*',
                                                p_boot < 0.1 ~ '.',
                                                T ~ ' '),
                    sig_approx = dplyr::case_when(p_approx < 0.001 ~ '***',
                                                  p_approx < 0.01 ~ '**',
                                                  p_approx < 0.05 ~ '*',
                                                  p_approx < 0.1 ~ '.',
                                                  T ~ ' ')) %>%
      dplyr::select(-p_positive)
    
    
    zi <- dplyr::left_join(zi,zi_boot) %>%
      dplyr::mutate(p_boot = dplyr::case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                                              sign(Estimate)==-1 ~ p_positive/2),
                    p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
                    sig_boot = dplyr::case_when(p_boot < 0.001 ~ '***',
                                                p_boot < 0.01 ~ '**',
                                                p_boot < 0.05 ~ '*',
                                                p_boot < 0.1 ~ '.',
                                                T ~ ' '),
                    sig_approx = dplyr::case_when(p_approx < 0.001 ~ '***',
                                                  p_approx < 0.01 ~ '**',
                                                  p_approx < 0.05 ~ '*',
                                                  p_approx < 0.1 ~ '.',
                                                  T ~ ' ')) %>%
      dplyr::select(-p_positive)
    
    evi <- dplyr::left_join(evi,evi_boot) %>%
      dplyr::mutate(p_boot = dplyr::case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                                              sign(Estimate)==-1 ~ p_positive/2),
                    p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
                    sig_boot = dplyr::case_when(p_boot < 0.001 ~ '***',
                                                p_boot < 0.01 ~ '**',
                                                p_boot < 0.05 ~ '*',
                                                p_boot < 0.1 ~ '.',
                                                T ~ ' '),
                    sig_approx = dplyr::case_when(p_approx < 0.001 ~ '***',
                                                  p_approx < 0.01 ~ '**',
                                                  p_approx < 0.05 ~ '*',
                                                  p_approx < 0.1 ~ '.',
                                                  T ~ ' ')) %>%
      dplyr::select(-p_positive)
    
    pareto <- dplyr::left_join(pareto,pareto_boot) %>%
      dplyr::mutate(p_boot = dplyr::case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                                              sign(Estimate)==-1 ~ p_positive/2),
                    p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
                    sig_boot = dplyr::case_when(p_boot < 0.001 ~ '***',
                                                p_boot < 0.01 ~ '**',
                                                p_boot < 0.05 ~ '*',
                                                p_boot < 0.1 ~ '.',
                                                T ~ ' '),
                    sig_approx = dplyr::case_when(p_approx < 0.001 ~ '***',
                                                  p_approx < 0.01 ~ '**',
                                                  p_approx < 0.05 ~ '*',
                                                  p_approx < 0.1 ~ '.',
                                                  T ~ ' ')) %>%
      dplyr::select(-p_positive)
    
    
    mean_props_boot <-  object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Props') %>% purrr::map(colMeans) %>% purrr::reduce(rbind) %>%
      dplyr::as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% tidyr::pivot_longer(everything(),names_to = 'state') %>%
      dplyr::group_by(state) %>%
      dplyr::summarize(boot_mean_prop = mean(value),boot_median_prop = median(value), boot_se_prop = sd(value))
    
    mean_props <- dplyr::left_join(mean_props,mean_props_boot)
    
    alpha_nb_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('Alpha.NB') %>% purrr::reduce(c)
    
    alpha_nb <- c(alpha_nb, mean(alpha_nb_boot),median(alpha_nb_boot),sd(alpha_nb_boot))
    names(alpha_nb) <- c('Alpha_NB','boot_mean_Alpha_NB','boot_median_Alpha_NB','boot_se_Alpha_NB')
    
    C_est_boot <- object$bootstraps %>% purrr::map('par.mat') %>% purrr::map('C')%>% purrr::reduce(c)
    C_est <- c(C_est, mean(C_est_boot),median(C_est_boot),sd(C_est_boot))
    names(C_est) <- c('C','boot_mean_C','boot_median_C','boot_se_C')
    
    if(!boot_mean){
      nb <- nb %>% dplyr::select(-boot_mean)
      zi <- zi %>% dplyr::select(-boot_mean)
      evi <- evi %>% dplyr::select(-boot_mean)
      pareto <- pareto %>% dplyr::select(-boot_mean)
      mean_props<- mean_props %>% dplyr::select(-boot_mean_prop)
    }
    
    if(!boot_se){
      nb <- nb %>% dplyr::select(-boot_se)
      zi <- zi %>% dplyr::select(-boot_se)
      evi <- evi %>% dplyr::select(-boot_se)
      pareto <- pareto %>% dplyr::select(-boot_se)
      mean_props<- mean_props %>% dplyr::select(-boot_se_prop)
    }
    if(!p_approx){
      nb <- nb %>% dplyr::select(-c(p_approx,sig_approx))
      zi <- zi %>% dplyr::select(-c(p_approx,sig_approx))
      evi <- evi %>% dplyr::select(-c(p_approx,sig_approx))
      pareto <- pareto %>% dplyr::select(-c(p_approx,sig_approx))
    }
    
  }
  
  
  
  if(!is.null(object$bootstraps)){
  res <- list(negative_binomial = list(summary = nb, alpha_nb = alpha_nb),
              pareto = list(summary = pareto, C = C_est),
              zero_inflation = list(summary = zi, props = mean_props),
              extreme_value_inflation = list(summary = evi, props = mean_props),
              n_failed_bootstraps = n_failed_bootstraps)
  }else{
    res <- list(negative_binomial = list(summary = nb, alpha_nb = alpha_nb),
                pareto = list(summary = pareto, C = C_est),
                zero_inflation = list(summary = zi, props = mean_props),
                extreme_value_inflation = list(summary = evi, props = mean_props))
  }
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


