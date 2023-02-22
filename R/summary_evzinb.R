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


  nb <- tibble(Variable = names(object$par.mat$Beta.NB),Estimate = object$par.mat$Beta.NB)
  zi <- tibble(Variable = names(object$par.mat$Beta.multinom.ZC),Estimate = object$par.mat$Beta.multinom.ZC)
  evi <- tibble(Variable = names(object$par.mat$Beta.multinom.PL),Estimate = object$par.mat$Beta.multinom.PL)
  pareto <- tibble(Variable = names(object$par.mat$Beta.PL),Estimate = object$par.mat$Beta.PL)

  nobs <- nrow(object$x.nb)
  npar <- length(object$par.all)
  nb_boot <- object$bootstraps %>% map('par.mat') %>% map('Beta.NB') %>% bind_rows() %>% pivot_longer(everything(),names_to = 'Variable') %>%
    group_by(Variable) %>%
    summarize(boot_mean = mean(value),
              boot_median = median(value),
              boot_se = sd(value),
              p_positive = mean(value>0))

  zi_boot <- object$bootstraps %>% map('par.mat') %>% map('Beta.multinom.ZC') %>% bind_rows() %>% pivot_longer(everything(),names_to = 'Variable') %>%
    group_by(Variable) %>%
    summarize(boot_mean = mean(value),
              boot_median = median(value),
              boot_se = sd(value),
              p_positive = mean(value>0))

  evi_boot <- object$bootstraps %>% map('par.mat') %>% map('Beta.multinom.PL') %>% bind_rows() %>% pivot_longer(everything(),names_to = 'Variable') %>%
    group_by(Variable) %>%
    summarize(boot_mean = mean(value),
              boot_median = median(value),
              boot_se = sd(value),
              p_positive = mean(value>0))

  pareto_boot <- object$bootstraps %>% map('par.mat') %>% map('Beta.PL') %>% bind_rows() %>% pivot_longer(everything(),names_to = 'Variable') %>%
    group_by(Variable) %>%
    summarize(boot_mean = mean(value),
              boot_median = median(value),
              boot_se = sd(value),
              p_positive = mean(value>0))

  nb <- left_join(nb,nb_boot) %>%
    mutate(p_boot = case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                              sign(Estimate)==-1 ~ p_positive/2),
           p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
           sig_boot = case_when(p_boot < 0.001 ~ '***',
                                p_boot < 0.01 ~ '**',
                                p_boot < 0.05 ~ '*',
                                p_boot < 0.1 ~ '.',
                                T ~ ' '),
           sig_approx = case_when(p_approx < 0.001 ~ '***',
                                  p_approx < 0.01 ~ '**',
                                  p_approx < 0.05 ~ '*',
                                  p_approx < 0.1 ~ '.',
                                  T ~ ' ')) %>%
    select(-p_positive)


  zi <- left_join(zi,zi_boot) %>%
    mutate(p_boot = case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                              sign(Estimate)==-1 ~ p_positive/2),
           p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
           sig_boot = case_when(p_boot < 0.001 ~ '***',
                                p_boot < 0.01 ~ '**',
                                p_boot < 0.05 ~ '*',
                                p_boot < 0.1 ~ '.',
                                T ~ ' '),
           sig_approx = case_when(p_approx < 0.001 ~ '***',
                                  p_approx < 0.01 ~ '**',
                                  p_approx < 0.05 ~ '*',
                                  p_approx < 0.1 ~ '.',
                                  T ~ ' ')) %>%
    select(-p_positive)

  evi <- left_join(evi,evi_boot) %>%
    mutate(p_boot = case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                              sign(Estimate)==-1 ~ p_positive/2),
           p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
           sig_boot = case_when(p_boot < 0.001 ~ '***',
                                p_boot < 0.01 ~ '**',
                                p_boot < 0.05 ~ '*',
                                p_boot < 0.1 ~ '.',
                                T ~ ' '),
           sig_approx = case_when(p_approx < 0.001 ~ '***',
                                  p_approx < 0.01 ~ '**',
                                  p_approx < 0.05 ~ '*',
                                  p_approx < 0.1 ~ '.',
                                  T ~ ' ')) %>%
    select(-p_positive)

  pareto <- left_join(pareto,pareto_boot) %>%
    mutate(p_boot = case_when(sign(Estimate)==1 ~ (1-p_positive)/2,
                              sign(Estimate)==-1 ~ p_positive/2),
           p_approx = (1-pt(abs(Estimate/boot_se),df=nobs-npar))*2,
           sig_boot = case_when(p_boot < 0.001 ~ '***',
                                p_boot < 0.01 ~ '**',
                                p_boot < 0.05 ~ '*',
                                p_boot < 0.1 ~ '.',
                                T ~ ' '),
           sig_approx = case_when(p_approx < 0.001 ~ '***',
                                  p_approx < 0.01 ~ '**',
                                  p_approx < 0.05 ~ '*',
                                  p_approx < 0.1 ~ '.',
                                  T ~ ' ')) %>%
    select(-p_positive)

  alpha_nb_boot <- object$bootstraps %>% map('par.mat') %>% map('Alpha.NB') %>% reduce(c)
  alpha_nb <- c(object$par.mat$Alpha.NB, mean(alpha_nb_boot),median(alpha_nb_boot),sd(alpha_nb_boot))
  names(alpha_nb) <- c('Alpha_NB','boot_mean_Alpha_NB','boot_median_Alpha_NB','boot_se_Alpha_NB')

  C_est_boot <- object$bootstraps %>% map('par.mat') %>% map('C')%>% reduce(c)
  C_est <- c(object$par.mat$C, mean(C_est_boot),median(C_est_boot),sd(C_est_boot))
  names(C_est) <- c('C','boot_mean_C','boot_median_C','boot_se_C')

  mean_props <- object$par.mat$Props %>% as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% pivot_longer(everything(),names_to = 'state') %>%
    group_by(state) %>%
    summarize(mean_prop = mean(value)) %>%
    slice(c(3,1,2))

  mean_props_boot <-  object$bootstraps %>% map('par.mat') %>% map('Props') %>% map(colMeans) %>% reduce(rbind) %>%
    as_tibble(.name_repair = ~c('zero','negative_binomial','pareto')) %>% pivot_longer(everything(),names_to = 'state') %>%
    group_by(state) %>%
    summarize(boot_mean_prop = mean(value),boot_median_prop = median(value), boot_se_prop = sd(value))

  mean_props <- left_join(mean_props,mean_props_boot)


  if(!boot_mean){
    nb <- nb %>% select(-boot_mean)
    zi <- zi %>% select(-boot_mean)
    evi <- evi %>% select(-boot_mean)
    pareto <- pareto %>% select(-boot_mean)
    mean_props<- mean_props %>% select(-boot_mean_prop)
  }

  if(!boot_se){
    nb <- nb %>% select(-boot_se)
    zi <- zi %>% select(-boot_se)
    evi <- evi %>% select(-boot_se)
    pareto <- pareto %>% select(-boot_se)
    mean_props<- mean_props %>% select(-boot_se_prop)
  }
  if(!p_approx){
    nb <- nb %>% select(-c(p_approx,sig_approx))
    zi <- zi %>% select(-c(p_approx,sig_approx))
    evi <- evi %>% select(-c(p_approx,sig_approx))
    pareto <- pareto %>% select(-c(p_approx,sig_approx))
  }

  res <- list(negative_binomial = list(summary = nb, alpha_nb = alpha_nb),
              pareto = list(summary = pareto, C = C_est),
              zero_inflation = list(summary = zi, props = mean_props),
              extreme_value_inflation = list(summary = evi, props = mean_props))
  class(res) <- 'summary.evzinb'

  return(res)


}
