harmonic_calc <- function(pr_count,count,pr_pareto,C,pareto_alpha){
  pr_count*count + pr_pareto*C*(1+pareto_alpha)/pareto_alpha
}

explog_calc <- function(pr_count,count,pr_pareto,C,pareto_alpha){
  pr_count*count + pr_pareto*exp(1/pareto_alpha)
}


#' Predictions from evzinb object
#'
#' @param object An evzinb object for which to produce predicted values
#' @param newdata Optional new data (tibble) to produce predicted values from
#' @param type Character string, 'harmonic' for the harmonic mean and 'explog' for exponentiated expected log, 'counts' for predicted count of the negative binomial component, 'pareto_alpha' for the predicted pareto alpha value, 'states' for the predicted component states (prior), 'count_state' for predicted probability of the count state, 'evinf' for predicted probability of the pareto state,'zi' for the predicted probability of the zero state, 'all' for all predicted values, and 'quantile' for quantile prediction.
#' @param ... Other arguments passed to predict function
#' @param quantile Quantile for which to produce quantile prediction
#' @param multicore Should multicore be used when calculating quantile prediction? Often it is enough to run quantile prediction on a single core, but in cases of large data or very skewed distributions it may be useful to run multicore
#' @param ncores Number of cores to be used for multicore.
#' @param pred Type of prediction to be used, defaults to the original prediction from the fitted model, with alternatives being the bootstrapped median or mean. Note that bootstrap mean may yield infinite values, especially when doing quantile prediction
#' @param confint Should confidence intervals be made for the predictions? Note: only available for vector type predictions and not 'states' and 'all'.
#' @param conf_level What confidence level should be used for confidence intervals
#'
#' @return A vector of predicted values for type 'harmonic', 'explog', 'counts', 'pareto_alpha','zi','evinf', 'count_state', and 'quantile' or a tibble of predicted values for type 'states' and 'all' or if confint=T
#' @export
#'
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' predict(model)
#' predict(model, type='all') # Getting all of the available predicted values
predict.evzinb <- function(object,newdata=NULL, type = c('harmonic','explog','counts','pareto_alpha','zi','evinf','count_state','states','all', 'quantile'), pred = c('original','bootstrap_median','bootstrap_mean'),quantile=NULL,confint=F, conf_level=0.9, multicore = F,ncores=NULL,...){
  
  pred <- match.arg(pred, c('original','bootstrap_median','bootstrap_mean'))
  
  type <- match.arg(type,c('harmonic','explog','counts','pareto_alpha','zi','evinf','count_state','states','all', 'quantile'))
  
  if(type %in% c('states','all') & confint){
    stop('Confidence interval prediction only available for vector outputs')
  }
  
  if(pred %in% c('bootstrap_median','bootstrap_mean') | confint){
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))    
    nboots <- length(object$bootstraps)
    if(is.null(newdata)){
      newdata <- object$data$data
    }
    prbs_boot <- foreach::foreach(i = 1:nboots) %do%
      dplyr::bind_cols(prob_from_evzinb(object$bootstraps[[i]],
                       newdata = newdata),tibble::tibble(id=1:nrow(newdata)))
    cnts_boot <- foreach::foreach(i = 1:nboots) %do%
      dplyr::bind_cols(counts_from_evzinb(object$bootstraps[[i]],
                       newdata = newdata),tibble::tibble(id=1:nrow(newdata)))
    alphs_boot <- foreach::foreach(i = 1:nboots) %do%
      dplyr::bind_cols(fitted_alpha_from_evzinb(object$bootstraps[[i]],
                       newdata = newdata),tibble::tibble(id=1:nrow(newdata)))
    C_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('C') %>% purrr::reduce(c)
    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
      q_boot <- foreach::foreach(i=1:nboots) %do%
        tibble::tibble(q=quantiles_from_evzinb(object$bootstraps[[i]],quantile,newdata = newdata,return_data = F,multicore = multicore,ncores = ncores), id=1:nrow(newdata))
      }else{
        q_boot <- NULL
      }
    }else{
      q_boot <- NULL
    }
    harmonic_boot <- foreach::foreach(i = 1:nboots) %do%
      tibble::tibble(harmonic = harmonic_calc(prbs_boot[[i]]$pr_count,cnts_boot[[i]]$count,pr_pareto = prbs_boot[[i]]$pr_pareto,C = C_boot[i],pareto_alpha = alphs_boot[[i]]$pareto_alpha), id = 1:nrow(newdata))
    
    explog_boot <- foreach::foreach(i = 1:nboots) %do%
      tibble::tibble(explog = explog_calc(prbs_boot[[i]]$pr_count,cnts_boot[[i]]$count,pr_pareto = prbs_boot[[i]]$pr_pareto,C = C_boot[i],pareto_alpha = alphs_boot[[i]]$pareto_alpha), id = 1:nrow(newdata))
  }
  
  
  if(pred == 'original'){
  if(type %in% c('quantile','all')){
    if(type == 'quantile' & is.null(quantile)){
      stop('quantile must be provided for quantile prediction')
    }else if(!is.null(quantile)){
    q <- quantiles_from_evzinb(object,quantile,newdata = newdata,return_data = F,multicore = multicore,ncores = ncores)
  }else{
    q <- NULL
  }
  }
  ## Estimate component probabilities for all individuals
  prbs <- prob_from_evzinb(object,
                           newdata = newdata)
  ## Estimate mu_nb for all individuals
  cnts <- counts_from_evzinb(object,
                             newdata = newdata)
  ## Estimate pareto alpha for all individuals
  alphs <- fitted_alpha_from_evzinb(object,
                                    newdata = newdata)
  C_est <- object$coef$C
  # if(min(alphs)<1e-02){
  #   warning('Fitted pareto alpha-values below 1e-02 detected. Setting those alphas to 1e-02.')
  #   alphs <- alphs %>% mutate(pred_alpha = case_when(pred_alpha<1e-02 ~ 1e-02,
  #                                               T~pred_alpha))
  # }
  harmonic <- harmonic_calc(pr_count = prbs$pr_count,count = cnts$count,pr_pareto = prbs$pr_pareto,C = C_est,pareto_alpha = alphs$pareto_alpha)
  
  explog <- explog_calc(pr_count = prbs$pr_count,
                        count = cnts$count,
                        pr_pareto = prbs$pr_pareto,
                        C = C_est,
                        pareto_alpha = alphs$pareto_alpha)
  
  }else if(pred=='bootstrap_median'){
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-id)
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-id)
    alphs <- alphs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-id)
    
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
        dplyr::summarize_all(median) %>% dplyr::select(-id) %>% dplyr::pull(q)
    }else{
      q <- NULL
    }
    
    harmonic <- harmonic_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% dplyr::summarize(harmonic = median(harmonic)) %>% dplyr::pull(harmonic)
    
    explog <- explog_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% dplyr::summarize(explog = median(explog)) %>% dplyr::pull(explog)
    
    C_est <- median(C_boot)
  }else if(pred=='bootstrap_mean'){
    warning('Bootstrapped mean predictions are experimental and may yield infinite values')
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-id)
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-id)
    alphs <- alphs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-id)
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
        dplyr::summarize_all(mean) %>% dplyr::select(-id) %>% dplyr::pull(q)
    }else{
      q <- NULL
    }
    C_est <- mean(C_boot)
    harmonic <- harmonic_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% dplyr::summarize(harmonic = mean(harmonic)) %>% dplyr::pull(harmonic)
    
    explog <- explog_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% dplyr::summarize(explog = mean(explog)) %>% dplyr::pull(explog)
  }


    if(confint){
      if(type %in% c('states','all')){
        stop("Confidence interval prediction only available for vector predictions (not 'states' or 'all')")
      }
      qs <- c((1-conf_level)/2,1-(1-conf_level)/2)
      
      if(type == 'harmonic'){
      ci <- harmonic_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
        dplyr::summarize(ci_lb = quantile(harmonic,qs[1]),
                         ci_ub = quantile(harmonic,qs[2])) %>%
        dplyr::select(-id)
      return(dplyr::bind_cols(tibble::tibble(harmonic=harmonic),ci))
      }
      if(type == 'explog'){
        ci <- explog_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
          dplyr::summarize(ci_lb = quantile(explog,qs[1]),
                           ci_ub = quantile(explog,qs[2])) %>% dplyr::select(-id)
        return(dplyr::bind_cols(tibble::tibble(explog=explog),ci))
      }
      if(type == 'counts'){
        ci <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
          dplyr::summarize(ci_lb = quantile(count,qs[1]),
                           ci_ub = quantile(count,qs[2])) %>% dplyr::select(-id)
        return(dplyr::bind_cols(tibble::tibble(count=cnts$count),ci))
        
      }
      if(type == 'pareto_alpha'){
        ci <- alphs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
          dplyr::summarize(ci_lb = quantile(pareto_alpha,qs[1]),
                           ci_ub = quantile(pareto_alpha,qs[2])) %>% dplyr::select(-id)
        return(dplyr::bind_cols(tibble::tibble(pareto_alpha=alphs$pareto_alpha),ci))
      }
      if(type == 'zi'){
        ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
          dplyr::summarize(ci_lb = quantile(pr_zc,qs[1]),
                           ci_ub = quantile(pr_zc,qs[2])) %>% dplyr::select(-id)
        return(dplyr::bind_cols(tibble::tibble(pr_zc=prbs$pr_zc),ci))
      }
      if(type == 'evinf'){
        ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
          dplyr::summarize(ci_lb = quantile(pr_pareto,qs[1]),
                           ci_ub = quantile(pr_pareto,qs[2])) %>% dplyr::select(-id)
        return(dplyr::bind_cols(tibble::tibble(pr_pareto=prbs$pr_pareto),ci))
      }
      if(type == 'count_state'){
        ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
          dplyr::summarize(ci_lb = quantile(pr_count,qs[1]),
                           ci_ub = quantile(pr_count,qs[2])) %>% dplyr::select(-id)
        return(dplyr::bind_cols(tibble::tibble(pr_count=prbs$pr_count),ci))
      }
      if(type == 'quantile'){
        warning('Confidence interval prediction with Quantiles may yield unstable results')
        ci <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
          dplyr::summarize(ci_lb = quantile(q,qs[1]),
                           ci_ub = quantile(q,qs[2])) %>% dplyr::select(-id)
        q_name <- paste0('q',100*quantile)
        return(dplyr::bind_cols(tibble::tibble(!!q_name:=q),ci))
      }

    }else{
      if(type == "harmonic"){
        return(harmonic)
      }
      if(type == "explog"){
        return(explog)
      }
      if(type == "counts"){
        return(cnts$count)
      }
      if(type == "pareto_alpha"){
        return(alphs$pareto_alpha)
      }
      if(type == "zi"){
        return(prbs$pr_zc)
      }
      if(type == "evinf"){
        return(prbs$pr_pareto)
      }
      if(type == "count_state"){
        return(prbs$pr_count)
      }
      if(type == 'states'){
        return(prbs)
      }
      if(type == 'quantile'){
        return(q)
      }
      if(type == 'all'){
        q_name <- paste0('q',100*quantile)
        return(dplyr::bind_cols(tibble::tibble(harmonic = harmonic,
                              explog = explog,
                              !!q_name := q,
                              prbs,cnts,alphs)))
      }
    }
  
}

#' Predictions from evinb object
#'
#' @param object An evinb object for which to produce predicted values
#' @param newdata Optional new data (tibble) to produce predicted values from
#' @param type Character string, 'harmonic' for the harmonic mean and 'explog' for exponentiated expected log, 'counts' for predicted count of the negative binomial component, 'pareto_alpha' for the predicted pareto alpha value, 'states' for the predicted component states (prior), 'count_state' for predicted probability of the count state, 'evinf' for predicted probability of the pareto state, 'all' for all predicted values, and 'quantile' for quantile prediction.
#' @param ... Other arguments passed to predict function
#' @param quantile Quantile for which to produce quantile prediction
#' @param multicore Should multicore be used when calculating quantile prediction? Often it is enough to run quantile prediction on a single core, but in cases of large data or very skewed distributions it may be useful to run multicore
#' @param ncores Number of cores to be used for multicore.
#' @param pred Type of prediction to be used, defaults to the original prediction from the fitted model, with alternatives being the bootstrapped median or mean. Note that bootstrap mean may yield infinite values, especially when doing quantile prediction
#' @param confint Should confidence intervals be made for the predictions? Note: only available for vector type predictions and not 'states' and 'all'.
#' @param conf_level What confidence level should be used for confidence intervals
#'
#' @return A vector of predicted values for type 'harmonic', 'explog', 'counts', 'pareto_alpha','evinf', 'count_state', and 'quantile' or a tibble of predicted values for type 'states' and 'all' or if confint=T
#' @export
#'
#' @examples data(genevzinb)
#' model <- evinb(y~x1+x2+x3,data=genevzinb)
#' predict(model)
#' predict(model, type='all') # Getting all of the available predicted values
predict.evinb <- function(object,newdata=NULL, type = c('harmonic','explog','counts','pareto_alpha','evinf','count_state','states','all', 'quantile'), pred = c('original','bootstrap_median','bootstrap_mean'),quantile=NULL,confint=F, conf_level=0.9, multicore = F,ncores=NULL,...){
  
  pred <- match.arg(pred, c('original','bootstrap_median','bootstrap_mean'))
  
  type <- match.arg(type,c('harmonic','explog','counts','pareto_alpha','evinf','count_state','states','all', 'quantile'))
  
  if(type %in% c('states','all') & confint){
    stop('Confidence interval prediction only available for vector outputs')
  }
  
  if(pred %in% c('bootstrap_median','bootstrap_mean') | confint){
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))    
    
    nboots <- length(object$bootstraps)
    if(is.null(newdata)){
      newdata <- object$data$data
    }
    prbs_boot <- foreach::foreach(i = 1:nboots) %do%
      dplyr::bind_cols(prob_from_evinb(object$bootstraps[[i]],
                                        newdata = newdata),tibble::tibble(id=1:nrow(newdata)))
    cnts_boot <- foreach::foreach(i = 1:nboots) %do%
      dplyr::bind_cols(counts_from_evzinb(object$bootstraps[[i]],
                                          newdata = newdata),tibble::tibble(id=1:nrow(newdata)))
    alphs_boot <- foreach::foreach(i = 1:nboots) %do%
      dplyr::bind_cols(fitted_alpha_from_evzinb(object$bootstraps[[i]],
                                                newdata = newdata),tibble::tibble(id=1:nrow(newdata)))
    C_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('C') %>% purrr::reduce(c)
    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q_boot <- foreach::foreach(i=1:nboots) %do%
          tibble::tibble(q=quantiles_from_evinb(object$bootstraps[[i]],quantile,newdata = newdata,return_data = F,multicore = multicore,ncores = ncores), id=1:nrow(newdata))
      }else{
        q_boot <- NULL
      }
    }else{
      q_boot <- NULL
    }
    harmonic_boot <- foreach::foreach(i = 1:nboots) %do%
      tibble::tibble(harmonic = harmonic_calc(prbs_boot[[i]]$pr_count,cnts_boot[[i]]$count,pr_pareto = prbs_boot[[i]]$pr_pareto,C = C_boot[i],pareto_alpha = alphs_boot[[i]]$pareto_alpha), id = 1:nrow(newdata))
    
    explog_boot <- foreach::foreach(i = 1:nboots) %do%
      tibble::tibble(explog = explog_calc(prbs_boot[[i]]$pr_count,cnts_boot[[i]]$count,pr_pareto = prbs_boot[[i]]$pr_pareto,C = C_boot[i],pareto_alpha = alphs_boot[[i]]$pareto_alpha), id = 1:nrow(newdata))
  }
  
  
  if(pred == 'original'){
    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q <- quantiles_from_evinb(object,quantile,newdata = newdata,return_data = F,multicore = multicore,ncores = ncores)
      }else{
        q <- NULL
      }
    }
    ## Estimate component probabilities for all individuals
    prbs <- prob_from_evinb(object,
                             newdata = newdata)
    ## Estimate mu_nb for all individuals
    cnts <- counts_from_evzinb(object,
                               newdata = newdata)
    ## Estimate pareto alpha for all individuals
    alphs <- fitted_alpha_from_evzinb(object,
                                      newdata = newdata)
    C_est <- object$coef$C
    # if(min(alphs)<1e-02){
    #   warning('Fitted pareto alpha-values below 1e-02 detected. Setting those alphas to 1e-02.')
    #   alphs <- alphs %>% mutate(pred_alpha = case_when(pred_alpha<1e-02 ~ 1e-02,
    #                                               T~pred_alpha))
    # }
    harmonic <- harmonic_calc(pr_count = prbs$pr_count,count = cnts$count,pr_pareto = prbs$pr_pareto,C = C_est,pareto_alpha = alphs$pareto_alpha)
    
    explog <- explog_calc(pr_count = prbs$pr_count,count = cnts$count,pr_pareto = prbs$pr_pareto,C = C_est,pareto_alpha = alphs$pareto_alpha)
    
  }else if(pred=='bootstrap_median'){
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-id)
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-id)
    alphs <- alphs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-id)
    
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
        dplyr::summarize_all(median) %>% dplyr::select(-id) %>% dplyr::pull(q)
    }else{
      q <- NULL
    }
    
    harmonic <- harmonic_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% dplyr::summarize(harmonic = median(harmonic)) %>% dplyr::pull(harmonic)
    
    explog <- explog_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% dplyr::summarize(explog = median(explog)) %>% dplyr::pull(explog)
    
    C_est <- median(C_boot)
  }else if(pred=='bootstrap_mean'){
    warning('Bootstrapped mean predictions are experimental and may yield infinite values')
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-id)
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-id)
    alphs <- alphs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-id)
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>%
        dplyr::summarize_all(mean) %>% dplyr::select(-id) %>% dplyr::pull(q)
    }else{
      q <- NULL
    }
    C_est <- mean(C_boot)
    harmonic <- harmonic_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% dplyr::summarize(harmonic = mean(harmonic)) %>% dplyr::pull(harmonic)
    
    explog <- explog_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% dplyr::summarize(explog = mean(explog)) %>% dplyr::pull(explog)
  }
  
  
  if(confint){
    if(type %in% c('states','all')){
      stop("Confidence interval prediction only available for vector predictions (not 'states' or 'all')")
    }
    qs <- c((1-conf_level)/2,1-(1-conf_level)/2)
    
    if(type == 'harmonic'){
      ci <- harmonic_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
        dplyr::summarize(ci_lb = quantile(harmonic,qs[1]),
                         ci_ub = quantile(harmonic,qs[2])) %>%
        dplyr::select(-id)
      return(dplyr::bind_cols(tibble::tibble(harmonic=harmonic),ci))
    }
    if(type == 'explog'){
      ci <- explog_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
        dplyr::summarize(ci_lb = quantile(explog,qs[1]),
                         ci_ub = quantile(explog,qs[2])) %>% dplyr::select(-id)
      return(dplyr::bind_cols(tibble::tibble(explog=explog),ci))
    }
    if(type == 'counts'){
      ci <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
        dplyr::summarize(ci_lb = quantile(count,qs[1]),
                         ci_ub = quantile(count,qs[2])) %>% dplyr::select(-id)
      return(dplyr::bind_cols(tibble::tibble(count=cnts$count),ci))
      
    }
    if(type == 'pareto_alpha'){
      ci <- alphs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
        dplyr::summarize(ci_lb = quantile(pareto_alpha,qs[1]),
                         ci_ub = quantile(pareto_alpha,qs[2])) %>% dplyr::select(-id)
      return(dplyr::bind_cols(tibble::tibble(pareto_alpha=alphs$pareto_alpha),ci))
    }
    
    if(type == 'evinf'){
      ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
        dplyr::summarize(ci_lb = quantile(pr_pareto,qs[1]),
                         ci_ub = quantile(pr_pareto,qs[2])) %>% dplyr::select(-id)
      return(dplyr::bind_cols(tibble::tibble(pr_pareto=prbs$pr_pareto),ci))
    }
    if(type == 'count_state'){
      ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
        dplyr::summarize(ci_lb = quantile(pr_count,qs[1]),
                         ci_ub = quantile(pr_count,qs[2])) %>% dplyr::select(-id)
      return(dplyr::bind_cols(tibble::tibble(pr_count=prbs$pr_count),ci))
    }
    if(type == 'quantile'){
      warning('Confidence interval prediction with Quantiles may yield unstable results')
      ci <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(id) %>% 
        dplyr::summarize(ci_lb = quantile(q,qs[1]),
                         ci_ub = quantile(q,qs[2])) %>% dplyr::select(-id)
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(!!q_name:=q),ci))
    }
    
  }else{
    if(type == "harmonic"){
      return(harmonic)
    }
    if(type == "explog"){
      return(explog)
    }
    if(type == "counts"){
      return(cnts$count)
    }
    if(type == "pareto_alpha"){
      return(alphs$pareto_alpha)
    }
    
    if(type == "evinf"){
      return(prbs$pr_pareto)
    }
    if(type == "count_state"){
      return(prbs$pr_count)
    }
    if(type == 'states'){
      return(prbs)
    }
    if(type == 'quantile'){
      return(q)
    }
    if(type == 'all'){
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(harmonic = harmonic,
                                             explog = explog,
                                             !!q_name := q),
                              prbs,cnts,alphs))
    }
  }
  
}
