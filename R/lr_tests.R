


#' Likelihood ratio test for individual variables of evzinb
#'
#' @param object EVZINB or EVZINB object to perform likelihood ratio test on
#' @param vars Either a list of character vectors with variable names which to be restricted in the LR test or a character vector of variable names. If a list, each character vector of the list will be run separately, allowing for multiple variables to be restricted as once. If a character vector, parameter 'single' can be used to determine whether all variables in the vector should be restricted at once (single = FALSE) or if the variables should be restricted one by one (single = TRUE)
#' @param single Logical. Determining whether variables in 'vars' should be restricted individually (single = TRUE) or all at once (single = FALSE)
#' @param bootstrap Should LR tests be conducted on each bootstrapped sample or only on the original sample.
#' @param multicore Logical. Should the function be run in parallel?
#' @param ncores Number of cores to use if multicore = TRUE
#' @param verbose Logical. Should the function be verbose?
#'
#' @return A tibble with one row per performed LR test
#' @export
#'
#' @examples 
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 10)
#'  lr_test(model,'x1')
#'  
lr_test <- function(object, vars, single = TRUE, bootstrap = FALSE, multicore = FALSE, ncores = NULL,verbose = FALSE){
  i <- 'temp_iter'
  j <- 'temp_iter'
  
  if(!is(vars, 'list')){
    if(single){
      formulas_dfs <- foreach::foreach(i = 1:length(vars)) %do%
        formula_var_remover(object$formulas,vars[i])
    }else{
      formulas_dfs <- foreach::foreach(i = 1:1) %do%
        formula_var_remover(object$formulas,vars)
      vars <- paste(vars,collapse ='_')
    }
  }else{
    formulas_dfs <- foreach::foreach(i = 1:length(vars)) %do%
      formula_var_remover(object$formulas,vars[[i]])
    vars <- vars  %>% purrr::map(~paste(.x,collapse = '_')) %>% purrr::reduce(c)
  }
    
  if(inherits(object,'evzinb')){
    reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %do%
      evzinb(formulas_dfs[[i]]$formulas$nb,
             formulas_dfs[[i]]$formulas$zi,
             formulas_dfs[[i]]$formulas$evinf,
             formulas_dfs[[i]]$formulas$pareto,
             data = object$data$data,
             bootstrap = FALSE,
             verbose = verbose)
  }else{
    reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %do%
      evinb(formulas_dfs[[i]]$formulas$nb,
             formulas_dfs[[i]]$formulas$evinf,
             formulas_dfs[[i]]$formulas$pareto,
             data = object$data,
             bootstrap = FALSE,
            verbose = verbose)
  }
    
    logliks <- reruns %>% purrr::map('log.lik') %>% purrr::reduce(c)
    dfs <- formulas_dfs %>% purrr::map('df') %>% purrr::reduce(c)
    
    res_full <- tibble::tibble(vars = vars,
                          loglik_full = object$log.lik,
                          loglik_restricted = logliks,
                          df = dfs) %>%
      dplyr::mutate(statistic = .data$loglik_full - .data$loglik_restricted,
             prob = 1-pchisq(.data$statistic,.data$df))
    
    if(bootstrap){
    if(inherits(object,'evzinb')){
      if(multicore){
        if(is.null(ncores)){
          ncores <- parallel::detectCores()-1
        }
        
        doParallel::registerDoParallel(ncores)
      }
      
    reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %:%
      foreach::foreach(j = 1:length(object$bootstraps)) %dopar%{
        if(verbose){
          message(paste('Running bootstrap',j,'for formula',i))
        }
        try(evzinb(formulas_dfs[[i]]$formulas$nb,
                   formulas_dfs[[i]]$formulas$zi,
                   formulas_dfs[[i]]$formulas$evinf,
                   formulas_dfs[[i]]$formulas$pareto,
                   data = object$data$data[object$bootstraps[[j]]$boot_id,],
                   bootstrap = FALSE))  
      }
      
    }else{
      reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %:%
        foreach::foreach(j = 1:length(object$bootstraps)) %dopar%{
          if(verbose){
            message(paste('Running bootstrap',j,'for formula',i))
          }
        try(evinb(formulas_dfs[[i]]$formulas$nb,
               formulas_dfs[[i]]$formulas$evinf,
               formulas_dfs[[i]]$formulas$pareto,
               data = object$data$data[object$bootstraps[[j]]$boot_id,],
               bootstrap = FALSE))
        }
    }
      

   logliks_boot_reduced <-  foreach::foreach(i = 1:length(formulas_dfs)) %:%
        foreach::foreach(j = 1:length(object$bootstraps),.final = unlist) %do%{
          if('try-error' %in% class(reruns[[i]][[j]])){
            NA
          }else{
            reruns[[i]][[j]]$log.lik
          }
        }
   
   logliks_boot <- object$bootstraps %>% purrr::map('log.lik') %>% purrr::reduce(c)
   
   statistics_boot <- foreach::foreach(i = 1:length(formulas_dfs)) %do%
     (logliks_boot - logliks_boot_reduced[[i]])
     
   
   p_vals <- tibble::tibble(vars = vars,
                    ks_p = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      stats::ks.test(na.omit(statistics_boot[[i]]),pchisq, dfs[i])$p.value,
                    chisq_mean = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      1-pchisq(mean(na.omit(statistics_boot[[i]])),dfs[i]),
                    chisq_median = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      1-pchisq(median(na.omit(statistics_boot[[i]])),dfs[i]),
                    prop_sig = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      mean(na.omit(statistics_boot[[i]]) > qchisq(0.95,dfs[i])),
                    n_failed_bootstraps = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      sum(is.na(statistics_boot[[i]])))
   
   res_boot <- foreach::foreach(i = 1:length(statistics_boot)) %do%
     tibble::tibble(ll_reduced = logliks_boot_reduced[[i]],
            ll_full = logliks_boot,
            statistic = statistics_boot[[i]])
   names(res_boot) <- vars
   
   res_full <- res_full %>% dplyr::left_join(p_vals,by = 'vars')
   
   out <- list(results = res_full,
               boot_results = res_boot)
   
   return(out)
   
  }else{
  return(res_full)
  }
}

  
formula_var_remover <- function(formulas,vars){
  i <- 'temp_iter'
  
  terms_nb <- attr(terms.formula(formulas$formula_nb),'term.labels')
  if(!is.null(formulas$formula_zi)){
  terms_zi <- attr(terms.formula(formulas$formula_zi),'term.labels')
  }else{
    terms_zi <- NULL
  }
  terms_evi <- attr(terms.formula(formulas$formula_evi),'term.labels')
  terms_pareto <- attr(terms.formula(formulas$formula_pareto),'term.labels')
  
  nb_rem_pos <- foreach::foreach(i = 1:length(vars)) %do%{
         stringi::stri_detect_regex(terms_nb,paste0("^",vars[i],"$"))+
      stringi::stri_detect_coll(terms_nb,paste0("(",vars[i],")"))
  }
  nb_rem_pos <- nb_rem_pos %>% purrr::reduce(`+`)
  
  zi_rem_pos <- foreach::foreach(i = 1:length(vars)) %do%{
    stringi::stri_detect_regex(terms_zi,paste0("^",vars[i],"$"))+
      stringi::stri_detect_coll(terms_zi,paste0("(",vars[i],")"))
  }
  zi_rem_pos <- zi_rem_pos %>% purrr::reduce(`+`)
  
  evi_rem_pos <- foreach::foreach(i = 1:length(vars)) %do%{
    stringi::stri_detect_regex(terms_evi,paste0("^",vars[i],"$"))+
      stringi::stri_detect_coll(terms_evi,paste0("(",vars[i],")"))
  }
  
  evi_rem_pos <- evi_rem_pos %>% purrr::reduce(`+`)
  pareto_rem_pos <- foreach::foreach(i = 1:length(vars)) %do%{
    stringi::stri_detect_regex(terms_pareto,paste0("^",vars[i],"$"))+
      stringi::stri_detect_coll(terms_pareto,paste0("(",vars[i],")"))
  }
  
  pareto_rem_pos <- pareto_rem_pos %>% purrr::reduce(`+`)
  
  df <- sum(c(nb_rem_pos,zi_rem_pos,evi_rem_pos,pareto_rem_pos))
  
  new_nb <- as.formula(paste0(formulas$formula_nb[[2]],'~',paste(terms_nb[!nb_rem_pos],collapse = '+')))
  
  if(!is.null(formulas$formula_zi)){
  new_zi <- as.formula(paste0(formulas$formula_zi[[2]],'~',paste(terms_zi[!zi_rem_pos],collapse = '+')))
  }else{
    new_zi = NULL
  }
  new_evi <- as.formula(paste0(formulas$formula_evi[[2]],'~',paste(terms_evi[!evi_rem_pos],collapse = '+')))
  new_pareto <- as.formula(paste0(formulas$formula_pareto[[2]],'~',paste(terms_pareto[!pareto_rem_pos],collapse = '+')))
  
  out <- list()
  out$df <- df
  
  if(is.null(formulas$formula_zi)){
    out$formulas <- list(nb = new_nb,
                         evinf = new_evi,
                         pareto = new_pareto)
    
  }else{
  out$formulas <- list(nb = new_nb,
                       zi = new_zi,
                       evinf = new_evi,
                       pareto = new_pareto)
  
  }
  return(out)
  }
  
#   
# ad_test_chisq <- function(x,df=NULL){
#   
#   DNAME <- deparse(substitute(x))
#   x <- sort(x[complete.cases(x)])
#   n <- length(x)
#   if (n < 8){
#     stop("sample size must be greater than 7")
#   }
#   if(is.null(df)){
#     logp1 <- pchisq(x, mean(x), log.p = TRUE)
#     logp2 <- 1-pchisq(x, mean(x), log.p = TRUE)
#   }else{
#     logp1 <- pchisq(x, df, log.p = TRUE)
#     logp2 <- 1-pchisq(x, df, log.p = TRUE)
#   }
#   
#   h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
#   A <- -n - mean(h)
#   AA <- (1 + 0.75/n + 2.25/n^2) * A
#   if(AA < 0.2){
#     pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
#   }
#   else if (AA < 0.34) {
#     pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
#   }
#   else if (AA < 0.6) {
#     pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
#   }
#   else if (AA < 10) {
#     pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
#   }
#   else{
#     pval <- 3.7e-24
#   }
#   RVAL <- list(statistic = c(A = A), p.value = pval, method = "Anderson-Darling test", 
#                data.name = DNAME)
#   class(RVAL) <- "htest"
#   return(RVAL)
# }
# 

