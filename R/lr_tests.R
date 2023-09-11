


#' Likelihood ratio test for individual variables of evzinb
#'
#' @param object EVZINB or EVZINB object to perform likelihood ratio test on
#' @param vars Either a list of character vectors with variable names which to be restricted in the LR test or a character vector of variable names. If a list, each character vector of the list will be run separately, allowing for multiple variables to be restricted as once. If a character vector, parameter 'single' can be used to determine whether all variables in the vector should be restricted at once (single = FALSE) or if the variables should be restricted one by one (single = TRUE)
#' @param single Logical. Determining whether variables in 'vars' should be restricted individually (single = TRUE) or all at once (single = FALSE)
#' @param bootstrap Should LR tests be conducted on each bootstrapped sample or only on the original sample. Not yet implemented.
#'
#' @return A tibble with one row per performed LR test
#' @export
#'
#' @examples 
#' \dontrun{
#' data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#'  lr_test_evzinb(model,'x1')
#'  }
lr_test_evzinb <- function(object, vars, single = TRUE, bootstrap = FALSE){
  i <- 'temp_iter'
  
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
    
    reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %do%
      evzinb(formulas_dfs[[i]]$formulas$nb,
             formulas_dfs[[i]]$formulas$zi,
             formulas_dfs[[i]]$formulas$evinf,
             formulas_dfs[[i]]$formulas$pareto,
             data = object$data$data,
             bootstrap = F)
    
    logliks <- reruns %>% purrr::map('log.lik') %>% purrr::reduce(c)
    dfs <- formulas_dfs %>% purrr::map('df') %>% purrr::reduce(c)
    
    out <- tibble::tibble(vars = vars,
                          loglik_full = object$log.lik,
                          loglik_restricted = logliks,
                          df = dfs) %>%
      dplyr::mutate(statistic = .data$loglik_full - .data$loglik_restricted,
             prob = 1-pchisq(.data$statistic,.data$df))
    return(out)
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
  out$formulas <- list(nb = new_nb,
                       zi = new_zi,
                       evinf = new_evi,
                       pareto = new_pareto)
  
  
  return(out)
  }
  
  
