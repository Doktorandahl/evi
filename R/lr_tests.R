


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
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#'  lr_test_evzinb(model,'x1')
lr_test_evzinb <- function(object, vars, single = TRUE, bootstrap = FALSE){
  if(class(vars) != 'list'){
    if(single){
      formulas_dfs <- foreach::foreach(i = 1:length(vars)) %do%
        formula_var_remover(object$formulas,vars[i])
      vars <- paste(vars,collapse ='_')
    }else{
      formulas_dfs <- foreach::foreach(i = 1:1) %do%
        formula_var_remover(object$formulas,vars)
    }
  }else{
    formulas_dfs <- foreach::foreach(i = 1:length(vars)) %do%
      formula_var_remover(object$formulas,vars[[i]])
    vars <- vars  %>% purrr::map(~paste(.x,collapse = '_')) %>% purrr::reduce(c)
  }
    
    reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %do%
      evzinb(formulas_dfs[[i]]$formulas$nb,
             formulas_dfs[[i]]$formulas$zi,
             formulas_dfs[[i]]$formulas$evi,
             formulas_dfs[[i]]$formulas$pareto,
             data = object$data$data,
             bootstrap = F)
    
    logliks <- reruns %>% purrr::map('log.lik') %>% purrr::reduce(c)
    dfs <- formulas_dfs %>% purrr::map('df') %>% purrr::reduce(c)
    
    out <- tibble::tibble(vars = vars,
                          loglik_full = object$log.lik,
                          loglik_restricted = logliks,
                          df = dfs) %>%
      mutate(statistic = loglik_full - loglik_restricted,
             prob = 1-pchisq(statistic,df))
    return(out)
  }

  
formula_var_remover <- function(formulas,vars){
    
  vars_nb <- all.vars(formulas$formula_nb)
  vars_zi <- all.vars(formulas$formula_zi)
  vars_evi <- all.vars(formulas$formula_evi)
  vars_pareto <- all.vars(formulas$formula_pareto)
  allvars <- c(vars_nb,vars_zi,vars_evi,vars_pareto)
  
  df <- foreach::foreach(i = 1:length(vars)) %do%
    sum(stringr::str_count(allvars,paste0("^",vars[i],"$")))
  
  df <- df %>% purrr::reduce(c) %>% sum()
  
  
  
  
  new_nb <- as.formula(paste0(vars_nb[1],'~',paste(vars_nb[!(vars_nb %in% vars)][-1],collapse = '+')))
  if(!is.null(formulas$formula_zi)){
  new_zi <- as.formula(paste0(vars_zi[1],'~',paste(vars_zi[!(vars_zi %in% vars)][-1],collapse = '+')))
  }else{
    new_zi = NULL
  }
  new_evi <- as.formula(paste0(vars_evi[1],'~',paste(vars_evi[!(vars_evi %in% vars)][-1],collapse = '+')))
  new_pareto <- as.formula(paste0(vars_pareto[1],'~',paste(vars_pareto[!(vars_pareto %in% vars)][-1],collapse = '+')))
  
  out <- list()
  out$df <- df
  out$formulas <- list(nb = new_nb,
                       zi = new_zi,
                       evi = new_evi,
                       pareto = new_pareto)
  
  
  return(out)
  }
  
  