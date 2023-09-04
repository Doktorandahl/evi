


#' Bootstrap coefficient extractor
#'
#' @param object a fitted model with bootstraps of class evzinb, evinb, nbboot, or zinbboot
#' @param ... Component to be extracted (not for nbboot). Alternatives are 'nb','zi','evi','pareto', and 'all'
#'
#' @return A tibble with coefficient values, one row per bootstrap and component
#' @export
#'
#' @examples
#' data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' coefficient_extractor(model, component = 'all')
coefficient_extractor <- function(object,...){
  UseMethod('coefficient_extractor')
}


#' Bootstrap coefficient extractor
#'
#' @param object A fitted evzinb model with bootstraps
#' @param component Which component should be extracted
#' @param ... Not in use
#'
#' @return A tibble with coefficient values, one row per bootstrap and component
#' @export
#'
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' coefficient_extractor(model, component = 'all')
coefficient_extractor.evzinb <- function(object,component = c('nb','zi','evi','pareto','all'),...){

  component <- match.arg(component,c('nb','zi','evi','pareto','all'))
  
  object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
  
  nb_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows()
  
  zi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.ZC') %>% dplyr::bind_rows()
    
  
  evi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() 
  
  pareto_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() 

if(component == 'nb'){
  return(nb_boot)
}else if(component == 'zi'){
  return(zi_boot)
}else if(component == 'evi'){
  return(evi_boot)
}else if(component == 'pareto'){
  return(pareto_boot)
}else{
  out <- dplyr::bind_rows(nb_boot %>% dplyr::mutate(.component = 'nb'),
                   zi_boot %>% dplyr::mutate(.component = 'zi'),
                   evi_boot %>% dplyr::mutate(.component = 'evi'),
                   pareto_boot %>% dplyr::mutate(.component = 'pareto'))
  return(out)
}
  
  
  }
  
#' Bootstrap coefficient extractor
#'
#' @param object A fitted evinb model with bootstraps
#' @param component Which component should be extracted
#' @param ... Not in use
#'
#' @return A tibble with coefficient values, one row per bootstrap and component
#' @export
#'
#' @examples data(genevzinb)
#' model <- evinb(y~x1+x2+x3,data=genevzinb)
#' coefficient_extractor(model, component = 'all')
coefficient_extractor.evinb <- function(object,component = c('nb','evi','pareto','all'),...){
  
  component <- match.arg(component,c('nb','evi','pareto','all'))
  
  object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
  
  nb_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows()
  evi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows() 
  
  pareto_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows() 

if(component == 'nb'){
  return(nb_boot)
}else if(component == 'evi'){
  return(evi_boot)
}else if(component == 'pareto'){
  return(pareto_boot)
}else{
  out <- dplyr::bind_rows(nb_boot %>% dplyr::mutate(.component = 'nb'),
                   evi_boot %>% dplyr::mutate(.component = 'evi'),
                   pareto_boot %>% dplyr::mutate(.component = 'pareto'))
  return(out)
}


}

#' Bootstrap coefficient extractor
#'
#' @param object A fitted evinb model with bootstraps
#' @param component Which component should be extracted
#' @param ... Not in use
#'
#' @return A tibble with coefficient values, one row per bootstrap and component
#' @export
#'
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' zinb_comp <- compare_models(model)
#' coefficient_extractor(zinb_comp$zinb)
coefficient_extractor.zinbboot <- function(object,component = c('nb','zi','all'),...){
  
  component <- match.arg(component,c('nb','zi','all'))
  
  object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
  
  
  nb_boot <- object$bootstraps %>% purrr::map('coefficients') %>% purrr::map('count') %>% dplyr::bind_rows() 
  
  zi_boot <- object$bootstraps %>% purrr::map('coefficients') %>% purrr::map('zero') %>% dplyr::bind_rows() 
  
  
  
  if(component == 'nb'){
    return(nb_boot)
  }else if(component == 'zi'){
    return(zi_boot)
  }else{
    out <- dplyr::bind_rows(nb_boot %>% dplyr::mutate(.component = 'nb'),
                     zi_boot %>% dplyr::mutate(.component = 'zi'))
    return(out)
  }
  
  
}

#' Bootstrap coefficient extractor
#'
#' @param object A fitted nbboot model with bootstraps
#' @param ... Not in use
#'
#' @return A tibble with coefficient value, one row per bootstrap
#' @export
#'
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' zinb_comp <- compare_models(model)
#' coefficient_extractor(zinb_comp$nb)
coefficient_extractor.nbboot <- function(object,...){

  object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
  
  nb_boot <- object$bootstraps %>% purrr::map('coefficients') %>% dplyr::bind_rows() 
  
    return(nb_boot)
  
  
  
}
