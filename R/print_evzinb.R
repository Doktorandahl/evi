



#' EVZINB print function
#'
#' @param x A fitted evzinb model
#' @param ... Not used
#' @return An evzinb print function
#' @export
#'
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' print(model)
print.evzinb <- function(x,...){
  cat('\n','Fitted EVZINB model with formulas:',
      '\n NB:', paste(x$formulas$formula_nb)[2],paste(x$formulas$formula_nb)[1],
      paste(x$formulas$formula_nb)[3],
      '\n ZI:', paste(x$formulas$formula_zi)[2],paste(x$formulas$formula_zi)[1],
      paste(x$formulas$formula_nb)[3],
      '\n EVI:', paste(x$formulas$formula_evi)[2],paste(x$formulas$formula_evi)[1],
      paste(x$formulas$formula_evi)[3],
      '\n Pareto:', paste(x$formulas$formula_pareto)[2],paste(x$formulas$formula_pareto)[1],paste(x$formulas$formula_pareto)[3],
      '\n ______',
      '\n Converged:', x$converge,
      '\n Number of bootstraps: ', length(x$bootstraps),
      '\n Number of failed bootstraps: ', length(x$bootstraps) - x$bootstraps %>% purrr::discard(~'try-error' %in% class(.x)) %>% length(),
      '\n Parameters: ', length(x$par.all))
}

#' EVINB print function
#'
#' @param x A fitted evinb model
#' @param ... Not used
#' @return An evinb print function
#' @export
#'
#' @examples data(genevzinb)
#' model <- evinb(y~x1+x2+x3,data=genevzinb)
#' print(model)
print.evinb <- function(x,...){
  cat('\n','Fitted EVZINB model with formulas:',
      '\n NB:', paste(x$formulas$formula_nb)[2],paste(x$formulas$formula_nb)[1],
      paste(x$formulas$formula_nb)[3],
      '\n EVI:', paste(x$formulas$formula_evi)[2],paste(x$formulas$formula_evi)[1],
      paste(x$formulas$formula_evi)[3],
      '\n Pareto:', paste(x$formulas$formula_pareto)[2],paste(x$formulas$formula_pareto)[1],paste(x$formulas$formula_pareto)[3],
      '\n ______',
      '\n Converged:', x$converge,
      '\n Number of bootstraps: ', length(x$bootstraps),
      '\n Number of failed bootstraps: ', length(x$bootstraps) - x$bootstraps %>% purrr::discard(~'try-error' %in% class(.x)) %>% length(),
      '\n Parameters: ', length(x$par.all))
}

#' Compare_models print function
#'
#' @param x A fitted evinb model
#' @param ... Not used
#' @return An evinb print function
#' @export
#'
#' @examples data(genevzinb)
#' model <- evinb(y~x1+x2+x3,data=genevzinb)
#' print(model)
print.evzinbcomp <- function(x,...){
  comp_class <- class(x[[1]])
  
  cat('\n','Model comparison of ', comp_class,
      '\n ', 'Number of compared models: ', length(x)-1,
      '\n Number of bootstraps:', length(x[[1]]$bootstraps))
}
