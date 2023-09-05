
#' Out of bag predictive performance of EVZINB and EVINB models
#'
#' @param object A fitted evzinb or evinb with bootstraps on which to conduct out-of-bag evaluation
#' @param predict_type What type of prediction should be made? Harmonic mean, or exp(log(prediction))?
#' @param metric What metric should be used for the out of bag evaluation? Default options include rmsle, rmse, mse, and mae. Can also take a user supplied function of the form function(y_pred,y_true) which returns a single value
#'
#' @return A vector of oob evaluation metrics of the length of the number of bootstraps in the evzinb/evinb object.
#' @export
#'
#' @examples data(genevzinb)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb)
#' oob_evaluation(model)
oob_evaluation <- function(object,predict_type = c('harmonic','explog'),
                           metric = c('rmsle','rmse','mse','mae')){
  
  if(is.character(metric)){
    metric <- match.arg(metric, c('rmsle','rmse','mse','mae'))
    if(metric=='rmsle'){
      ev_metric <- MLmetrics::RMSLE
    }else{
  ev_metric <- getFromNamespace(toupper(metric),"MLmetrics")
    }
  }else{
    ev_metric <- metric
  }
  predict_type <- match.arg(predict_type, c('harmonic','explog'))
  
  evals <- foreach::foreach(i = 1:length(object$bootstraps)) %do%
    try(oob_inner(object$bootstraps[[i]],object$data,predict_type,ev_metric,model_type = class(object)))
  
  evals <- purrr::reduce(evals,c)
  
  return(evals)
}


oob_inner <- function(bootstrap,data,predict_type,ev_metric,model_type = c('evzinb','evinb')){
  
  oob_data <- data$data[-bootstrap$boot_id,]
if(model_type=="evzinb"){
  predictions <- predict.evzinb(bootstrap,newdata=oob_data,type=predict_type)
}else{
  predictions <- predict.evinb(bootstrap,newdata=oob_data,type=predict_type)
}
  ev_metric(predictions,data$y[-bootstrap$boot_id])
  
}


# oob_inner_old <- function(bootstrap,data,predict_type){
#   data <- data[-bootstrap$boot_id,]
#   
#   predictions <- predict.evzinb(bootstrap,newdata=data,type='both',include_index = model.response(model.frame(bootstrap$formulas$formula_nb,data)))
#   predictions <- predictions %>% mutate(sq_logerror_harmonic = (log1p(index)-log1p(harmonic_pred))^2,
#                                         sq_error_harmonic = (index-harmonic_pred)^2,
#                                         sq_logerror_explog = (log1p(index)-log1p(explog_pred))^2,
#                                         sq_error_explog = (index-explog_pred)^2)
#   out <- predictions %>% summarize(rmsle_harmonic = sqrt(mean(sq_logerror_harmonic)),
#                                    rmse_harmonic = sqrt(mean(sq_error_harmonic)),
#                                    rmsle_explog = sqrt(mean(sq_logerror_explog)),
#                                    rmse_explog = sqrt(mean(sq_error_explog)))
#   if(predict_type=='harmonic'){
#     out <- out %>% select(-contains('explog'))
#   }else if(predict_type == 'explog'){
#     out <- out %>% select(-contains('harmonic'))
#   }
#   return(out)
# }
# 
# 
# 
  
  
