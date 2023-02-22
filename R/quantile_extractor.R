#' Extracting full mixture quantiles from an evzinb object
#'
#' @param evzinb  An evzinb object for which to produce quantiles
#' @param quantile The quantile for which to produce predictions
#' @param newdata  Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#' @param multicore_ok Logical: Should parallel processing be used to obtain results
#' @param ncores Number of cores if multicore is used
#'
#' @return A vector of predicted quantiles, or if return_data=T, a tibble with the predicted quantile attached last
#' @export
#'
#' @examples quantiles_from_evzinb(evzinb,0.95)
quantiles_from_evzinb <- function(evzinb, quantile,
                                  newdata = NULL,
                                  return_data =F,multicore_ok = T,ncores=NULL){



  ## Estimate component probabilities for all individuals
  prbs <- prob_from_evzinb(evzinb,
                            newdata = newdata)
  ## Estimate mu_nb for all individuals
  cnts <- counts_from_evzinb(evzinb,
                              newdata = newdata)
  ## Estimate pareto alpha for all individuals
  alphs <- fitted_alpha_from_evzinb(evzinb,
                                     newdata = newdata)
  if(min(alphs)<1e-02){
    warning('Fitted pareto alpha-values below 1e-02 detected. Setting those alphas to 1e-02.')
    alphs <- alphs %>% mutate(alpha = case_when(alpha<1e-02 ~ 1e-02,
                                                T~alpha))
  }

  all_pars <- bind_cols(prbs,cnts,alphs) %>% mutate(C = evzinb$par.mat$C,
                                                    alpha_nb = evzinb$par.mat$Alpha.NB)



  individual_dists <- all_pars %>% transpose() %>% map(~mistr::mixdist(mistr::binomdist(1,0),
                                                                nbinomdist2(mu=.x$count,size=.x$alpha_nb),
                                                                mistr::paretodist(scale = .x$C,shape=.x$alpha),
                                                                weights = c(.x$pr_zc,.x$pr_count,.x$pr_pl)))
  if(!multicore_ok){
    q <- individual_dists %>% map(~mistr::mistr_q(.x,quantile)) %>% reduce(c) %>% round()
  }else{
    if(is.null(ncores)){
      registerDoParallel(cores=detectCores()-1)
    }else{
      registerDoParallel(cores=ncores)
    }
    q <- foreach(i = 1:length(individual_dists),.final = unlist)%dopar%
      mistr::mistr_q(individual_dists[[i]],quantile)
}

  if(return_data){
    return(newdata %>% mutate(q = q))
  }else{
    return(q)
  }
}


#' Extracting state probabilities from an evzinb object
#'
#' @param evzinb An evzinb object for which to produce probabilities
#' @param newdata Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#'
#' @return A tibble with the predicted state probabilities. If return_data=T this is appended to the data or newdata
#' @export
#'
#' @examples
prob_from_evzinb <- function(evzinb,newdata=NULL, return_data = F){

  if(is.null(newdata)){
    x.multinom.zc <- evzinb$x.multinom.zc
    x.multinom.pl <- evzinb$x.multinom.pl
  }else{
    x.multinom.zc <- as.matrix(model.frame(evzinb$formulas$formula_zi[-2],newdata))
    x.multinom.pl <- as.matrix(model.frame(evzinb$formulas$formula_evi[-2],newdata))
   }

  pr_zc <- exp(cbind(1,x.multinom.zc)%*%evzinb$par.mat$Beta.multinom.ZC)/
    (1+exp(cbind(1,x.multinom.zc)%*%evzinb$par.mat$Beta.multinom.ZC)+
       exp(cbind(1,x.multinom.pl)%*%evzinb$par.mat$Beta.multinom.PL))


  pr_pl <- exp(cbind(1,x.multinom.pl)%*%evzinb$par.mat$Beta.multinom.PL)/
    (1+exp(cbind(1,x.multinom.zc)%*%evzinb$par.mat$Beta.multinom.ZC)+
       exp(cbind(1,x.multinom.pl)%*%evzinb$par.mat$Beta.multinom.PL))

  pr_count <- 1-pr_zc-pr_pl

  out <- tibble(pr_zc=as.numeric(pr_zc),
                pr_count = as.numeric(pr_count),
                pr_pl = as.numeric(pr_pl))

  if(return_data){
    if(is.null(newdata)){
      out <- bind_cols(evzinb$data,out)

    }else{
    out <- bind_cols(newdata,out)
    }
  }

  return(out)


}

#' Extracting fitted count values of the NB component of an evzinb object
#'
#' @param evzinb An evzinb object for which to produce counts
#' @param newdata Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#'
#' @return A tibble with the predicted nb counts. If return_data=T this is appended to the data or newdata
#' @export
#'
#' @examples
counts_from_evzinb <- function(evzinb,newdata=NULL, return_data = F){

  if(is.null(newdata)){
    x.nb <- evzinb$x.nb
  }else{
    x.nb <- as.matrix(model.frame(evzinb$formulas$formula_nb[-2],newdata))
  }

  count <- exp(cbind(1,x.nb)%*%evzinb$par.mat$Beta.NB)


  out <- tibble(count=as.numeric(count))

  if(return_data){
    if(is.null(newdata)){
      out <- bind_cols(evzinb$data,out)

    }else{
      out <- bind_cols(newdata,out)
    }
  }
  return(out)


}

#' Extracting fitted alpha values of the pareto component of an evzinb object
#'
#' @param evzinb An evzinb object for which to produce quantiles
#' @param newdata Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#'
#' @return A tibble with the predicted pareto alpha values. If return_data=T this is appended to the data or newdata
#' @export
#'
#' @examples
fitted_alpha_from_evzinb <- function(evzinb,newdata=NULL, return_data = F){
  if(is.null(newdata)){
    x.pl <- evzinb$x.pl
  }else{
    x.pl <- as.matrix(model.frame(evzinb$formulas$formula_pareto[-2],newdata))
  }

  alpha <- exp(cbind(1,x.pl)%*%evzinb$par.mat$Beta.PL)


  out <- tibble(alpha=as.numeric(alpha))

  if(return_data){
    if(is.null(newdata)){
      out <- bind_cols(evzinb$data,out)

    }else{
      out <- bind_cols(newdata,out)
    }
  }
  return(out)


}
