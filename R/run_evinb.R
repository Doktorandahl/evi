
#' Running an extreme value and zero inflated negative binomial model
#'
#' @param formula_nb Formula for the negative binomial (count) component of the model
#' @param formula_zi Formula for the zero-inflation component of the model
#' @param formula_evi Formula for the extreme-value inflation component of the model
#' @param formula_pareto Formula for the pareto (extreme value) component of the model
#' @param data Data to run the model on
#'
#' @return An object of class 'evi' containing XX
#' @export
#'
#' @examples run_evinb(y~x1+x2,y~x1+x3,y~x1,y~x4, data=data_test)
run_evinb <- function(formula_nb,
                       formula_zi,
                       formula_evi,
                       formula_pareto,
                       data,
                       max.diff.par = 1e-3,
                       max.no.em.steps = 200,
                       max.no.em.steps.warmup = 5,
                       c.lim=c(70,300),
                       max.upd.par.zc.multinomial=0.5,
                       max.upd.par.pl.multinomial=0.5,
                       max.upd.par.nb = 0.5,
                       max.upd.par.pl = 0.5,
                       no.m.bfgs.steps.multinomial=3,
                       no.m.bfgs.steps.nb = 3,
                       no.m.bfgs.steps.pl = 3,
                       pdf.pl.type="approx",
                       eta.int = c(-1,1),
                       init.Beta.multinom.ZC = NULL,
                       init.Beta.multinom.PL = NULL,
                       init.Beta.NB = NULL,
                       init.Beta.PL = NULL,
                       init.Alpha.NB = 0.01,
                       init.C = 200){

mf_nb <- model.frame(formula_nb,data)
mf_zi <- model.frame(formula_zi,data)
mf_evi <- model.frame(formula_evi,data)
mf_pareto <- model.frame(formula_pareto,data)

OBS.Y <- as.matrix(model.response(mf_nb))

OBS.X.obj <- list()
OBS.X.obj$X.multinom.ZC <- as.matrix(mf_zi[,-1])
OBS.X.obj$X.multinom.PL <- as.matrix(mf_evi[,-1])
OBS.X.obj$X.NB <- as.matrix(mf_nb[,-1])
OBS.X.obj$X.PL <- as.matrix(mf_pareto[,-1])
Control <- list(max.diff.par = max.diff.par,
                max.no.em.steps = max.no.em.steps,
                max.no.em.steps.warmup = max.no.em.steps.warmup,
                c.lim=c.lim,
                max.upd.par.zc.multinomial=max.upd.par.zc.multinomial,
                max.upd.par.pl.multinomial=max.upd.par.pl.multinomial,
                max.upd.par.nb = max.upd.par.nb,
                max.upd.par.pl = max.upd.par.pl,
                no.m.bfgs.steps.multinomial=no.m.bfgs.steps.multinomial,
                no.m.bfgs.steps.nb = no.m.bfgs.steps.nb,
                no.m.bfgs.steps.pl = no.m.bfgs.steps.pl,
                pdf.pl.type=pdf.pl.type,
                eta.int = eta.int)

Ini.Val <- list()
if(is.null(init.Beta.multinom.ZC) | length(init.Beta.multinom.ZC) != ncol(mf_zi)){
Ini.Val$Beta.multinom.ZC <- rep(0,ncol(mf_zi))
Ini.Val$Beta.multinom.ZC[1] <- -100
}else{
  Ini.Val$Beta.multinom.ZC <- init.Beta.multinom.ZC
  Ini.Val$Beta.multinom.ZC[1] <- -100
}

if(is.null(init.Beta.multinom.PL) | length(init.Beta.multinom.PL) != ncol(mf_evi)){
  Ini.Val$Beta.multinom.PL <- rep(0,ncol(mf_evi))
}else{
  Ini.Val$Beta.multinom.PL <- init.Beta.multinom.PL
}

if(is.null(init.Beta.NB) | length(init.Beta.NB) != ncol(mf_evi)){
  Ini.Val$Beta.NB <- rep(0,ncol(mf_evi))
}else{
  Ini.Val$Beta.NB <- init.Beta.NB
}

if(is.null(init.Beta.PL) | length(init.Beta.PL) != ncol(mf_pareto)){
  Ini.Val$Beta.PL <- rep(0,ncol(mf_pareto))
}else{
  Ini.Val$Beta.PL <- init.Beta.PL
}

Ini.Val$Alpha.NB <- init.Alpha.NB
Ini.Val$C <- init.C


evinb <- zerinfl.nb.pl.regression.fun(OBS.Y,OBS.X.obj,Ini.Val,Control)
evinb$par.mat$Beta.multinom.ZC <- as.numeric(evinb$par.mat$Beta.multinom.ZC)
evinb$par.mat$Beta.multinom.PL <- as.numeric(evinb$par.mat$Beta.multinom.PL)
evinb$par.mat$Beta.NB <- as.numeric(evinb$par.mat$Beta.NB)
evinb$par.mat$Beta.PL <- as.numeric(evinb$par.mat$Beta.PL)

names(evinb$par.mat$Beta.NB) <- c('(Intercept)',all.vars(formula_nb)[-1])
names(evinb$par.mat$Beta.multinom.ZC) <- c('(Intercept)',all.vars(formula_zi)[-1])
names(evinb$par.mat$Beta.multinom.PL) <- c('(Intercept)',all.vars(formula_evi)[-1])
names(evinb$par.mat$Beta.PL) <- c('(Intercept)',all.vars(formula_pareto)[-1])



evinb$formulas <- list(formula_nb = formula_nb,
                        formula_zi = formula_zi,
                        formula_evi = formula_evi,
                        formula_pareto = formula_pareto)
evinb$data <-  data %>% dplyr::select(all_of(unique(c(all.vars(formula_nb),
                                                       all.vars(formula_zi),
                                                       all.vars(formula_evi),
                                                       all.vars(formula_pareto))))) %>%
  na.omit()


class(evinb) <- 'evinb'
return(evinb)

}


#' Running an extreme value inflated negative binomial model with bootstrapping
#'
#' @param formula_nb Formula for the negative binomial (count) component of the model
#' @param formula_zi Formula for the zero-inflation component of the model
#' @param formula_evi Formula for the extreme-value inflation component of the model
#' @param formula_pareto Formula for the pareto (extreme value) component of the model
#' @param data Data to run the model on
#'
#' @return An object of class 'evi' containing XX
#' @export
#'
#' @examples run_evinb(y~x1+x2,y~x1+x3,y~x1,y~x4, data=data_test)
bootrun_evinb <- function(evinb,block = NULL, timing=T){
  
  tim <- Sys.time()
  if(is.null(block)){
    boot_id <- sample(1:nrow(evinb$x.nb),nrow(evinb$x.nb),replace = T)
  }else{
    uniques <- unique(block)
    boot_block_id <- sample(uniques,length(uniques),replace=T)
    boot_id <- boot_block_id %>% purrr::map(~which(block == .x)) %>% purrr::reduce(c)
  }
  OBS.Y <- evinb$y[boot_id,]
  
  OBS.X.obj <- list()
  OBS.X.obj$X.multinom.ZC <- evinb$x.multinom.zc[boot_id,]
  OBS.X.obj$X.multinom.PL <- evinb$x.multinom.pl[boot_id,]
  OBS.X.obj$X.NB <- evinb$x.nb[boot_id,]
  OBS.X.obj$X.PL <- evinb$x.pl[boot_id,]
  Control <- evinb$control
  
  Ini.Val <- list()
  Ini.Val$Beta.multinom.ZC <- as.numeric(evinb$par.mat$Beta.multinom.ZC)
  Ini.Val$Beta.multinom.PL <- as.numeric(evinb$par.mat$Beta.multinom.PL)
  Ini.Val$Beta.NB <- as.numeric(evinb$par.mat$Beta.NB)
  Ini.Val$Beta.PL <- as.numeric(evinb$par.mat$Beta.PL)
  Ini.Val$Alpha.NB <- evinb$par.mat$Alpha.NB
  Ini.Val$C <- evinb$par.mat$C
  evinb_boot <- zerinfl.nb.pl.regression.fun(OBS.Y,OBS.X.obj,Ini.Val,Control)
  evinb_boot$par.mat$Beta.multinom.ZC <- as.numeric(evinb_boot$par.mat$Beta.multinom.ZC)
  evinb_boot$par.mat$Beta.multinom.PL <- as.numeric(evinb_boot$par.mat$Beta.multinom.PL)
  evinb_boot$par.mat$Beta.NB <- as.numeric(evinb_boot$par.mat$Beta.NB)
  evinb_boot$par.mat$Beta.PL <- as.numeric(evinb_boot$par.mat$Beta.PL)
  
  names(evinb_boot$par.mat$Beta.NB) <- c('(Intercept)',all.vars(evinb$formulas$formula_nb)[-1])
  names(evinb_boot$par.mat$Beta.multinom.ZC) <- c('(Intercept)',all.vars(evinb$formulas$formula_zi)[-1])
  names(evinb_boot$par.mat$Beta.multinom.PL) <- c('(Intercept)',all.vars(evinb$formulas$formula_evi)[-1])
  names(evinb_boot$par.mat$Beta.PL) <- c('(Intercept)',all.vars(evinb$formulas$formula_pareto)[-1])
  
  evinb_boot[['x.nb']] <- NULL
  evinb_boot[['x.pl']] <- NULL
  evinb_boot[['y']] <- NULL
  evinb_boot[['x.multinom.pl']] <- NULL
  evinb_boot[['x.multinom.zc']] <- NULL
  evinb_boot$boot_id <- boot_id
  if(timing){
    evinb_boot$time <- difftime(Sys.time(),tim,units='secs')
  }
  return(evinb_boot)
}

#' Running an extreme value inflated negative binomial model with bootstrapping
#'
#' @param formula_nb Formula for the negative binomial (count) component of the model
#' @param formula_zi Formula for the zero-inflation component of the model
#' @param formula_evi Formula for the extreme-value inflation component of the model
#' @param formula_pareto Formula for the pareto (extreme value) component of the model
#' @param data Data to run the model on
#' 
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#'
#' @return An object of class 'evi' containing XX
#' @export
#'
#' @examples run_evinb(y~x1+x2,y~x1+x3,y~x1,y~x4, data=data_test)
run_evinb_boot <- function(formula_nb,
                            formula_zi,
                            formula_evi,
                            formula_pareto,
                            data,
                            n_bootstraps = 100,
                            multicore = T,
                            ncores = 6,
                            block = NULL,
                           boot_seed = NULL,
                           max.diff.par = 1e-3,
                            max.no.em.steps = 200,
                            max.no.em.steps.warmup = 5,
                            c.lim=c(50,1000),
                            max.upd.par.zc.multinomial=0.5,
                            max.upd.par.pl.multinomial=0.5,
                            max.upd.par.nb = 0.5,
                            max.upd.par.pl = 0.5,
                            no.m.bfgs.steps.multinomial=3,
                            no.m.bfgs.steps.nb = 3,
                            no.m.bfgs.steps.pl = 3,
                            pdf.pl.type="approx",
                            eta.int = c(-1,1),
                            init.Beta.multinom.ZC = NULL,
                            init.Beta.multinom.PL = NULL,
                            init.Beta.NB = NULL,
                            init.Beta.PL = NULL,
                            init.Alpha.NB = 0.01,
                            init.C = 200){
  
  full_run <- run_evinb(formula_nb = formula_nb,
                         formula_zi = formula_zi,
                         formula_evi = formula_evi,
                         formula_pareto = formula_pareto,
                         data = data,
                         max.diff.par = max.diff.par,
                         max.no.em.steps = max.no.em.steps,
                         max.no.em.steps.warmup = max.no.em.steps.warmup,
                         c.lim = c.lim,
                         max.upd.par.zc.multinomial = max.upd.par.zc.multinomial,
                         max.upd.par.pl.multinomial = max.upd.par.pl.multinomial,
                         max.upd.par.nb = max.upd.par.nb,
                         max.upd.par.pl = max.upd.par.pl,
                         no.m.bfgs.steps.multinomial,
                         no.m.bfgs.steps.nb,
                         no.m.bfgs.steps.pl,
                         pdf.pl.type,
                         eta.int,
                         init.Beta.multinom.ZC,
                         init.Beta.multinom.PL,
                         init.Beta.NB,
                         init.Beta.PL,
                         init.Alpha.NB,
                         init.C)
  
  full_run$block <- block
  if(!is.null(block)){
    if(is.character(block)){
      block2 <- data %>% dplyr::select(all_of(block)) %>% dplyr::pull()
      full_run$data <- dplyr::bind_cols(data %>% dplyr::select(all_of(block)),full_run$data)
    }
  }else{
    block2 <- NULL
  }
  
  
  doParallel::registerDoParallel(cores = ncores)
  boots <- foreach::foreach(i=1:n_bootstraps,.options.RNG = boot_seed) %dorng%
    try(bootrun_evinb(full_run,block2))
  doParallel::stopImplicitCluster()
  
  names(boots) <- paste('bootstrap_',1:length(boots),sep="")
  out <- c(full_run,
           list(bootstraps = boots))
  class(out) <- 'evinb'
  return(out)
  
}
