
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
#' @examples run_evzinb(y~x1+x2,y~x1+x3,y~x1,y~x4, data=data_test)
run_evzinb <- function(formula_nb,
                       formula_zi,
                       formula_evi,
                       formula_pareto,
                       data,
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
}else{
  Ini.Val$Beta.multinom.ZC <- init.Beta.multinom.ZC
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

          
evzinb <- zerinfl.nb.pl.regression.fun(OBS.Y,OBS.X.obj,Ini.Val,Control)
evzinb$par.mat$Beta.multinom.ZC <- as.numeric(evzinb$par.mat$Beta.multinom.ZC)
evzinb$par.mat$Beta.multinom.PL <- as.numeric(evzinb$par.mat$Beta.multinom.PL)
evzinb$par.mat$Beta.NB <- as.numeric(evzinb$par.mat$Beta.NB)
evzinb$par.mat$Beta.PL <- as.numeric(evzinb$par.mat$Beta.PL)

names(evzinb$par.mat$Beta.NB) <- c('(Intercept)',all.vars(formula_nb)[-1])
names(evzinb$par.mat$Beta.multinom.ZC) <- c('(Intercept)',all.vars(formula_zi)[-1])
names(evzinb$par.mat$Beta.multinom.PL) <- c('(Intercept)',all.vars(formula_evi)[-1])
names(evzinb$par.mat$Beta.PL) <- c('(Intercept)',all.vars(formula_pareto)[-1])



evzinb$formulas <- list(formula_nb = formula_nb,
                        formula_zi = formula_zi,
                        formula_evi = formula_evi,
                        formula_pareto = formula_pareto)
evzinb$data <-  data %>% dplyr::select(all_of(unique(c(all.vars(formula_nb),
                                                all.vars(formula_zi),
                                                all.vars(formula_evi),
                                                all.vars(formula_pareto))))) %>%
  na.omit()


class(evzinb) <- 'evzinb'
return(evzinb)

}



#' Running an extreme value and zero inflated negative binomial model with bootstrapping
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
#' @examples run_evzinb(y~x1+x2,y~x1+x3,y~x1,y~x4, data=data_test)
run_evzinb_boot <- function(formula_nb,
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

  full_run <- run_evzinb(formula_nb = formula_nb,
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
    try(bootrun_evzinb(full_run,block2))
  doParallel::stopImplicitCluster()

  names(boots) <- paste('bootstrap_',1:length(boots),sep="")
out <- c(full_run,
         list(bootstraps = boots))
class(out) <- 'evzinb'
  return(out)

}

#' Running an extreme value and zero inflated negative binomial model with bootstrapping
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
#' @examples run_evzinb(y~x1+x2,y~x1+x3,y~x1,y~x4, data=data_test)
bootrun_evzinb <- function(evzinb,block = NULL, timing=T){

  tim <- Sys.time()
  if(is.null(block)){
  boot_id <- sample(1:nrow(evzinb$x.nb),nrow(evzinb$x.nb),replace = T)
  }else{
    uniques <- unique(block)
    boot_block_id <- sample(uniques,length(uniques),replace=T)
    boot_id <- boot_block_id %>% purrr::map(~which(block == .x)) %>% purrr::reduce(c)
}
  OBS.Y <- evzinb$y[boot_id,]

  OBS.X.obj <- list()
  OBS.X.obj$X.multinom.ZC <- evzinb$x.multinom.zc[boot_id,]
  OBS.X.obj$X.multinom.PL <- evzinb$x.multinom.pl[boot_id,]
  OBS.X.obj$X.NB <- evzinb$x.nb[boot_id,]
  OBS.X.obj$X.PL <- evzinb$x.pl[boot_id,]
  Control <- evzinb$control

  Ini.Val <- list()
    Ini.Val$Beta.multinom.ZC <- as.numeric(evzinb$par.mat$Beta.multinom.ZC)
    Ini.Val$Beta.multinom.PL <- as.numeric(evzinb$par.mat$Beta.multinom.PL)
    Ini.Val$Beta.NB <- as.numeric(evzinb$par.mat$Beta.NB)
    Ini.Val$Beta.PL <- as.numeric(evzinb$par.mat$Beta.PL)
    Ini.Val$Alpha.NB <- evzinb$par.mat$Alpha.NB
  Ini.Val$C <- evzinb$par.mat$C
  evzinb_boot <- zerinfl.nb.pl.regression.fun(OBS.Y,OBS.X.obj,Ini.Val,Control)
  evzinb_boot$par.mat$Beta.multinom.ZC <- as.numeric(evzinb_boot$par.mat$Beta.multinom.ZC)
  evzinb_boot$par.mat$Beta.multinom.PL <- as.numeric(evzinb_boot$par.mat$Beta.multinom.PL)
  evzinb_boot$par.mat$Beta.NB <- as.numeric(evzinb_boot$par.mat$Beta.NB)
  evzinb_boot$par.mat$Beta.PL <- as.numeric(evzinb_boot$par.mat$Beta.PL)

  names(evzinb_boot$par.mat$Beta.NB) <- c('(Intercept)',all.vars(evzinb$formulas$formula_nb)[-1])
  names(evzinb_boot$par.mat$Beta.multinom.ZC) <- c('(Intercept)',all.vars(evzinb$formulas$formula_zi)[-1])
  names(evzinb_boot$par.mat$Beta.multinom.PL) <- c('(Intercept)',all.vars(evzinb$formulas$formula_evi)[-1])
  names(evzinb_boot$par.mat$Beta.PL) <- c('(Intercept)',all.vars(evzinb$formulas$formula_pareto)[-1])

  evzinb_boot[['x.nb']] <- NULL
  evzinb_boot[['x.pl']] <- NULL
  evzinb_boot[['y']] <- NULL
  evzinb_boot[['x.multinom.pl']] <- NULL
  evzinb_boot[['x.multinom.zc']] <- NULL
  evzinb_boot$boot_id <- boot_id
  if(timing){
    evzinb_boot$time <- difftime(Sys.time(),tim,units='secs')
  }
  return(evzinb_boot)
}


