
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
                       data){

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
Control <- list(max.diff.par = 1e-3, max.no.em.steps = 200, max.no.em.steps.warmup = 5,c.lim=c(70,300),
                max.upd.par.zc.multinomial=0.5, max.upd.par.pl.multinomial=0.5, max.upd.par.nb = 0.5,
                max.upd.par.pl = 0.5,no.m.bfgs.steps.multinomial=3,no.m.bfgs.steps.nb = 3,
                no.m.bfgs.steps.pl = 3,pdf.pl.type="approx", eta.int = c(-1,1))
Control$c.lim <- c(20,1000)

Ini.Val <- list()
Ini.Val$Beta.multinom.ZC <- rep(0,ncol(mf_zi))
Ini.Val$Beta.multinom.PL <- rep(0,ncol(mf_evi))
Ini.Val$Beta.NB <- rep(0,ncol(mf_nb))
Ini.Val$Beta.PL <- rep(0,ncol(mf_pareto))
Ini.Val$Alpha.NB <- 0.01
Ini.Val$C <- 200


zerinfl.nb.pl.regression.fun(OBS.Y,OBS.X.obj,Ini.Val,Control)

}
