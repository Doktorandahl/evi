library(mvtnorm)
library(maxLik)

#x.obj <- X.obj
#y <- OBS.Y
#ini.val <- Ini.Val
#control <- Control
# x.obj <- OBS.X.obj
# y <- OBS.Y
# control <- Control
# ini.val <- Ini.Val
# ini.val <- TruePar
# 
# 
# x.obj <- OBS.X.obj.r
# y <- OBS.Y.r
# control <- Control.r
# ini.val <- Ini.Val.r

# x.obj <- OBS.X.obj.ins
# y <- OBS.Y.ins
# ini.val <- Ini.Val
# control <- Control
zerinfl.nb.pl.reg.cond.c.fun <- function(y,x.obj,ini.val,control){
  
  #Functions
  ell.nb.i.fun <- function(theta.nb){
    
    beta.nb <- theta.nb[1:n.beta.nb]
    alpha.nb <- theta.nb[n.beta.nb + 1]
    
    ell.nb.i <- 0
    xtb.nb.i <- x.nb.extended[i,]%*%beta.nb
    mu.i <- exp(xtb.nb.i)
    
    #sum ln(j + a^-1)
    if(y[i]>0){
      for(j in 0:(y[i]-1)){
        ell.nb.i <- ell.nb.i + log(j + 1/alpha.nb)
      }
    }
    
    #gamma(y+1)
    if(y[i]>0){
      for(j in 1:y[i]){
        ell.nb.i <- ell.nb.i - log(j)
      }
    }
    
    ell.nb.i <- ell.nb.i - (1/alpha.nb)*log(1 + alpha.nb*mu.i) - y[i]*log(1 + alpha.nb*mu.i) + y[i]*log(alpha.nb) + y[i]*log(mu.i)
    
  
    return(ell.nb.i)
  }
  
  delldtheta.nb.i.fun <- function(theta.nb){
    beta.nb <- theta.nb[1:n.beta.nb]
    alpha.nb <- theta.nb[n.beta.nb + 1]
    xtb.nb.i <- x.nb.extended[i,]%*%beta.nb
    mu.i <- exp(xtb.nb.i)
    
    if(y[i]==0){
      delldalpha.nb.i <- log(1 + alpha.nb*mu.i)/alpha.nb^2 - mu.i/(alpha.nb*(1+alpha.nb*mu.i))
    }else{
    
      delldalpha.nb.i <- log(1 + alpha.nb*mu.i)
      for(j in 0:(y[i]-1)){
        delldalpha.nb.i <- delldalpha.nb.i - 1/(j+1/alpha.nb)
      }
      delldalpha.nb.i <- delldalpha.nb.i/alpha.nb^2
      
      delldalpha.nb.i <- delldalpha.nb.i + (y[i]-mu.i)/(alpha.nb*(1+alpha.nb*mu.i))
    }
    delldbeta.nb.i <- as.numeric((y[i]-mu.i)/(1+alpha.nb*mu.i))*x.nb.extended[i,]
    
    delldtheta.nb.i <- c(delldbeta.nb.i,delldalpha.nb.i)
    return(delldtheta.nb.i)
  }
  
  d2elldtheta2.nb.i.fun <- function(theta.nb){
    beta.nb <- theta.nb[1:n.beta.nb]
    alpha.nb <- theta.nb[n.beta.nb + 1]
    xtb.nb.i <- x.nb.extended[i,]%*%beta.nb
    mu.i <- exp(xtb.nb.i)
    
    hessian.nb.i <- matrix(NA,nrow=length(theta.nb),ncol=length(theta.nb))
    
    d2elldbeta2.nb.i <- as.numeric(-1.0*mu.i*(1+alpha.nb.old*y[i])/(1+alpha.nb.old*mu.i)^2)*x.nb.extended[i,]%*%t(x.nb.extended[i,])
    
    d2elldalpha2.nb.i <- 0
    
    if(y[i]>0){
      for(j in 0:(y[i]-1)){
        d2elldalpha2.nb.i <- d2elldalpha2.nb.i - (j/(1+alpha.nb.old*j))^2
      }
    }
    
    d2elldalpha2.nb.i <- d2elldalpha2.nb.i - 2/alpha.nb.old^3*log(1 + alpha.nb.old*mu.i) + (2*(1/alpha.nb.old)^2*mu.i)/(1+alpha.nb.old*mu.i) + (y[i]+1/alpha.nb.old)*mu.i^2/(1+alpha.nb.old*mu.i)^2
    d2elldalphadbeta.nb.i <- as.numeric(-1.0*mu.i*(y[i]-mu.i)/(1+alpha.nb.old*mu.i)^2)*x.nb.extended[i,]
    
    hessian.nb.i[1:n.beta.nb,1:n.beta.nb] <- d2elldbeta2.nb.i
    hessian.nb.i[n.beta.nb+1,n.beta.nb+1] <- d2elldalpha2.nb.i
    hessian.nb.i[1:n.beta.nb,n.beta.nb+1] <- d2elldalphadbeta.nb.i
    hessian.nb.i[n.beta.nb+1,1:n.beta.nb] <- d2elldalphadbeta.nb.i
    return(hessian.nb.i)
  }
  
  ell.pl.i.fun <- function(beta.pl){
    xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    cdivy <- c.pl/y[i]
    cdivyp1 <- c.pl/(y[i]+1)
    ##OBS!! The restriction that y.i>c.pl will be taken care of outside
   
    ell.pl.i <- log(cdivy^exp.xtb.pl.i - cdivyp1^exp.xtb.pl.i)
    return(ell.pl.i)
  }
  
  delldbeta.pl.i.fun.exact <- function(beta.pl){
    xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    cdivy <- c.pl/y[i]
    cdivyp1 <- c.pl/(y[i]+1)
  
    numerator <- log(cdivy)*exp.xtb.pl.i*cdivy^exp.xtb.pl.i - log(cdivyp1)*exp.xtb.pl.i*cdivyp1^exp.xtb.pl.i
    denominator <- cdivy^exp.xtb.pl.i - cdivyp1^exp.xtb.pl.i
  
    delldbeta.pl.i <- x.pl.extended[i,]*as.numeric(numerator/denominator)
    return(delldbeta.pl.i)
  }
  
  delldbeta.pl.i.fun.approx <- function(beta.pl){
    xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    delldbeta.pl.i <- x.pl.extended[i,]%*%(1 + log(c.pl)*exp.xtb.pl.i - log(y[i])*exp.xtb.pl.i)
    return(delldbeta.pl.i)
  }
  
  d2elldbeta2.pl.i.fun.exact <- function(beta.pl){
    xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    cdivy <- c.pl/y[i]
    cdivyp1 <- c.pl/(y[i]+1)
    a.1 <- cdivy^exp.xtb.pl.i - cdivyp1^exp.xtb.pl.i
    a.2 <- log(cdivy)*exp.xtb.pl.i*cdivy^exp.xtb.pl.i - log(cdivyp1)*exp.xtb.pl.i*cdivyp1^exp.xtb.pl.i
    grad.1.scalar <- log(cdivy)*(1+log(cdivy)*exp.xtb.pl.i)*exp.xtb.pl.i*cdivy^exp.xtb.pl.i   -  log(cdivyp1)*(1 + log(cdivyp1)*exp.xtb.pl.i)*exp.xtb.pl.i*cdivyp1^exp.xtb.pl.i
    grad.2.scalar <- log(cdivy)*exp.xtb.pl.i*cdivy^exp.xtb.pl.i - log(cdivyp1)*exp.xtb.pl.i*cdivyp1^exp.xtb.pl.i 
    denominator.scalar <- (cdivy^exp.xtb.pl.i-cdivyp1^exp.xtb.pl.i)^2
    #hessian.scalar <- a.1*grad.1.scalar-a.2*grad.2.scalar/denominator.scalar
    hessian.scalar <- (a.1*grad.1.scalar-a.2*grad.2.scalar)/denominator.scalar   #Should be like this. Check later. Weird that it went so well anyway
    hessian.pl.i <- x.pl.extended[i,]%*%t(x.pl.extended[i,])*as.numeric(hessian.scalar)
    return(hessian.pl.i)
  }
  
  d2elldbeta2.pl.i.fun.approx <- function(beta.pl){
    xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    hessian.pl.i <- x.pl.extended[i,]%*%t(x.pl.extended[i,])*as.numeric(log(c.pl)*exp.xtb.pl.i - log(y[i])*exp.xtb.pl.i)
    return(hessian.pl.i)
  }
  
  log.lik.fun <- function(beta.zc.multinomial,beta.pl.multinomial,theta.nb,beta.pl,c.pl){
    
    #Model-implied proportions based on covariates only
    props <- matrix(NA,nrow=n,ncol=3)
    for(i in 1:n){
      
      denominator <- 1 + exp(t(beta.zc.multinomial)%*%x.multinom.zc.extended[i,]) + exp(t(beta.pl.multinomial)%*%x.multinom.pl.extended[i,])
      props[i,1] <- exp(t(beta.zc.multinomial)%*%x.multinom.zc.extended[i,])/denominator
      props[i,2] <- 1/denominator
      props[i,3] <-exp(t(beta.pl.multinomial)%*%x.multinom.pl.extended[i,])/denominator
    }
    
    func.val <- 0
    for(i in 1:n){
      beta.nb <- theta.nb[1:n.beta.nb]
      alpha.nb <- theta.nb[n.beta.nb + 1]
      ell.nb.i <- 0
      xtb.nb.i <- x.nb.extended[i,]%*%beta.nb
      mu.i <- exp(xtb.nb.i)
      #sum ln(j + a^-1)
      if(y[i]>0){
        for(j in 0:(y[i]-1)){
          ell.nb.i <- ell.nb.i + log(j + 1/alpha.nb)
        }
      }
      #gamma(y+1)
      if(y[i]>0){
        for(j in 1:y[i]){
          ell.nb.i <- ell.nb.i - log(j)
        }
      }
      ell.nb.i <- ell.nb.i - (1/alpha.nb)*log(1 + alpha.nb*mu.i) - y[i]*log(1 + alpha.nb*mu.i) + y[i]*log(alpha.nb) + y[i]*log(mu.i)
      
      if(y[i]==0){
        func.val <- func.val + log(props[i,1] + props[i,2]*exp(ell.nb.i))
      }else if(y[i]>0 & y[i]<c.pl){
        func.val <- func.val + log(props[i,2]*exp(ell.nb.i))
      }else{
        xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
        exp.xtb.pl.i <- exp(xtb.pl.i)
        cdivy <- c.pl/y[i]
        cdivyp1 <- c.pl/(y[i]+1)
        ell.pl.i <- log(cdivy^exp.xtb.pl.i - cdivyp1^exp.xtb.pl.i)
        
        func.val <- func.val + log(props[i,2]*exp(ell.nb.i) + props[i,3]*exp(ell.pl.i))
      }
    }
    
    return(func.val)
  }
 
  
  #Initialize parameters and data
  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL

  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL
  
  n <- max(nrow(x.nb),nrow(x.pl),nrow(x.multinom.zc),nrow(x.multinom.pl))
  
  if(is.null(dim(x.multinom.zc))){
    x.multinom.zc.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.multinom.zc.extended <- cbind(1,x.multinom.zc)
  }
  if(is.null(dim(x.multinom.pl))){
    x.multinom.pl.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.multinom.pl.extended <- cbind(1,x.multinom.pl)
  }
  if(is.null(dim(x.nb))){
    x.nb.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.nb.extended <- cbind(1,x.nb)
  }
  if(is.null(dim(x.pl))){
    x.pl.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.pl.extended <- cbind(1,x.pl)
  }
  
  n.beta.zc.multinomial <- dim(x.multinom.zc.extended)[2]
  n.beta.pl.multinomial <- dim(x.multinom.pl.extended)[2]
  n.theta.multinomial <- n.beta.zc.multinomial + n.beta.pl.multinomial
  n.beta.nb <- dim(x.nb.extended)[2]
  n.beta.pl <- dim(x.pl.extended)[2]
  
  beta.zc.multinomial.old <- ini.val$Beta.multinom.ZC
  beta.pl.multinomial.old <- ini.val$Beta.multinom.PL
  beta.nb.old <- ini.val$Beta.NB
  alpha.nb.old <- ini.val$Alpha.NB
  theta.nb.old <- c(beta.nb.old,alpha.nb.old)
  beta.pl.old <- ini.val$Beta.PL
  c.pl <- ini.val$C
  
  #Initial log likelihood value
  func.val.initial <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
  
  
  #Model-implied proportions based on covariates only
  props.old <- matrix(NA,nrow=n,ncol=3)
  for(i in 1:n){
    denominator <- 1 + exp(t(beta.zc.multinomial.old)%*%x.multinom.zc.extended[i,]) + exp(t(beta.pl.multinomial.old)%*%x.multinom.pl.extended[i,])
    props.old[i,1] <- exp(t(beta.zc.multinomial.old)%*%x.multinom.zc.extended[i,])/denominator
    props.old[i,2] <- 1/denominator
    props.old[i,3] <-exp(t(beta.pl.multinomial.old)%*%x.multinom.pl.extended[i,])/denominator
    
    
    
  }
  
  
  # The maximum change after one EM step needs to be small in order for the algorithm to stop
  max.abs.par.diff <- 100
  
  
  ################
  #Start EM algorithm
  ###################
  i.em <- 1
  func.val.vec <- c()
  func.val.old <- -1e50
  while(i.em<control$max.no.em.steps & max.abs.par.diff>control$max.diff.par){
    
    props.old <- matrix(NA,nrow=n,ncol=3)
    for(i in 1:n){
      denominator <- 1 + exp(t(beta.zc.multinomial.old)%*%x.multinom.zc.extended[i,]) + exp(t(beta.pl.multinomial.old)%*%x.multinom.pl.extended[i,])
      props.old[i,1] <- exp(t(beta.zc.multinomial.old)%*%x.multinom.zc.extended[i,])/denominator
      props.old[i,2] <- 1/denominator
      props.old[i,3] <-exp(t(beta.pl.multinomial.old)%*%x.multinom.pl.extended[i,])/denominator
    }
  
    ##Set parameters for comparison after one full EM step 
    props.start.em <- props.old
    beta.zc.multinomial.start.em <- beta.zc.multinomial.old
    beta.pl.multinomial.start.em <- beta.pl.multinomial.old
    beta.nb.start.em <- beta.nb.old
    alpha.nb.start.em <- abs(alpha.nb.old)
    beta.pl.start.em <- beta.pl.old
    
    par.start.em <- c(beta.zc.multinomial.start.em,beta.pl.multinomial.start.em,beta.nb.start.em,alpha.nb.start.em,beta.pl.start.em)
    
    
    #Calculate responsibilities (posterior component probabilities, i.e., after observing y)
    resp <- matrix(NA,nrow=n,ncol=3)
    
    for (i in 1:n){
      #print(i)
      #Marginal probability mass function of y and x (sum over components) - marginal.yx.i
      if(y[i]==0){
        marginal.yx.i <- props.old[i,1] + props.old[i,2]*exp(ell.nb.i.fun(theta.nb.old)) 
        resp[i,1] <- props.old[i,1]/marginal.yx.i
        resp[i,2] <- props.old[i,2]*exp(ell.nb.i.fun(theta.nb.old))/marginal.yx.i
        resp[i,3] <- 0
      }else if(y[i]>0 & y[i]<c.pl){
        marginal.yx.i <- props.old[i,2]*exp(ell.nb.i.fun(theta.nb.old)) 
        resp[i,1] <- 0
        resp[i,2] <- 1
        resp[i,3] <- 0
      }else if(y[i]>=c.pl){
        marginal.yx.i <- props.old[i,2]*exp(ell.nb.i.fun(theta.nb.old))  + props.old[i,3]*exp(ell.pl.i.fun(beta.pl.old))
        resp[i,1] <- 0
        resp[i,2] <- props.old[i,2]*exp(ell.nb.i.fun(theta.nb.old))/marginal.yx.i
        resp[i,3] <- props.old[i,3]*exp(ell.pl.i.fun(beta.pl.old))/marginal.yx.i
      }
    }
    
    
    func.val.before.nb <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
    
    #Update theta.nb, i.e., beta.nb and alpha.nb
    theta.nb.old <- c(beta.nb.old,alpha.nb.old)
    
    for(i.bfgs.bn in 1:control$no.m.bfgs.steps.nb){
      theta.nb.old[n.beta.nb+1] <- abs(theta.nb.old[n.beta.nb+1])
      d2Qdtheta2.nb <- matrix(0,nrow=n.beta.nb+1,ncol=n.beta.nb+1)
      dQdtheta.nb <- matrix(0,nrow=n.beta.nb+1,ncol=1)
      
      for(i in 1:n){

        dQdtheta.nb <- dQdtheta.nb + delldtheta.nb.i.fun(theta.nb.old)*resp[i,2]
        d2Qdtheta2.nb <- d2Qdtheta2.nb + d2elldtheta2.nb.i.fun(theta.nb.old)*resp[i,2]
      }
      
      change.nb.bfgs <- -solve(d2Qdtheta2.nb)%*%dQdtheta.nb
      
      if(sum(is.na(change.nb.bfgs))==0){
      
        if(max(abs(change.nb.bfgs))>control$max.upd.par.nb){
          change.nb.bfgs <- control$max.upd.par.nb/max(abs(change.nb.bfgs))*change.nb.bfgs
        }
        
        theta.nb.old <- theta.nb.old + change.nb.bfgs
      }
    }
    beta.nb.old <- theta.nb.old[1:n.beta.nb]
    alpha.nb.old <- abs(theta.nb.old[n.beta.nb+1])
  

    func.val.after.nb <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
    func.val.old <- func.val.after.nb
    
    ################################################
    ##################Optimize theta.nb.old
    nb.log.lik.optim <- function(eta){
      #Go back eta times the step that was already made
      theta.nb.old <- theta.nb.old + eta*change.nb.bfgs
      return(-1.0*log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl))
    }
    
    if(func.val.after.nb<func.val.before.nb){
      
      if(sum(is.na(change.nb.bfgs))==0){
      
      ###########################
        change.nb.obj <- optimise(f=nb.log.lik.optim,interval=control$eta.int)
        eta.nb <- change.nb.obj$minimum
        
        theta.nb.after.optim <- theta.nb.old + eta.nb*change.nb.bfgs
        
        func.val.after.nb.optim <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.after.optim,beta.pl.old,c.pl)
        
        theta.nb.old <- theta.nb.after.optim
        func.val.old <- func.val.after.nb.optim
      }
    }
    ##################################
    ##############################################
    
 
    
    func.val.before.pl <- func.val.old
    #Update pl parameters
    for(i.bfgs.pl in 1:control$no.m.bfgs.steps.pl){
      
      d2Qdbeta2.pl <- matrix(0,nrow=n.beta.pl,ncol=n.beta.pl)
      dQdbeta.pl <- matrix(0,nrow=n.beta.pl,ncol=1)
      
      if(control$pdf.pl.type=="exact"){
      
        for(i in 1:n){
          if(y[i]>=c.pl){
            print(d2elldbeta2.pl.i.fun.exact(beta.pl.old)*resp[i,3])
            dQdbeta.pl <- dQdbeta.pl + delldbeta.pl.i.fun.exact(beta.pl.old)*resp[i,3]
            d2Qdbeta2.pl <- d2Qdbeta2.pl + d2elldbeta2.pl.i.fun.exact(beta.pl.old)*resp[i,3]
          }
        }
      }else{
        for(i in 1:n){
          #print(i)
          if(y[i]>=c.pl){
            dQdbeta.pl <- dQdbeta.pl + delldbeta.pl.i.fun.approx(beta.pl.old)*resp[i,3]
            d2Qdbeta2.pl <- d2Qdbeta2.pl + d2elldbeta2.pl.i.fun.approx(beta.pl.old)*resp[i,3]
          }
        }
      }
      
      
        
      change.pl.bfgs <- -solve(d2Qdbeta2.pl)%*%dQdbeta.pl
      if(max(abs(change.pl.bfgs))>control$max.upd.par.nb){
        change.pl.bfgs <- control$max.upd.par.pl/max(abs(change.pl.bfgs))*change.pl.bfgs
      }
        
      beta.pl.old <- beta.pl.old + change.pl.bfgs
      
      
      
    }
    
    func.val.after.pl <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
    func.val.old <- func.val.after.pl
    ################################################
    ##################Optimize theta.pl.old
    pl.log.lik.optim <- function(eta){
      #Go back eta times the step that was already made
      beta.pl.old <- beta.pl.old + eta*change.pl.bfgs
      return(-1.0*log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl))
    }
    
    if(func.val.after.pl<func.val.before.pl){
    ###########################
      if(det(d2Qdbeta2.pl)>0){
        change.pl.obj <- optimise(f=pl.log.lik.optim,interval=control$eta.int)
        eta.pl <- change.pl.obj$minimum
        
        beta.pl.after.optim <- beta.pl.old + eta.pl*change.pl.bfgs
        
        func.val.after.pl.optim <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.after.optim,c.pl)
        
        beta.pl.old <- beta.pl.after.optim
        func.val.old <- func.val.after.pl.optim
      }
    }
    
    ##################################
    ##############################################
    

    func.val.before.zc.multinomial <- func.val.old
    
    #Update beta.multinomials
    
    for(i.bfgs.multinomial in 1:control$no.m.bfgs.steps.multinomial){
      dQdbeta.zc.multinomial <- matrix(0,nrow=n.beta.zc.multinomial,ncol=1)
      d2Qdbeta2.zc.multinomial <- matrix(0,nrow=n.beta.zc.multinomial,ncol=n.beta.zc.multinomial)
      for(i in 1:n){
        denominator <- as.numeric(1 + exp(t(beta.zc.multinomial.old)%*%x.multinom.zc.extended[i,]) + exp(t(beta.pl.multinomial.old)%*%x.multinom.pl.extended[i,]))
        dQdbeta.zc.multinomial.i <- as.matrix(x.multinom.zc.extended[i,])*as.numeric(resp[i,1] - exp(t(beta.zc.multinomial.old)%*%x.multinom.zc.extended[i,])/denominator)
        
        d2Qdbeta2.zc.multinomial.i <- -t(x.multinom.zc.extended[i,,drop=FALSE])%*%x.multinom.zc.extended[i,,drop=FALSE]/denominator
        
        dQdbeta.zc.multinomial <- dQdbeta.zc.multinomial + dQdbeta.zc.multinomial.i
        d2Qdbeta2.zc.multinomial <- d2Qdbeta2.zc.multinomial + d2Qdbeta2.zc.multinomial.i
      }
      
      change.zc.multinomial.bfgs <- -solve(d2Qdbeta2.zc.multinomial)%*%dQdbeta.zc.multinomial
      
      if(max(abs(change.zc.multinomial.bfgs))>control$max.upd.par.zc.multinomial){
        change.zc.multinomial.bfgs <- control$max.upd.par.zc.multinomial/max(abs(change.zc.multinomial.bfgs))*change.zc.multinomial.bfgs
      }
      
  
      beta.zc.multinomial.old <- beta.zc.multinomial.old + change.zc.multinomial.bfgs
    
    }
    
    func.val.after.zc.multinomial <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
    func.val.old <- func.val.after.zc.multinomial
    ################################################
    ##################Optimize theta.pl.old
    zc.multinomial.log.lik.optim <- function(eta){
      #Go back eta times the step that was already made
      beta.zc.multinomial.old <- beta.zc.multinomial.old + eta*change.zc.multinomial.bfgs
      return(-1.0*log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl))
    }
    
    
    if(func.val.after.zc.multinomial<func.val.before.zc.multinomial){
      ###########################
      change.zc.multinomial.obj <- optimise(f=zc.multinomial.log.lik.optim,interval=control$eta.int)
      eta.zc.multinomial <- change.zc.multinomial.obj$minimum
      
      beta.zc.multinomial.after.optim <- beta.zc.multinomial.old + eta.zc.multinomial*change.zc.multinomial.bfgs
      
      func.val.after.zc.multinomial.optim <- log.lik.fun(beta.zc.multinomial.after.optim,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
      
      beta.zc.multinomial.old <- beta.zc.multinomial.after.optim
      func.val.old <- func.val.after.zc.multinomial.optim
    }
    
    ##################################
    ##############################################
    

    
    func.val.before.pl.multinomial <- func.val.old
    
    for(i.bfgs.multinomial in 1:control$no.m.bfgs.steps.multinomial){
      dQdbeta.pl.multinomial <- matrix(0,nrow=n.beta.pl.multinomial,ncol=1)
      d2Qdbeta2.pl.multinomial <- matrix(0,nrow=n.beta.pl.multinomial,ncol=n.beta.pl.multinomial)
      for(i in 1:n){
        denominator <- as.numeric(1 + exp(t(beta.zc.multinomial.old)%*%x.multinom.zc.extended[i,]) + exp(t(beta.pl.multinomial.old)%*%x.multinom.pl.extended[i,]))
        dQdbeta.pl.multinomial.i <- as.matrix(x.multinom.pl.extended[i,])*as.numeric(resp[i,3] - exp(t(beta.pl.multinomial.old)%*%x.multinom.pl.extended[i,])/denominator)
        
        d2Qdbeta2.pl.multinomial.i <- -t(x.multinom.pl.extended[i,,drop=FALSE])%*%x.multinom.pl.extended[i,,drop=FALSE]/denominator
        
        dQdbeta.pl.multinomial <- dQdbeta.pl.multinomial + dQdbeta.pl.multinomial.i
        d2Qdbeta2.pl.multinomial <- d2Qdbeta2.pl.multinomial + d2Qdbeta2.pl.multinomial.i
      }
      

      change.pl.multinomial.bfgs <- -solve(d2Qdbeta2.pl.multinomial)%*%dQdbeta.pl.multinomial
      
      if(max(abs(change.pl.multinomial.bfgs))>control$max.upd.par.pl.multinomial){
        change.pl.multinomial.bfgs <- control$max.upd.par.pl.multinomial/max(abs(change.pl.multinomial.bfgs))*change.pl.multinomial.bfgs
      }
      
      beta.pl.multinomial.old <- beta.pl.multinomial.old + change.pl.multinomial.bfgs
      
      
    }
    
    
    func.val.after.pl.multinomial <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
    
    ################################################
    ##################Optimize theta.pl.old
    pl.multinomial.log.lik.optim <- function(eta){
      #Go back eta times the step that was already made
      beta.pl.multinomial.old <- beta.pl.multinomial.old + eta*change.pl.multinomial.bfgs
      return(-1.0*log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl))
    }
    
  
    ###########################
    if(func.val.after.pl.multinomial<func.val.before.pl.multinomial){
      change.pl.multinomial.obj <- optimise(f=pl.multinomial.log.lik.optim,interval=control$eta.int)
      eta.pl.multinomial <- change.pl.multinomial.obj$minimum
      
      beta.pl.multinomial.after.optim <- beta.pl.multinomial.old + eta.pl.multinomial*change.pl.multinomial.bfgs
      
      func.val.after.pl.multinomial.optim <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.after.optim,theta.nb.old,beta.pl.old,c.pl)
      
      beta.pl.multinomial.old <- beta.pl.multinomial.after.optim
      func.val.old <- func.val.after.pl.multinomial.optim
    }
    ##################################
    ##############################################
    
    
    
    ########################
    ############################## Finished updating parameters
    
    par.end.em <- c(beta.zc.multinomial.old,beta.pl.multinomial.old,beta.nb.old,alpha.nb.old,beta.pl.old)

    par.diff <- par.end.em - par.start.em 
    max.abs.par.diff <- max(abs(par.diff))
    
    
    
 
    # #Calculate the log lik after initialization of parameters
    # func.val <- 0
    # for(i in 1:n){
    #   if(y[i]==0){
    #     func.val <- func.val + log(props.old[i,1] + props.old[i,2]*exp(ell.nb.i.fun(c(beta.nb.old,alpha.nb.old))))
    #   }else if(y[i]>0 & y[i]<c.pl){
    #     func.val <- func.val + log(props.old[i,2]*exp(ell.nb.i.fun(c(beta.nb.old,alpha.nb.old))))
    #   }else{
    #     func.val <- func.val + log(props.old[i,2]*exp(ell.nb.i.fun(c(beta.nb.old,alpha.nb.old))) + props.old[i,3]*exp(ell.pl.i.fun(beta.pl.old)))
    #   }
    # }

    #par.all.tmp <- c(beta.nb.multinomial.old,beta.pl.multinomial.old,beta.nb.old,alpha.nb.old,beta.pl.old)
    
    #if(func.val<func.val.old){
    #  print("The function decreased")
      # par.diff.mod <- 0.01*par.diff
      # par.all.tmp <- par.start.em + par.diff.mod
      # beta.nb.multinomial.old <- par.all.tmp[1:n.beta.nb.multinomial]
      # par.all.tmp <- par.all.tmp[-c(1:n.beta.nb.multinomial)]
      # beta.pl.multinomial.old <- par.all.tmp[1:n.beta.pl.multinomial]
      # par.all.tmp <- par.all.tmp[-c(1:n.beta.pl.multinomial)]
      # beta.nb.old <- par.all.tmp[1:n.beta.nb]
      # par.all.tmp <- par.all.tmp[-c(1:n.beta.nb)]
      # alpha.nb.old <- par.all.tmp[1]
      # par.all.tmp <- par.all.tmp[-1]
      # beta.pl.old <- par.all.tmp
      # rm(par.all.tmp)
    #}    

    func.val.old <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
        
    #func.val.old <- func.val
    
    cat('The max abs diff in parameters is ',max.abs.par.diff , '. The function value is ', round(func.val.old,4), '\n', sep='')
    
    func.val.vec[i.em] <- func.val.old
    i.em <- i.em + 1
    # print(beta.nb.multinomial.old)
    # print(beta.pl.multinomial.old)
    # print(beta.nb.old)
    # print(beta.pl.old)
    # print(alpha.nb.old)
  }
  
  
  #Gather results after finishing the estimation
  res <- list()
  par.mat <- list()
  par.mat$Props <- props.old
  par.mat$Beta.multinom.ZC <- beta.zc.multinomial.old
  par.mat$Beta.multinom.PL <- beta.pl.multinomial.old
  par.mat$Beta.NB <- beta.nb.old
  par.mat$Alpha.NB <- alpha.nb.old
  par.mat$Beta.PL <- as.numeric(beta.pl.old)
  par.mat$C <- c.pl
  
  res$par.mat <- par.mat
  
  res$par.all <- par.end.em
  
  res$resp <- resp

  res$log.lik.vec <- func.val.vec
  res$log.lik <- func.val.old
  
  if(i.em<control$max.no.em.steps){
    res$converge <- TRUE
  }else{
    res$converge <- FALSE
  }
  return(res)
  
}


# x.obj <- OBS.X.obj
# y <- OBS.Y
# control <- Control
# ini.val <- Ini.Val
# # ini.val <- TruePar

zerinfl.nb.pl.regression.fun <- function(y,x.obj,ini.val,control){

  #Functions
  ell.nb.i.fun <- function(theta.nb){
    
    beta.nb <- theta.nb[1:n.beta.nb]
    alpha.nb <- theta.nb[n.beta.nb + 1]
    
    ell.nb.i <- 0
    xtb.nb.i <- x.nb.extended[i,]%*%beta.nb
    mu.i <- exp(xtb.nb.i)
    
    #sum ln(j + a^-1)
    if(y[i]>0){
      for(j in 0:(y[i]-1)){
        ell.nb.i <- ell.nb.i + log(j + 1/alpha.nb)
      }
    }
    
    #gamma(y+1)
    if(y[i]>0){
      for(j in 1:y[i]){
        ell.nb.i <- ell.nb.i - log(j)
      }
    }
    
    ell.nb.i <- ell.nb.i - (1/alpha.nb)*log(1 + alpha.nb*mu.i) - y[i]*log(1 + alpha.nb*mu.i) + y[i]*log(alpha.nb) + y[i]*log(mu.i)
    
    
    return(ell.nb.i)
  }
  
  delldtheta.nb.i.fun <- function(theta.nb){
    beta.nb <- theta.nb[1:n.beta.nb]
    alpha.nb <- theta.nb[n.beta.nb + 1]
    xtb.nb.i <- x.nb.extended[i,]%*%beta.nb
    mu.i <- exp(xtb.nb.i)
    
    if(y[i]==0){
      delldalpha.nb.i <- log(1 + alpha.nb*mu.i)/alpha.nb^2 - mu.i/(alpha.nb*(1+alpha.nb*mu.i))
    }else{
      
      delldalpha.nb.i <- log(1 + alpha.nb*mu.i)
      for(j in 0:(y[i]-1)){
        delldalpha.nb.i <- delldalpha.nb.i - 1/(j+1/alpha.nb)
      }
      delldalpha.nb.i <- delldalpha.nb.i/alpha.nb^2
      
      delldalpha.nb.i <- delldalpha.nb.i + (y[i]-mu.i)/(alpha.nb*(1+alpha.nb*mu.i))
    }
    delldbeta.nb.i <- as.numeric((y[i]-mu.i)/(1+alpha.nb*mu.i))*x.nb.extended[i,]
    
    delldtheta.nb.i <- c(delldbeta.nb.i,delldalpha.nb.i)
    return(delldtheta.nb.i)
  }
  
  d2elldtheta2.nb.i.fun <- function(theta.nb){
    beta.nb <- theta.nb[1:n.beta.nb]
    alpha.nb <- theta.nb[n.beta.nb + 1]
    xtb.nb.i <- x.nb.extended[i,]%*%beta.nb
    mu.i <- exp(xtb.nb.i)
    
    hessian.nb.i <- matrix(NA,nrow=length(theta.nb),ncol=length(theta.nb))
    
    d2elldbeta2.nb.i <- as.numeric(-1.0*mu.i*(1+alpha.nb.old*y[i])/(1+alpha.nb.old*mu.i)^2)*x.nb.extended[i,]%*%t(x.nb.extended[i,])
    
    d2elldalpha2.nb.i <- 0
    
    if(y[i]>0){
      for(j in 0:(y[i]-1)){
        d2elldalpha2.nb.i <- d2elldalpha2.nb.i - (j/(1+alpha.nb.old*j))^2
      }
    }
    
    d2elldalpha2.nb.i <- d2elldalpha2.nb.i - 2/alpha.nb.old^3*log(1 + alpha.nb.old*mu.i) + (2*(1/alpha.nb.old)^2*mu.i)/(1+alpha.nb.old*mu.i) + (y[i]+1/alpha.nb.old)*mu.i^2/(1+alpha.nb.old*mu.i)^2
    d2elldalphadbeta.nb.i <- as.numeric(-1.0*mu.i*(y[i]-mu.i)/(1+alpha.nb.old*mu.i)^2)*x.nb.extended[i,]
    
    hessian.nb.i[1:n.beta.nb,1:n.beta.nb] <- d2elldbeta2.nb.i
    hessian.nb.i[n.beta.nb+1,n.beta.nb+1] <- d2elldalpha2.nb.i
    hessian.nb.i[1:n.beta.nb,n.beta.nb+1] <- d2elldalphadbeta.nb.i
    hessian.nb.i[n.beta.nb+1,1:n.beta.nb] <- d2elldalphadbeta.nb.i
    return(hessian.nb.i)
  }
  
  ell.pl.i.fun <- function(beta.pl){
    xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    cdivy <- c.pl/y[i]
    cdivyp1 <- c.pl/(y[i]+1)
    ##OBS!! The restriction that y.i>c.pl will be taken care of outside
    ell.pl.i <- log(cdivy^exp.xtb.pl.i - cdivyp1^exp.xtb.pl.i)
    return(ell.pl.i)
  }
  
  delldbeta.pl.i.fun.exact <- function(beta.pl){
    xtb.pl.i <- x.pl.extendqed[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    cdivy <- c.pl/y[i]
    cdivyp1 <- c.pl/(y[i]+1)
    
    numerator <- log(cdivy)*exp.xtb.pl.i*cdivy^exp.xtb.pl.i - log(cdivyp1)*exp.xtb.pl.i*cdivyp1^exp.xtb.pl.i
    denominator <- cdivy^exp.xtb.pl.i - cdivyp1^exp.xtb.pl.i
    
    delldbeta.pl.i <- x.pl.extended[i,]*as.numeric(numerator/denominator)
    return(delldbeta.pl.i)
  }
  
  delldbeta.pl.i.fun.approx <- function(beta.pl){
    xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    delldbeta.pl.i <- x.pl.extended[i,]%*%(1 + log(c.pl)*exp.xtb.pl.i - log(y[i])*exp.xtb.pl.i)
    return(delldbeta.pl.i)
  }
  
  d2elldbeta2.pl.i.fun.exact <- function(beta.pl){
    xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    cdivy <- c.pl/y[i]
    cdivyp1 <- c.pl/(y[i]+1)
    a.1 <- cdivy^exp.xtb.pl.i - cdivyp1^exp.xtb.pl.i
    a.2 <- log(cdivy)*exp.xtb.pl.i*cdivy^exp.xtb.pl.i - log(cdivyp1)*exp.xtb.pl.i*cdivyp1^exp.xtb.pl.i
    grad.1.scalar <- log(cdivy)*(1+log(cdivy)*exp.xtb.pl.i)*exp.xtb.pl.i*cdivy^exp.xtb.pl.i   -  log(cdivyp1)*(1 + log(cdivyp1)*exp.xtb.pl.i)*exp.xtb.pl.i*cdivyp1^exp.xtb.pl.i
    grad.2.scalar <- log(cdivy)*exp.xtb.pl.i*cdivy^exp.xtb.pl.i - log(cdivyp1)*exp.xtb.pl.i*cdivyp1^exp.xtb.pl.i 
    denominator.scalar <- (cdivy^exp.xtb.pl.i-cdivyp1^exp.xtb.pl.i)^2
    #hessian.scalar <- a.1*grad.1.scalar-a.2*grad.2.scalar/denominator.scalar
    hessian.scalar <- (a.1*grad.1.scalar-a.2*grad.2.scalar)/denominator.scalar   #Should be like this. Check later. Weird that it went so well anyway
    hessian.pl.i <- x.pl.extended[i,]%*%t(x.pl.extended[i,])*as.numeric(hessian.scalar)
    return(hessian.pl.i)
  }
  
  d2elldbeta2.pl.i.fun.approx <- function(beta.pl){
    xtb.pl.i <- x.pl.extended[i,]%*%beta.pl
    exp.xtb.pl.i <- exp(xtb.pl.i)
    hessian.pl.i <- x.pl.extended[i,]%*%t(x.pl.extended[i,])*as.numeric(log(c.pl)*exp.xtb.pl.i - log(y[i])*exp.xtb.pl.i)
    return(hessian.pl.i)
  }
  
  #Initialize data and parameters
  
  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL
  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL
  
  n <- max(nrow(x.nb),nrow(x.pl),nrow(x.multinom.zc),nrow(x.multinom.pl))
  
  if(is.null(dim(x.multinom.zc))){
    x.multinom.zc.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.multinom.zc.extended <- cbind(1,x.multinom.zc)
  }
  if(is.null(dim(x.multinom.pl))){
    x.multinom.pl.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.multinom.pl.extended <- cbind(1,x.multinom.pl)
  }
  if(is.null(dim(x.nb))){
    x.nb.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.nb.extended <- cbind(1,x.nb)
  }
  if(is.null(dim(x.pl))){
    x.pl.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.pl.extended <- cbind(1,x.pl)
  }
  
  n.beta.zc.multinomial <- dim(x.multinom.zc.extended)[2]
  n.beta.pl.multinomial <- dim(x.multinom.pl.extended)[2]
  n.theta.multinomial <- n.beta.zc.multinomial + n.beta.pl.multinomial
  n.beta.nb <- dim(x.nb.extended)[2]
  n.beta.pl <- dim(x.pl.extended)[2]
  
  prel.val <- ini.val

  n <- max(nrow(x.nb),nrow(x.pl),nrow(x.multinom.zc),nrow(x.multinom.pl))
  
  #Determine which values of c to be investigated
  c.range <- unique(sort(y))[unique(sort(y))>=control$c.lim[1] & unique(sort(y))<=control$c.lim[2]]  
  
  ###When c is not well known we only run short versions of zerinfl.nb.pl.reg.cond.c.fun() in order to save time
  control.warmup <- control
  control.warmup$max.no.em.steps <- control$max.no.em.steps.warmup
  
  ###Initializing the warm-up phase when the update stops after control$max.no.em.steps.warmup EM steps
  #If the update of c is less than 1, stop the warmup
  c.abs.diff <- 100
  log.lik.vec.all <- NULL
  cat('Begin warm-up', '\n', sep='')
  while(c.abs.diff>0){
    log.lik.vec <- c()
    est.obj <- zerinfl.nb.pl.reg.cond.c.fun(y,x.obj,prel.val,control.warmup)
    log.lik.vec.all <- c(log.lik.vec.all,est.obj$log.lik.vec)
    prel.val <- est.obj$par.mat
    props.old <- prel.val$Props
    beta.nb.old <- prel.val$Beta.NB
    alpha.nb.old <- prel.val$Alpha.NB
    beta.pl.old <- prel.val$Beta.PL
    c.pl <- prel.val$C
  
    
    for(c.index in 1:length(c.range)){
      c.pl <- c.range[c.index]
      
      
      
      
      func.val <- 0
      for(i in 1:n){
        if(y[i]==0){
          func.val <- func.val + log(props.old[i,1] + props.old[i,2]*exp(ell.nb.i.fun(c(beta.nb.old,alpha.nb.old))))
        }else if(y[i]>0 & y[i]<c.pl){
          func.val <- func.val + log(props.old[i,2]*exp(ell.nb.i.fun(c(beta.nb.old,alpha.nb.old))))
        }else{
          func.val <- func.val + log(props.old[i,2]*exp(ell.nb.i.fun(c(beta.nb.old,alpha.nb.old))) + props.old[i,3]*exp(ell.pl.i.fun(beta.pl.old)))
        }
      }
      log.lik.vec[c.index] <- func.val
    }
    
    c.pl.new <- c.range[which(log.lik.vec==max(log.lik.vec))]
    
    c.abs.diff <- abs(c.pl.new - prel.val$C)
    
    prel.val$C <- c.pl.new
    func.val <- max(log.lik.vec)
    
    log.lik.vec.all <- c(log.lik.vec.all,func.val)
    
    cat('The new c is ',c.pl.new , '. The function value is ', round(func.val,4), '\n', sep='')
    prel.val$C <- c.pl.new
  }
  
  
  
  ###The end phase, when c is presumably rather accurate
  c.abs.diff <- 100
  cat('End warm-up. Run until convergence', '\n', sep='')
  while(c.abs.diff>0){
    log.lik.vec <- c()
    
    est.obj <- zerinfl.nb.pl.reg.cond.c.fun(y,x.obj,prel.val,control)
    log.lik.vec.all <- c(log.lik.vec.all,est.obj$log.lik.vec)
    prel.val <- est.obj$par.mat
    props.old <- prel.val$Props
    beta.nb.old <- prel.val$Beta.NB
    alpha.nb.old <- prel.val$Alpha.NB
    beta.pl.old <- prel.val$Beta.PL
    c.pl <- prel.val$C
    
    
    for(c.index in 1:length(c.range)){
      c.pl <- c.range[c.index]
      func.val <- 0
      for(i in 1:n){
        if(y[i]==0){
          func.val <- func.val + log(props.old[i,1] + props.old[i,2]*exp(ell.nb.i.fun(c(beta.nb.old,alpha.nb.old))))
        }else if(y[i]>0 & y[i]<c.pl){
          func.val <- func.val + log(props.old[i,2]*exp(ell.nb.i.fun(c(beta.nb.old,alpha.nb.old))))
        }else{
          func.val <- func.val + log(props.old[i,2]*exp(ell.nb.i.fun(c(beta.nb.old,alpha.nb.old))) + props.old[i,3]*exp(ell.pl.i.fun(beta.pl.old)))
        }
      }
      log.lik.vec[c.index] <- func.val
    }
    
    c.pl.new <- c.range[which(log.lik.vec==max(log.lik.vec))]
    
    c.abs.diff <- abs(c.pl.new - prel.val$C)
    
    prel.val$C <- c.pl.new
    func.val <- max(log.lik.vec)
    log.lik.vec.all <- c(log.lik.vec.all,func.val)
    
    cat('The new c is ',c.pl.new , '. The function value is ', round(func.val,4), '\n', sep='')
    prel.val$C <- c.pl.new
  }
  
  final.val <- prel.val
  
  #Mean conditional on negative binomial component 
  #Pareto shape parameter 
  mu.nb.vec <- c()
  alpha.pl.vec <- c()
  mean.pl.vec <- c()
  exp.E.log.y <- c()  #exp(E(log(y))), a third alternative for prediction
  median.pl.vec <- c()
  E.inv.y <- c()
  y.hat.plmedian <- c()
  y.hat.plmean <- c()
  y.hat.plexpElogy <- c()
  y.hat.pl.E.inv.y <- c()
  for(i in 1:n){
    mu.nb.vec[i] <- exp(x.nb.extended[i,]%*%prel.val$Beta.NB)
    alpha.pl.vec[i] <- exp(x.pl.extended[i,]%*%prel.val$Beta.PL)
    exp.E.log.y[i] <- c.pl.new*exp(1/alpha.pl.vec[i])
    E.inv.y[i] <- alpha.pl.vec[i]/(c.pl.new*(alpha.pl.vec[i]+1))
    if(alpha.pl.vec[i]>1){
      mean.pl.vec[i] <- alpha.pl.vec[i]*c.pl.new/(alpha.pl.vec[i]-1)  
    }else{
      mean.pl.vec[i] <- NA
    }
    median.pl.vec[i] <- median.pl.vec[i] <- c.pl.new*(2)^(1/alpha.pl.vec[i])
    y.hat.plmedian[i] <- props.old[i,2]*mu.nb.vec[i] + props.old[i,3]*median.pl.vec[i]
    y.hat.plmean[i] <- props.old[i,2]*mu.nb.vec[i] + props.old[i,3]*mean.pl.vec[i]
    y.hat.plexpElogy[i] <- props.old[i,2]*mu.nb.vec[i] + props.old[i,3]*exp.E.log.y[i]
    y.hat.pl.E.inv.y[i] <- props.old[i,2]*mu.nb.vec[i] + props.old[i,3]/E.inv.y[i]
  }

  par.all <- c(final.val$Beta.multinom.ZC,final.val$Beta.multinom.PL,final.val$Beta.NB,final.val$Alpha.NB,final.val$Beta.PL,final.val$C)
  n.par <- length(par.all)  
  BIC <-  log(n)*n.par -2*func.val
  AIC <- 2*n.par - 2*func.val
  
  #Gather results from estimation
  out <- list()
  
  out$control <- control
  out$par.mat <- final.val
  out$log.lik.vec.all <- log.lik.vec.all
  out$log.lik <- log.lik.vec.all[length(log.lik.vec.all)]
  out$resp <- est.obj$resp
  out$converge <- est.obj$converge
  out$ini.val <- ini.val
  out$x.nb <- x.nb
  out$x.pl <- x.pl
  out$x.multinom.zc <- x.multinom.zc
  out$x.multinom.pl <- x.multinom.pl
  out$median.pl.vec <- median.pl.vec
  out$mean.pl.vec <- mean.pl.vec
  out$y <- y
  out$y.hat.plmedian <- y.hat.plmedian
  out$y.hat.plmean <- y.hat.plmean
  out$y.hat.plexpElogy <- y.hat.plexpElogy
  out$mu.nb.vec <- mu.nb.vec
  out$alpha.pl.vec <- alpha.pl.vec
  out$mean.pl.vec <- mean.pl.vec
  out$median.pl.vec <- median.pl.vec
  out$exp.E.log.y <- exp.E.log.y
  out$y.hat.pl.E.inv.y <- y.hat.pl.E.inv.y
  out$par.all <- par.all
  out$BIC <- BIC
  out$AIC <- AIC
  
  return(out)
}

#est.par <- Est.Obj$par.mat
#x.obj <- OBS.X.obj

#Takes parameters and x as input and predicts y (using mean or median of the pl component)
prediction.znbpl.fun <- function(x.obj,est.par){
  
  

  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL
  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL
  
  n <- max(nrow(x.multinom.zc),nrow(x.multinom.pl),nrow(x.nb),nrow(x.pl))
  
  if(is.null(x.multinom.zc)){
    x.multinom.zc.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.multinom.zc.extended <- cbind(1,x.multinom.zc)
  }
  if(is.null(x.multinom.pl)){
    x.multinom.pl.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.multinom.pl.extended <- cbind(1,x.multinom.pl)
  }
  if(is.null(x.nb)){
    x.nb.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.nb.extended <- cbind(1,x.nb)
  }
  if(is.null(x.pl)){
    x.pl.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.pl.extended <- cbind(1,x.pl)
  }
  
  n.beta.zc.multinomial <- dim(x.multinom.zc.extended)[2]
  n.beta.pl.multinomial <- dim(x.multinom.pl.extended)[2]
  n.theta.multinomial <- n.beta.zc.multinomial + n.beta.pl.multinomial
  n.beta.nb <- dim(x.nb.extended)[2]
  n.beta.pl <- dim(x.pl.extended)[2]
  
  beta.zc.multinomial.old <- est.par$Beta.multinom.ZC
  beta.pl.multinomial.old <- est.par$Beta.multinom.PL
  beta.nb.old <- est.par$Beta.NB
  alpha.nb.old <- est.par$Alpha.NB
  theta.nb.old <- c(beta.nb.old,alpha.nb.old)
  beta.pl.old <- est.par$Beta.PL
  c.pl <- est.par$C
  
  props.old <- matrix(NA,nrow=n,ncol=3)
  for(i in 1:n){
    denominator <- 1 + exp(t(beta.zc.multinomial.old)%*%x.multinom.zc.extended[i,]) + exp(t(beta.pl.multinomial.old)%*%x.multinom.pl.extended[i,])
    props.old[i,1] <- exp(t(beta.zc.multinomial.old)%*%x.multinom.zc.extended[i,])/denominator
    props.old[i,2] <- 1/denominator
    props.old[i,3] <-exp(t(beta.pl.multinomial.old)%*%x.multinom.pl.extended[i,])/denominator
  }
  
  #Mean conditional on negative binomial component 
  #Pareto shape parameter 
  mu.nb.vec <- c()
  alpha.pl.vec <- c()
  mean.pl.vec <- c()
  median.pl.vec <- c()
  exp.E.log.y <- c()
  E.inv.y <- c()
  y.hat.plmedian <- c()
  y.hat.plmean <- c()
  y.hat.plexpElogy <- c()
  y.hat.pl.E.inv.y <- c()
  for(i in 1:n){
    mu.nb.vec[i] <- exp(x.nb.extended[i,]%*%est.par$Beta.NB)
    alpha.pl.vec[i] <- exp(x.pl.extended[i,]%*%est.par$Beta.PL)
    exp.E.log.y[i] <- c.pl*exp(1/alpha.pl.vec[i])
    E.inv.y[i] <- alpha.pl.vec[i]/(c.pl*(alpha.pl.vec[i]+1))
    if(alpha.pl.vec[i]>1){
      mean.pl.vec[i] <- alpha.pl.vec[i]*c.pl/(alpha.pl.vec[i]-1)  
    }else{
      mean.pl.vec[i] <- NA
    }
    median.pl.vec[i] <- c.pl*(2)^(1/alpha.pl.vec[i])
    y.hat.plmedian[i] <- props.old[i,2]*mu.nb.vec[i] + props.old[i,3]*median.pl.vec[i]
    y.hat.plmean[i] <- props.old[i,2]*mu.nb.vec[i] + props.old[i,3]*mean.pl.vec[i]
    y.hat.plexpElogy[i] <- props.old[i,2]*mu.nb.vec[i] + props.old[i,3]*exp.E.log.y[i]
    y.hat.pl.E.inv.y[i] <- props.old[i,2]*mu.nb.vec[i] + props.old[i,3]/E.inv.y[i]
  }
  out <- list()
  out$y.hat.plmean <- y.hat.plmean
  out$y.hat.plmedian <- y.hat.plmedian
  out$y.hat.plexpElogy <- y.hat.plexpElogy
  out$y.hat.pl.E.inv.y <- y.hat.pl.E.inv.y
  return(out)
}


#znb <- Znb.Obj
#x <- OBS.X.obj$X.multinom.ZC
prediction.znb.fun <- function(x,znb){
  coefficients <- znb$coefficients
  theta <- znb$theta
  n <- dim(x)[1]
  
  if(is.null(x)){
    x.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.extended <- cbind(1,x)
  }
  
  n.beta <- dim(x.extended)[2]
 
  beta.zc.multinomial.old <- coefficients$zero
  beta.nb.old <- coefficients$count
  alpha.nb.old <- 1/theta   #Not sure. Could be 1/theta
  theta.nb.old <- c(beta.nb.old,alpha.nb.old)
  
  
  props.old <- matrix(NA,nrow=n,ncol=2)
  for(i in 1:n){
    denominator <- 1 + exp(t(beta.zc.multinomial.old)%*%x.extended[i,])
    props.old[i,1] <- exp(t(beta.zc.multinomial.old)%*%x.extended[i,])/denominator
    props.old[i,2] <- 1/denominator
  }
  
  #Mean conditional on negative binomial component 
  #Pareto shape parameter 
  mu.nb.vec <- c()
  y.hat.mean <- c()
  
  
  for(i in 1:n){
    mu.nb.vec[i] <- exp(x.extended[i,]%*%beta.nb.old)
    y.hat.mean[i] <- props.old[i,2]*mu.nb.vec[i]

  }
  return(y.hat.mean)
}


#nb <- nb.obj
#x <- OBS.X.ins
prediction.nb.fun <- function(x,nb){
  coefficients <- nb$coefficients
  theta <- nb$theta
  n <- dim(x)[1]
  
  if(is.null(x)){
    x.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.extended <- cbind(1,x)
  }
  
  n.beta <- dim(x.extended)[2]
  
  beta.nb.old <- coefficients
  alpha.nb.old <- 1/theta   #Not sure. Could be 1/theta
  theta.nb.old <- c(beta.nb.old,alpha.nb.old)
  
  
  #Mean conditional on negative binomial component 
  #Pareto shape parameter 
  mu.nb.vec <- c()
  y.hat.mean <- c()
  
  
  for(i in 1:n){
    mu.nb.vec[i] <- exp(x.extended[i,]%*%beta.nb.old)
    y.hat.mean[i] <- mu.nb.vec[i]
    
  }
  return(y.hat.mean)
}



marginal.effect.nb.fun <- function(sign.level,j,x,x.lim,dx,beta.nb){
  n <- dim(x)[1]
  
  beta.nb <- as.matrix(beta.nb)
  
  if(is.null(x)){
    x.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.extended <- cbind(1,x)
  }
  
  j <- j + 1   #To compensate for the intercept
  
  beta.nb.j <- beta.nb[j]  
  
  #i <- 1  #individual by individual
  x.j.all <- seq(x.lim[1],x.lim[2],by=dx)
  nx.j <- length(x.j.all)
  
  margEff.nb <- matrix(NA,nrow=n,ncol=nx.j)
  for(l in 1:nx.j){
    for(i in 1:n){
      x.il.extended <- as.matrix(x.extended[i,])
      x.il.extended[j] <- x.j.all[l]    
      exp.beta.nb.x <- exp(t(beta.nb)%*%x.il.extended)
      margEff.nb[i,l] <- exp.beta.nb.x
    }  
    
  }
  mean.margEff.nb <- colMeans(margEff.nb)
  sort.margEff.nb <- matrix(NA,nrow=n,ncol=nx.j)
  upper.conf.nb <- c()
  lower.conf.nb <- c()
  lower.index <- round(0.5*sign.level*n)
  upper.index <- round((1-0.5*sign.level)*n)
  for(l in 1:nx.j){
    sort.margEff.nb[,l] <- sort(margEff.nb[,l])
    lower.conf.nb[l] <- sort.margEff.nb[lower.index,l]
    upper.conf.nb[l] <- sort.margEff.nb[upper.index,l]
  }
  out <- list()
  out$x <- x.j.all
  out$mean.margEff.nb <- mean.margEff.nb
  out$lower.conf.nb <- lower.conf.nb
  out$upper.conf.nb <- upper.conf.nb
  return(out)
}

#j is an integer representing the chosen covariate
#dx is the interval in the chosen covariate
#x.lim is the min and max of the chosen covariate
marginal.effect.znb.fun <- function(sign.level,j,x,x.lim,dx,gamma.zc,beta.nb){
  n <- dim(x)[1]
  
  beta.nb <- as.matrix(beta.nb)
  gamma.zc <- as.matrix(gamma.zc)

  if(is.null(x)){
    x.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.extended <- cbind(1,x)
  }
  
  j <- j + 1   #To compensate for the intercept
  
  beta.nb.j <- beta.nb[j]  
  gamma.zc.j <- gamma.zc[j]
  
    #i <- 1  #individual by individual
  x.j.all <- seq(x.lim[1],x.lim[2],by=dx)
  nx.j <- length(x.j.all)
    
  margEff.znb <- matrix(NA,nrow=n,ncol=nx.j)
  prop.z <- matrix(NA,nrow=n,ncol=nx.j)
  prop.nb <- matrix(NA,nrow=n,ncol=nx.j)
  for(l in 1:nx.j){
    for(i in 1:n){
      x.il.extended <- as.matrix(x.extended[i,])
      x.il.extended[j] <- x.j.all[l]    
      exp.gamma.zc.x <- exp(t(gamma.zc)%*%x.il.extended)
      exp.beta.nb.x <- exp(t(beta.nb)%*%x.il.extended)
      margEff.znb[i,l] <- exp.beta.nb.x/(1+exp.gamma.zc.x)
      #margEff.znb[i,l] <- beta.nb.j - gamma.zc.j*exp.gamma.zc.x/(1+exp.gamma.zc.x)
      prop.nb[i,l] <- 1/(exp.gamma.zc.x+1)
      prop.z[i,l] <- 1-prop.nb[i,l]
    }  
    
  }
  mean.margEff.znb <- colMeans(margEff.znb)
  sort.margEff.znb <- matrix(NA,nrow=n,ncol=nx.j)
  mean.prop.z <- colMeans(prop.z)
  mean.prop.nb <- colMeans(prop.nb)
  upper.conf.znb <- c()
  lower.conf.znb <- c()
  upper.conf.prop.z <- c()
  lower.conf.prop.z <- c()
  upper.conf.prop.nb <- c()
  lower.conf.prop.nb <- c()
  sort.prop.z <- matrix(NA,nrow=n,ncol=nx.j)
  sort.prop.nb <- matrix(NA,nrow=n,ncol=nx.j)
  lower.index <- round(0.5*sign.level*n)
  upper.index <- round((1-0.5*sign.level)*n)
  for(l in 1:nx.j){
    sort.margEff.znb[,l] <- sort(margEff.znb[,l])
    sort.prop.z[,l] <- sort(prop.z[,l])
    sort.prop.nb[,l] <- sort(prop.nb[,l])
    lower.conf.znb[l] <- sort.margEff.znb[lower.index,l]
    lower.conf.prop.z[l] <- sort.prop.z[lower.index,l]
    lower.conf.prop.nb[l] <- sort.prop.nb[lower.index,l]
    upper.conf.znb[l] <- sort.margEff.znb[upper.index,l]
    upper.conf.prop.z[l] <- sort.prop.z[upper.index,l]
    upper.conf.prop.nb[l] <- sort.prop.nb[upper.index,l]
  }
  out <- list()
  out$x <- x.j.all
  out$mean.margEff.znb <- mean.margEff.znb
  out$lower.conf.znb <- lower.conf.znb
  out$upper.conf.znb <- upper.conf.znb
  out$mean.prop.z <- mean.prop.z
  out$mean.prop.nb <- mean.prop.nb
  out$lower.conf.prop.z <- lower.conf.prop.z
  out$upper.conf.prop.z <- upper.conf.prop.z
  out$lower.conf.prop.nb <- lower.conf.prop.nb
  out$upper.conf.prop.nb <- upper.conf.prop.nb
  return(out)
}


# sign.level <- 0.05
# j <- 9 #is ilustrative
# #j <- 11
# dx <- 0.2
# x <- OBS.X
# x.lim <- c(min(x[,j]),max(x[,j]))
# #x.lim <- c(0,75)
# #x.lim <- c(min(x[,j]),50) #for j=7
# gamma.zc <- Est.Obj$par.mat$Beta.multinom.ZC
# gamma.pl <- Est.Obj$par.mat$Beta.multinom.PL[1:5]
# gamma.pl <- c(gamma.pl,0,Est.Obj$par.mat$Beta.multinom.PL[6:12])
# beta.nb <- Est.Obj$par.mat$Beta.NB
# beta.pl <- Est.Obj$par.mat$Beta.PL
# beta.pl <- c(beta.pl,rep(0,12))
# c.pl <- Est.Obj$par.mat$C

marginal.effect.znbpl.fun <- function(sign.level,j,x,x.lim,dx,gamma.zc,gamma.pl,beta.nb,beta.pl,c.pl){
  n <- dim(x)[1]

  gamma.zc <- as.matrix(gamma.zc)
  gamma.pl <- as.matrix(gamma.pl)
  beta.nb <- as.matrix(beta.nb)
  beta.pl <- as.matrix(beta.pl)

  
  
  if(is.null(x)){
    x.extended <- matrix(1,nrow=n,ncol=1)
  }else{
    x.extended <- cbind(1,x)
  }
  
  j <- j + 1   #To compensate for the intercept
  
  gamma.zc.j <- gamma.zc[j]
  gamma.pl.j <- gamma.pl[j]
  beta.nb.j <- beta.nb[j]  
  beta.pl.j <- beta.pl[j]
  
  
  x.j.all <- seq(x.lim[1],x.lim[2],by=dx)
  nx.j <- length(x.j.all)
  
  margEff.znbpl <- matrix(NA,nrow=n,ncol=nx.j)
  prop.z <- matrix(NA,nrow=n,ncol=nx.j)
  prop.nb <- matrix(NA,nrow=n,ncol=nx.j)
  prop.p <- matrix(NA,nrow=n,ncol=nx.j)
  for(l in 1:nx.j){
    for(i in 1:n){
      x.il.extended <- as.matrix(x.extended[i,])
      x.il.extended[j] <- x.j.all[l]    
      gamma.plTx <- t(gamma.pl)%*%x.il.extended
      exp.gamma.zc.x <- exp(t(gamma.zc)%*%x.il.extended)
      exp.gamma.pl.x <- exp(t(gamma.pl)%*%x.il.extended)
      exp.beta.nb.x <- exp(t(beta.nb)%*%x.il.extended)
      exp.beta.pl.x <- exp(t(beta.pl)%*%x.il.extended)
      
      numerator <- exp.beta.nb.x + c.pl*exp.gamma.pl.x/exp.beta.pl.x + c.pl*exp.gamma.pl.x
      denominator <- 1 + exp.gamma.zc.x + exp.gamma.pl.x
       
      margEff.znbpl[i,l] <- numerator/denominator
      prop.nb[i,l] <- 1/(exp.beta.pl.x + exp.gamma.zc.x+1)
      prop.p[i,l] <- exp.beta.pl.x/(exp.beta.pl.x + exp.gamma.zc.x+1)
      prop.z[i,l] <- exp.gamma.zc.x/(exp.beta.pl.x + exp.gamma.zc.x+1)
    }  
  
  }
  mean.margEff.znbpl <- colMeans(margEff.znbpl)
  upper.conf.znbpl <- c()
  lower.conf.znbpl <- c()
  sort.margEff.znbpl <- matrix(NA,nrow=n,ncol=nx.j)
  mean.prop.z <- colMeans(prop.z)
  mean.prop.nb <- colMeans(prop.nb)
  mean.prop.p <- colMeans(prop.p)
  upper.conf.prop.z <- c()
  lower.conf.prop.z <- c()
  upper.conf.prop.nb <- c()
  lower.conf.prop.nb <- c()
  upper.conf.prop.p <- c()
  lower.conf.prop.p <- c()
  sort.prop.z <- matrix(NA,nrow=n,ncol=nx.j)
  sort.prop.nb <- matrix(NA,nrow=n,ncol=nx.j)
  sort.prop.p <- matrix(NA,nrow=n,ncol=nx.j)
  lower.index <- round(0.5*sign.level*n)
  upper.index <- round((1-0.5*sign.level)*n)
  for(l in 1:nx.j){
    sort.margEff.znbpl[,l] <- sort(margEff.znbpl[,l])
    lower.conf.znbpl[l] <- sort.margEff.znbpl[lower.index,l]
    upper.conf.znbpl[l] <- sort.margEff.znbpl[upper.index,l]
    sort.prop.z[,l] <- sort(prop.z[,l])
    sort.prop.nb[,l] <- sort(prop.nb[,l])
    sort.prop.p[,l] <- sort(prop.p[,l])
    lower.conf.prop.z[l] <- sort.prop.z[lower.index,l]
    lower.conf.prop.nb[l] <- sort.prop.nb[lower.index,l]
    lower.conf.prop.p[l] <- sort.prop.p[lower.index,l]
    upper.conf.prop.z[l] <- sort.prop.z[upper.index,l]
    upper.conf.prop.nb[l] <- sort.prop.nb[upper.index,l]
    upper.conf.prop.p[l] <- sort.prop.p[upper.index,l]
  }
  out <- list()
  out$x <- x.j.all
  out$mean.margEff.znbpl <- mean.margEff.znbpl
  out$lower.conf.znbpl <- lower.conf.znbpl
  out$upper.conf.znbpl <- upper.conf.znbpl
  out$mean.prop.z <- mean.prop.z
  out$mean.prop.nb <- mean.prop.nb
  out$mean.prop.p <- mean.prop.p
  out$lower.conf.prop.z <- lower.conf.prop.z
  out$upper.conf.prop.z <- upper.conf.prop.z
  out$lower.conf.prop.nb <- lower.conf.prop.nb
  out$upper.conf.prop.nb <- upper.conf.prop.nb
  out$lower.conf.prop.p <- lower.conf.prop.p
  out$upper.conf.prop.p <- upper.conf.prop.p
  return(out)
  #plot(x.j.all,mean.margEff.znbpl)
}




