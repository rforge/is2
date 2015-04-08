## IS2 algorithm functions

## define the is2 class
setClass(
         'is2',
         contains='pfilterd2.pomp',
         slots=c(
           transform = "logical",
           ivps = 'character',
           pars = 'character',
           Nis = 'integer',
           particles = 'function',
           var.factor='numeric',
           ic.lag='integer',
           cooling.type='character',
           cooling.fraction='numeric',
           method='character',
           random.walk.sd = 'numeric',
           conv.rec = 'matrix'
           )
         )


default.pomp.particles.fun <- function (Np, center, sd, ...) {
  matrix(
         data=rnorm(
           n=Np*length(center),
           mean=center,
           sd=sd
           ),
         nrow=length(center),
         ncol=Np,
         dimnames=list(
           names(center),
           NULL
           )
         )
}

cooling.function <- function (type, perobs, fraction, ntimes) {
  switch(
         type,
         geometric={
           factor <- fraction^(1/50)
           if (perobs) {
             function (nt, m) {
               alpha <- factor^(nt/ntimes+m-1)
               list(alpha=alpha,gamma=alpha^2)
             }
           } else {
             function (nt, m) {
               alpha <- factor^(m-1)
               list(alpha=alpha,gamma=alpha^2)
             }
           }
         },
         hyperbolic={
           if (fraction>=1)
             stop(
                  "iterated smoothing error: ",sQuote("cooling.fraction"),
                  " must be < 1 when cooling.type = ",
                  sQuote("hyperbolic"),
                  call.=FALSE
                  )
           if (perobs) {
             scal <- (50*ntimes*fraction-1)/(1-fraction)
             function (nt, m) {
               alpha <- (1+scal)/(scal+nt+ntimes*(m-1))
               list(alpha=alpha,gamma=alpha^2)
             }
           } else {
             scal <- (50*fraction-1)/(1-fraction)
             function (nt, m) {
               alpha <- (1+scal)/(scal+m-1)
               list(alpha=alpha,gamma=alpha^2)
             }

           }
         },
         stop("unrecognized cooling schedule type ",sQuote(type))
         )
}

is.cooling <- function (factor, n) {   # default geometric cooling schedule
  alpha <- factor^(n-1)
  list(alpha=alpha,gamma=alpha^2)
}

is2.cooling <- function (frac, nt, m, n) {   # cooling schedule for is2
  ## frac is the fraction of cooling after 50 iterations
  scal <- (50*n*frac-1)/(1-frac)
  alpha <- (1+scal)/(scal+nt+n*(m-1))
  list(alpha=alpha)
}

powerlaw.cooling <- function (init = 1, delta = 0.1, eps = (1-delta)/2, n) {
  m <- init
  if (n <= m) {                         # linear cooling regime
    alpha <- (m-n+1)/m
    gamma <- alpha
  } else {                              # power-law cooling regime
    alpha <- ((n/m)^(delta+eps))/n
    gamma <- ((n/m)^(delta+1))/n/n
  }
  list(alpha=alpha,gamma=gamma)
}

is2.internal <- function (object, Nis,
                          start, pars, ivps,
                          particles,
                          rw.sd,
                          Np, var.factor, ic.lag, lag,
                          cooling.type, cooling.fraction, cooling.factor,
                          method,
                          tol, max.fail,
                          verbose, transform, .ndone = 0L,
                          paramMatrix = NULL,
                          .getnativesymbolinfo = TRUE,
                          ...) {
  
  gnsi <- as.logical(.getnativesymbolinfo)

  transform <- as.logical(transform)
  
  if (length(start)==0)
    stop(
         "iterated smoothing error: ",sQuote("start")," must be specified if ",
         sQuote("coef(object)")," is NULL",
         call.=FALSE
         )
  
  if (transform)
    start <- partrans(object,start,dir="inverse")
  
  start.names <- names(start)
  if (is.null(start.names))
    stop("iterated smoothing error: ",sQuote("start")," must be a named vector",call.=FALSE)
  
  rw.names <- names(rw.sd)
  rwsdMat <-rw.sd
  if (is.null(rw.names) || any(rw.sd<0))
    stop("iterated smoothing error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
  if (!all(rw.names%in%start.names))
    stop("iterated smoothing error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
  rw.names <- names(rw.sd[rw.sd>0])
  if (length(rw.names) == 0)
    stop("iterated smoothing error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)
  
  if (
      !is.character(pars) ||
      !is.character(ivps) ||
      !all(pars%in%start.names) ||
      !all(ivps%in%start.names) ||
      any(pars%in%ivps) ||
      any(ivps%in%pars) ||
      !all(pars%in%rw.names) ||
      !all(ivps%in%rw.names)
      )
    stop(
         "iterated smoothing error: ",
         sQuote("pars")," and ",sQuote("ivps"),
         " must be mutually disjoint subsets of ",
         sQuote("names(start)"),
         " and must have a positive random-walk SDs specified in ",
         sQuote("rw.sd"),
         call.=FALSE
         )
  
  if (!all(rw.names%in%c(pars,ivps))) {
    extra.rws <- rw.names[!(rw.names%in%c(pars,ivps))]
    warning(
            ngettext(length(extra.rws),"iterated smoothing warning: the variable ",
                     "iterated smoothing warning: the variables "),
            paste(sQuote(extra.rws),collapse=", "),
            ngettext(length(extra.rws)," has positive random-walk SD specified, but is included in neither ",
                     " have positive random-walk SDs specified, but are included in neither "),
            sQuote("pars")," nor ",sQuote("ivps"),
            ngettext(length(extra.rws),". This random walk SD will be ignored.",
                     ". These random walk SDs will be ignored."),
            call.=FALSE
            )
  }
  rw.sd <- rw.sd[c(pars,ivps)]
  rw.names <- names(rw.sd)
  
  ntimes <- length(time(object))
  if (is.null(Np)) stop("iterated smoothing error: ",sQuote("Np")," must be specified",call.=FALSE)
  if (is.function(Np)) {
    Np <- try(
              vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
              silent=FALSE
              )
    if (inherits(Np,"try-error"))
      stop("if ",sQuote("Np")," is a function, it must return a single positive integer")
  }
  if (length(Np)==1)
    Np <- rep(Np,times=ntimes+1)
  if (any(Np<=0))
    stop("number of particles, ",sQuote("Np"),", must always be positive")
  if (!is.numeric(Np))
    stop(sQuote("Np")," must be a number, a vector of numbers, or a function")
  Np <- as.integer(Np)
  
  ic.lag <- as.integer(ic.lag)
  if ((length(ic.lag)!=1)||(ic.lag<1))
    stop("iterated smoothing error: ",sQuote("ic.lag")," must be a positive integer",call.=FALSE)
  if (ic.lag>ntimes) {
    warning(
            "iterated smoothing warning: ",sQuote("ic.lag")," = ",ic.lag," > ",ntimes,
            " = length(time(",sQuote("object"),"))",
            " is nonsensical.  Setting ",sQuote("ic.lag")," = ",ntimes,".",
            call.=FALSE
            )
    ic.lag <- length(time(object))
  }
  if ((length(pars)==0)&&(ic.lag<length(time(object)))) {
    warning(
            "iterated smoothing warning: only IVPs are to be estimated, yet ",sQuote("ic.lag")," = ",ic.lag,
            " < ",ntimes," = length(time(",sQuote("object"),")),",
            " so unnecessary work is to be done.",
            call.=FALSE
            )
  }
  
  ## the following deals with the deprecated option 'cooling.factor'
  if (!missing(cooling.factor)) {
    warning(sQuote("cooling.factor")," is deprecated.\n",
            "See ",sQuote("?is2")," for instructions on specifying the cooling schedule.",
            call.=FALSE)
    cooling.factor <- as.numeric(cooling.factor)
    if ((length(cooling.factor)!=1)||(cooling.factor<0)||(cooling.factor>1))
      stop("iterated smoothing error: ",sQuote("cooling.factor")," must be a number between 0 and 1",call.=FALSE)
    if (missing(cooling.fraction)) {
      cooling.fraction <- cooling.factor^50
    } else {
      warning("specification of ",sQuote("cooling.factor"),
              " is overridden by that of ",sQuote("cooling.fraction"),
              call.=FALSE)
    }
  }

  if (missing(cooling.fraction))
    stop("iterated smoothing error: ",sQuote("cooling.fraction")," must be specified",call.=FALSE)
  cooling.fraction <- as.numeric(cooling.fraction)
  if ((length(cooling.fraction)!=1)||(cooling.fraction<0)||(cooling.fraction>1))
    stop("iterated smoothing error: ",sQuote("cooling.fraction")," must be a number between 0 and 1",call.=FALSE)
  
  cooling <- cooling.function(
                              type=cooling.type,
                              perobs=(method=="is2")||(method=="ris1"),
                              fraction=cooling.fraction,
                              ntimes=ntimes
                              )

  if ((length(var.factor)!=1)||(var.factor < 0))
    stop("iterated smoothing error: ",sQuote("var.factor")," must be a positive number",call.=FALSE)
  
  Nis <- as.integer(Nis)
  if (Nis<0)
    stop("iterated smoothing error: ",sQuote("Nis")," must be a positive integer",call.=FALSE)

  theta <- start
  
  sigma <- rep(0,length(start))
  names(sigma) <- start.names
  
  rw.sd <- rw.sd[c(pars,ivps)]
  rw.names <- names(rw.sd)
  
  sigma[rw.names] <- rw.sd
  
  conv.rec <- matrix(
                     data=NA,
                     nrow=Nis+1,
                     ncol=length(theta)+2,
                     dimnames=list(
                       seq(.ndone,.ndone+Nis),
                       c('loglik','nfail',names(theta))
                       )
                     )
  conv.rec[1L,] <- c(NA,NA,theta)
  
  if (!all(is.finite(theta[c(pars,ivps)]))) {
    stop(
         sQuote("is2")," cannot estimate non-finite parameters.\n",
         "The following ",if (transform) "transformed ", "parameters are non-finite: ",
         paste(
               sQuote(c(pars,ivps)[!is.finite(theta[c(pars,ivps)])]),
               collapse=","
               ),
         call.=FALSE
         )
  }
  
  obj <- as(object,"pomp")
  
  if (Nis>0) {
    tmp.is2 <- new("is2",object,particles=particles,Np=Np[1L])
  } else {
    pfp <- obj
  }
  
  have.parmat <- !(is.null(paramMatrix) || length(paramMatrix)==0)
  newtheta<-theta
  
  for (n in seq_len(Nis)) { ## iterate the smoothing

    ## get the intensity of artificial noise from the cooling schedule
    cool.sched <- cooling(nt=1,m=.ndone+n)
    sigma.n <- sigma*cool.sched$alpha
    names(sigma.n)<-names(theta)
    
    ## initialize the parameter portions of the particles
    P <- try(
             particles(
                       tmp.is2,
                       Np=Np[1L],
                       center=theta,
                       sd=sigma.n*var.factor
                       ),
             silent = FALSE
             )
    if (inherits(P,"try-error")) 
      stop("iterated smoothing error: error in ",sQuote("particles"),call.=FALSE)

    pfp <- try(
               pfilter2.internal(
                                object=obj,
                                params=P, 
                                Np=Np,
                                tol=tol,
                                max.fail=max.fail,
                                pred.mean=(n==Nis),
                                pred.var=TRUE,
                                filter.mean=TRUE,
                                cooling=cooling,
                                cooling.m=.ndone+n,
                                .corr=(method=="is1"),
                                .wn =(method=="ris1")||(method=="is1"),
                                .rw.sd=sigma.n[pars],
                                .transform=transform,
                                save.states=FALSE, 
                                save.params=FALSE,
                                lag=lag,
                                verbose=verbose,
                                .getnativesymbolinfo=gnsi
                                ),
               silent=FALSE
               )
    if (inherits(pfp,"try-error")) 
      stop("iterated smoothing error: error in ",sQuote("pfilter2"),call.=FALSE)

    gnsi <- FALSE
    
    ## update parameters
    switch(
           method,
           is2={
             paramMatrix <- pfp@paramMatrix
               oldtheta1<-theta
               oldtheta<-theta
               npars<-length(pars)
               names(oldtheta)<-names(theta)
               newtheta[c(pars,ivps)]<-0
               npars<-length(theta)
               Hessian<-array(0,dim=c(npars,npars))
               colnames(Hessian)<-names(theta)
               rownames(Hessian)<-names(theta)
               if (lag>0){
                 
                 phat<-pfp@phats
                 names(phat)<-names(theta)
                 covhat<-pfp@covhats
                 colnames(covhat)<-names(theta)
                 rownames(covhat)<-names(theta)
                 Hessian[c(pars,ivps),c(pars,ivps)]<-covhat[c(pars,ivps),c(pars,ivps)]-ntimes*diag(sigma.n[c(pars,ivps)]^2)
                 newtheta[c(pars,ivps)] <- phat[c(pars,ivps)]-ntimes*oldtheta[c(pars,ivps)]  
                 Hessian[c(pars,ivps),c(pars,ivps)]<-0.5*(Hessian[c(pars,ivps),c(pars,ivps)]+t(Hessian[c(pars,ivps),c(pars,ivps)]))
                 v1 <- cool.sched$gamma*rwsdMat[c(pars,ivps)]^2  
                 newtheta[c(pars,ivps)]<- solve(Hessian[c(pars,ivps),c(pars,ivps)])%*%newtheta[c(pars,ivps)]*v1  
                 theta[c(pars,ivps)]  <-  oldtheta1[c(pars,ivps)]-newtheta[c(pars,ivps)]
               }
           },
           ris1={
             paramMatrix <- pfp@paramMatrix
             oldtheta1<-theta
             oldtheta<-theta
             npars<-length(pars)
             names(oldtheta)<-names(theta)
             newtheta[c(pars,ivps)]<-0
             npars<-length(theta)
             Hessian<-array(0,dim=c(npars,npars))
             colnames(Hessian)<-names(theta)
             rownames(Hessian)<-names(theta)
             if (lag>0){
               
               phat<-pfp@phats
               names(phat)<-names(theta)
               covhat<-pfp@covhats
               colnames(covhat)<-names(theta)
               rownames(covhat)<-names(theta)
               Hessian[c(pars,ivps),c(pars,ivps)]<-covhat[c(pars,ivps),c(pars,ivps)]-ntimes*diag(sigma.n[c(pars,ivps)]^2)
               newtheta[c(pars,ivps)] <- phat[c(pars,ivps)]-ntimes*oldtheta[c(pars,ivps)]  
               Hessian[c(pars,ivps),c(pars,ivps)]<-0.5*(Hessian[c(pars,ivps),c(pars,ivps)]+t(Hessian[c(pars,ivps),c(pars,ivps)]))
               v1 <- cool.sched$gamma*rwsdMat[c(pars,ivps)]^2  
               newtheta[c(pars,ivps)]<- solve(Hessian[c(pars,ivps),c(pars,ivps)])%*%newtheta[c(pars,ivps)]*v1  
               theta[c(pars,ivps)]  <-  oldtheta1[c(pars,ivps)]-newtheta[c(pars,ivps)]
             }
           },
           
           is1={
             paramMatrix <- pfp@paramMatrix
             oldtheta1<-theta
             oldtheta<-theta
             npars<-length(pars)
             names(oldtheta)<-names(theta)
             newtheta[c(pars,ivps)]<-0
             npars<-length(theta)
             Hessian<-array(0,dim=c(npars,npars))
             colnames(Hessian)<-names(theta)
             rownames(Hessian)<-names(theta)
             if (lag>0){
               
               phat<-pfp@phats
               names(phat)<-names(theta)
               covhat<-pfp@covhats
               colnames(covhat)<-names(theta)
               rownames(covhat)<-names(theta)
               Hessian[c(pars,ivps),c(pars,ivps)]<-covhat[c(pars,ivps),c(pars,ivps)]-ntimes*diag(sigma.n[c(pars,ivps)]^2)
               newtheta[c(pars,ivps)] <- phat[c(pars,ivps)]-ntimes*oldtheta[c(pars,ivps)]  
               for(i in 1:ntimes){
                 for(j in 1:lag){
                   covhat<-pfp@pcovhats[,,j,i]
                   colnames(covhat)<-names(theta)
                   rownames(covhat)<-names(theta)
                   Hessian[c(pars,ivps),c(pars,ivps)]<-Hessian[c(pars,ivps),c(pars,ivps)]+2*covhat[c(pars,ivps),c(pars,ivps)]
                   
                 }
               }
               Hessian[c(pars,ivps),c(pars,ivps)]<-0.5*(Hessian[c(pars,ivps),c(pars,ivps)]+t(Hessian[c(pars,ivps),c(pars,ivps)]))
               v1 <- cool.sched$gamma*rwsdMat[c(pars,ivps)]^2  
               newtheta[c(pars,ivps)]<- solve(Hessian[c(pars,ivps),c(pars,ivps)])%*%newtheta[c(pars,ivps)]*v1  
               theta[c(pars,ivps)]  <-  oldtheta1[c(pars,ivps)]-newtheta[c(pars,ivps)]
             }
           },
           
           
           
           stop("unrecognized method ",sQuote(method))
           )
    theta[ivps] <- pfp@filter.mean[ivps,ic.lag]
    conv.rec[n+1,-c(1,2)] <- theta
    conv.rec[n,c(1,2)] <- c(pfp@loglik,pfp@nfail)
    
    if (verbose) cat("IS2 iteration ",n," of ",Nis," completed\n")
    
  } ### end of main loop

  ## back transform the parameter estimate if necessary
  if (transform) theta <- partrans(pfp,theta,dir="forward")
  
  new(
      "is2",
      pfp,
      transform=transform,
      params=theta,
      ivps=ivps,
      pars=pars,
      Nis=Nis,
      particles=particles,
      var.factor=var.factor,
      ic.lag=ic.lag,
      random.walk.sd=sigma[rw.names],
      tol=tol,
      conv.rec=conv.rec,
      method=method,
      cooling.type=cooling.type,
      cooling.fraction=cooling.fraction,
      paramMatrix= array(data=numeric(0),dim=c(0,0)),
      lag=lag
      )
}

setMethod(
          "is2",
          signature=signature(object="pomp"),
          function (object, Nis = 1,
                    start,
                    pars, ivps = character(0),
                    particles, rw.sd,
                    Np, ic.lag, var.factor,lag,
                    cooling.type = c("geometric","hyperbolic"),
                    cooling.fraction, cooling.factor,
                    method = c("is2","ris1","is1"),
                    tol = 1e-17, max.fail = Inf,
                    verbose = getOption("verbose"),
                    transform = FALSE,
                    ...) {
            
            transform <- as.logical(transform)
            method <- match.arg(method)
            
            if (missing(start)) start <- coef(object)
            if (missing(rw.sd))
              stop("iterated smoothing error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
            if (missing(ic.lag)) {
              if (length(ivps)>0) {
                stop("iterated smoothing error: ",sQuote("ic.lag"),
                     " must be specified if ",sQuote("ivps"),
                     " are",call.=FALSE)
              } else {
                ic.lag <- length(time(object))
              }
            }
            if (missing(lag)){
              lag=0
            }
            if (missing(pars)) {
              rw.names <- names(rw.sd)[rw.sd>0]
              pars <- rw.names[!(rw.names%in%ivps)]
            }
            if (missing(Np))
              stop("iterated smoothing error: ",sQuote("Np")," must be specified",call.=FALSE)
            if (missing(var.factor))
              stop("iterated smoothing error: ",sQuote("var.factor")," must be specified",call.=FALSE)

            cooling.type <- match.arg(cooling.type)
            
            if (missing(particles)) { # use default: normal distribution
              particles <- default.pomp.particles.fun
            } else {
              particles <- match.fun(particles)
              if (!all(c('Np','center','sd','...')%in%names(formals(particles))))
                stop(
                     "iterated smoothing error: ",
                     sQuote("particles"),
                     " must be a function of prototype ",
                     sQuote("particles(Np,center,sd,...)"),
                     call.=FALSE
                     )
            }
            
            is2.internal(
                         object=object,
                         Nis=Nis,
                         start=start,
                         pars=pars,
                         ivps=ivps,
                         particles=particles,
                         rw.sd=rw.sd,
                         Np=Np,
                         cooling.type=cooling.type,
                         cooling.factor=cooling.factor,
                         cooling.fraction=cooling.fraction,
                         var.factor=var.factor,
                         ic.lag=ic.lag,
                         lag=lag,
                         method=method,
                         tol=tol,
                         max.fail=max.fail,
                         verbose=verbose,
                         transform=transform,
                         ...
                         )
            
          }
          )


setMethod(
          "is2",
          signature=signature(object="pfilterd2.pomp"),
          function (object, Nis = 1, Np, tol,
                    ...) {
            
            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol
            
            is2(
                object=as(object,"pomp"),
                Nis=Nis,
                Np=Np,
                tol=tol,
                ...
                )
          }
          )

setMethod(
          "is2",
          signature=signature(object="is2"),
          function (object, Nis,
                    start,
                    pars, ivps,
                    particles, rw.sd,
                    Np, ic.lag, lag, var.factor,
                    cooling.type, cooling.fraction,
                    method,
                    tol,
                    transform,
                    ...) {
            
            if (missing(Nis)) Nis <- object@Nis
            if (missing(start)) start <- coef(object)
            if (missing(pars)) pars <- object@pars
            if (missing(ivps)) ivps <- object@ivps
            if (missing(particles)) particles <- object@particles
            if (missing(rw.sd)) rw.sd <- object@random.walk.sd
            if (missing(ic.lag)) ic.lag <- object@ic.lag
            if (missing(var.factor)) var.factor <- object@var.factor
            if (missing(cooling.type)) cooling.type <- object@cooling.type
            if (missing(cooling.fraction)) cooling.fraction <- object@cooling.fraction
            if (missing(method)) method <- object@method
            if (missing(transform)) transform <- object@transform
            transform <- as.logical(transform)

            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol
            if (missing(lag)) lag <- object@lag
            is2(
                object=as(object,"pomp"),
                Nis=Nis,
                start=start,
                pars=pars,
                ivps=ivps,
                particles=particles,
                rw.sd=rw.sd,
                Np=Np,
                cooling.type=cooling.type,
                cooling.fraction=cooling.fraction,
                var.factor=var.factor,
                ic.lag=ic.lag,
                lag=lag,
                method=method,
                tol=tol,
                transform=transform,
                ...
                )
          }
          )

setMethod(
          'continue',
          signature=signature(object='is2'),
          function (object, Nis = 1,
                    ...) {
            
            ndone <- object@Nis
            
            obj <- is2(
                       object=object,
                       Nis=Nis,
                       .ndone=ndone,
                       paramMatrix=object@paramMatrix,
                       ...
                       )
            
            object@conv.rec[ndone+1,c('loglik','nfail')] <- obj@conv.rec[1L,c('loglik','nfail')]
            obj@conv.rec <- rbind(
                                  object@conv.rec,
                                  obj@conv.rec[-1L,colnames(object@conv.rec)]
                                  )
            obj@Nis <- as.integer(ndone+Nis)
            
            obj
          }
          )

is2.profileDesign <- function (object, profile, lower, upper, nprof, ivps, 
                               rw.sd, Np, ic.lag,lag, var.factor, cooling.factor,option, cooling.fraction, paramMatrix, ...)
{
  if (missing(profile)) profile <- list()
  if (missing(lower)) lower <- numeric(0)
  if (missing(upper)) upper <- lower
  if (length(lower)!=length(upper))
    stop(sQuote("lower")," and ",sQuote("upper")," must be of the same length")
  pars <- names(lower)
  if (missing(ivps)) ivps <- character(0)
  Np <- as.integer(Np)
  
  pd <- do.call(profileDesign,c(profile,list(lower=lower,upper=upper,nprof=nprof)))
  
  object <- as(object,"pomp")
  
  pp <- coef(object)
  idx <- !(names(pp)%in%names(pd))
  if (any(idx)) pd <- cbind(pd,as.list(pp[idx]))
  
  ans <- vector(mode="list",length=nrow(pd))
  for (k in seq_len(nrow(pd))) {
    ans[[k]] <- list(
                     mf=is2(
                       object,
                       Nis=0,
                       start=unlist(pd[k,]),
                       pars=pars,
                       ivps=ivps,
                       rw.sd=rw.sd,
                       Np=Np,
                       ic.lag=ic.lag,
                       lag=lag,
                       var.factor=var.factor,
                       cooling.factor=cooling.factor,
                       option=option,
                       cooling.fraction=cooling.fraction,
                       paramMatrix=paramMatrix,
                       ...
                       )
                     )
  }
  
  ans
}
