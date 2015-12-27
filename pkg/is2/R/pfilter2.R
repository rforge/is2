## particle fixed lag smoothing for is2 codes

setClass(
    "pfilterd2.pomp",
    contains="pomp",
    slots=c(
        pred.mean="array",
        pred.var="array",
        filter.mean="array",
        paramMatrix="array",
        eff.sample.size="numeric",
        cond.loglik="numeric",
        saved.states="list",
        saved.params="list",
        seed="integer",
        Np="integer",
        tol="numeric",
        nfail="integer",
        loglik="numeric",
        phats="numeric",
        covhats="array",
        pcovhats="array",
        lag = "numeric"
    ),
    prototype=prototype(
        pred.mean=array(data=numeric(0),dim=c(0,0)),
        pred.var=array(data=numeric(0),dim=c(0,0)),
        filter.mean=array(data=numeric(0),dim=c(0,0)),
        paramMatrix=array(data=numeric(0),dim=c(0,0)),
        eff.sample.size=numeric(0),
        cond.loglik=numeric(0),
        saved.states=list(),
        saved.params=list(),
        seed=as.integer(NA),
        Np=as.integer(NA),
        tol=as.double(NA),
        nfail=as.integer(NA),
        loglik=as.double(NA),
        phats=numeric(0),
        covhats=array(data=numeric(0),dim=c(0,0)),
        pcovhats=array(data=numeric(0),dim=c(0,0,0,0)),
        lag = as.integer(NA)
    )
)

ancestor<-function(plist, t, lag, index){
    for (i in 0:(lag-1)){
        index=plist[[t-i]][index]
    }
    return(index)
}

smoothing<-function(aparticles, xparticles, pparticles, wparticles, nt, ntimes, lag, nvars, npars, rw, Np){
    at=rep(0,Np)        #current parent index
    bt=rep(0,Np)
    nlength<-length(rw)
    phat<-rep(0,npars)  #smoothed par
    pcovhat<-matrix(0,npars,lag) #covariance
    lcovhat<-array(0,dim=c(npars,npars,lag)) 
  
    if(lag>0){
        kk=nt+lag
        if(nt==ntimes){
        }
        else{
            if(kk<=ntimes){
                at=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=lag ))
                for ( jj in 1: npars ){
                    phat[jj]=sum(pparticles[[nt]][jj,at]*wparticles[[kk]])
                }
                for ( jj in 1: npars ){
                    for(nn in 1:(kk-nt)){
                        bt=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=lag+1-nn )) 
                        pcovhat[jj,nn] =sum(pparticles[[nt+nn]][jj,bt]*wparticles[[kk]])
                    }
                }
                for ( jj in 1: npars ){
                    for (ll in 1: npars){
                        for(nn in 1:lag){
                            bt=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=lag+1-nn ))
                            lcovhat[jj,ll,nn] = sum((pparticles[[nt+nn]][ll,bt]-pcovhat[ll,nn])*(pparticles[[nt]][jj,at]-phat[jj])*wparticles[[kk]])/(Np-1)
                        }
                    }
                }
            }  
            else{
                at=unlist(lapply(aparticles[[ntimes]],ancestor,  plist=aparticles, t=ntimes,lag=(ntimes-nt) ))
                kk=ntimes
                for ( jj in 1: npars ){
                    phat[jj]=sum(pparticles[[nt]][jj,at]*wparticles[[ntimes]])
                }
                for ( jj in 1: npars ){
                    for(nn in 1:(kk-nt)){
                        bt=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=nn )) 
                        pcovhat[jj,nn] =sum(pparticles[[nt]][jj,bt]*wparticles[[kk]])
                    }
                }
                for ( jj in 1: npars ){
                    for (ll in 1: npars){
                        for(nn in 1:(kk-nt)){
                            bt=unlist(lapply(aparticles[[kk]],ancestor,  plist=aparticles, t=kk,lag=nn ))
                            lcovhat[jj,ll,nn] = sum((pparticles[[nt]][ll,bt]-pcovhat[ll,nn])*(pparticles[[ntimes]][jj,at]-phat[jj])*wparticles[[kk]])/(Np-1)
                        }
                    }
                }
            }
        }
    }
    return(lcovhat)
}

pfilter2.internal <- function (object, params, Np, tol, max.fail,
                                            pred.mean, pred.var, filter.mean,
                                            cooling, cooling.m, .mif2 = FALSE, .wn=FALSE,.corr=FALSE,
                                            .rw.sd, seed, verbose,
                                            save.states, save.params,lag,
                                            .transform, .getnativesymbolinfo = TRUE){
    
    ptsi.inv <- ptsi.for <- gnsi.rproc <- gnsi.dmeas <- as.logical(.getnativesymbolinfo)
    mif2 <- as.logical(.mif2)
    corr <- as.logical(.corr)
    wn <- as.logical(.wn)
    Sumsigma2 <- 100	#asymptotic sum of sigma square
    transform <- as.logical(.transform)
    
    if (missing(seed)) seed <- NULL
    if (missing(lag)) lag <- 0
    if (!is.null(seed)){
        if (!exists(".Random.seed",where=.GlobalEnv)){ # need to initialize the RNG
            runif(1)
        }
        save.seed <- get(".Random.seed",pos=.GlobalEnv)
        set.seed(seed)
    }
  
    if (length(params)==0)
        stop(sQuote("pfilter2")," error: ",sQuote("params")," must be specified",call.=FALSE)
  
    if (missing(tol))
        stop(sQuote("pfilter2")," error: ",sQuote("tol")," must be specified",call.=FALSE)
  
    one.par <- FALSE
    times <- time(object,t0=TRUE)
    ntimes <- length(times)-1
  
    if (missing(Np))
        Np <- NCOL(params)
    if (is.function(Np)){
        Np <- try(
            vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
            silent=FALSE
        )
        if (inherits(Np,"try-error"))
            stop("if ",sQuote("Np")," is a function, it must return a single positive integer",call.=FALSE)
    }
    if (length(Np)==1)
        Np <- rep(Np,times=ntimes+1)
    if (any(Np<=0))
        stop("number of particles, ",sQuote("Np"),", must always be positive",call.=FALSE)
    if (!is.numeric(Np))
        stop(sQuote("Np")," must be a number, a vector of numbers, or a function",call.=FALSE)
    Np <- as.integer(Np)
  
    if (is.null(dim(params))){
        one.par <- TRUE                  # there is only one parameter vector
        coef(object) <- params           # set params slot to the parameters
        params <- matrix(
            params,
            nrow=length(params),
            ncol=Np[1L],
            dimnames=list(
                names(params),
                NULL
            )
        )
    }
    paramnames <- rownames(params)
    if (is.null(paramnames))
        stop(sQuote("pfilter2")," error: ",sQuote("params")," must have rownames",call.=FALSE)
  
    x <- init.state(
        object,
        params=if (transform){
            partrans(object,params,dir="forward", .getnativesymbolinfo=ptsi.for)
        }
        else{
            params
        }
    )
    statenames <- rownames(x)
    nvars <- nrow(x)
    ptsi.for <- FALSE
  
    ## set up storage for saving samples from filtering distributions
    if (save.states)
        xparticles <- vector(mode="list",length=ntimes)
    else
        xparticles <- list()
    if (save.params)
        pparticles <- vector(mode="list",length=ntimes)
    else
        pparticles <- list()
    random.walk <- !missing(.rw.sd)
    if (random.walk){
        rw.names <- names(.rw.sd)
        if (is.null(rw.names)||!is.numeric(.rw.sd))
            stop(sQuote("pfilter2")," error: ",sQuote(".rw.sd")," must be a named vector",call.=FALSE)
        if (any(!(rw.names%in%paramnames)))
            stop(
                sQuote("pfilter2")," error: the rownames of ",
                sQuote("params")," must include all of the names of ",
                sQuote(".rw.sd"),"",call.=FALSE
           )
        sigma <- .rw.sd
    } 
    else{
        rw.names <- character(0)
        sigma <- NULL
    }
  
    loglik <- rep(NA,ntimes)
    eff.sample.size <- numeric(ntimes)
    nfail <- 0
    npars <- length(rw.names)
  
    ## set up storage for prediction means, variances, etc.
    if (pred.mean)
        pred.m <- matrix(
            data=0,
            nrow=nvars+npars,
            ncol=ntimes,
            dimnames=list(c(statenames,rw.names),NULL)
        )
    else
        pred.m <- array(data=numeric(0),dim=c(0,0))
  
    if (pred.var)
        pred.v <- matrix(
            data=0,
            nrow=nvars+npars,
            ncol=ntimes,
            dimnames=list(c(statenames,rw.names),NULL)
        )
    else
        pred.v <- array(data=numeric(0),dim=c(0,0))
  
    if (filter.mean)
        if (random.walk)
            filt.m <- matrix(
                data=0,
                nrow=nvars+length(paramnames),
                ncol=ntimes,
                dimnames=list(c(statenames,paramnames),NULL)
            )
        else
            filt.m <- matrix(
                data=0,
                nrow=nvars,
                ncol=ntimes,
                dimnames=list(statenames,NULL)
            )
    else
        filt.m <- array(data=numeric(0),dim=c(0,0))
  
    ##########################################
    # Fixed-lag Smoothing 
    ##########################################
    if ((lag<0)||(lag>ntimes))
        stop("Lag, ",sQuote("lag"),", must greater than 0 and less than ntimes",call.=FALSE)
    
    npars<- length(paramnames)
    phats<-rep(0,npars)
    names(phats)<-paramnames
  
    covhats <- array(
        0,
        dim=c(npars,npars)
    )
    pcovhats <- array(
        0,
        dim=c(0,0,0, 0)
    )
    if (lag>0 && !corr){
        asparticles <- vector(mode="list",length=(lag+1))
        xsparticles <- vector(mode="list",length=lag)
        psparticles <- vector(mode="list",length=lag)
    }
    
    if (lag>0 && corr){
        asparticles <- vector(mode="list",length=2*lag+1)
		xsparticles <- vector(mode="list",length=2*lag)        
		wsparticles <- vector(mode="list",length=2*lag)
		psparticles <- vector(mode="list",length=2*lag)
		Cindex <- matrix(0, nrow=lag, ncol = Np[1])
    }
  
    ##########################################
  

    for (nt in seq_len(ntimes)) {
        if (mif2) {	  
          cool.sched <- cooling(nt=nt,m=cooling.m)
          sigma1 <- sigma*cool.sched$alpha
        } 
        else {
          sigma1 <- sigma
        }
        
        ## transform the parameters if necessary
        if (transform) tparams <- partrans(object,params,dir="forward", .getnativesymbolinfo=ptsi.for)
            ptsi.for <- FALSE
    
        ## advance the state variables according to the process model
        X <- try(
            rprocess(
                object,
                xstart=x,
                times=times[c(nt,nt+1)],
                params=if (transform) tparams else params,
                offset=1,
                .getnativesymbolinfo=gnsi.rproc
            ),
            silent=FALSE
        )
        if (inherits(X,'try-error'))
            stop(sQuote("pfilter2")," error: process simulation error",call.=FALSE)
        gnsi.rproc <- FALSE
    
        if (pred.var){ ## check for nonfinite state variables and parameters
            problem.indices <- unique(which(!is.finite(X),arr.ind=TRUE)[,1L])
            if (length(problem.indices)>0){  # state variables
                stop(
                    sQuote("pfilter2")," error: non-finite state variable(s): ",
                    paste(rownames(X)[problem.indices],collapse=', '),
                    call.=FALSE
                )
            }
            if (random.walk){ # parameters (need to be checked only if 'random.walk=TRUE')
                problem.indices <- unique(which(!is.finite(params[rw.names,,drop=FALSE]),arr.ind=TRUE)[,1L])
                if (length(problem.indices)>0){
                    stop(
                        sQuote("pfilter2")," error: non-finite parameter(s): ",
                        paste(rw.names[problem.indices],collapse=', '),
                        call.=FALSE
                    )
                }
            }
        }
    
        ## determine the weights
        weights <- try(
            dmeasure(
                object,
                y=object@data[,nt,drop=FALSE],
                x=X,
                times=times[nt+1],
                params=if (transform) tparams else params,
                log=FALSE,
                .getnativesymbolinfo=gnsi.dmeas
            ),
            silent=FALSE
        )
        if (inherits(weights,'try-error'))
            stop(sQuote("pfilter2")," error: error in calculation of weights",call.=FALSE)
        if (any(!is.finite(weights))){
            stop(sQuote("pfilter2")," error: ",sQuote("dmeasure")," returns non-finite value",call.=FALSE)
        }
        gnsi.dmeas <- FALSE
    
        ## compute prediction mean, prediction variance, filtering mean,
        ## effective sample size, log-likelihood
        ## also do resampling if filtering has not failed
        xx <- try(
            .Call(
                pfilter2_computations,
                X,params,Np[nt+1],
                random.walk,
                sigma1,
                pred.mean,pred.var,
                filter.mean,one.par,
                weights,tol
            ),
            silent=FALSE
        )
        if (inherits(xx,'try-error')){
            stop(sQuote("pfilter2")," error",call.=FALSE)
        }
        all.fail <- xx$fail
        loglik[nt] <- xx$loglik
        eff.sample.size[nt] <- xx$ess
    
        x <- xx$states
        #random walk change to white noise for lag>0
        if(lag>0 && wn){
            params <- params
        }
        else{
            params <- xx$params
        }
    
        if (pred.mean)
            pred.m[,nt] <- xx$pm
        if (pred.var)
            pred.v[,nt] <- xx$pv
        if (filter.mean)
            filt.m[,nt] <- xx$fm
    
        if (all.fail){ ## all particles are lost
            nfail <- nfail+1
            if (verbose)
                message("filtering failure at time t = ",times[nt+1])
            if (nfail>max.fail)
                stop(sQuote("pfilter2")," error: too many filtering failures",call.=FALSE)
        }   
    
        if (save.states){
            xparticles[[nt]] <- x
        }
        
        if (save.params){
            pparticles[[nt]] <- params
        }
    
        if (verbose && (nt%%5==0))
            cat("pfilter2 timestep",nt,"of",ntimes,"finished\n")
    
        ##########################################
        if (lag>0 && !corr){
            if(nt<(lag+1)){
                xsparticles[[nt]] <- x
                psparticles[[nt]] <- params
                asparticles[[nt+1]] <- xx$pa+1  #offset 1 from C
                
            }
            if(nt>lag && nt<=ntimes){
		index<-unlist(ancestor(plist=asparticles,t=lag+1, lag=lag, 1:(Np[1])))
                psparticles[[1]][!is.finite(psparticles[[1]])] <- 0 
                
                C<-cov.wt(t(psparticles[[1]][,index]),wt=xx$weight)
                phats<-phats+C$center
                covhats<-covhats+C$cov/ntimes
				if (lag>1){                
					for (i in 1:(lag-1)){
                		psparticles[[i]]<-psparticles[[i+1]]
						asparticles[[i+1]]<-asparticles[[i+2]]  
              		}
				}
            	psparticles[[lag]]<-params
				asparticles[[lag+1]]<-xx$pa+1
            }
            if(nt==ntimes){
				index<-unlist(ancestor(plist=asparticles,t=lag+1, lag=lag, 1:(Np[1])))
                psparticles[[1]][!is.finite(psparticles[[1]])] <- 0 
                C<-cov.wt(t(psparticles[[1]][,index]),wt=xx$weight)
                phats<-phats+C$center
                covhats<-covhats+C$cov/ntimes
				if (lag>1){                
					for (i in 1:(lag-1)){
                    	psparticles[[i]]<-psparticles[[i+1]]
						asparticles[[i+1]]<-asparticles[[i+2]]  
                	}
				}
                
            }
        }
        if(lag>0 && corr){
	    if(nt<(2*lag+1)){
                psparticles[[nt]] <- params
		wsparticles[[nt]] <- xx$weight
                asparticles[[nt+1]] <- xx$pa+1  #offset 1 from C
                
            }
            if(nt>(2*lag) && nt<=ntimes){  
				index<-unlist(ancestor(plist=asparticles,t=lag+1, lag=lag, 1:(Np[1])))
				psparticles[[1]][!is.finite(psparticles[[1]])] <- 0
				C<-cov.wt(t(psparticles[[1]][,index]),wt=wsparticles[[lag+1]])
				phats<-phats+C$center
				covhats<-covhats+C$cov/Np[1]
				if(lag>1){
					for(i in 2:lag){
						index1<-unlist(ancestor(plist=asparticles,t=lag+1, lag=lag+1-i, 1:(Np[1])))
						psparticles[[i]][!is.finite(psparticles[[i]])] <- 0
						C<-cov.wt(cbind(t(psparticles[[i]][,index1]),t(psparticles[[1]][,index])),wt=wsparticles[[lag+i]])
						covhats<-covhats+C$cov/Np[1]
 					}

				}
				C<-cov.wt(cbind(t(psparticles[[lag+1]][,1:(Np[1])]),t(psparticles[[1]][,index])),wt=xx$weight)
				covhats<-covhats+C$cov[1:nrow(params),(nrow(params)+1):(2*nrow(params))]/Np[1]
				for (i in 1:(2*lag-1)){
                			psparticles[[i]]<-psparticles[[i+1]]
					wsparticles[[i]]<-wsparticles[[i+1]]
					asparticles[[i+1]]<-asparticles[[i+2]]  
              			}
				
            			psparticles[[2*lag]] <- params
				wsparticles[[2*lag]] <- xx$weight
                		asparticles[[2*lag+1]] <- xx$pa+1
			}

			if(nt==ntimes){
            			for(j in 1:2*lag){
					index<-unlist(ancestor(plist=asparticles,t=lag+1, lag=lag, 1:(Np[1])))
					psparticles[[1]][!is.finite(psparticles[[1]])] <- 0
					C<-cov.wt(t(psparticles[[1]][,index]),wt=wsparticles[[lag+1]])
					phats<-phats+C$center
					covhats<-covhats+C$cov/Np[1]
				
					if(lag>1){
						for(i in 2:lag){
							index1<-unlist(ancestor(plist=asparticles,t=lag+1, lag=lag+1-i, 1:(Np[1])))
							psparticles[[i]][!is.finite(psparticles[[i]])] <- 0
							C<-cov.wt(cbind(t(psparticles[[i]][,index1]),t(psparticles[[1]][,index])),wt=wsparticles[[lag+i]])
							covhats<-covhats+C$cov/Np[1]
						}
					}
					C<-cov.wt(cbind(t(psparticles[[lag+1]][,1:(Np[1])]),t(psparticles[[1]][,index])),wt=xx$weight)
					covhats<-covhats+C$cov[1:nrow(params),(nrow(params)+1):(2*nrow(params))]/Np[1]
				
				
					for (i in 1:(2*lag-1)){
                				psparticles[[i]]<-psparticles[[i+1]]
						wsparticles[[i]]<-wsparticles[[i+1]]
						asparticles[[i+1]]<-asparticles[[i+2]]  
              				}
				
            				psparticles[[2*lag]] <- params
					wsparticles[[2*lag]] <- xx$weight
                			asparticles[[2*lag+1]] <- xx$pa+1
				}
      	    		}
            	
        }
    }
    
    if (!is.null(seed)){
        assign(".Random.seed",save.seed,pos=.GlobalEnv)
        seed <- save.seed
    }

    if (nfail>0)
        warning(sprintf(ngettext(nfail,msg1="%d filtering failure occurred in ",
            msg2="%d filtering failures occurred in "),nfail),
            sQuote("pfilter2"),call.=FALSE
        )
    new(
        "pfilterd2.pomp",
        object,
        pred.mean=pred.m,
        pred.var=pred.v,
        filter.mean=filt.m,
        paramMatrix= params,
        eff.sample.size=eff.sample.size,
        cond.loglik=loglik,
        saved.states=xparticles,
        saved.params=pparticles,
        seed=as.integer(seed),
        Np=as.integer(Np),
        tol=tol,
        nfail=as.integer(nfail),
        loglik=sum(loglik),
        phats=phats,
        covhats=covhats,
        pcovhats=pcovhats,
        lag=lag
    )
}

setMethod(
    "pfilter2",
    signature=signature(object="pomp"),
    function(
        object, params, Np,
        tol = 1e-17,
        max.fail = Inf,
        pred.mean = FALSE,
        pred.var = FALSE,
        filter.mean = FALSE,
        save.states = FALSE,
        save.params = FALSE,
        lag=0,
        seed = NULL,
        verbose = getOption("verbose"),
        ...
    ){
        if (missing(params)) 
            params <- coef(object)
        pfilter2.internal(
            object=object,
            params=params,
            Np=Np,
            tol=tol,
            max.fail=max.fail,
            pred.mean=pred.mean,
            pred.var=pred.var,
            filter.mean=filter.mean,
            save.states=save.states,
            save.params=save.params,
            lag=lag,
            seed=seed,
            verbose=verbose,
            .transform=FALSE,
            ...
        )
    }
)

setMethod(
    "pfilter2",
    signature=signature(object="pfilterd2.pomp"),
    function (object, params, Np, tol, ...){
        if (missing(params)) 
            params <- coef(object)
        if (missing(Np)) 
            Np <- object@Np
        if (missing(tol)) 
            tol <- object@tol
        pfilter2(
            object=as(object,"pomp"),
            params=params,
            Np=Np,
            tol=tol,
            ...
        )
    }
)
