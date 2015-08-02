\name{is2}
\docType{methods}
\alias{is2}
\alias{is2,is2-method}
\alias{is2-is2}
\alias{is2,pfilterd2.pomp-method}
\alias{is2-pfilterd2.pomp}
\alias{is2,pomp-method}
\alias{is2-pomp}
\alias{continue}
\alias{continue,is2-method}
\alias{continue-is2}
\alias{is2-class}
\title{Iterated Smoothing}
\description{
  Iterated smoothing algorithms for estimating the parameters of a partially-observed Markov process.
}
\usage{
\S4method{is2}{pomp}(object, Nis = 1, start, pars, ivps = character(0),
    particles, rw.sd, Np, ic.lag, var.factor, lag,
    cooling.type = c("geometric","hyperbolic"), 
    cooling.fraction, cooling.factor,
    method = c("is2","ris1","is1"),
    tol = 1e-17, max.fail = Inf,
    verbose = getOption("verbose"), transform = FALSE, \dots)
\S4method{is2}{pfilterd2.pomp}(object, Nis = 1, Np, tol, \dots)
\S4method{is2}{is2}(object, Nis, start, pars, ivps,
    particles, rw.sd, Np, ic.lag, lag, var.factor,
    cooling.type, cooling.fraction,
    method, tol, transform, \dots)
\S4method{continue}{is2}(object, Nis = 1, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nis}{
    The number of filtering iterations to perform.
  }
  \item{start}{
    named numerical vector;
    the starting guess of the parameters.
  }
  \item{pars}{
    optional character vector naming the ordinary parameters to be estimated.
    Every parameter named in \code{pars} must have a positive random-walk standard deviation specified in \code{rw.sd}.
    Leaving \code{pars} unspecified is equivalent to setting it equal to the names of all parameters with a positive value of \code{rw.sd} that are not \code{ivps}.
  }
  \item{ivps}{
    optional character vector naming the initial-value parameters (IVPs) to be estimated.
    Every parameter named in \code{ivps} must have a positive random-walk standard deviation specified in \code{rw.sd}.
    If \code{pars} is empty, i.e., only IVPs are to be estimated, see below \dQuote{Using is2 to estimate initial-value parameters only}.
  }
  \item{particles}{
    Function of prototype \code{particles(Np,center,sd,...)} which sets up the starting particle matrix by drawing a sample of size \code{Np} from the starting particle distribution centered at \code{center} and of width \code{sd}.
    If \code{particles} is not supplied by the user, the default behavior is to draw the particles from a multivariate normal distribution with mean \code{center} and standard deviation \code{sd}.
  }
  \item{rw.sd}{
    numeric vector with names; the intensity of the random walk to be applied to parameters.
    The random walk is only applied to parameters named in \code{pars} (i.e., not to those named in \code{ivps}).
    The algorithm requires that the random walk be nontrivial, so each element in \code{rw.sd[pars]} must be positive.
    \code{rw.sd} is also used to scale the initial-value parameters (via the \code{particles} function).
    Therefore, each element of \code{rw.sd[ivps]} must be positive.
    The following must be satisfied:
    \code{names(rw.sd)} must be a subset of \code{names(start)},
    \code{rw.sd} must be non-negative (zeros are simply ignored),
    the name of every positive element of \code{rw.sd} must be in either \code{pars} or \code{ivps}.
  }
  \item{Np}{
    the number of particles to use in filtering.
    This may be specified as a single positive integer, in which case the same number of particles will be used at each timestep.
    Alternatively, if one wishes the number of particles to vary across timestep, one may specify \code{Np} either as a vector of positive integers (of length \code{length(time(object,t0=TRUE))}) or as a function taking a positive integer argument.
    In the latter case, \code{Np(k)} must be a single positive integer, representing the number of particles to be used at the \code{k}-th timestep:
    \code{Np(0)} is the number of particles to use going from \code{timezero(object)} to \code{time(object)[1]},
    \code{Np(1)}, from \code{timezero(object)} to \code{time(object)[1]},
    and so on, while when \code{T=length(time(object,t0=TRUE))},
    \code{Np(T)} is the number of particles to sample at the end of the time-series.
  }
  \item{ic.lag}{
    a positive integer;
    the timepoint for fixed-lag smoothing of initial-value parameters.
    The \code{is2} update for initial-value parameters consists of replacing them by their filtering mean at time \code{times[ic.lag]}, where \code{times=time(object)}.
    It makes no sense to set \code{ic.lag>length(times)};
    if it is so set, \code{ic.lag} is set to \code{length(times)} with a warning.
  }
  \item{var.factor}{
    a positive number;
    the scaling coefficient relating the width of the starting particle distribution to \code{rw.sd}.
    In particular, the width of the distribution of particles at the start of the first is2 iteration will be \code{random.walk.sd*var.factor}.
  }
  \item{lag}{
    a positive integer;
    the timepoint for fixed-lag smoothing parameters estimation.
  }
  \item{cooling.type, cooling.fraction, cooling.factor}{
    specifications for the cooling schedule, i.e., the manner in which the intensity of the parameter perturbations is reduced with successive filtering iterations.
    \code{cooling.type} specifies the nature of the cooling schedule.
    When \code{cooling.type}-\code{="geometric"}, on the n-th is2 iteration, the relative perturbation intensity is \code{cooling.fraction}-\code{^(n/50)}.
    When \code{cooling.type="hyperbolic"}, on the n-th is2 iteration, the relative perturbation intensity is \code{(s+1)/(s+n)}, where \code{(s+1)/}-\code{(s+50)=cooling.fraction}.
    \code{cooling.fraction} is the relative magnitude of the parameter perturbations after 50 is2 iterations.
    \code{cooling.factor} is now deprecated:
    to achieve the old behavior, use \code{cooling.type="geometric"} and \code{cooling.fraction}-
    \code{=(cooling.factor)^50}.
  }
  \item{method}{
    \code{method} sets the update rule used in the algorithm.
    \code{method="is2"} uses the iterated smoothing update rule (Ionides 2015 submitted);
    \code{method="ris1"} uses the reduced iterated smoothing update rule (Doucet 2013, Ionides 2015 submitted);
    \code{method="is1"} uses the reduced iterated smoothing update rule (Doucet 2013);
  }
  \item{tol}{
    See the description under \code{\link{pfilter2}}.
  }
  \item{max.fail}{
    See the description under \code{\link{pfilter2}}.
  }
  \item{verbose}{
    logical; if TRUE, print progress reports.
  }
  \item{transform}{
    logical;
    if \code{TRUE}, optimization is performed on the transformed scale.
  }
  \item{\dots}{
    additional arguments that override the defaults.
  }
}
\section{Re-running is2 Iterations}{
  To re-run a sequence of is2 iterations, one can use the \code{is2} method on a \code{is2} object.
  By default, the same parameters used for the original is2 run are re-used (except for \code{weighted}, \code{tol}, \code{max.fail}, and \code{verbose}, the defaults of which are shown above).
  If one does specify additional arguments, these will override the defaults.
}
\section{Continuing is2 Iterations}{
  One can resume a series of is2 iterations from where one left off using the \code{continue} method.
  A call to \code{is2} to perform \code{Nis=m} iterations followed by a call to \code{continue} to perform \code{Nis=n} iterations will produce precisely the same effect as a single call to \code{is2} to perform \code{Nis=m+n} iterations.
  By default, all the algorithmic parameters are the same as used in the original call to \code{is2}.
  Additional arguments will override the defaults.
}
\section{Using is2 to estimate initial-value parameters only}{
  One can use is2 fixed-lag smoothing to estimate only initial value parameters (IVPs).
  In this case, \code{pars} is left empty and the IVPs to be estimated are named in \code{ivps}.
  If \code{theta} is the current parameter vector, then at each is2 iteration, \code{Np} particles are drawn from a distribution centered at \code{theta} and with width proportional to \code{var.factor*rw.sd}, a particle filtering operation is performed, and \code{theta} is replaced by the filtering mean at \code{time(object)[ic.lag]}.
  Note the implication that, when \code{is2} is used in this way on a time series any longer than \code{ic.lag}, unnecessary work is done.
  If the time series in \code{object} is longer than \code{ic.lag}, consider replacing \code{object} with \code{window(object,end=ic.lag)}. 
}
\section{Details}{
  If \code{particles} is not specified, the default behavior is to draw the particles from a multivariate normal distribution.
  It is the user's responsibility to ensure that, if the optional \code{particles} argument is given, that the \code{particles} function satisfies the following conditions:

  \code{particles} has at least the following arguments:
  \code{Np}, \code{center}, \code{sd}, and \code{\dots}.
  \code{Np} may be assumed to be a positive integer;
  \code{center} and \code{sd} will be named vectors of the same length.
  Additional arguments may be specified;
  these will be filled with the elements of the \code{userdata} slot of the underlying \code{pomp} object (see \code{\link{pomp}}).

  \code{particles} returns a \code{length(center)} x \code{Np} matrix with rownames matching the names of \code{center} and \code{sd}.
  Each column represents a distinct particle.

  The center of the particle distribution returned by \code{particles} should be \code{center}.
  The width of the particle distribution should vary monotonically with \code{sd}.
  In particular, when \code{sd=0}, the \code{particles} should return matrices with \code{Np} identical columns, each given by the parameters specified in \code{center}.
}
\references{
  D. Nguyen, E. L. Ionides, 
  A second-order iterated smoothing,
  Journal of Statistics and Computing, 2015 (submitted).

  A., Doucet, P. E. Jacob,  and S. Rubenthaler,
  Derivative-free estimation of the score vector and observed information matrix with application to state-space models, 
  Preprint arXiv:1304.5768.

  A. A. King, D. Nguyen, and E. L. Ionides, 
  Statistical inference for partially observed Markov processes via the R package pomp, 
  Journal of Statistical Software, 2015 (just-accepted).
}

\author{{Dao Nguyen \email{nguyenxd at umich dot edu}}, {Edward L. Ionides \email{ionides at umich dot edu}}}

\seealso{
  \code{\link{is2-methods}}, \code{\link{pfilter2}}.
  See the \dQuote{intro_to_is2} vignette for examples.
}
\keyword{ts}
