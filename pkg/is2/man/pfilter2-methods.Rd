\name{pfilter2-methods}
\alias{logLik,pfilterd2.pomp-method}
\alias{logLik-pfilterd2.pomp}
\alias{$,pfilterd2.pomp-method}
\alias{$-pfilterd2.pomp}
\alias{pred.mean}
\alias{pred.mean,pfilterd2.pomp-method}
\alias{pred.mean-pfilterd2.pomp}
\alias{pred.var}
\alias{pred.var,pfilterd2.pomp-method}
\alias{pred.var-pfilterd2.pomp}
\alias{filter.mean}
\alias{filter.mean,pfilterd2.pomp-method}
\alias{filter.mean-pfilterd2.pomp}
\alias{eff.sample.size}
\alias{eff.sample.size,pfilterd2.pomp-method}
\alias{eff.sample.size-pfilterd2.pomp}
\alias{cond.logLik}
\alias{cond.logLik,pfilterd2.pomp-method}
\alias{cond.logLik-pfilterd2.pomp}
\alias{as,pfilterd2.pomp-method}
\alias{as.data.frame.pfilterd2.pomp}
\alias{coerce,pfilterd2.pomp,data.frame-method}
\title{Methods of the "pfilterd2.pomp" class}
\description{Methods of the "pfilterd2.pomp" class.}
\usage{
\S4method{logLik}{pfilterd2.pomp}(object, \dots)
\S4method{pred.mean}{pfilterd2.pomp}(object, pars, \dots)
\S4method{pred.var}{pfilterd2.pomp}(object, pars, \dots)
\S4method{filter.mean}{pfilterd2.pomp}(object, pars, \dots)
\S4method{eff.sample.size}{pfilterd2.pomp}(object, \dots)
\S4method{cond.logLik}{pfilterd2.pomp}(object, \dots)
\S4method{as}{pfilterd2.pomp}(object, class)
\S4method{coerce}{pfilterd2.pomp,data.frame}(from, to = "data.frame", strict = TRUE)
\S3method{as.data.frame}{pfilterd2.pomp}(x, row.names, optional, \dots)
}
\arguments{
  \item{object, x}{
    An object of class \code{pfilterd2.pomp} or inheriting class \code{pfilterd2.pomp}.
  }
  \item{pars}{Names of parameters.}
  \item{class}{
    character;
    name of the class to which \code{object} should be coerced.
  }
  \item{from, to}{
    the classes between which coercion should be performed.
  }
  \item{strict}{
    ignored.
  }
  \item{row.names, optional}{
    ignored.
  }
  \item{\dots}{
    Additional arguments unused at present.
  }
}
\author{{Dao Nguyen \email{nguyenxd at umich dot edu}}, {Edward L. Ionides \email{ionides at umich dot edu}}}
\seealso{\code{\link{pfilter2}}}
\keyword{models}
\keyword{ts}
