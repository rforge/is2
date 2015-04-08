\name{is2-methods}
\docType{methods}
\alias{is2-methods}
\alias{logLik,is2-method}
\alias{logLik-is2}
\alias{conv.rec}
\alias{conv.rec,is2-method}
\alias{conv.rec-is2}
\alias{conv.rec,is2List-method}
\alias{conv.rec-is2List}
\alias{plot-is2}
\alias{plot,is2-method}
\alias{plot-is2List}
\alias{plot,is2List-method}
\alias{is2List-class}
\alias{c-is2}
\alias{c,is2-method}
\alias{c-is2List}
\alias{c,is2List-method}
\alias{[-is2List}
\alias{[,is2List-method}
\alias{compare.is2}
\title{Methods of the "is2" class}
\description{Methods of the \code{is2} class.}
\usage{
\S4method{logLik}{is2}(object, \dots)
\S4method{conv.rec}{is2}(object, pars, transform = FALSE, \dots)
\S4method{conv.rec}{is2List}(object, \dots)
\S4method{plot}{is2}(x, y, \dots)
\S4method{plot}{is2List}(x, y, \dots)
\S4method{c}{is2}(x, \dots, recursive = FALSE)
\S4method{c}{is2List}(x, \dots, recursive = FALSE)
compare.is2(z)
}
\arguments{
  \item{object}{The \code{is2} object.}
  \item{pars}{Names of parameters.}
  \item{x}{The \code{is2} object.}
  \item{y, recursive}{Ignored.}
  \item{z}{A \code{is2} object or list of \code{is2} objects.}
  \item{transform}{
    optional logical;
    should the parameter transformations be applied?
    See \code{\link[=coef-pomp]{coef}} for details.
  }
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
  }
}
\section{Methods}{
  \describe{
    \item{conv.rec}{
      \code{conv.rec(object, pars = NULL)} returns the columns of the convergence-record matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{logLik}{
      Returns the value in the \code{loglik} slot.
    }
    \item{c}{
      Concatenates \code{is2} objects into an \code{is2List}.
    }
    \item{plot}{
      Plots a series of diagnostic plots.
    }
    \item{compare.is2}{
      Deprecated: use \code{plot} instead.
    }
  }
}
\author{{Dao Nguyen \email{nguyenxd at umich dot edu}}, {Edward L. Ionides \email{ionides at umich dot edu}}}
\seealso{\code{\link{is2}}, \code{\link{pfilter2}}}
\keyword{models}
\keyword{ts}
