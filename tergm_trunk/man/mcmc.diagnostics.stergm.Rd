\name{mcmc.diagnostics}
\alias{mcmc.diagnostics.stergm}
\title{Conduct MCMC diagnostics on an ergm or stergm fit}
\description{
  This function prints diagnistic information and creates simple diagnostic plots for the MCMC sampled
  statistics produced from a \code{\link{stergm}} fit.
}
\usage{
\method{mcmc.diagnostics}{stergm} (object,  
                                  center=TRUE,
                                  curved=TRUE,
                                  vars.per.page=3,
                                  \dots) 
}
\arguments{
  \item{object}{A stergm object.  See documentation for
    \code{\link{stergm}}.
  }
  \item{center}{Logical:  If TRUE, ; center the samples
    on the observed statistics.
  }
  \item{curved}{Logical:  If TRUE, summarize the curved statistics 
    (evaluated at the MLE of any non-linear parameters), 
    rather than the raw components of the curved statistics.
  }
  \item{vars.per.page}{Number of rows (one variable per row) 
    per plotting page.  Ignored
    if \code{latticeExtra} package is not installed.
  }
  \item{\dots}{Additional arguments, to be passed to plotting functions.
  }
}

\details{ The plots produced are a trace of the sampled output and a
  density estimate for each variable in the chain.  The diagnostics
  printed include correlations and convergence diagnostics.
  
  In fact, an \code{object} contains the matrix of
  statistics from the MCMC run as component \code{$sample}.
  This matrix is actually an object of class \code{mcmc} and
  can be used directly in the \code{coda} package to assess MCMC
  convergence. \emph{Hence all MCMC diagnostic methods available
  in \code{coda} are available directly.} See the examples and
  \url{http://www.mrc-bsu.cam.ac.uk/bugs/classic/coda04/readme.shtml}.

  More information can be found by looking at the documentation of
  \code{\link{stergm}}.
  
}
\references{
%  Warnes, G.W. (2000).  Multi-Chain and Parallel Algorithms for Markov
%  Chain Monte Carlo. Dissertation, Department of Biostatistics,
%  University of Washington,
%
  Raftery, A.E. and Lewis, S.M. (1992).  One long run with diagnostics:
  Implementation strategies for Markov chain Monte Carlo. Statistical
  Science, 7, 493-497.

  Raftery, A.E. and Lewis, S.M. (1995).  The number of iterations,
  convergence diagnostics and generic Metropolis algorithms.  In
  Practical Markov Chain Monte Carlo (W.R. Gilks, D.J. Spiegelhalter and
  S. Richardson, eds.). London, U.K.: Chapman and Hall.

  This function is based on the \code{coda} package 
  It is based on the the
  R function \code{raftery.diag} in \code{coda}.  \code{raftery.diag},
  in turn, is based on the FORTRAN program \code{gibbsit} written by
  Steven Lewis which is available from the Statlib archive.
}

\value{\code{\link{mcmc.diagnostics.ergm}} returns 
  some degeneracy information, if it is included in the original
  object.  The function is mainly used for its side effect, which is
  to produce plots and summary output based on those plots.
}
\seealso{\code{\link{ergm}}, \code{\link{stergm}},\code{network} package,
         \code{coda} package, 
         \code{\link{summary.ergm}}
        }
\keyword{models}
