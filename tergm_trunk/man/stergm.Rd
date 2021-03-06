\name{stergm}
\alias{stergm}
\alias{stergm.object}
\title{Separable Temporal Exponential Family Random Graph Models
}
\description{
  \code{\link{stergm}} is used for finding Separable Temporal ERGMs'
  (STERGMs) Conditional MLE (CMLE) (Krivitsky and Handcock, 2010) and
  Equilibrium Generalized Method of Moments Estimator (EGMME) (Krivitsky, 2009).
}
\usage{
  stergm (nw, 
          formation, 
          dissolution,
          constraints = ~.,
          estimate, 
          times=NULL, 
          offset.coef.form=NULL, 
          offset.coef.diss=NULL,
          targets=NULL, 
          target.stats=NULL,
          eval.loglik=TRUE,
          control=control.stergm(),
          verbose=FALSE, 
          \dots)
}
\arguments{
  \item{nw}{A \code{\link[network]{network}} object (for EGMME); or
    \code{\link[networkDynamic]{networkDynamic}} object, a
    \code{\link{network.list}} object, or a \code{\link{list}}
    containing networks (for CMLE and CMPLE).

    \code{simulate.network} understands the \code{\link{lasttoggle}}
    "API".
  }
  \item{formation, dissolution}{
    One-sided \code{\link{ergm}}-style formulas for the formation and
    dissolution models, respectively.
  }

  \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being modeled,
    using syntax similar to the \code{formula} argument. Multiple constraints
    may be given, separated by \dQuote{+} operators.
    Together with the model terms in the formula and the reference measure, the constraints
    define the distribution of networks being modeled.

    It is also possible to specify a proposal function directly
    by passing a string with the function's name. In that case,
    arguments to the proposal should be specified through the
    \code{prop.args} argument to \code{\link{control.ergm}}.

    The default is \code{~.}, for an unconstrained model.

    See the \link[=ergm-constraints]{ERGM constraints} documentation for
    the constraints implemented in the \code{\link[=ergm-package]{ergm}}
    package. Other packages may add their own constraints.

    For STERGMs in particular, the constraints apply to the
    post-formation and the post-dissolution network, rather than the
    final network. This means, for example, that if the degree of all
    vertices is constrained to be less than or equal to three, and a vertex
    begins a time step with three edges, then, even if one of its edges
    is dissolved during its time step, it won't be able to form another
    edge until the next time step. This behavior may change in the future.
    
    Note that not all possible combinations of constraints are supported.
  }
  
  \item{estimate}{One of "EGMME" for Equilibrium Generalized Method
    of Moments Estimation, based on a single network with some temporal
    information and making an assumption that it is a product of a
    STERGM process running to its stationary (equilibrium) distribution; "CMLE" 
    for Conditional Maximum Likelihood Estimation, modeling a transition
    between two networks, or
    "CMPLE" for Conditional Maximum PseudoLikelihood Estimation, using
    MPLE instead of MLE.
    CMPLE is extremely inaccurate at this time.
  }
  \item{times}{For CMLE and CMPLE estimation, times or indexes at which
    the networks whose transition is to be modeled are
    observed. Default to \code{c(0,1)} if \code{nw} is a
    \code{\link[networkDynamic]{networkDynamic}} and to \code{c(1:2)}
    (first two elements) if \code{nw} is a \code{\link{network.list}} or
    a \code{\link{list}}. Unused for EGMME.
  }
  \item{offset.coef.form}{Numeric vector to specify offset formation
    parameters.
  }
  \item{offset.coef.diss}{Numeric vector to specify offset dissolution
    parameters.
  }
  \item{targets}{One-sided \code{\link{ergm}}-style formula specifying statistics whose
    moments are used for the EGMME. Unused for CMLE and CMPLE.
  }

  \item{target.stats}{
    A vector specifying the values of the \code{targets} statistics that
    EGMME will try to match. 
    Defaults to the
    statistics of \code{nw}. Unused for CMLE and CMPLE.
  }
  \item{eval.loglik}{Whether or not to calculate the log-likelihood of
    a CMLE STERGM fit. See \code{\link{ergm}} for details.}
  \item{control}{
    A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.stergm}}. 
  }
  \item{verbose}{
    logical or integer; if TRUE or positive, the program will print out progress
    information. Higher values result in more output.
  }
  \item{\dots}{Additional arguments, to be passed to lower-level 
    functions.
  }  
}
\section{Model Terms}{
  See \code{\link{ergm}} and  \code{\link{ergm-terms}} for details. At
  this time, only linear ERGM terms are allowed.
}
\value{\code{\link{stergm}} returns an object of class \code{\link{stergm}} that is a list
  consisting of the following elements:
  \item{formation, dissolution}{Formation and dissolution formulas,
    respectively.}
  \item{targets}{The targets formula.}
  \item{target.stats}{The target statistics.}
  \item{estimate}{The type of estimate.}
  \item{opt.history}{A matrix containing the full trace of the EGMME
    optimization process: coefficients tried and target statistics
    simulated.}
  \item{sample}{An \code{\link{mcmc}} object containing target
    statistics sampled at the estimate.}
  \item{covar}{The full estimated variance-covariance matrix of the
    parameter estimates for EGMME. (Note that although the CMLE
    formation parameter estimates are independent of the dissolution
    parameter estimates due to the separability assumption, this is not
    necessarily the case for EGMME.) }
  \item{formation.fit, dissolution.fit}{For CMLE and CMPLE,
    \code{\link{ergm}} objects from fitting formation and dissolution,
    respectively. For EGMME, stripped down \code{\link{ergm}}-like lists.}
  \item{network}{For \code{estimate=="EGMME"}, the original network; for
    \code{estimate=="CMLE"} or \code{estimate=="CMPLE"}, a
    \code{\link{network.list}} (a discrete series of networks) to which the
    model was fit.}
  \item{control}{The control parameters used to fit the model.}
  
  See the method \code{\link{print.stergm}} for details on how
  an \code{\link{stergm}} object is printed.  Note that the
  method \code{\link{summary.stergm}} returns a summary of the
  relevant parts of the \code{\link{stergm}} object in concise summary
  format.
}

\examples{
\donttest{
# EGMME Example
par(ask=FALSE)
n<-30
g0<-network.initialize(n,dir=FALSE)

#                     edges, degree(1), mean.age
target.stats<-c(      n*1/2,    n*0.6,        20)

dynfit<-stergm(g0,formation = ~edges+degree(1), dissolution = ~edges,
               targets = ~edges+degree(1)+mean.age,
               target.stats=target.stats, estimate="EGMME",
               control=control.stergm(SA.plot.progress=TRUE))

par(ask=TRUE)
mcmc.diagnostics(dynfit)
summary(dynfit)

# CMLE Example
data(samplk)

# Fit a transition from Time 1 to Time 2
samplk12 <- stergm(list(samplk1, samplk2),
                   formation=~edges+mutual+transitiveties+cyclicalties,
                   dissolution=~edges+mutual+transitiveties+cyclicalties,
                   estimate="CMLE")

mcmc.diagnostics(samplk12)
summary(samplk12)

# Fit a transition from Time 1 to Time 2 and from Time 2 to Time 3 jointly
samplk123 <- stergm(list(samplk1, samplk2, samplk3),
                    formation=~edges+mutual+transitiveties+cyclicalties,
                    dissolution=~edges+mutual+transitiveties+cyclicalties,
                    estimate="CMLE")

mcmc.diagnostics(samplk123)
summary(samplk123)
}
}

\references{
  \itemize{
    \item Krivitsky PN, Handcock MS (2010). A Separable Model for Dynamic
  Networks. \url{http://arxiv.org/abs/1011.1937}
  
  \item Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
  Egocentric Data with Durational Information. \emph{Pennsylvania State
    University Department of Statistics Technical Report},
  2012(2012-01). \url{http://stat.psu.edu/Research/2012-technical-reports}
}
}

\seealso{ergm, network, \%v\%, \%n\%, \code{\link{ergm-terms}}}