\name{simulate.ergm}
\alias{simulate.ergm}
\alias{simulate.formula}
\alias{simulate.formula.ergm}

\title{Draw from the distribution of an Exponential Family Random Graph Model}
\usage{
\method{simulate}{formula}(object, nsim=1, seed=NULL,
                           coef, 
                           response=NULL, reference=~Bernoulli,
                           constraints=~.,
                           monitor=NULL,
                           basis=NULL,
                           statsonly=FALSE,
                           sequential=TRUE,
                           control=control.simulate.formula(),
                           verbose=FALSE, 
                           \dots)
\method{simulate}{ergm}(object, nsim=1, seed=NULL, 
                        coef=object$coef,
                        response=object$response, reference=object$reference,
                        constraints=object$constraints,
                        monitor=NULL,
                        statsonly=FALSE,
                        sequential=TRUE,
                        control=control.simulate.ergm(),
                        verbose=FALSE, 
                        \dots)
}
\description{
  \code{\link[stats]{simulate}} is used to draw from exponential family 
  random network models in their natural parameterizations.
  See \code{\link{ergm}} for more information on these models.
}
\arguments{
  \item{object}{an \R object. Either
    a \code{\link{formula}} or an \code{\link{ergm}} object.
    The \code{\link{formula}} should be of the form \code{y ~ <model terms>},
    where \code{y} is a network object or a matrix that can be
    coerced to a \code{\link[network]{network}} object.  For the details 
    on the possible
    \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
    \code{\link[network]{network}} object in \R, use the 
    \code{network()} function,
    then add nodal attributes to it using the \code{\%v\%}
    operator if necessary.
  }
  \item{nsim}{Number of networks to be randomly drawn
    from the given distribution on the set of all networks,
    returned by the Metropolis-Hastings algorithm.
  }
  \item{seed}{Random number integer seed.  
    See \code{\link[base]{set.seed}}.
  }
  \item{coef}{Vector of parameter values for the model from which the 
    sample is to be drawn.  If \code{object} is of class \code{ergm}, the 
    default value is the vector of estimated coefficients.
  }
  \item{response}{\emph{EXPERIMENTAL.} Name of the edge attribute whose value is to be
    modeled. Defaults to \code{NULL} for simple presence or
    absence, modeled via binary ERGM terms. Passing
    anything but \code{NULL} uses valued ERGM terms.}
  
  \item{reference}{\emph{EXPERIMENTAL.} A one-sided formula specifying
    the reference measure (\eqn{h(y)}) to be used. (Defaults to \code{~Bernoulli}.)
    See help for \link[=ergm-references]{ERGM reference measures} implemented in the
    \code{\link[=ergm-package]{ergm}} package.}
  
  \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being
    simulated. See the documentation for a similar argument for
    \code{\link{ergm}} and see \link[=ergm-constraints]{list of
    implemented constraints} for more information. For
    \code{simulate.formula}, defaults to no constraints. For
    \code{simulate.ergm}, defaults to using the same constraints as
    those with which \code{object} was fitted.
  }
  \item{monitor}{A one-sided formula specifying one or more terms whose
    value is to be monitored. These terms are appeneded to the model, along
    with a coefficient of 0, so their statistics are returned.
  }
  \item{basis}{An optional \code{\link[network]{network}} object to start
    the Markov chain.  If omitted, the default is the left-hand-side of 
    the \code{formula}.
    If neither a left-hand-side nor a \code{basis} is present, an error 
    results because the characteristics of the network (e.g., size and
    directedness) must be specified.
  }
  \item{statsonly}{Logical:  If TRUE, return only the network statistics, 
    not the network(s) themselves.
  }
  \item{sequential}{Logical:  If FALSE, each of the \code{nsim} simulated
    Markov chains begins at the initial network.  If TRUE, the end of one 
    simulation is used as the start of the next.
    Irrelevant when \code{nsim=1}.
  }
  \item{control}{A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.simulate.ergm}} or
    \code{\link{control.simulate.formula}}, which have different defaults.
  }
  \item{verbose}{Logical:  If TRUE, extra information is printed as
    the Markov chain progresses.
  }
  \item{\dots}{Further arguments passed to or used by methods.
  }
}

\details{
  A sample of networks is randomly drawn from the specified model.  The
  model is specified by the first argument of the function.  If
  the first argument is a \code{\link{formula}} then this defines 
  the model.  If the
  first argument is the output of a call to \code{\link{ergm}} then the model
  used for that call is the one fit - and unless \code{coef} is specified,
  the sample is from the MLE of the parameters.  
  If neither of those are given as the
  first argument then a Bernoulli network is generated with the probability
  of ties defined by \code{prob} or \code{coef}.

  Note that the first network is sampled after \code{burnin} + \code{interval}
  steps, and any subsequent networks are sampled each
  \code{interval} steps after the first.

  More information can be found by looking at the documentation of
  \code{\link{ergm}}.
}

\value{
  If \code{statsonly==TRUE} a matrix containing the simulated network
  statistics. If \code{control$parallel>0}, the statistics from each
  Markov chain are stacked.


  Otherwise, if \code{nsim==1}, an object of class \code{network}.  If
  \code{nsim>1}, it returns an object of class
  \code{\link{network.list}}: a list of networks with the following
  \code{\link{attr}}-style attributes on the list: \item{formula}{The
  \code{\link{formula}} used to generate the sample.}  \item{stats}{The
  \eqn{\code{nsim}\times p} matrix of network statistics, where \eqn{p}
  is the number of network statistics specified in the model.}
  \item{control}{Control parameters used to generate the sample.}
  \item{constraints}{Constraints used to generate the sample.}
  \item{reference}{The reference measure for the sample.}
  \item{monitor}{The monitoring formula.}  \item{response}{The edge
    attribute used as a response.}

  If \code{statsonly==FALSE && control$parallel>0} the returned networks
  are "interleaved", in the sense that for \code{y[i,j]} is the
  \code{j}th network from MCMC chain \code{i}, the sequence returned if
  \code{control$parallel==2} is \code{list(y[1,1], y[2,1], y[1,2],
  y[2,2], y[1,3], y[2,3], ...)}. This is different from the behavior
  when \code{statsonly==TRUE}. This detail may change in the future.
  
  This object has summary and print methods.
}
\seealso{\code{\link{ergm}}, \code{\link[network]{network}}}
\examples{
#
# Let's draw from a Bernoulli model with 16 nodes
# and density 0.5 (i.e., coef = c(0,0))
#
g.sim <- simulate(network(16) ~ edges + mutual, coef=c(0, 0))
#
# What are the statistics like?
#
summary(g.sim ~ edges + mutual)
#
# Now simulate a network with higher mutuality
#
g.sim <- simulate(network(16) ~ edges + mutual, coef=c(0,2))
#
# How do the statistics look?
#
summary(g.sim ~ edges + mutual)
#
# Let's draw from a Bernoulli model with 16 nodes
# and tie probability 0.1
#
g.use <- network(16,density=0.1,directed=FALSE)
#
# Starting from this network let's draw 3 realizations
# of a edges and 2-star network
#
g.sim <- simulate(~edges+kstar(2), nsim=3, coef=c(-1.8,0.03),
               basis=g.use, control=control.simulate(
                 MCMC.burnin=1000,
                 MCMC.interval=100))
g.sim
summary(g.sim)
#
# attach the Florentine Marriage data
#
data(florentine)
#
# fit an edges and 2-star model using the ergm function
#
gest <- ergm(flomarriage ~ edges + kstar(2))
summary(gest)
#
# Draw from the fitted model (satatistics only), and observe the number
# of triangles as well.
#
g.sim <- simulate(gest, nsim=10, 
            monitor=~triangles, statsonly=TRUE,
            control=control.simulate.ergm(MCMC.burnin=1000, MCMC.interval=100))
g.sim
}
\keyword{models}
