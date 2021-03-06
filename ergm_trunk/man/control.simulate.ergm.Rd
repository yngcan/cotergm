\name{control.simulate}
\alias{control.simulate}
\alias{control.simulate.formula}
\alias{control.simulate.formula.ergm}
\alias{control.simulate.ergm}
\title{ Auxiliary for Controlling ERGM Simulation }
\description{
  Auxiliary function as user interface for fine-tuning ERGM simulation.
}
\usage{
control.simulate(MCMC.burnin=1000,
                 MCMC.interval=1000,
                 MCMC.prop.weights="default",
                 MCMC.prop.args=list(),
                 MCMC.init.maxedges=20000,
                 MCMC.packagenames=c(),
                 MCMC.runtime.traceplot=FALSE,  

                 network.output="network",
                                                     
                 parallel=0,
                 parallel.type=NULL,
                 parallel.version.check=TRUE,
                 \dots)

control.simulate.formula(MCMC.burnin=1000,
                         MCMC.interval=1000,
                         MCMC.prop.weights="default",
                         MCMC.prop.args=list(),
                         MCMC.init.maxedges=20000,
                         MCMC.packagenames=c(),
                         MCMC.runtime.traceplot=FALSE,  

                         network.output="network",
                                                     
                         parallel=0,
                         parallel.type=NULL,
                         parallel.version.check=TRUE,
                         \dots)

control.simulate.formula.ergm(MCMC.burnin=1000,
                              MCMC.interval=1000,
                              MCMC.prop.weights="default",
                              MCMC.prop.args=list(),
                              MCMC.init.maxedges=20000,
                              MCMC.packagenames=c(),
                              MCMC.runtime.traceplot=FALSE,  

                              network.output="network",
                                                     
                              parallel=0,
                              parallel.type=NULL,
                              parallel.version.check=TRUE,
                              \dots)

control.simulate.ergm(MCMC.burnin=NULL,
                      MCMC.interval=NULL,
                      MCMC.prop.weights=NULL,
                      MCMC.prop.args=NULL,
                      MCMC.init.maxedges=NULL,
                      MCMC.packagenames=NULL,
                      MCMC.runtime.traceplot=FALSE,

                      network.output="network",

                      parallel=0,
                      parallel.type=NULL,
                      parallel.version.check=TRUE,
                      \dots)
}
\details{
    This function is only used within a call to the \code{\link{simulate}} function.
    See the \code{usage} section in \code{\link{simulate.ergm}} for details.
}
\arguments{
  \item{MCMC.prop.weights}{Specifies the proposal distribution used in the MCMC
    Metropolis-Hastings algorithm.  Possible choices are \code{"TNT"} or 
    \code{"random"}; the \code{"default"} is one of these two, depending on the
    constraints in place (as defined by the \code{constraints}
    argument of the \code{\link{ergm}} function), though not all weights
    may be used with all constraints.  
    The \code{TNT} (tie / no tie) option puts roughly equal weight on selecting a
    dyad with or without a tie as a candidate for toggling, whereas the
    \code{random} option puts equal weight on all possible dyads, though the
    interpretation of \code{random} may change according to the constraints in
    place.  When no constraints are in place, the default is TNT, which 
    appears to improve Markov chain mixing particularly for networks with a
    low edge density, as is typical of many realistic social networks.
  }
  \item{MCMC.prop.args}{An alternative, direct way of specifying additional arguments to proposal. 
  }
  \item{MCMC.burnin}{Number of proposals before any MCMC sampling
    is done. It typically is set to a fairly large number.
  }
  \item{MCMC.interval}{Number of proposals between sampled statistics.
  }
  \item{MCMC.init.maxedges}{Maximum number of edges expected in network.
  }
  \item{MCMC.packagenames}{Names of packages in which to look for change statistic
    functions in addition to those autodetected. This argument should not be
    needed outside of very strange setups.
  }
  \item{MCMC.runtime.traceplot}{Logical:  If TRUE, plot traceplots of 
    the MCMC sample after every MCMC MLE iteration.
  }
  \item{network.output}{R class with which to output 
    networks. The options are "network" (default) and "edgelist.compressed" 
    (which saves space but only supports networks without vertex attributes)
  }
  \item{parallel}{Number of threads in which to run the
    sampling. Defaults to 0 (no parallelism). See the entry on
    \link[=ergm-parallel]{parallel processing} for
    details and troubleshooting.
  }
  \item{parallel.type}{API to use for parallel
    processing. Supported values are \code{"MPI"} and
    \code{"SOCK"}. Defaults to using the \code{snow} package
    default.
  }
  \item{parallel.version.check}{Logical:  If TRUE, check that the version of
    \code{\link[=ergm-package]{ergm}} running on the slave nodes is the
    same as that running on the master node.
  }
  \item{\dots}{Additional arguments, passed to other functions    
    This argument is helpful because it collects any control 
    parameters that have been deprecated; a warning message is
    printed in case of deprecated arguments.
  }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{simulate.ergm}}, \code{\link{simulate.formula}}. 
\code{\link{control.ergm}} performs a 
similar function for
\code{\link{ergm}};
\code{\link{control.gof}} performs a
similar function for \code{\link{gof}}.   }
%\keyword{simulation}
\keyword{models}
