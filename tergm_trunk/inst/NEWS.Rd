\name{NEWS}
\title{NEWS file for the \code{\link[=tergm-package]{tergm}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.
\section{Changes in version 3.1.3}{

  \subsection{BUG FIXES}{
    \itemize{
      \item Some Metropolis-Hastings proposal functions would sometimes
      return incorrect acceptance probabilities when combined with the
      \code{\link[=ergm-constraints]{bd}} constraint. This has been
      fixed.

      \item In \code{\link{simulate.networkDynamic}}, an error in which
      vertices were queried and updated in the simulation has been
      fixed.
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{\link{simulate.networkDynamic}} and
      \code{\link{simulate.network}} functions now take an additional
      argument, \code{time.offset}. See help for those functions for
      details.
    }
  }

  % \subsection{NEW FEATURES}{
  %   \itemize{
  %   }
  % }
  
  % \subsection{DEPRECATED & DEFUNCT}{
  %   \itemize{
  %   }
  % }
  
  % \subsection{C-LEVEL FACILITIES}{
  %   \itemize{
  %     \item It is not clear that we need a subsection called "C-LEVEL FACILITIES".
  %   }
  % }

  % \subsection{UTILITIES}{
  %   \itemize{
  %     \item Likewise, "UTILITIES" is probably unnecessary.
  %   }
  % }

  % \subsection{INSTALLATION}{
  %   \itemize{
  %     \item ... as is "INSTALLATION".
  %   }
  % }

  % \subsection{PACKAGE INSTALLATION}{
  %   \itemize{
  %     \item ditto.
  %   }
  % }

}


\section{Changes in version 3.1.2}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item The Tie-NonTie (TNT) proposal has been implemented for
      dissolution phase models. This should improve mixing and inference
      for these models.

      \item \code{stergm} and \code{simulate.stergm} now determine the
      number of Metropolis-Hastings steps per time step adaptively,
      stopping when the formation/dissolution process appears to have
      converged.
    }
  }
  
  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item The previously deprecated start and end attr values attached 
      to networks returned by \code{\link{simulate.network}} and 
      \code{\link{simulate.networkDynamic}} have been removed and 
      replaced with observation spells recorded as a 
      \code{\link{net.obs.period}} network attribute.

      \item \code{MCMC.burnin} and \code{MCMC.interval} arguments to
      \code{control.stergm}, \code{control.simulate.network}, and
      \code{control.simulate.stergm} have been replaced by a different
      mechanism. See the help for the respective control functions for
      more information.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      
      \item No longer generates deprecation warnings about start and end
       attrs when various internal functions use as.data.frame. Resolved
      by replacing start and end with a \code{\link{net.obs.period}}
      object.  

      \item A bug in the check for whether terms were amenable to being
      fit using STERGM CMLE for multiple transitions has been fixed.

      \item Some minor documentation typos have been fixed.
    }
  }
}

\section{Changes in version 3.1.1}{
  \itemize{
    \item Version 3.1.1 has been skipped to ease upgrading for those using a
    preview release.
  }
}

\section{Changes in version 3.1}{

  \subsection{THIS IS AN INITIAL RELEASE}{
    \itemize{
      \item This package consists of the dynamic network modeling code that has been
      split out of the \code{\link[=ergm-package]{ergm}} package.
      
      \item Changes listed in the following sections are relative to
      \code{\link[=ergm-package]{ergm}} version 3.0.
      }
  }
  
  \subsection{NEW FEATURES}{
    \itemize{
      
      \item Although fitting the EGMME for dissolution was possible
      before, it was impractical due to nonidentifiability: for example,
      one cannot fit both \code{\link{edges}} formation and
      \code{\link{edges}} dissolution with only one \code{\link{edges}}
      target statistic. Four new statistics have been added and
      documented that focus on targeting observed and hazard:
      \code{\link{mean.age}}, \code{\link{edgecov.mean.age}},
      \code{\link{degree.mean.age}}, \code{\link{degrange.mean.age}},
      \code{\link{edge.ages}}, \code{\link{edgecov.ages}}, and
      \code{\link{edges.ageinterval}}. This allow jointly fitting
      formation and dissolution.

      \item In addition to the progress plot, EGMME fitting routines can
      now plot the estimated gradient matrix and the matrix of
      correlations among the target statistics.

      \item EGMME fitting is now more adaptive when determining when to
      stop the optimization and return the result.

      \item EGMME fitting can now take advantage of multiple CPUs,
      cores, or cluster nodes for faster and more robust fitting.

      \item EGMME and CMLE fitting can now accommodate a
      \code{constraints} argument. However, note that the constraints
      apply to post-formation (\eqn{y^+=y^0\cup\y^1}) and
      post-dissolution (\eqn{y^-=y^0\cap\y^1}) networks, not to the
      final network (\eqn{y^1}). This may change in the future.

      \item CMLE can now be fit to more than two networks. (Not all ERGM
      terms and constraints can be used in this mode, however.)

      \item A new (sort of) function, \code{\link{tergm.godfather}} has
      been fixed and documented; it can be used to apply a specific set
      of changes to a network, returning statistics of interest as it
      evolves. In particular, it can be used to \dQuote{retrace} the evolution
      of a \code{\link{networkDynamic}}, calculating statistics at
      discrete time points along the way.

      \item \code{\link[=tergm-package]{tergm}} now implements a
      formula-based \code{\link[=summary.networkDynamic]{summary}}
      method for \code{\link{networkDynamic}} LHS, to compute dynamic
      network statistics at specified time points.

      \item Fitting CMLE to network series with transitioned-from
      networks having missing dyads is now possible, using automatic
      imputation. See \code{\link{impute.network.list}} and
      \code{\link{control.tergm}}.

      STERGM \code{\link[=simulate.stergm]{simulate}} can now be used
      for CMLE fits, takes a number of new arguments, and can be used to
      \dQuote{resume} a simulation from a \code{\link{networkDynamic}}
      object.

      \code{\link{gof}} methods have been implemented for CMLE fits.      
      
    }
  }
  
  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{\link{stergm}} EGMME initial fitting code has been
      vastly improved.      
      
      \item The EGMME fitting algorithm has been vastly improved and is
      now a lot more adaptive and able to recover from problems. Many
      bugs have also been fixed.

    }
  }

  % \subsection{C-LEVEL FACILITIES}{
  %   \itemize{
  %     \item It is not clear that we need a subsection called "C-LEVEL FACILITIES".
  %   }
  % }

  % \subsection{UTILITIES}{
  %   \itemize{
  %     \item Likewise, "UTILITIES" is probably unnecessary.
  %   }
  % }

  % \subsection{INSTALLATION}{
  %   \itemize{
  %     \item ... as is "INSTALLATION".
  %   }
  % }

  % \subsection{PACKAGE INSTALLATION}{
  %   \itemize{
  %     \item ditto.
  %   }
  % }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item A number of \code{\link{control.stergm}} parameters had been
      renamed and otherwise chaged.

      \item Argument \code{statsonly} to the
      \code{\link{simulate.stergm}} and related functions has been
      deprecated in favor of \code{output}.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      
      \item The conditional MPLE (CMPLE) for formation and dissolution
      is now fit correctly. This also means that the starting values for
      the CMLE are much better.

      \item Bugs in EGMME code related to handling of bipartite networks
      have been fixed. 
    }
  }
}
