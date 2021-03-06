#==============================================================================
# This file contains the 2 following functions for getting an MCMC sample
#      <ergm.getMCMCsample>
#      <ergm.mcmcslave>
#==============================================================================




#########################################################################################
# The <ergm.getMCMCsample> function samples networks using an MCMC algorithm via
# <MCMC_wrapper.C>. Unlike its <ergm.getMCMCsample> counterpart, this function is
# caple of running in multiple threads.  Note that the returned stats will be relative to
# the original network, i.e., the calling function must shift the statistics if required. 
# The calling function must also attach column names to the statistics matrix if required.
#
# --PARAMETERS--
#   nw        :  a network object
#   model     :  a model for the given 'nw' as returned by <ergm.getmodel>
#   MHproposal:  a list of the parameters needed for Metropolis-Hastings proposals and
#                the result of calling <MHproposal>
#   eta0      :  the initial eta coefficients
#   verbose   :  whether the C functions should be verbose; default=FALSE
#   control:  list of MCMC tuning parameters; those recognized include
#       parallel    : the number of threads in which to run the sampling
#       packagenames: names of packages; this is only relevant if "ergm" is given
#       samplesize  : the number of networks to be sampled
#       interval    : the number of proposals to ignore between sampled networks
#       burnin      : the number of proposals to initially ignore for the burn-in
#                     period
#
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal, control, verbose should certainly
# be part of Clist.  But this is a project for another day!
#
# --RETURNED--
#   the sample as a list containing:
#     statsmatrix:  the stats matrix for the sampled networks, RELATIVE TO THE ORIGINAL
#                   NETWORK!
#     newnetwork :  the edgelist of the final sampled network
#     nedges     :  the number of edges in the 'newnetwork' ??
#
#########################################################################################

ergm.getMCMCsample <- function(nw, model, MHproposal, eta0, control, 
                                        verbose, response=NULL, ...) {
  
  if(is.network(nw[[1]])){ # I.e., we are dealing with a list of initial networks.
    nnw <- length(nw)
    # FIXME: This doesn't have to be the case:
    if(control$parallel!=nnw)
      stop("Number of initial networks passed to ergm.getMCMCsample must equal control$parallel (for now).")
  }else nnw <- 1
     
  Clist <- if(nnw>1) lapply(nw, ergm.Cprepare, model,
response=response) else ergm.Cprepare(nw, model, response=response)

  if(control$parallel==0){
    flush.console()
    z <- ergm.mcmcslave(Clist,MHproposal,eta0,control,verbose,...)
    
	
    if(z$status == 1){ # MCMC_TOO_MANY_EDGES, exceeding even control$MCMC.max.maxedges
      return(list(status=z$status))
    }
    
    if(z$status == 2){ # MCMC_MH_FAILED
      # MH proposal failed somewhere. Throw an error.
      stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
    }
    
    if(!is.null(z$burnin.failed) && z$burnin.failed) warning("Burn-in failed to converge after retries.")
    
    statsmatrix <- matrix(z$s, nrow=control$MCMC.samplesize,
                          ncol=Clist$nstats,
                          byrow = TRUE)
    newnetwork <- newnw.extract(nw,z,response=response,
                                output=control$network.output)
    newnetworks <- list(newnetwork)

    burnin.total <- z$burnin.total
  }else{
    control.parallel <- control
    control.parallel$MCMC.samplesize <- ceiling(control$MCMC.samplesize / control$parallel)
    
    cl <- ergm.getCluster(control, verbose)
    #
    #   Run the jobs with rpvm or Rmpi
    #
    flush.console()
    outlist <- {
      if(nnw>1) clusterMap(cl,ergm.mcmcslave,
                           Clist, MoreArgs=list(MHproposal,eta0,control.parallel,verbose,...))
      else clusterCall(cl,ergm.mcmcslave,
                       Clist,MHproposal,eta0,control.parallel,verbose,...)
    }
    #
    #   Process the results
    #
    statsmatrix <- NULL
    newnetworks <- list()
    burnin.total <- c()
    for(i in (1:control$parallel)){
      z <- outlist[[i]]

      if(z$status == 1){ # MCMC_TOO_MANY_EDGES, exceeding even control$MCMC.max.maxedges
        return(list(status=z$status))
      }
      
      if(z$status == 2){ # MCMC_MH_FAILED
        # MH proposal failed somewhere. Throw an error.
        stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
      }
      
      if(!is.null(z$burnin.failed) && z$burnin.failed) warning("Burn-in failed to converge after retries.")
      
      statsmatrix <- rbind(statsmatrix,
                           matrix(z$s, nrow=control.parallel$MCMC.samplesize,
                                  ncol=(if(nnw==1) Clist else Clist[[i]])$nstats,
                                  byrow = TRUE))
      newnetworks[[i]]<-newnw.extract(if(nnw==1) nw else nw[[i]],z,
                                  response=response,output=control$network.output)
      burnin.total <- c(burnin.total, z$burnin.total)
    }
    newnetwork<-newnetworks[[1]]

    if(verbose){cat("parallel samplesize=",nrow(statsmatrix),"by",
                    control.parallel$MCMC.samplesize,"\n")}
    
    ergm.stopCluster(cl)
  }

  if(verbose){cat("new nodalstatus = ",z$newnodalstatus,"\n")}
  
  colnames(statsmatrix) <- model$coef.names

  statsmatrix[is.na(statsmatrix)] <- 0
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, newnetworks=newnetworks, status=0, burnin.total=burnin.total)
}





###############################################################################
# The <ergm.mcmcslave> function is that which the slaves will call to perform
# a validation on the mcmc equal to their slave number. It also returns an
# MCMC sample.
#
# --PARAMETERS--
#   Clist     : the list of parameters returned by <ergm.Cprepare>
#   MHproposal: the MHproposal list as returned by <getMHproposal>
#   eta0      : the canonical eta parameters
#   control: a list of parameters for controlling the MCMC algorithm;
#               recognized components include:
#       samplesize  : the number of networks to be sampled
#       interval    : the number of proposals to ignore between sampled networks
#       burnin      : the number of proposals to initially ignore for the burn-in
#                     period
#   verbose   : whether the C code should be verbose (T or F)
#   ...       : optional arguments
#
# --RETURNED--
#   the MCMC sample as a list of the following:
#     s         : the statsmatrix
#     newnwtails: the vector of tails for the new network- is this the final
#                 network sampled? - is this the original nw if 'maxedges' is 0
#     newnwheads: the vector of heads for the new network - same q's
#
###############################################################################

ergm.mcmcslave <- function(Clist,MHproposal,eta0,control,verbose,...) {
  # A subroutine to allow caller to override some settings or resume
  # from a pervious run.
  dorun <- function(prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL, maxedges=NULL,maxnstatus=NULL,MHproposal=NULL){
    
    numnetworks <- 0
    if(is.null(prev.run)){ # Start from Clist
      nedges <- c(Clist$nedges,0,0)
      tails <- Clist$tails
      heads <- Clist$heads
	  nodalstatus <- Clist$nodalstatus
      weights <- Clist$weights
      stats <- rep(0, Clist$nstats)
    }else{ # Pick up where we left off
      nedges <- prev.run$newnwtails[1]
      tails <- prev.run$newnwtails[2:(nedges+1)]
      heads <- prev.run$newnwheads[2:(nedges+1)]
	  nodalstatus <- prev.run$newnodalstatus
      weights <- prev.run$newnwweights[2:(nedges+1)]
      nedges <- c(nedges,0,0)
      stats <- matrix(prev.run$s,
                      ncol=Clist$nstats,
                      byrow = TRUE)
      stats <- stats[nrow(stats),]
    }
 
    if(is.null(burnin)) burnin <- control$MCMC.burnin
    if(is.null(samplesize)) samplesize <- control$MCMC.samplesize
    if(is.null(interval)) interval <- control$MCMC.interval
    if(is.null(maxedges)) maxedges <- control$MCMC.init.maxedges

	if(all(is.na(eta0))) eta0=rep(0.5,length(eta0))

	eta0[!sapply(eta0,is.finite)]<-sign(eta0[!sapply(eta0,is.finite)])*1000
	
	#######
	y0nw <- MHproposal$arguments$constraints[[1]]$nw
	y0e <- as.edgelist(y0nw) 
	y0nedges<-dim(y0e)[1]
	y0tails<-y0e[,1]
	y0heads<-y0e[,2]
	y0nodalstatus <- rep(0,Clist$n)
	if(!any(is.null(get.vertex.attribute(y0,"status"))) & length(get.vertex.attribute(y0,"status"))==Clist$n){
		if(all(get.vertex.attribute(y0,"status") %in% c("-","+")))
			y0nodalstatus <- match(get.vertex.attribute(y0,"status"),c("-","+")) 
		else if (all(get.vertex.attribute(y0,"status") %in% c(1,2))){
			y0nodalstatus <- get.vertex.attribute(y0,"status")
		}}
	#######
    repeat{
      if(is.null(Clist$weights)){
        z <- .C("MCMC_wrapper",
                as.integer(numnetworks), as.integer(nedges),
                as.integer(tails), as.integer(heads),
                as.integer(Clist$n),
				as.integer(y0tails),
				as.integer(y0heads),
				as.integer(Clist$n),
				as.integer(y0nedges),
				as.double(c(y0nodalstatus)),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms),
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring),
                as.character(MHproposal$name), as.character(MHproposal$pkgname),
                as.double(c(Clist$inputs,MHproposal$inputs)), 
				# nodalstatus updates each run
				as.double(c(nodalstatus)),
				as.double(.deinf(eta0)),
                as.integer(samplesize),
                s = as.double(rep(stats, samplesize)),
                as.integer(burnin), 
                as.integer(interval),
                newnwtails = integer(maxedges),
                newnwheads = integer(maxedges),
				newnodalstatus = integer(Clist$n),
                as.integer(verbose), as.integer(MHproposal$arguments$constraints$bd$attribs),
                as.integer(MHproposal$arguments$constraints$bd$maxout), as.integer(MHproposal$arguments$constraints$bd$maxin),
                as.integer(MHproposal$arguments$constraints$bd$minout), as.integer(MHproposal$arguments$constraints$bd$minin),
                as.integer(MHproposal$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal$arguments$constraints$bd$attribs)),
                as.integer(maxedges),
                status = integer(1),
                PACKAGE="ergm")
        
        # save the results
        z<-list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads, newnodalstatus=z$newnodalstatus, status=z$status, maxedges=maxedges)
		
      }else{
        z <- .C("WtMCMC_wrapper",
                as.integer(length(nedges)), as.integer(nedges),
                as.integer(tails), as.integer(heads), as.double(weights),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms),
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring),
                as.character(MHproposal$name), as.character(MHproposal$pkgname),
                as.double(c(Clist$inputs,MHproposal$inputs)), 
				as.double(c(Clist$nodalstatus)),
				as.double(.deinf(eta0)),
                as.integer(samplesize),
                s = as.double(rep(stats, samplesize)),
                as.integer(burnin), 
                as.integer(interval),
                newnwtails = integer(maxedges),
                newnwheads = integer(maxedges),
                newnwweights = double(maxedges),
                as.integer(verbose), 
                as.integer(maxedges),
                status = integer(1),
                PACKAGE="ergm") 
        # save the results
        z<-list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads, newnwweights=z$newnwweights, status=z$status, maxedges=maxedges)
      }
      if(z$status!=1) return(z) # Handle everything except for MCMC_TOO_MANY_EDGES elsewhere.

      # The following is only executed (and the loop continued) if too many edges.
      maxedges <- maxedges * 10
      if(!is.null(control$MCMC.max.maxedges)){
        if(maxedges == control$MCMC.max.maxedges*10) # True iff the previous maxedges exactly equaled control$MCMC.max.maxedges and that was too small.
          return(z) # This will kick the too many edges problem upstairs, so to speak.
        maxedges <- min(maxedges, control$MCMC.max.maxedges)
      }

    }
  }
  
  if(control$MCMC.burnin>0 && NVL(control$MCMC.burnin.retries,0)>0){
    out <- NULL
    burnin.stats <- NULL
    burnin.total <- 0
    samplesize <- min(control$MCMC.samplesize,control$MCMC.burnin)
    burnin <- 0
    interval <- ceiling(control$MCMC.burnin/samplesize)
    
    for(try in seq_len(control$MCMC.burnin.retries+1)){
      out<-dorun(prev.run=out,
                 burnin = burnin,
                 samplesize = samplesize,
                 interval = interval,
                 maxedges = out$maxedges, # note that the first time through, maxedges=NULL$maxedges, which is NULL.
		 MHproposal=MHproposal)
      # Stop if something went wrong.
      if(out$status!=0) return(out)

      burnin.total <- burnin.total + burnin + (samplesize-1)*interval
      
      # Get the array of the burnin draws. Note that all draws get stored.
      burnin.stats <- rbind(burnin.stats,
                            matrix(out$s, nrow=samplesize,
                                   ncol=Clist$nstats,
                                   byrow = TRUE)
                            )
      colnames(burnin.stats) <- names(Clist$diagnosable)

      if(nrow(burnin.stats)>=samplesize*8){
        burnin.stats <- burnin.stats[seq_len(floor(nrow(burnin.stats)/2))*2,,drop=FALSE]
        interval <- interval*2
        if(verbose) cat("Increasing thinning to",interval,".\n")
      }
      
      # Extract the last draws for diagnostics.
      burnin.stats.last <- burnin.stats[-seq_len((1-control$MCMC.burnin.check.last)*nrow(burnin.stats)),,drop=FALSE]
      
      burnin.esteq.last <-
        if(all(c("theta","etamap") %in% names(list(...)))) .ergm.esteq(list(...)$theta, list(etamap=list(...)$etamap), burnin.stats.last)
        else burnin.stats.last[,Clist$diagnosable,drop=FALSE]

      if(control$MCMC.runtime.traceplot) plot(window(mcmc(burnin.esteq.last,thin=max(1,floor(nrow(burnin.esteq.last)/1000)))),ask=FALSE,smooth=TRUE,density=FALSE)
      
      burnin.test <- suppressWarnings(try(geweke.diag.mv(burnin.esteq.last,.3,.3))) # Note that the first and the last are different. More similar sample sizes give more power, and the gap between the samples is the same as before (0.4).
      if(inherits(burnin.test,"try-error")){
        if(verbose) cat("Burn-in convergence test failed. Rerunning.\n")
        if(try == control$MCMC.burnin.retries+1) burnin.failed <- TRUE
      }else if(burnin.test$parameter["df"]<control$MCMC.burnin.min.df){
        if(verbose) cat("Insufficient burn-in sample size (",burnin.test$parameter["df"],"<",control$MCMC.burnin.min.df,") to test convergence. Rerunning.\n")
        if(try == control$MCMC.burnin.retries+1) burnin.failed <- TRUE
      }else if(burnin.test$p.value < control$MCMC.burnin.check.alpha){
        failed <- effectiveSize(burnin.esteq.last)
        failed <- order(failed)
        if(verbose) cat("Burn-in failed to converge or mixed very poorly, with p-value =",burnin.test$p.value,". Ranking of statistics from worst-mixing to best-mixing:", paste.and(names(Clist$diagnosable[Clist$diagnosable])[failed]), ". Rerunning.\n")
        if(try == control$MCMC.burnin.retries+1) burnin.failed <- TRUE
      }else{
        if(verbose){
          print(burnin.test)
          cat("Burn-in converged. Proceeding to the sampling run.\n")
        }
        burnin.failed <- FALSE
        break
      }
    }

    # Do the actual sampling run. Note that we've already done the burnin.
    out <- dorun(prev.run=out, burnin=0, maxedges=out$maxedges,MHproposal=MHproposal)
    if(control$MCMC.runtime.traceplot) {
      stats <- matrix(out$s, nrow=control$MCMC.samplesize,
                      ncol=Clist$nstats,
                      byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]
      
      plot(mcmc(stats,start=1,end=control$MCMC.samplesize*control$MCMC.interval,thin=control$MCMC.interval),ask=FALSE,smooth=TRUE,density=FALSE)
    }
    out$burnin.failed <- burnin.failed
    out$burnin.total <- burnin.total
    out
  } else {
    out <- dorun(MHproposal=MHproposal)
	out$burnin.total <- 0
    if(control$MCMC.runtime.traceplot) {
      stats <- matrix(out$s, nrow=control$MCMC.samplesize,
                      ncol=Clist$nstats,
                      byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]
      
      plot(mcmc(stats,start=control$MCMC.burnin+1,control$MCMC.burnin+control$MCMC.samplesize*control$MCMC.interval,thin=control$MCMC.interval),ask=FALSE,smooth=TRUE,density=FALSE)
    }
    out
  }
}
