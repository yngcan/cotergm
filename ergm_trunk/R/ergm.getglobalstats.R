###########################################################################
# The <ergm.getglobalstats> function calculates and returns the global
# statistics for a given network via <network_stats_wrapper.C>
#
# --PARAMETERS--
#   nw:  a network object
#   m :  the model in use with network nw, as returned by <ergm.getmodel>
#
#
# --RETURNED--
#   gs:  a vector of the global statistics
#
#############################################################################

ergm.getglobalstats <- function(nw, m, response=NULL,MHproposal=NULL) {
  Clist <- ergm.Cprepare(nw, m, response=response)
  # Adjust to global values. This needs to happen before the C call,
  # so that an s_function, if exists could override.
                                                                
  # New method:  Use $emptynwstats added to m$terms by the InitErgmTerm function
  # Read the comments at the top of InitErgm.R or InitErgmTerm.R for 
  # an explanation of the $emptynwstats mechanism
  gs <- rep(0, Clist$nstats)
  i <- 1
  for (j in 1:length(m$terms)) {
    tmp <- m$terms[[j]]
    k <- tmp$inputs[2] # Number of statistics for this model term
    if (!is.null(tmp$emptynwstats)) {
      gs[i:(i+k-1)] <- gs[i:(i+k-1)] + tmp$emptynwstats
    }
    i <- i + k
  }

  
  
#  if(is.null(MHproposal))
  #######
	if(!is.null(MHproposal)){
  y0nw <- MHproposal$arguments$constraints[[1]]$nw
  y0e <- as.edgelist(y0nw) 
  y0nedges<-dim(y0e)[1]
  y0tails<-y0e[,1]
  y0heads<-y0e[,2]
  y0nodalstatus <- rep(0,Clist$n)
  if(!any(is.null(get.vertex.attribute(y0nw,"status"))) & length(get.vertex.attribute(y0nw,"status"))==Clist$n){
	  if(all(get.vertex.attribute(y0nw,"status") %in% c("-","+")))
		  y0nodalstatus <- match(get.vertex.attribute(y0nw,"status"),c("-","+")) 
	  else if (all(get.vertex.attribute(y0nw,"status") %in% c(1,2))){
		  y0nodalstatus <- get.vertex.attribute(y0nw,"status")
	  }}
}
  #######
  
  
  # Note that the empty network statistics are passed to the C
  # code. The reason is that if an s_??? function exists, it can
  # overwrite them, since it can compute the whole thing, while if
  # only the d_??? function exists, it needs to add on to empty
  # network statistics.
  
  # *** don't forget, tails are passes in first now, notheads  
  gs <- if(is.null(response))
         .C("network_stats_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads), as.integer(!is.null(Clist$time)), as.integer(Clist$time), as.integer(NVL(Clist$lasttoggle,0)),
            as.integer(Clist$nedges),
            as.integer(Clist$n),
			as.integer(y0tails),as.integer(y0heads),as.integer(Clist$n),as.integer(y0nedges),as.double(y0nodalstatus),
            as.integer(Clist$dir), as.integer(Clist$bipartite), 
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring), 
            as.double(Clist$inputs),
			as.double(Clist$nodalstatus),
            gs = as.double(gs),
            PACKAGE="ergm"
            )$gs
         else
         .C("wt_network_stats_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads), as.double(Clist$weights), as.integer(!is.null(Clist$time)), as.integer(Clist$time), as.integer(NVL(Clist$lasttoggle,0)),
            as.integer(Clist$nedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite), 
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring), 
            as.double(Clist$inputs),
			as.double(Clist$nodalstatus),
            gs = as.double(gs),
            PACKAGE="ergm"
            )$gs
  names(gs) <- m$coef.names

  gs
}


