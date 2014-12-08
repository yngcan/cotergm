

#################################

InitConstraint.atleastnonminus<-function(conlist, lhs.nw, nw=NULL, ...){
	if(is.null(nw)) stop("Formation constraint ``nonMinus'' requires a baseline network.",call.=FALSE)
	if(network.naedgecount(nw)) stop("Baseline network passed to formation constraint ``nonMinus'' may not have missing dyads.")
	conlist$atleastnonminus<-list(nw=nw)
	conlist$atleastnonminus$free.dyads <- function(){
		standardize.network(invert.network(nw)) & standardize.network(invert.network(get.MM.dyads(nw)))
	}
	
	conlist
}


InitConstraint.atmostnonminus<-function(conlist, lhs.nw, nw=NULL, ...){
	if(is.null(nw)) stop("Formation constraint ``nonMinus'' requires a baseline network.",call.=FALSE)
	if(network.naedgecount(nw)) stop("Baseline network passed to formation constraint ``nonMinus'' may not have missing dyads.")
	conlist$atmostnonminus<-list(nw=nw)
	conlist$atmostnonminus$free.dyads <- function(){
		standardize.network(nw) & standardize.network(invert.network(get.MM.dyads(nw)))
	}
	
	conlist
}


InitConstraint.atleastnonplus<-function(conlist, lhs.nw, nw=NULL, ...){
	if(is.null(nw)) stop("Formation constraint ``nonPlus'' requires a baseline network.",call.=FALSE)
	if(network.naedgecount(nw)) stop("Baseline network passed to formation constraint ``nonPlus'' may not have missing dyads.")
	conlist$atleastnonplus<-list(nw=nw)
	conlist$atleastnonplus$free.dyads <- function(){
		standardize.network(invert.network(nw)) & standardize.network(invert.network(get.PP.dyads(nw)))
	}
	
	conlist
}


InitConstraint.atmostnonplus<-function(conlist, lhs.nw, nw=NULL, ...){
	if(is.null(nw)) stop("Formation constraint ``nonPlus'' requires a baseline network.",call.=FALSE)
	if(network.naedgecount(nw)) stop("Baseline network passed to formation constraint ``nonPlus'' may not have missing dyads.")
	conlist$atmostnonplus<-list(nw=nw)
	conlist$atmostnonplus$free.dyads <- function(){
		standardize.network(nw) & standardize.network(invert.network(get.PP.dyads(nw)))
	}
	
	conlist
}