# TODO: Add comment
# 
# Author: kirk
###############################################################################
InitErgmTerm.plus<-function(nw, arglist, ...) {
	a <- check.ErgmTerm(nw, arglist,
			varnames = NULL,
			vartypes = NULL,
			defaultvalues = list(),
			required = NULL)
	
	list(name="plus", coef.names="plus", dependence=TRUE,
			minval = 0, maxval = network.size(nw))
}




################################################################################
InitErgmTerm.nodematch_cotergm<-function (nw, arglist, ...) {
	### Check the network and arguments to make sure they are appropriate.
	a <- check.ErgmTerm(nw, arglist, 
			varnames = c("attrname", "diff", "keep","y0edge","y0nodal"),
			vartypes = c("character", "logical", "numeric", "numeric","numeric"),
			defaultvalues = list(NULL, FALSE, NULL,0,0),
			required = c(TRUE, FALSE, FALSE,FALSE,FALSE))
	### Process the arguments
	nodecov <-
			if(length(a$attrname)==1)
				get.node.attr(nw, a$attrname)
			else{
				do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
			}
	u <- sort(c(1,2))  # deal with when nodes are all of the same status
	if (!is.null(a$keep)) {
		u <- u[a$keep]
	}
	#   Recode to numeric
#	nodecov <- match(nodecov,u,nomatch=length(u)+1)
	# All of the "nomatch" should be given unique IDs so they never match:
#	dontmatch <- nodecov==(length(u)+1)
#	nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
#	ui <- seq(along=u)
	### Construct the list to return
	
	
	if (a$diff) {
		coef.names <- paste("nodematch", paste(a$attrname,collapse="."), u, sep=".")
		inputs <- c(u,a$y0edge,a$y0nodal,nodecov)
	} else {
#		coef.names <- paste("nodematch", paste(a$attrname,collapse="."), sep=".")
#		inputs <- nodecov
		stop("diff is FALSE is not valid for cotergm")
	}
	
	list(name="nodematch_cotergm",                                 #name: required
			coef.names = coef.names,                          #coef.names: required
			inputs =  inputs,
			dependence = TRUE, # So we don't use MCMC if not necessary
			minval = 0
	)
}





################################################################################
InitErgmTerm.gwesp_cotergm<-function(nw, arglist, initialfit=FALSE, ...) {
# the following line was commented out in <InitErgm.gwesp>:
#   ergm.checkdirected("gwesp", is.directed(nw), requirement=FALSE)
# so, I've not passed 'directed=FALSE' to <check.ErgmTerm>  
	a <- check.ErgmTerm(nw, arglist,
			varnames = c("alpha","fixed","cutoff"),
			vartypes = c("numeric","logical","numeric"),
			defaultvalues = list(0, FALSE, 30),
			required = c(FALSE, FALSE, FALSE))
	alpha<-a$alpha;fixed<-a$fixed
	cutoff<-a$cutoff
	alpha=alpha[1] # Not sure why anyone would enter a vector here, but...
	if(!initialfit && !fixed){ # This is a curved exponential family model
#   d <- 1:(network.size(nw)-2)
		maxesp <- min(cutoff,network.size(nw)-2)
		d <- 1:maxesp
		ld<-length(d)
		if(ld==0){return(NULL)}
		map <- function(x,n,...){
			i <- 1:n
			x[1]*exp(x[2])*(1-(1-exp(-x[2]))^i)
		}
		gradient <- function(x,n,...){
			i <- 1:n
			a <- 1-exp(-x[2])
			exp(x[2]) * rbind(1-a^i, x[1] * (1 - a^i - i*a^(i-1) ) )
		}
		if(is.directed(nw)){dname <- "tesp"}else{dname <- "esp"}
		list(name=dname, coef.names=paste("esp#",d,sep=""), 
				inputs=c(d), params=list(gwesp=NULL,gwesp.alpha=alpha),
				map=map, gradient=gradient)
	}else{
		if (initialfit && !fixed)  # First pass to get MPLE coefficient
			coef.names <- "gwesp"
		else # fixed == TRUE
			coef.names <- paste("gwesp.fixed",alpha,sep="")
		if(is.directed(nw)){dname <- "gwtesp_cotergm"}else{dname <- "gwesp_cotergm"}
		list(name=dname, coef.names=coef.names, inputs=c(alpha))
	}
}


