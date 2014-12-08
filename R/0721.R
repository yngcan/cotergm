
# TODO: Add comment
# 
# Author: kirk
###############################################################################
source("/home/kirk/Dropbox/hiv_network/workspace/kk_utils/utils.R", echo=FALSE, encoding="UTF-8")

.libPaths(cotergm_path_linux)
library("network")
library("ergm",lib=cotergm_path_linux)
library("tergm",lib=cotergm_path_linux)
#?ergm
#?stergm

set.seed(1234)
control <- control.ergm()

y0 <- network.initialize(10,directed=FALSE)

y0 <- san(y0~edges,target.stats=4)

y0 %v% "status" <- c(rep(c("+","-"),each=5))

summary(y0~edges)

nw <-  y0

#as.edgelist(y0)

#debug(san.formula)


#nw <- san(nw~edges,target.stats=10,constraints=~atleast(y0))
#
nw <- san(nw~edges,target.stats=10,constraints=~atleastnonminus(y0))

#par(mfrow=c(1,2))
#p1 <-  plot(y0,label=get.vertex.attribute(y0,"status"),label.cex=2)
#plot(nw,label=get.vertex.attribute(nw,"status"),coord=p1,label.cex=2)

edge.add <- as.edgelist(nw)[setdiff.mat(as.edgelist(y0),as.edgelist(nw)),]

edge.add

nw <- y0

nw <- san(nw~edges,target.stats=10,constraints=~atleastnonplus(y0))

debug(san.formula)

#par(mfrow=c(1,2))
#p1 <-  plot(y0,label=get.vertex.attribute(y0,"status"),label.cex=2)
#plot(nw,label=get.vertex.attribute(nw,"status"),coord=p1,label.cex=2)


edge.add <- as.edgelist(nw)[setdiff.mat(as.edgelist(y0),as.edgelist(nw)),]

edge.add




nw <-y0

nw <- san(nw~edges,target.stats=2,constraints=~atmostnonminus(y0))


as.edgelist(nw)

edge.delete <- as.edgelist(y0)[setdiff.mat(as.edgelist(nw),as.edgelist(y0)),]

edge.delete



nw <-y0

nw <- san(nw~edges,target.stats=2,constraints=~atmostnonplus(y0))


as.edgelist(nw)

edge.delete <- as.edgelist(y0)[setdiff.mat(as.edgelist(nw),as.edgelist(y0)),]

edge.delete


nwsan <- san(nw~nodematch("status",diff=FALSE),target.stats=c(5))


summary(nwsan~nodematch("status",diff=FALSE))

#debug(InitErgmTerm.nodematch)





#
#MHproposal$arguments$constraints$atleastnonminus$free.dyads()
#
#as.edgelist(nw)
#
#setdiff.mat(as.edgelist(y0),as.edgelist(nw))
#
#ind <- intersect.mat(get.MM.dyads(nw,is.edgelist=T),as.edgelist(nw))
#
#
#delete.edges(nw,ind)
#summary(nw~edges)
#
#ergm(nw~edges,constraints=~atleastnonminus(y0))
#
#
#nw <- san(nw~edges,target.stats=10,constraints=~atleastnonminus(y0))
#
#
#
#
#
#
#as.edgelist(y0)
#as.edgelist(nw)
#get.MM.dyads(nw,is.edgelist=T)
#
#
#nw <- san(nw~edges,target.stats=18,constraints=~atleast(y0))
#
#
#
#
#summary(nw~edges)
#
#
#
#
#
#ergm(nw~edges,constraint=~atleast(y0))
#
#log(14/41/(1-14/41))
#
#nw %v% "status" <- sample(c("+","-"),10,rep=TRUE)
#
#as.edgelist(nw)
#
#
#get.MM.dyads(nw)
#
#get.MM.dyads(nw,is.edgelist=T)
#
#debug(InitMHP.formationPlusMLE)
#
#MHproposal <- MHproposal(~atleastnonminus(y0), weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class="c",reference=~Bernoulli,response=NULL)
#
#ergm(nw~edges,constraint=atleastnonminus(y0))
#
#
#ergm(nw~edges)
#debug(ergm.pl)
#
#data(samplk)
## Fit a transition from Time 1 to Time 2
#samplk12 <- stergm(list(samplk1, samplk2),
#		formation=~edges+mutual+transitiveties+cyclicalties,
#		dissolution=~edges+mutual+transitiveties+cyclicalties,
#		estimate="CMLE")
#
#debug(InitMHP.formationMLETNT)
#
#
#
#
#
#z <- .C("SAN_wrapper", as.integer(length(nedges)), as.integer(nedges), 
#		as.integer(tails), as.integer(heads), as.integer(Clist$n), 
#		as.integer(Clist$dir), as.integer(Clist$bipartite), as.integer(Clist$nterms), 
#		as.character(Clist$fnamestring), as.character(Clist$snamestring), 
#		as.character(MHproposal$name), as.character(MHproposal$pkgname), 
#		as.double(c(Clist$inputs, MHproposal$inputs)), as.double(.deinf(eta0)), 
#		as.double(.deinf(tau)), as.integer(1), s = as.double(stats), 
#		as.integer(if (i == 1 | !sequential) control$SAN.burnin else control$SAN.interval), 
#		as.integer(control$SAN.interval), newnwtails = integer(maxedges), 
#		newnwheads = integer(maxedges), as.double(control$invcov), 
#		as.integer(verb), as.integer(MHproposal$arguments$constraints$bd$attribs), 
#		as.integer(MHproposal$arguments$constraints$bd$maxout), 
#		as.integer(MHproposal$arguments$constraints$bd$maxin), 
#		as.integer(MHproposal$arguments$constraints$bd$minout), 
#		as.integer(MHproposal$arguments$constraints$bd$minin), 
#		as.integer(MHproposal$arguments$constraints$bd$condAllDegExact), 
#		as.integer(length(MHproposal$arguments$constraints$bd$attribs)), 
#		
#		
#		
		