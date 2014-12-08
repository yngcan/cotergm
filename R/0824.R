# TODO: Add comment
# 
# Author: kirk
###############################################################################

source("/home/kirk/Dropbox/hiv_network/workspace/kk_utils/utils.R", echo=FALSE, encoding="UTF-8")
source("/home/kirk/Dropbox/hiv_network/workspace/cotergm/R/functions.R",echo=FALSE, encoding="UTF-8")


.libPaths(cotergm_path_linux)
library("network")
library("ergm",lib=cotergm_path_linux)
library("tergm",lib=cotergm_path_linux)

#
#.libPaths(cran_path_linux)
#library("network")
#library("ergm",lib=cran_path_linux)
#library("tergm",lib=cran_path_linux)
#

#?ergm
#?stergm

set.seed(1234)


n <- 20
y0 <- network.initialize(n,directed=FALSE)

y0 %v% "status" <- c(rep(c("+","-"),each=n/2))


#"+" is 2, "-" is 1
add.edges(y0,c(1,2,3,4,1,2,3,4,5,6,7,8),c(2,3,4,5,6,7,8,9,6,7,8,9))

formula <- nw~nodematch_cotergm("status",diff=FALSE)
#y0 <- (y0~edges,target.stats=4)

#	Rdev -d gdb -f  /home/kirk/Dropbox/hiv_network/workspace/cotergm/R/0824.R
#
#	dir /home/kirk/Dropbox/hiv_network/workspace/cotergm/ergm_trunk/src
#   dir /home/kirk/Dropbox/hiv_network/workspace/cotergm/tergm_trunk/src
##	
#	cd  /home/kirk/Dropbox/hiv_network/workspace/cotergm 
#	Rdev CMD build ergm_trunk
#	Rdev CMD INSTALL -l /home/kirk/Project/R-3.1.1/library_cotergm/ ergm_3.2-12576.1-2014.02.06-00.25.15.tar.gz 
#
#	Rdev CMD build tergm_trunk
#	Rdev CMD INSTALL -l /home/kirk/Project/R-3.1.1/library_cotergm/ tergm_3.2-12531.1-2014.01.24-11.40.49.tar.gz
#	
##
summary(y0~edges)

#nw <- san(nw~edges,target.stats=10,constraints=~atleast(y0))

#debug(san.formula)

nw <- network.initialize(10,directed=FALSE)
nw %v% "status" <- c(rep(c("+","-"),each=5))


nw <- san(nw~nodematch_cotergm("status",diff=FALSE),target.stats=10)

as.edgelist(nw)

plot.network(nw,label=get.vertex.attribute(nw,"status"))



summary(formula)

#debug(ergm)



FP <- FP.nw(y0,nw,"status")

FM <- FM.nw(y0,nw,"status")

DP <- DP.nw(y0,nw,"status")

DM <- DM.nw(y0,nw,"status")

#as.edgelist(y0)
#as.edgelist(DM)
#

#vertex.col=as.factor(get.vertex.attribute(nw,"status"))
#
#layout_mat <- matrix(c(1,2,6,1,3,6,1,4,6,1,5,6),4,3,byrow=T)
#par(mar=c(1,1,3,1))
#layout(layout_mat)
#
#p1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col=vertex.col,label.cex=2,main="y0",cex.main=2)
#
#p1.1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col=vertex.col,label.cex=2,coord=p1,main="F+",cex.main=2)
#plot(FP,coord=p1,new=FALSE,vertex.cex=2,vertex.col=vertex.col,label.cex=2,edge.col=3)
#
#
#p2 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col=vertex.col,label.cex=2,coord=p1,main="F-",cex.main=2)
#plot(FM,coord=p1,new=FALSE,vertex.cex=2,vertex.col=vertex.col,label.cex=2,edge.col=3)
#
#p3 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col=vertex.col,label.cex=2,coord=p1,main="D+",cex.main=2)
#plot(DP,coord=p1,new=FALSE,vertex.cex=2,vertex.col=vertex.col,label.cex=2,edge.col=4)
#
#p4 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col=vertex.col,label.cex=2,coord=p1,main="D-",cex.main=2)
#plot(DM,coord=p1,new=FALSE,vertex.cex=2,vertex.col=vertex.col,label.cex=2,edge.col=4)
#
#p5 <- plot(nw,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col=vertex.col,label.cex=2,coord=p1,main="y1",cex.main=2)


debug(ergm.MCMLE)

as.edgelist(FP)

summary(FP~nodematch_cotergm("status",diff=TRUE))
#as.edgelist(FP)

fit1 <- ergm(FP~nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=0,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=30), constraints=~atleastnonminus(y0),verbose=0)
names(fit1)



mcmc.diagnostics(fit1)


#
#names(fit1)
#
#
#fit1$sample
#
#fit0 <- ergm(network(10,directed=FALSE,density=0.2)~edges+degree(1))
#
#debug(mcmc.diagnostics)
#mcmc.diagnostics(fit0)





#
#fit2 <- ergm(FM~nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=0,MCMC.samplesize=100,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=10), constraints=~atleastnonplus(y0),verbose=5)
#
#
#
#summary(FP~nodematch_cotergm("status",diff=TRUE))
#
#
#
#
#
#
#
#
#
#
#sum(apply(as.edgelist(FP),1,function(x)x[1]<=5 & x[2]<=5))
#
#
#apply(as.edgelist(FP),1,function(x)x[1]>5 & x[2]>5)
#
#
#fit1
#


#
#
#fit2 <- ergm(FM~edges+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=0,MCMC.samplesize=100,MCMC.interval=10,MCMLE.maxit=20,MCMLE.min.effectiveSize=3), constraints=~atleastnonplus(y0),verbose=F)
#
#fit2
#
#
#
#
#fit3 <- ergm(DP~edges+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=0,MCMC.samplesize=100,MCMC.interval=10,MCMLE.maxit=20,MCMLE.min.effectiveSize=3), constraints=~atmostnonminus(y0),verbose=F)
#
#fit3
#
#
#fit4 <- ergm(DM~edges+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=0,MCMC.samplesize=100,MCMC.interval=10,MCMLE.maxit=20,MCMLE.min.effectiveSize=3), constraints=~atmostnonplus(y0),verbose=F)
#
#fit4
#
#
#
#
#
#
#
#
#
#
#?simulate.ergm
#
#S1 <- simulate(formula,coef=0, constraints=~atleastnonminus(y0), nsim=1000)
#S2 <- simulate(formula,coef=0, constraints=~atleastnonplus(y0), nsim=1000)
#S3 <- simulate(formula,coef=0, constraints=~atmostnonminus(y0), nsim=1000)
#S4 <- simulate(formula,coef=0, constraints=~atmostnonplus(y0), nsim=1000)
#
#
#L1 <- lapply(S1,as.edgelist)
#L2 <- lapply(S2,as.edgelist)
#L3 <- lapply(S3,as.edgelist)
#L4 <- lapply(S4,as.edgelist)
#
#
#plot(nw)
#
#
#
#net1
#
#
#net2

#
#fit1 <- ergm(nw~nodematch_cotergm("status",diff=FALSE), constraints=~atleastnonminus(y0),verbose=6)
#




#
#fit1 <- ergm(nw~degree(1),constraints=~atleast(y0),control=control.ergm(MCMC.burnin=0),verbose=F)


#
#nw <- san(nw~nodematch("status",diff=FALSE),target.stats=5,constraints=~atleast(y0),verbose=5)


#debug(san.formula)
#
#
#debug(ergm.Cprepare)
#
#
#debug(ergm.getmodel)
#debug(san.formula)


#
##par(mfrow=c(1,2))
##p1 <-  plot(y0,label=get.vertex.attribute(y0,"status"),label.cex=2)
##plot(nw,label=get.vertex.attribute(nw,"status"),coord=p1,label.cex=2)
#
#edge.add <- as.edgelist(nw)[setdiff.mat(as.edgelist(y0),as.edgelist(nw)),]
#
#edge.add
#
#nw <- y0
#
#nw <- san(nw~edges,target.stats=10,constraints=~atleastnonplus(y0))
#
#debug(san.formula)
#
##par(mfrow=c(1,2))
##p1 <-  plot(y0,label=get.vertex.attribute(y0,"status"),label.cex=2)
##plot(nw,label=get.vertex.attribute(nw,"status"),coord=p1,label.cex=2)
#
#
#edge.add <- as.edgelist(nw)[setdiff.mat(as.edgelist(y0),as.edgelist(nw)),]
#
#edge.add
#
#
#
#
#nw <-y0
#
#nw <- san(nw~edges,target.stats=2,constraints=~atmostnonminus(y0))
#
#
#as.edgelist(nw)
#
#edge.delete <- as.edgelist(y0)[setdiff.mat(as.edgelist(nw),as.edgelist(y0)),]
#
#edge.delete
#
#
#
#nw <-y0
#
#nw <- san(nw~edges,target.stats=2,constraints=~atmostnonplus(y0))
#
#
#as.edgelist(nw)
#
#edge.delete <- as.edgelist(y0)[setdiff.mat(as.edgelist(nw),as.edgelist(y0)),]
#
#edge.delete
#
#
#nwsan <- san(nw~nodematch("status",diff=FALSE),target.stats=c(5))
#
#
#summary(nwsan~nodematch("status",diff=FALSE))

#debug(InitErgmTerm.nodematch)

