# TODO: Add comment
# 
# Author: ki
###############################################################################
source("/home/kirk/Dropbox/hiv_network/workspace/kk_utils/utils.R", echo=FALSE, encoding="UTF-8")
source("/home/kirk/Dropbox/hiv_network/workspace/cotergm/R/functions.R",echo=FALSE, encoding="UTF-8")

fig.path = "figures/"
.libPaths(cotergm_path_linux)
library("network",lib=cotergm_path_linux)
library("ergm",lib=cotergm_path_linux)
library("tergm",lib=cotergm_path_linux)

#.libPaths(cran_path_linux)
#library("network")
#library("ergm",lib=cran_path_linux)
#library("tergm",lib=cran_path_linux)

#?ergm
#?stergm

?control.ergm

n <- 20
y0 <- network.initialize(n,directed=FALSE)

y0 %v% "status" <- c(rep(c(2,1),each=n/2))

#"+" is 2, "-" is 1
y0 <- san(y0 ~ edges,target.stats=choose(n,2)*0.2)


formula <- nw~edges+plus+nodematch_cotergm("status",diff=TRUE)
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
#summary(y0~edges)

#nw <- san(nw~edges,target.stats=10,constraints=~atleast(y0))

#debug(san.formula)

nw <- network.initialize(n,directed=FALSE)

nw %v% "status" <- c(sample(c(2,1),n,replace=T))

nw <- san(nw~edges,target.stats=choose(n,2)*0.2)


#plot.network(nw,label=get.vertex.attribute(nw,"status"))


FP <- FP.nw(y0,nw,"status")

FM <- FM.nw(y0,nw,"status")

DP <- DP.nw(y0,nw,"status")

DM <- DM.nw(y0,nw,"status")

#as.edgelist(y0)
#as.edgelist(DM)
#

#
#pdf(paste0(fig.path,"0922.pdf"),height=10,width=10)
#layout_mat <- matrix(c(7,2,8,1,3,6,1,4,6,0,5,0),4,3,byrow=T)
#par(mar=c(2,2,3,2))
#layout(layout_mat)
##
#p1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,main="y0",cex.main=2,displaylabels=T,cex.lab=1)
#
#
#p1.1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="F+",cex.main=2,displaylabels=T,cex.lab=1)
#plot(FP,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=3)
##
##
#p2 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="F-",cex.main=2,displaylabels=T,cex.lab=1)
#plot(FM,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=3)
#
#p3 <- plot(y0,edge.col="yellow",edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="D+",cex.main=2,displaylabels=T,cex.lab=1)
#plot(DP,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=1)
#
#p4 <- plot(y0,edge.col="yellow",edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="D-",cex.main=2,displaylabels=T,cex.lab=1)
#plot(DM,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=1)
#
#p5 <- plot(nw,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="y1",cex.main=2,displaylabels=T,cex.lab=1,vertex.sides=ifelse(y0 %v% "status" == nw %v% "status",50,3))
##
#
#p7 <- plot.new()
#legend("topleft",legend=c("+","-"),pch=19,col=c(2,1),cex=2,bty="n")
#
#p8 <- plot.new()
#legend("topleft",legend=c("tie formed","tie dissolved","nodal status changed"),lty=c(1,1,NA),pch=c(NA,NA,24), col=c(3,"yellow",1),cex=2,bty="n")
#
#
#dev.off()




summary(nw~edges+plus+nodematch_cotergm("status",diff=TRUE)+gwesp(-0.2,fixed=T))

summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE))

summary(y0~edges+plus+nodematch_cotergm("status",diff=TRUE))


table(FP%v%"status")

table(y0%v%"status")

#debug(ergm.mcmcslave)
#debug(dorun)

#debug(ergm)

plot(FP)

summary(FP~edges+gwesp_cotergm(-0.5,fixed=T))


options(warn=2)



fit1 <- ergm(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100,MCMC.runtime.traceplot=T), constraints=~atleastnonminus(y0),verbose=0)

mcmc.diagnostics(fit1, vars.per.page=4)



fit1





