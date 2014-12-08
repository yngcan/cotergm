
## ------------------------------------------------------------------------
source("/home/kirk/Dropbox/hiv_network/workspace/kk_utils/utils.R", echo=FALSE, encoding="UTF-8")
source("/home/kirk/Dropbox/hiv_network/workspace/cotergm/R/functions.R",echo=FALSE, encoding="UTF-8")
set.seed(123)
fig.path = "figures/"
# .libPaths(cotergm_path_linux)
library("statnet.common",lib=cotergm)
library("network",lib=cotergm)
library("ergm",lib=cotergm)
library("tergm",lib=cotergm)



## ------------------------------------------------------------------------
n <- 40
y0 <- network(n,density=0.2,directed=FALSE)
y0 %v% "status" <- c(rep(c(2,1),each=n/2))
#"+" is 2, "-" is 1
y0%v%"status"


## ------------------------------------------------------------------------
nw <-  network(n,density=0.2,directed=FALSE)
nw %v% "status" <- c(sample(c(2,1),n,replace=T))
nw%v%"status"


## ------------------------------------------------------------------------
# Formation Plus Network
FP <- FP.nw(y0,nw,"status")
# Formation Minus Network
FM <- FM.nw(y0,nw,"status")
# Dissolution Plus Network
DP <- DP.nw(y0,nw,"status")
# Dissolution Minus Network
DM <- DM.nw(y0,nw,"status")


## ------------------------------------------------------------------------
# pdf(paste0(fig.path,"0922.pdf"),height=10,width=10)
layout_mat <- matrix(c(7,2,8,1,3,6,1,4,6,0,5,0),4,3,byrow=T)
par(mar=c(2,2,3,2))
layout(layout_mat)
#
p1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,main="y0",cex.main=2,displaylabels=T,cex.lab=1)

p1.1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="F+",cex.main=2,displaylabels=T,cex.lab=1)
plot(FP,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=3)
#
#
p2 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="F-",cex.main=2,displaylabels=T,cex.lab=1)
plot(FM,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=3)

p3 <- plot(y0,edge.col="yellow",edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="D+",cex.main=2,displaylabels=T,cex.lab=1)
plot(DP,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=1)

p4 <- plot(y0,edge.col="yellow",edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="D-",cex.main=2,displaylabels=T,cex.lab=1)
plot(DM,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=1)

p5 <- plot(nw,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="y1",cex.main=2,displaylabels=T,cex.lab=1,vertex.sides=ifelse(y0 %v% "status" == nw %v% "status",50,3))
#
p7 <- plot.new()
legend("topleft",legend=c("+","-"),pch=19,col=c(2,1),cex=2,bty="n")

p8 <- plot.new()
legend("topleft",legend=c("tie formed","tie dissolved","nodal status changed"),lty=c(1,1,NA),pch=c(NA,NA,24), col=c(3,"yellow",1),cex=2,bty="n")
#
#dev.off()


## ------------------------------------------------------------------------
summary(y0~edges+plus+nodematch_cotergm("status",diff=TRUE))
summary(nw~edges+plus+nodematch_cotergm("status",diff=TRUE))
summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE))
summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE))
summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE))
summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE))


## ----, warning=FALSE-----------------------------------------------------
fit1 <- ergm(FP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0)
fit1
mcmc.diagnostics(fit1, vars.per.page=5)


fit1.1 <- ergm(FP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0)
fit1.1
mcmc.diagnostics(fit1.1, vars.per.page=5)



fit1.2 <- ergm(FP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0)
fit1.2
mcmc.diagnostics(fit1.2, vars.per.page=5)


## ----, warning=FALSE-----------------------------------------------------
fit2 <- ergm(FM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0)
fit2
mcmc.diagnostics(fit2, vars.per.page=5)

fit2.1 <- ergm(FM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0)
fit2.1
mcmc.diagnostics(fit2.1, vars.per.page=5)

fit2.2 <- ergm(FM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0)
fit2.2
mcmc.diagnostics(fit2.2, vars.per.page=5)


## ----, warning=FALSE-----------------------------------------------------
fit3 <- ergm(DP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0)
fit3
mcmc.diagnostics(fit3, vars.per.page=5)

fit3.1 <- ergm(DP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0)
fit3.1
mcmc.diagnostics(fit3.1, vars.per.page=5)

fit3.2 <- ergm(DP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0)
fit3.2
mcmc.diagnostics(fit3.2, vars.per.page=5)


## ----, warning=FALSE-----------------------------------------------------
fit4 <- ergm(DM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0)
fit4
mcmc.diagnostics(fit4, vars.per.page=5)

fit4.1 <- ergm(DM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0)
fit4.1
mcmc.diagnostics(fit4.1, vars.per.page=5)

fit4.2 <- ergm(DM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0)
fit4.2
mcmc.diagnostics(fit4.2, vars.per.page=5)


