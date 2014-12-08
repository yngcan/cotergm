# TODO: Add comment
# 
# Author: likl
###############################################################################



## ------------------------------------------------------------------------

fls <-list.files("../../",recursive=T)
source(paste0("../../",fls[grep("kk_utils/utils.R",fls)]))
source(paste0("../../",fls[grep("cotergm/R/functions.R",fls)]))
fig.path = "figures/"
# .libPaths(cotergm)
library("statnet.common",lib=cotergm)
library("network",lib=cotergm)
library("ergm",lib=cotergm)
library("tergm",lib=cotergm)
set.seed(1)





## ------------------------------------------------------------------------
n <- 10
y0 <- network(n,density=0.2,directed=FALSE)
y0 %v% "status" <- c(rep(c(2,1),each=n/2))
#"+" is 2, "-" is 1
y0%v%"status"


## ------------------------------------------------------------------------
nw <-  network(n,density=0.2,directed=FALSE)
nw %v% "status" <- c(sample(c(2,1),n,replace=T))



## ------------------------------------------------------------------------
# Formation Plus Network
FP <- FP.nw(y0,nw,"status")
# Formation Minus Network
FM <- FM.nw(y0,nw,"status")
# Dissolution Plus Network
DP <- DP.nw(y0,nw,"status")
# Dissolution Minus Network
DM <- DM.nw(y0,nw,"status")

## ----, warning=FALSE-----------------------------------------------------
fit2 <- ergm(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0)
fit2


mcmc.diagnostics(fit2,vars.per.page=5)
