# TODO: Add comment
# 
# Author: kirk
###############################################################################
## ------------------------------------------------------------------------
#source(paste0(addr,"kk_utils/utils.R"), echo=FALSE, encoding="UTF-8")
#
try(source("/home/kirk/Dropbox/hiv_network/workspace/kk_utils/utils.R", echo=FALSE, encoding="UTF-8"))
try(source("/home/kirk/Dropbox/hiv_network/workspace/cotergm/R/functions.R", echo=FALSE, encoding="UTF-8"))

#
#try(source("/home/local/ANT/likl/Dropbox/hiv_network/workspace/kk_utils/utils.R", echo=FALSE, encoding="UTF-8"))
#try(source("/home/local/ANT/likl/Dropbox/hiv_network/workspace/cotergm/R/functions.R", echo=FALSE, encoding="UTF-8"))

#source(paste0(addr,"cotergm/R/functions.R"),echo=FALSE, encoding="UTF-8")
set.seed(123)
fig.path = "figures/"
# .libPaths(cotergm_path_linux)

#install.packages("ergm",dependencies=TRUE,lib=cotergm)

library("statnet.common",lib=cotergm)
library("network",lib=cotergm)
library("ergm",lib=cotergm)
library("tergm",lib=cotergm)

## ------------------------------------------------------------------------
net1 <-as.network(read.table("./s50_data/s50-network1.dat"))
net2 <-as.network(read.table("./s50_data/s50-network2.dat"))
net3 <-as.network(read.table("./s50_data/s50-network3.dat"))
alcohol <- read.table("./s50_data/s50-alcohol.dat")
drugs <- read.table("./s50_data/s50-drugs.dat")
net <- net1

d2ud <- function(net){
	net <- as.matrix(net)
	for (i in 1:nrow(net)){
		for(j in 1:nrow(net)){
			net[i,j] = max(net[i,j],net[j,i])
		}}
	as.network(net,directed=FALSE)
}

net1 <- d2ud(net1)
net2 <- d2ud(net2)
net3 <- d2ud(net3)

al1 <- as.numeric(alcohol[,1]>3)+1
al2 <- as.numeric(alcohol[,2]>3)+1
net1 %v% "status" <- al1
net2 %v% "status" <- al2

y0 <- net1
nw <- net2

## ------------------------------------------------------------------------
# Formation Plus Network
FP <- FP.nw(y0,nw,"status")
# Formation Minus Network
FM <- FM.nw(y0,nw,"status")
# Dissolution Plus Network
DP <- DP.nw(y0,nw,"status")
# Dissolution Minus Network
DM <- DM.nw(y0,nw,"status")

#debug(ergm.mcmcslave)

## ------------------------------------------------------------------------
## ----, warning=FALSE-----------------------------------------------------
fit1 <- ergm(FP~edges+nodematch_cotergm("status",diff=TRUE,y0edge=2),control=control.ergm(MCMC.burnin=10,MCMC.samplesize=10,MCMC.interval=1,MCMLE.maxit=1,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=1)




#
#fit1 <- ergm(FP~nodematch_cotergm("status",diff=TRUE,y0edge=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0)
#
mcmc.diagnostics.ergm(fit1,vars.per.page=5)
#
#
#



# Rdev -d gdb -f  /home/kirk/Dropbox/hiv_network/workspace/cotergm/R/1027.R

# Rdev CMD build /home/local/ANT/likl/Dropbox/hiv_network/workspace/cotergm/ergm_trunk/ && Rdev CMD INSTALL -l ~/Project/R-3.1.1/library_cotergm/   /home/local/ANT/likl/Dropbox/hiv_network/workspace/cotergm/ergm_3.3.tar.gz 


#  	
#  dir /home/kirk/Dropbox/hiv_network/workspace/cotergm/ergm_trunk/src
#  dir /home/local/ANT/likl/Dropbox/hiv_network/workspace/cotergm/ergm_trunk/src
#Rdev -d "valgrind --tool=callgrind" --vanilla < /home/local/ANT/likl/Dropbox/hiv_network/workspace/cotergm/R/1027.R a.out > log.txt 2>&1



#mcmc.diagnostics(fit1, vars.per.page=5)
#
#
#fit1.1 <- ergm(FP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0)
#fit1.1
#mcmc.diagnostics(fit1.1, vars.per.page=5)
#
#
#
#fit1.2 <- ergm(FP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0)
#fit1.2
#mcmc.diagnostics(fit1.2, vars.per.page=5)
#
#
### ----, warning=FALSE-----------------------------------------------------
#fit2 <- ergm(FM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0)
#fit2
#mcmc.diagnostics(fit2, vars.per.page=5)
#
#fit2.1 <- ergm(FM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0)
#fit2.1
#mcmc.diagnostics(fit2.1, vars.per.page=5)
#
#fit2.2 <- ergm(FM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0)
#fit2.2
#mcmc.diagnostics(fit2.2, vars.per.page=5)
#
#
### ----, warning=FALSE-----------------------------------------------------
#fit3 <- ergm(DP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0)
#fit3
#mcmc.diagnostics(fit3, vars.per.page=5)
#
#fit3.1 <- ergm(DP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0)
#fit3.1
#mcmc.diagnostics(fit3.1, vars.per.page=5)
#
#fit3.2 <- ergm(DP~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0)
#fit3.2
#mcmc.diagnostics(fit3.2, vars.per.page=5)
#
#
### ----, warning=FALSE-----------------------------------------------------
#fit4 <- ergm(DM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0)
#fit4
#mcmc.diagnostics(fit4, vars.per.page=5)
#
#fit4.1 <- ergm(DM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0)
#fit4.1
#mcmc.diagnostics(fit4.1, vars.per.page=5)
#
#fit4.2 <- ergm(DM~edges+plus+gwesp_cotergm(-0.5, fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0)
#fit4.2
#mcmc.diagnostics(fit4.2, vars.per.page=5)



