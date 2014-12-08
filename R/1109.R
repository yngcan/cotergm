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

n <- 20
y0 <- network(n,directed=FALSE,density=0.2)
y0 %v% "status" <- c(rep(c(2,1),each=n/2))

#"+" is 2, "-" is 1
#y0 <- san(y0 ~ edges,target.stats=choose(n,2)*0.2)


#formula <- nw~edges+plus+nodematch_cotergm("status",diff=TRUE)
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

nw <-network(n,directed=FALSE,density=0.2)

nw %v% "status" <- c(sample(c(2,1),n,replace=T))




#plot.network(nw,label=get.vertex.attribute(nw,"status"))


FP <- FP.nw(y0,nw,"status")

FM <- FM.nw(y0,nw,"status")

DP <- DP.nw(y0,nw,"status")

DM <- DM.nw(y0,nw,"status")

#debug(ergm.mcmcslave)



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




## ------------------------------------------------------------------------
## ----, warning=FALSE-----------------------------------------------------
fit1 <- ergm(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=2),control=control.ergm(MCMC.burnin=1000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=10,MCMLE.min.effectiveSize=10,MCMC.max.interval=100, MCMC.runtime.traceplot=T), constraints=~atleastnonminus(y0),verbose=0)

#  dir /home/kirk/Dropbox/hiv_network/workspace/cotergm/ergm_trunk/src

mcmc.diagnostics(fit1)



debug(ergm.MCMLE)





