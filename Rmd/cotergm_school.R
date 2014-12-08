
## ------------------------------------------------------------------------
options(markdown.HTML.stylesheet = system.file('doc', 'style.css', package='EpiModel'), tidy=FALSE)
opts_chunk$set(cache=TRUE, comment=NA, fig.width=9, fig.height=5.4,fig.cap="",dpi=80)
op <- par()
par(mfrow=c(1,1),mar=c(3,3,1,1), mgp=c(2,1,0))
np <- par()

fls <-list.files("../../",recursive=T)
source(paste0("../../",fls[grep("kk_utils/utils.R",fls)]))
source(paste0("../../",fls[grep("cotergm/R/functions.R",fls)]))
fig.path = "figures/"
# .libPaths(cotergm)
library("statnet.common",lib=cotergm)
library("network",lib=cotergm)
library("ergm",lib=cotergm)
library("tergm",lib=cotergm)
set.seed(123)


## ------------------------------------------------------------------------
net1 <-as.network(read.table("../s50_data/s50-network1.dat"))
net2 <-as.network(read.table("../s50_data/s50-network2.dat"))
net3 <-as.network(read.table("../s50_data/s50-network3.dat"))
alcohol <- read.table("../s50_data/s50-alcohol.dat")
drugs <- read.table("../s50_data/s50-drugs.dat")
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
constraints=~atleastnonminus(y0)

MHproposal <- MHproposal(constraints,nw=y0,arguments=list())

summary(y0~edges+plus+nodematch_cotergm("status",diff=TRUE),MHproposal=MHproposal)
summary(nw~edges+plus+nodematch_cotergm("status",diff=TRUE),MHproposal=MHproposal)
summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=0),MHproposal=MHproposal)
summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=1),MHproposal=MHproposal)
summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2),MHproposal=MHproposal)
summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=0),MHproposal=MHproposal)
summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=1),MHproposal=MHproposal)
summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2),MHproposal=MHproposal)
summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=0),MHproposal=MHproposal)
summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=1),MHproposal=MHproposal)
summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2),MHproposal=MHproposal)
summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=0),MHproposal=MHproposal)
summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=1),MHproposal=MHproposal)
summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2),MHproposal=MHproposal)


## ----, warning=FALSE-----------------------------------------------------
#fit1 <- ergm(FP~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0edge=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0, eval.loglik=F)
#fit1
#mcmc.diagnostics(fit1, vars.per.page=5)

#fit1.1 <- ergm(FP~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0, eval.loglik=F)
#fit1.1
#mcmc.diagnostics(fit1.1, vars.per.page=5)

fit1.2 <- ergm(FP~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0, eval.loglik=F)
fit1.2
mcmc.diagnostics(fit1.2, vars.per.page=5)


## ----, warning=FALSE-----------------------------------------------------
#fit2 <- ergm(FM~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0edge=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0, eval.loglik=F)
#fit2
#mcmc.diagnostics(fit2, vars.per.page=5)

fit2.1 <- ergm(FM~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0, eval.loglik=F)
fit2.1
mcmc.diagnostics(fit2.1, vars.per.page=5)

#fit2.2 <- ergm(FM~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0, eval.loglik=F)
#fit2.2
#mcmc.diagnostics(fit2.2, vars.per.page=5)


## ----, warning=FALSE-----------------------------------------------------
#fit3 <- ergm(DP~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0nodal=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0, eval.loglik=F)
#fit3
#mcmc.diagnostics(fit3, vars.per.page=5)

fit3.1 <- ergm(DP~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0, eval.loglik=F)
fit3.1
mcmc.diagnostics(fit3.1, vars.per.page=5)

#fit3.2 <- ergm(DP~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0, eval.loglik=F)
#fit3.2
#mcmc.diagnostics(fit3.2, vars.per.page=5)


## ----, warning=FALSE-----------------------------------------------------
#fit4 <- ergm(DM~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0nodal=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0, eval.loglik=F)
#fit4
#mcmc.diagnostics(fit4, vars.per.page=5)

#fit4.1 <- ergm(DM~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0, eval.loglik=F)
#fit4.1
#mcmc.diagnostics(fit4.1, vars.per.page=5)

fit4.2 <- ergm(DM~edges+plus+gwesp_cotergm(-0.5,fixed=TRUE)+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0, eval.loglik=F)
fit4.2
mcmc.diagnostics(fit4.2, vars.per.page=5)


## ----,warning=FALSE------------------------------------------------------
effstring <- c("edges","plus","gwesp.fixed-0.5","nodematch.status.1","nodematch.status.2")
mat <- cbind(coef(fit1.2)[effstring],coef(fit2.1)[effstring],coef(fit3.1)[effstring],coef(fit4.2)[effstring])

va1 <- round(coef(fit1.2)[effstring],2)
va2 <- round(coef(fit2.1)[effstring],2)
va3 <- round(coef(fit3.1)[effstring],2)
va4 <- round(coef(fit4.2)[effstring],2)

sd1 <- round(summary(fit1.2)$coefs[,4],2)
names(sd1) <- rownames(summary(fit1.2)$coefs)
sd1 <- sd1[effstring]
sd2 <- round(summary(fit2.1)$coefs[,4],2)
names(sd2) <- rownames(summary(fit2.1)$coefs)
sd2 <- sd2[effstring]
sd3 <- round(summary(fit3.1)$coefs[,4],2)
names(sd3) <- rownames(summary(fit3.1)$coefs)
sd3 <- sd3[effstring]
sd4 <- round(summary(fit4.2)$coefs[,4],2)
names(sd4) <- rownames(summary(fit4.2)$coefs)
sd4 <- sd4[effstring]
mat <- cbind(paste0(va1," (",sd1,")"),
paste0(va2," (",sd2,")"),
paste0(va3," (",sd3,")"),
paste0(va4," (",sd4,")"))


rownames(mat) <- effstring
colnames(mat) <- c("FP,SS+","FM,SS-","DP,SI+","DP,SI-")
mat


## ----,waning=FALSE-------------------------------------------------------




