---
title: "coevolution STERGM"
author: "Kirk Li"
date: "09/23/2014"
output:
  html_document:
    number_sections: yes
    toc: yes
---
# Coevolution STERGM
###############################################################################

```r
opts_chunk$set(cache=TRUE)
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
```


```r
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
```





```r
# Formation Plus Network
FP <- FP.nw(y0,nw,"status")
# Formation Minus Network
FM <- FM.nw(y0,nw,"status")
# Dissolution Plus Network
DP <- DP.nw(y0,nw,"status")
# Dissolution Minus Network
DM <- DM.nw(y0,nw,"status")
```



```r
constraints=~atleastnonminus(y0)

MHproposal <- MHproposal(constraints,nw=y0,arguments=list())

summary(y0~edges+plus+nodematch_cotergm("status",diff=TRUE),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 74                 17                 32 
## nodematch.status.2 
##                 16
```

```r
summary(nw~edges+plus+nodematch_cotergm("status",diff=TRUE),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 81                 19                 29 
## nodematch.status.2 
##                 22
```

```r
summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=0),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 94                 14                 43 
## nodematch.status.2 
##                 23
```

```r
summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=1),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 94                 14                 40 
## nodematch.status.2 
##                 14
```

```r
summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 94                 14                  3 
## nodematch.status.2 
##                  9
```

```r
summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=0),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                103                 22                 32 
## nodematch.status.2 
##                 30
```

```r
summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=1),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                103                 22                 22 
## nodematch.status.2 
##                 26
```

```r
summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                103                 22                 10 
## nodematch.status.2 
##                  4
```

```r
summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=0),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 54                 14                 36 
## nodematch.status.2 
##                  7
```

```r
summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=1),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 54                 14                 32 
## nodematch.status.2 
##                  7
```

```r
summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 54                 14                  4 
## nodematch.status.2 
##                  0
```

```r
summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=0),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 50                 22                 15 
## nodematch.status.2 
##                 21
```

```r
summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=1),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 50                 22                 15 
## nodematch.status.2 
##                 16
```

```r
summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2),MHproposal=MHproposal)
```

```
##              edges               plus nodematch.status.1 
##                 50                 22                  0 
## nodematch.status.2 
##                  5
```


```r
fit1.edges <- ergm(FP~edges,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0, eval.loglik=F)

fit1.plus <- ergm(FP~plus,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 0.3914 
## Iteration 2 of at most 20: 
## Convergence test P-value: 4.2e-03 
## The log-likelihood improved by 0.009541 
## Iteration 3 of at most 20: 
## Convergence test P-value: 9.2e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by < 0.0001 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
fit1.edges_plus <- ergm(FP~edges+plus,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 3.5e-269 
## The log-likelihood improved by 0.4097 
## Iteration 2 of at most 20: 
## Convergence test P-value: 1.1e-10 
## The log-likelihood improved by 0.01311 
## Iteration 3 of at most 20: 
## Convergence test P-value: 1.4e-01 
## The log-likelihood improved by 0.007404 
## Iteration 4 of at most 20: 
## Convergence test P-value: 2.8e-02 
## The log-likelihood improved by 0.05826 
## Iteration 5 of at most 20: 
## Convergence test P-value: 7.8e-02 
## The log-likelihood improved by 0.02853 
## Iteration 6 of at most 20: 
## Convergence test P-value: 3.9e-01 
## The log-likelihood improved by 0.007823 
## Iteration 7 of at most 20: 
## Convergence test P-value: 1.4e-01 
## The log-likelihood improved by 0.007779 
## Iteration 8 of at most 20: 
## Convergence test P-value: 3.2e-02 
## The log-likelihood improved by 0.005597 
## Iteration 9 of at most 20: 
## Convergence test P-value: 3e-03 
## The log-likelihood improved by 0.02465 
## Iteration 10 of at most 20: 
## Convergence test P-value: 2e-04 
## The log-likelihood improved by 0.1051 
## Iteration 11 of at most 20: 
## Convergence test P-value: 2.1e-02 
## The log-likelihood improved by 0.04564 
## Iteration 12 of at most 20: 
## Convergence test P-value: 2.7e-01 
## The log-likelihood improved by 0.005615 
## Iteration 13 of at most 20: 
## Convergence test P-value: 3.5e-01 
## The log-likelihood improved by 0.003722 
## Iteration 14 of at most 20: 
## Convergence test P-value: 2.3e-03 
## The log-likelihood improved by 0.0349 
## Iteration 15 of at most 20: 
## Convergence test P-value: 1e-01 
## The log-likelihood improved by 0.0251 
## Iteration 16 of at most 20: 
## Convergence test P-value: 7.7e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.0001951 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
fit1.edges_plus_nodematch <- ergm(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonminus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 4.9e-234 
## The log-likelihood improved by 3.563 
## Iteration 2 of at most 20: 
## Convergence test P-value: 2.4e-46 
## The log-likelihood improved by 0.547 
## Iteration 3 of at most 20: 
## Convergence test P-value: 2.3e-31 
## The log-likelihood improved by 0.1308 
## Iteration 4 of at most 20: 
## Convergence test P-value: 2.3e-03 
## The log-likelihood improved by 0.06414 
## Iteration 5 of at most 20: 
## Convergence test P-value: 2.7e-01 
## The log-likelihood improved by 0.01536 
## Iteration 6 of at most 20: 
## Convergence test P-value: 4.8e-01 
## The log-likelihood improved by 0.01059 
## Iteration 7 of at most 20: 
## Convergence test P-value: 2.7e-01 
## The log-likelihood improved by 0.02147 
## Iteration 8 of at most 20: 
## Convergence test P-value: 6.6e-04 
## The log-likelihood improved by 0.0923 
## Iteration 9 of at most 20: 
## Convergence test P-value: 1.8e-03 
## The log-likelihood improved by 0.003759 
## Iteration 10 of at most 20: 
## Convergence test P-value: 3.7e-02 
## The log-likelihood improved by 0.01467 
## Iteration 11 of at most 20: 
## Convergence test P-value: 6.5e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.007387 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```


```r
fit2.edges <- ergm(FM~edges,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0, eval.loglik=F)

fit2.plus <- ergm(FM~plus,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 5.725 
## Iteration 2 of at most 20: 
## Convergence test P-value: 2e-120 
## The log-likelihood improved by 4.606 
## Iteration 3 of at most 20: 
## Convergence test P-value: 9.9e-14 
## The log-likelihood improved by 0.3637 
## Iteration 4 of at most 20: 
## Convergence test P-value: 9.4e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by < 0.0001 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
fit2.edges_plus <- ergm(FM~edges+plus,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 6.546 
## Iteration 2 of at most 20: 
## Convergence test P-value: 1.2e-84 
## The log-likelihood improved by 3.866 
## Iteration 3 of at most 20: 
## Convergence test P-value: 7.8e-32 
## The log-likelihood improved by 0.2075 
## Iteration 4 of at most 20: 
## Convergence test P-value: 2.3e-06 
## The log-likelihood improved by 0.02126 
## Iteration 5 of at most 20: 
## Convergence test P-value: 7.9e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.00345 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
fit2.edges_plus_nodematch <- ergm(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atleastnonplus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 7.703 
## Iteration 2 of at most 20: 
## Convergence test P-value: 5.1e-101 
## The log-likelihood improved by 1.784 
## Iteration 3 of at most 20: 
## Convergence test P-value: 1.2e-22 
## The log-likelihood improved by 0.09794 
## Iteration 4 of at most 20: 
## Convergence test P-value: 9.1e-03 
## The log-likelihood improved by 0.03835 
## Iteration 5 of at most 20: 
## Convergence test P-value: 3.2e-01 
## The log-likelihood improved by 0.02945 
## Iteration 6 of at most 20: 
## Convergence test P-value: 5.8e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.006642 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```


```r
fit3.edges <- ergm(DP~edges,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0, eval.loglik=F)

fit3.plus <- ergm(DP~plus,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 0.3849 
## Iteration 2 of at most 20: 
## Convergence test P-value: 1.4e-01 
## The log-likelihood improved by 0.00227 
## Iteration 3 of at most 20: 
## Convergence test P-value: 3.2e-01 
## The log-likelihood improved by 0.002628 
## Iteration 4 of at most 20: 
## Convergence test P-value: 4.5e-01 
## The log-likelihood improved by 0.001378 
## Iteration 5 of at most 20: 
## Convergence test P-value: 7.2e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.0003217 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
fit3.edges_plus <- ergm(DP~edges+plus,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 0.3902 
## Iteration 2 of at most 20: 
## Convergence test P-value: 8.4e-05 
## The log-likelihood improved by 0.02541 
## Iteration 3 of at most 20: 
## Convergence test P-value: 9.1e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.0006517 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
fit3.edges_plus_nodematch <- ergm(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=1),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonminus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 6.131 
## Iteration 2 of at most 20: 
## Convergence test P-value: 7.4e-10 
## The log-likelihood improved by 0.07171 
## Iteration 3 of at most 20: 
## Convergence test P-value: 3.1e-01 
## The log-likelihood improved by 0.01637 
## Iteration 4 of at most 20: 
## Convergence test P-value: 2.8e-03 
## The log-likelihood improved by 0.06527 
## Iteration 5 of at most 20: 
## Convergence test P-value: 4.2e-02 
## The log-likelihood improved by 0.03248 
## Iteration 6 of at most 20: 
## Convergence test P-value: 8.8e-03 
## The log-likelihood improved by 0.0259 
## Iteration 7 of at most 20: 
## Convergence test P-value: 5.4e-03 
## The log-likelihood improved by 0.0477 
## Iteration 8 of at most 20: 
## Convergence test P-value: 7.8e-05 
## The log-likelihood improved by 0.07253 
## Iteration 9 of at most 20: 
## Convergence test P-value: 3e-02 
## The log-likelihood improved by 0.02481 
## Iteration 10 of at most 20: 
## Convergence test P-value: 7.8e-02 
## The log-likelihood improved by 0.01764 
## Iteration 11 of at most 20: 
## Convergence test P-value: 8.7e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.002645 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
fit4.edges <- ergm(DM~edges,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0, eval.loglik=F)

fit4.plus <- ergm(DM~plus,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 10.48 
## Iteration 2 of at most 20: 
## Convergence test P-value: 3.7e-71 
## The log-likelihood improved by 1.767 
## Iteration 3 of at most 20: 
## Convergence test P-value: 1.6e-05 
## The log-likelihood improved by 0.1224 
## Iteration 4 of at most 20: 
## Convergence test P-value: 2.7e-01 
## The log-likelihood improved by 0.002831 
## Iteration 5 of at most 20: 
## Convergence test P-value: 5.1e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.002731 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
fit4.edges_plus <- ergm(DM~edges+plus,control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 5.605 
## Iteration 2 of at most 20: 
## Convergence test P-value: 5.6e-98 
## The log-likelihood improved by 4.085 
## Iteration 3 of at most 20: 
## Convergence test P-value: 2.9e-10 
## The log-likelihood improved by 0.1427 
## Iteration 4 of at most 20: 
## Convergence test P-value: 3.6e-01 
## The log-likelihood improved by 0.006937 
## Iteration 5 of at most 20: 
## Convergence test P-value: 1.3e-02 
## The log-likelihood improved by 0.02203 
## Iteration 6 of at most 20: 
## Convergence test P-value: 5.8e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.003472 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```

```r
fit4.edges_plus_nodematch <- ergm(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=2),control=control.ergm(MCMC.burnin=10000,MCMC.samplesize=10000,MCMC.interval=100,MCMLE.maxit=20,MCMLE.min.effectiveSize=20,MCMC.max.interval=100), constraints=~atmostnonplus(y0),verbose=0, eval.loglik=F)
```

```
## Iteration 1 of at most 20: 
## Convergence test P-value: 0e+00 
## The log-likelihood improved by 7.939 
## Iteration 2 of at most 20: 
## Convergence test P-value: 7e-105 
## The log-likelihood improved by 5.184 
## Iteration 3 of at most 20: 
## Convergence test P-value: 2e-25 
## The log-likelihood improved by 0.5053 
## Iteration 4 of at most 20: 
## Convergence test P-value: 3.4e-01 
## The log-likelihood improved by 0.01001 
## Iteration 5 of at most 20: 
## Convergence test P-value: 2.6e-02 
## The log-likelihood improved by 0.02306 
## Iteration 6 of at most 20: 
## Convergence test P-value: 7.6e-01 
## Convergence detected. Stopping.
## The log-likelihood improved by 0.004009 
## 
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
```



```r
fit1.list <- lapply(ls()[grep("fit1",ls())][c(1,4,2,3)],function(x)eval(parse(text=paste0('summary(',x,')$coef'))))
```

```
## Error in summary(fit1.list)$coef: $ operator is invalid for atomic vectors
```

```r
fit2.list <- lapply(ls()[grep("fit2",ls())][c(1,4,2,3)],function(x)eval(parse(text=paste0('summary(',x,')$coef'))))
```

```
## Error in summary(fit2.list)$coef: $ operator is invalid for atomic vectors
```

```r
fit3.list <- lapply(ls()[grep("fit3",ls())][c(1,4,2,3)],function(x)eval(parse(text=paste0('summary(',x,')$coef'))))
```

```
## Error in summary(fit3.list)$coef: $ operator is invalid for atomic vectors
```

```r
fit4.list <- lapply(ls()[grep("fit4",ls())][c(1,4,2,3)],function(x)eval(parse(text=paste0('summary(',x,')$coef'))))
```

```
## Error in summary(fit4.list)$coef: $ operator is invalid for atomic vectors
```

```r
fit1.list
```

```
## [[1]]
##        Estimate Std. Error MCMC %      p-value
## edges -3.457893  0.2270941     NA 4.897876e-45
## 
## [[2]]
##      Estimate Std. Error MCMC %    p-value
## plus 1.587755  0.7131211      0 0.02632147
## 
## [[3]]
##        Estimate Std. Error MCMC %      p-value
## edges -3.458554  0.2604048      0 8.341438e-36
## plus   1.610638  0.6978884      0 2.131790e-02
## 
## [[4]]
##                      Estimate Std. Error MCMC %      p-value
## edges              -3.9104041  0.3219941      0 9.713269e-31
## plus                0.4790889  0.7693421      0 5.336827e-01
## nodematch.status.2  1.8431451  0.4732825      0 1.085947e-04
```

```r
fit2.list
```

```
## [[1]]
##        Estimate Std. Error MCMC %      p-value
## edges -3.542457  0.1883518     NA 4.575177e-68
## 
## [[2]]
##       Estimate Std. Error MCMC %     p-value
## plus -1.735732  0.5453915      0 0.001503589
## 
## [[3]]
##        Estimate Std. Error MCMC %      p-value
## edges -3.554596  0.1924233      0 5.028375e-66
## plus  -1.764495  0.5253348      0 8.115170e-04
## 
## [[4]]
##                       Estimate Std. Error MCMC %      p-value
## edges              -3.55480252  0.2288698      0 4.978876e-49
## plus               -1.77812653  0.5964075      0 2.936985e-03
## nodematch.status.1 -0.02660311  0.3948599      0 9.462974e-01
```

```r
fit3.list
```

```
## [[1]]
##         Estimate Std. Error MCMC %   p-value
## edges 0.09531018   0.308957     NA 0.7592715
## 
## [[2]]
##      Estimate Std. Error MCMC %    p-value
## plus  1.63531   0.692017      0 0.02294643
## 
## [[3]]
##        Estimate Std. Error MCMC %    p-value
## edges 0.1256747  0.3466263      0 0.71883845
## plus  1.6378520  0.6792311      0 0.02057776
## 
## [[4]]
##                      Estimate Std. Error MCMC %   p-value
## edges              0.05552453  0.3481871      0 0.8741236
## plus               1.91469855  0.9195429      0 0.0439296
## nodematch.status.1 0.35364347  0.4625858      0 0.4491748
```

```r
fit4.list
```

```
## [[1]]
##        Estimate Std. Error MCMC %   p-value
## edges 0.3483067  0.2666054     NA 0.1966466
## 
## [[2]]
##       Estimate Std. Error MCMC %      p-value
## plus -1.728007  0.4657376      0 0.0004715494
## 
## [[3]]
##        Estimate Std. Error MCMC %     p-value
## edges  0.340181  0.2771656      0 0.224823315
## plus  -1.764267  0.5336227      0 0.001654617
## 
## [[4]]
##                     Estimate Std. Error MCMC %     p-value
## edges               0.237703  0.2771837      0 0.394854251
## plus               -2.498611  0.7683464      0 0.001960458
## nodematch.status.2  1.041502  0.5866911      0 0.081395997
```


