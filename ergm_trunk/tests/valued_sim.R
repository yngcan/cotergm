library(statnet.common)
library(ergm)

## StdNormal-reference
cat("Standard-normal-reference ERGM with mutuality by correlation\n")
load("testnet3d.RData")

# Joint distribution of (i,j) and (j,i)
mu<-1
sig<-2
rho<-.3

# Naural parameters of bivariate normal
# We ought to be able to specify them in a more interpretable way.
denom<--2*(1-rho^2)*sig^2
xx.coef<-1/denom+1/2 # 1/2 is from the reference measure
x.coef<-2*(-1+rho)*mu/denom
xy.coef<--2*rho/denom

theta<-c(x.coef,xy.coef,xx.coef)

cat("mean=",mu,", var=",sig^2,", corr=",rho,"\neta=(",paste(theta,collapse=","),")\n",sep="")

s<-simulate(testnet3d~sum+mutual("product")+sum(pow=2), nsim=1000, reference=~StdNormal, response="w", coef=theta,
            statsonly=TRUE, control=control.simulate(MCMC.burnin=10000))

cat("Simulated mean (statsonly):",mean(s[,1])/6,"\n",sep="")

s.full<-simulate(testnet3d~sum+mutual("product")+sum(pow=2), nsim=1000, reference=~StdNormal, response="w",
                 coef=theta, statsonly=FALSE, control=control.simulate(MCMC.burnin=10000))

s.cells<-sapply(s.full, function(x) as.matrix(x,m="a",a="w"),simplify=FALSE)
cat("Simulated means (target=",mu,"):\n",sep="")
print(matrix(c(NA,
  mean(sapply(s.cells,"[",1,2)),
  mean(sapply(s.cells,"[",1,3)),
  mean(sapply(s.cells,"[",2,1)),
  NA,
  mean(sapply(s.cells,"[",2,3)),
  mean(sapply(s.cells,"[",3,1)),
  mean(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))

cat("Simulated vars (target=",sig^2,"):\n",sep="")
print(matrix(c(NA,
  var(sapply(s.cells,"[",1,2)),
  var(sapply(s.cells,"[",1,3)),
  var(sapply(s.cells,"[",2,1)),
  NA,
  var(sapply(s.cells,"[",2,3)),
  var(sapply(s.cells,"[",3,1)),
  var(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))

cat("Simulated correlations (1,2) (1,3) (2,3) (target=",rho,"):\n",sep="")
print(c(cor(sapply(s.cells,"[",1,2),sapply(s.cells,"[",2,1)),
        cor(sapply(s.cells,"[",1,3),sapply(s.cells,"[",3,1)),
        cor(sapply(s.cells,"[",2,3),sapply(s.cells,"[",3,2))))

## StdNormal-reference with rank constraint
cat("Standard-normal-reference ERGM with rank constraint\n")
load("testrank3d.RData")

s.full<-simulate(testrank3d~sum, nsim=1000, reference=~StdNormal, response="w", coef=0, statsonly=FALSE,
                 constraints=~ranks, control=control.simulate(MCMC.burnin=10000))
s.cells<-sapply(s.full, function(x) as.matrix(x,m="a",a="w"), simplify=FALSE)
ref.sample<-pmax(rnorm(10000),rnorm(10000))

cat("Simulated means (target[1:2,]=+-",mean(ref.sample),";target[3,]=0):\n",sep="")
print(matrix(c(NA,
  mean(sapply(s.cells,"[",1,2)),
  mean(sapply(s.cells,"[",1,3)),
  mean(sapply(s.cells,"[",2,1)),
  NA,
  mean(sapply(s.cells,"[",2,3)),
  mean(sapply(s.cells,"[",3,1)),
  mean(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))

cat("Simulated vars (target[1:2,]=+-",var(ref.sample),";target[3,]=0):\n",sep="")
print(matrix(c(NA,
  var(sapply(s.cells,"[",1,2)),
  var(sapply(s.cells,"[",1,3)),
  var(sapply(s.cells,"[",2,1)),
  NA,
  var(sapply(s.cells,"[",2,3)),
  var(sapply(s.cells,"[",3,1)),
  var(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))

opttest({

set.seed(0)

## DiscUnif-reference
cat("======== Discrete-uniform-reference ERGM with minimum of -1 and maxium of 5\n")

f <- function(x, q, a, b) exp(q*x)/sum(exp(q*(a:b)))

E <- function(q, a, b) sum((a:b)*exp(q*(a:b))/sum(exp(q*(a:b))))

V <- function(q, a, b) sum(((a:b)-E(q,a,b))^2*exp(q*(a:b))/sum(exp(q*(a:b))))

a <- -1; b <- 5; coefs <- c(runif(1,-3,3),0)

for(coef in coefs){
  cat("==== statsonly=TRUE, coef=",coef,"\n",sep="")
  s <- simulate(testnet3d~sum, nsim=1000, reference=~DiscUnif(a,b), response="w", coef=coef, statsonly=TRUE, control=control.simulate(MCMC.burnin=100))
  test <- approx.hotelling.diff.test(s/6,mu0=E(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (mean).")}
  test <- approx.hotelling.diff.test((s-E(coef,a,b)*6)^2/6,mu0=V(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (variance).")}

  
  cat("==== statsonly=FALSE, coef=",coef,"\n",sep="")
  s.full<-simulate(testnet3d~sum, nsim=1000, reference=~DiscUnif(a,b), response="w", coef=coef, statsonly=FALSE, control=control.simulate(MCMC.burnin=100))
  s <- sapply(s.full, function(x) sum(as.matrix(x,m="a",a="w")),simplify=TRUE)
  test <- approx.hotelling.diff.test(s/6,mu0=E(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (mean).")}
  test <- approx.hotelling.diff.test((s-E(coef,a,b)*6)^2/6,mu0=V(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (variance).")}
}
  
## Unif-reference
cat("======== Continuous-uniform-reference ERGM with minimum of -1 and maxium of 5\n")

E <- function(q, a, b) if(isTRUE(all.equal(q,0))) (b+a)/2 else ((b*q-1)*exp(b*q)-(a*q-1)*exp(a*q))/(q*exp(b*q)-q*exp(a*q))
V <- function(q, a, b) if(isTRUE(all.equal(q,0))) (b-a)^2/12 else (((-b^2+2*a*b-a^2)*q^2-2)*exp(b*q+a*q)+exp(2*b*q)+exp(2*a*q))/(-2*q^2*exp(b*q+a*q)+q^2*exp(2*b*q)+q^2*exp(2*a*q))

a <- -1; b <- 5; coefs <- c(runif(1,-3,3),0)

for(coef in coefs){
  cat("==== statsonly=TRUE, coef=",coef,"\n",sep="")
  s <- simulate(testnet3d~sum, nsim=1000, reference=~Unif(a,b), response="w", coef=coef, statsonly=TRUE, control=control.simulate(MCMC.burnin=100))
  test <- approx.hotelling.diff.test(s/6,mu0=E(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (mean).")}
  test <- approx.hotelling.diff.test((s-E(coef,a,b)*6)^2/6,mu0=V(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (variance).")}

  
  cat("==== statsonly=FALSE, coef=",coef,"\n",sep="")
  s.full<-simulate(testnet3d~sum, nsim=1000, reference=~Unif(a,b), response="w", coef=coef, statsonly=FALSE, control=control.simulate(MCMC.burnin=100))
  s <- sapply(s.full, function(x) sum(as.matrix(x,m="a",a="w")),simplify=TRUE)
  test <- approx.hotelling.diff.test(s/6,mu0=E(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (mean).")}
  test <- approx.hotelling.diff.test((s-E(coef,a,b)*6)^2/6,mu0=V(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (variance).")}
}
}, "continuous uniform reference")
