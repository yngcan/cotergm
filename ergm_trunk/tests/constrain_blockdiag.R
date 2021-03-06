library(statnet.common)
opttest({
library(ergm)
n <- 10
a <- rep(1:4,1:4)

M <- matrix(0,n,n)

for(i in unique(a)){
  M[a==i,a==i]<-1
}
diag(M)<-0

#### Directed ####

y0 <- network.initialize(n, directed=TRUE)
y0 %v% "b" <- a

y <- simulate(y0~edges, coef=100, constraints=~blockdiag("b"), control=control.simulate.formula(MCMC.burnin=10000))

stopifnot(all(as.matrix(y)==M))

#### Undirected ####

y0 <- network.initialize(n, directed=FALSE)
y0 %v% "b" <- a

y <- simulate(y0~edges, coef=100, constraints=~blockdiag("b"), control=control.simulate.formula(MCMC.burnin=10000))

stopifnot(all(as.matrix(y)==M))

#### Unobserved ####

y0 <- network.initialize(n, directed=TRUE)
y0 %v% "b" <- a
y0[2,3]<-NA
y0[2,10]<-NA

y <- simulate(y0~edges, coef=100, constraints=~blockdiag("b")+observed, control=control.simulate.formula(MCMC.burnin=10000))

M[]<-0
M[2,3]<-1

stopifnot(all(as.matrix(y)==M))
}, "block diagonal constraint")
