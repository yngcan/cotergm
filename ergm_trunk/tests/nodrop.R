### Tests to make sure drop=FALSE works.
library(statnet.common)
opttest({
library(ergm)
data(sampson)

## Shouldn't need to drop.
# MPLE
summary(ergm(samplike~edges, control=control.ergm(drop=FALSE)))
# MCMC
summary(ergm(samplike~edges, control=control.ergm(drop=FALSE, force.main=TRUE)))

## Empty network.
y0 <- network.initialize(10)
# MPLE
summary(ergm(y0~edges, control=control.ergm(drop=FALSE)))
# MCMC
summary(ergm(y0~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0)))

## Full network.
y1 <- as.network(matrix(1,10,10))
# MPLE
summary(ergm(y1~edges, control=control.ergm(drop=FALSE)))
# MCMC
summary(ergm(y1~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0)))
}, "drop disabled")
