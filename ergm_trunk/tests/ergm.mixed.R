library(ergm)
opttest({
library(statnet.common)
# import synthetic network that looks like a molecule
data(molecule)
set.vertex.attribute(molecule,"atomic type",c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3))
#
# create a plot of the social network
# colored by atomic type
#
plot(molecule, vertex.col=molecule %v% "atomic type",vertex.cex=3)

# measure tendency to match within each atomic type
gest <- ergm(molecule ~ edges + kstar(2) + triangle + nodematch("atomic type"),
# mixed=TRUE,
             control=control.ergm(MCMC.samplesize=10000))
summary(gest)

# compare it to differential homophily
gest <- ergm(molecule ~ edges + kstar(2) + triangle + nodematch("atomic type",diff=TRUE),
# mixed=TRUE,
             control=control.ergm(MCMC.samplesize=10000))
summary(gest)
}, "ergm mixed")
