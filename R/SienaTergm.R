# TODO: Add comment
# Aug 7, 2014
# Author: Kirk Li
# Email: kirli@stat.washington.edu
###############################################################################

library(tergm)
library(RSiena)
data(samplk)

net1 <- network.initialize(100,directed=FALSE)

net1 <- san(net1~edges,target.stats=choose(100,2)/2)

summary(net1~density)

net2 <- invert.network(net1)

summary(net2~density)



net12 <- stergm(list(net1, net2),    
		formation=~edges,
		dissolution=~edges,
		estimate="CMLE")

net12


net12.s <- sienaDependent(
				array( c(as.matrix(net1), as.matrix(net2)),dim=c(100,100,2)))

mydata <- sienaDataCreate(net12.s)
mydata

myeff <- getEffects( mydata )
#myeff <- includeEffects( myeff, outPop, transTies, balance,transTrip, cycle3 )
myeff
myalgorithm <- sienaAlgorithmCreate(useStdInits = FALSE, projname = 'samp')
ans <- siena07( myalgorithm, data = mydata, effects = myeff)








samplk12 <- stergm(list(samplk1, samplk2),    formation=~edges+mutual+odegreepopularity+ttriple+transitiveties+balance+cycle(3),
    dissolution=~edges+mutual+odegreepopularity+ttriple+transitiveties+balance+cycle(3),
    estimate="CMLE")


samplk12.s <- sienaDependent(
    array( c(as.matrix(samplk1), as.matrix(samplk2)),dim=c(18,18,2)))

mydata <- sienaDataCreate(samplk12.s)
mydata

myeff <- getEffects( mydata )
myeff <- includeEffects( myeff, outPop, transTies, balance,transTrip, cycle3 )
myeff
myalgorithm <- sienaAlgorithmCreate(useStdInits = FALSE, projname = 'samp')
ans <- siena07( myalgorithm, data = mydata, effects = myeff)


summary(samplk12)