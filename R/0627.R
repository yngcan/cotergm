# TODO: Add comment
# 
# Author: kirk
###############################################################################



library(tergm,lib="/home/kirk/Project/R-3.1.0/library_trunk")



data(samplk)


# Fit a transition from Time 1 to Time 2
samplk12 <- stergm(list(samplk1, samplk2),
		formation=~edges+mutual+transitiveties+cyclicalties,
		dissolution=~edges+mutual+transitiveties+cyclicalties,
		estimate="CMLE")






debug(MHproposal.character)



library(statnet)

net1 <- network.initialize(10,directed=FALSE)

