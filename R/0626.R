# TODO: Add comment
# 
# Author: kecoli
###############################################################################

library(statnet)
data(samplk)
debug(stergm.CMLE)

# Fit a transition from Time 1 to Time 2
samplk12 <- stergm(list(samplk1, samplk2),
        formation=~edges+mutual+transitiveties+cyclicalties,
        dissolution=~edges+mutual+transitiveties+cyclicalties,
        estimate="CMLE")
#
#mcmc.diagnostics(samplk12)
#summary(samplk12)




?ergm