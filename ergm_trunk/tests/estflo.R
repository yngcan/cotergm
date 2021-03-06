library(statnet.common)
opttest({
library(ergm)
data(florentine)
# a markov graph fit to the Florentine data
gest <- ergm(flomarriage ~ edges + kstar(2), control=control.ergm(seed=16124))
gest
summary(gest)
#anova(gest)

#Newton-Raphson iterations:  4
#MCMC sample of size 1000 based on:
#   edges     star2
#-1.66463   0.01181
#
#Monte Carlo MLE Coefficients:
#    edges      star2
#-1.622292   0.006467

# While we are at it, test the constrainted version.
# (The edges term will be ignored because the constraint makes it irrelevant.)
# XXX uncomment:
#gest <- ergm(flomarriage ~ edges + kstar(2), constraints=~edges, control=control.ergm(seed=16124))
gest <- ergm(flomarriage ~ kstar(2), constraints=~edges, control=control.ergm(seed=16124))
gest
summary(gest)
}, "Florentine")
