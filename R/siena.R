# TODO: Add comment
# Aug 6, 2014
# Author: Kirk Li
# Email: kirli@stat.washington.edu
###############################################################################


######################################################
#              Lab_coevolution.R                     #
#                                                    #
# R script for the assignment on co-evolution        #
# of alcohol consumption and friendship.             #
# Written by Christian Steglich,                     #
# with some additions by Tom Snijders.               #
# Version March 26, 2014                             #
#                                                    #
######################################################


# load RSiena commands:
library(RSiena)

# This analysis uses the s50 data,
# an excerpt of which is included in RSiena.
# For a description, see
# http://www.stats.ox.ac.uk/~snijders/siena/s50_data.htm

# See short description on the help page:
?s50

# Identify dependent network variable:
friendship <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
# Identify dependent behavior variable:
drinking <- sienaDependent(s50a, type="behavior")

# Bind data together for Siena analysis:
CoEvolutionData <- sienaDataCreate(friendship,drinking)
# What do we have:
CoEvolutionData

# Create effects object for model specification:
CoEvolutionEffects <- getEffects(CoEvolutionData)

# First write descriptive results to file:
print01Report(CoEvolutionData, CoEvolutionEffects, modelname = 'CoEvol_init')
# Look at this file!

# Specify the model according to the assignment.
# First the structural effects that are not yet included by default:
CoEvolutionEffects <- includeEffects(CoEvolutionEffects,transTrip,transRecTrip,cycle3)
CoEvolutionEffects
# Now sender, receiver and homophily effects for friendship formation:
CoEvolutionEffects <- includeEffects(CoEvolutionEffects,
    egoX,altX,simX,interaction1="drinking")
CoEvolutionEffects
# Now the assimilation effect for drinking:
CoEvolutionEffects <- includeEffects(CoEvolutionEffects,
    name="drinking",avSim,interaction1="friendship")
CoEvolutionEffects
# Note that you need to additionally specify 'name="drinking"' because
# the default for 'name' is the network variable (here "friendship").

# Now create an algorithm object:
CoEvolutionAlgo <- sienaAlgorithmCreate(projname = 'CoEvol_results')

# Estimate the model:
CoEvolutionResults <- siena07(CoEvolutionAlgo, data=CoEvolutionData,
    effects=CoEvolutionEffects)
CoEvolutionResults
# If the t-ratios for convergence for the non-fixed effects
# are not satisfactorily small (all less than 0.1 in absolute value),
# run the estimation again:
CoEvolutionResults <- siena07(CoEvolutionAlgo, data=CoEvolutionData,
    effects=CoEvolutionEffects, prevAns=CoEvolutionResults)
# Test the specific parameter of interest:
?Wald.RSiena
Multipar.RSiena(CoEvolutionResults, 15)

# Replace Average Similarity by Total Similarity:
CoEvolutionEffects <- includeEffects(CoEvolutionEffects,include=FALSE,
    name="drinking",avSim,interaction1="friendship")
CoEvolutionEffects <- includeEffects(CoEvolutionEffects,
    name="drinking",totSim,interaction1="friendship")
CoEvolutionEffects
CoEvolutionResults <- siena07(CoEvolutionAlgo, data=CoEvolutionData,
    effects=CoEvolutionEffects)
CoEvolutionResults
# If the t-ratios for convergence for the non-fixed effects
# are not satisfactorily small (all less than 0.1 in absolute value),
# run the estimation again:
CoEvolutionResults <- siena07(CoEvolutionAlgo, data=CoEvolutionData,
    effects=CoEvolutionEffects, prevAns=CoEvolutionResults)
Multipar.RSiena(CoEvolutionResults, 15)

# Replace Total Similarity by Average Alter:
CoEvolutionEffects <- includeEffects(CoEvolutionEffects,include=FALSE,
    name="drinking",totSim,interaction1="friendship")
CoEvolutionEffects <- includeEffects(CoEvolutionEffects,
    name="drinking",avAlt,interaction1="friendship")
CoEvolutionEffects
CoEvolutionResults <- siena07(CoEvolutionAlgo, data=CoEvolutionData,
    effects=CoEvolutionEffects)
CoEvolutionResults
# If the t-ratios for convergence for the non-fixed effects
# are not satisfactorily small (all less than 0.1 in absolute value),
# run the estimation again:
CoEvolutionResults <- siena07(CoEvolutionAlgo, data=CoEvolutionData,
    effects=CoEvolutionEffects, prevAns=CoEvolutionResults)
Multipar.RSiena(CoEvolutionResults, 15)

# What could be called "Total Alter" is not implemented directly,
# but can be obtained as the interaction between Average Alter and Outdegree.
# We need to include the main effects in the model, but keep their
# parameters fixed at 0.

?includeInteraction
CoEvolutionEffects <- setEffect(CoEvolutionEffects,avAlt,fix=TRUE,test=TRUE,
    name="drinking",interaction1="friendship")
CoEvolutionEffects <- setEffect(CoEvolutionEffects,outdeg,fix=TRUE,test=TRUE,
    name="drinking",interaction1="friendship")
CoEvolutionEffects <- includeInteraction(CoEvolutionEffects,
    name="drinking",avAlt,outdeg,
    interaction1=c("friendship","friendship"))
CoEvolutionEffects
CoEvolutionResults <- siena07(CoEvolutionAlgo, data=CoEvolutionData,
    effects=CoEvolutionEffects)
CoEvolutionResults
# If the t-ratios for convergence for the non-fixed effects
# are not satisfactorily small (all less than 0.1 in absolute value),
# run the estimation again:
CoEvolutionResults <- siena07(CoEvolutionAlgo, data=CoEvolutionData,
    effects=CoEvolutionEffects, prevAns=CoEvolutionResults)
CoEvolutionResults
Multipar.RSiena(CoEvolutionResults, 17)

# To see the results of the score-type tests for
# the indegree and outdegree effects on drinking:
CoEvolutionEffects <- setEffect(CoEvolutionEffects,indeg,fix=TRUE,test=TRUE,
    name="drinking",interaction1="friendship")
CoEvolutionEffects
CoEvolutionResults <- siena07(CoEvolutionAlgo, data=CoEvolutionData,
    effects=CoEvolutionEffects, prevAns=CoEvolutionResults)
summary(CoEvolutionResults)
# The test for the outdegree effect comes close to significance,
# the other two are clearly not significant.

# We can conclude at this point that for all four specifications of
# social influence there appears a positive estimated effect,
# but the significance depends on the specification.
# This data set is too small to allow simultaneous estimation
# of all four influence effects and do a backward elimination.
# Choosing the most significant result is methodologically not permitted,
# because it amounts to chance capitalization:
# the level of significance would not be respected.
# The best procedure is to make a choice before seeing the data,
# based on theory or prior experience with similar data sets.
# The choice can also be made based on the fit of the model,
# if results from sienaGOF() differentiate between the models.
# This is not elaborated in the current script.
# Since the p-values in all four models are between 0.05 and 0.10
# (this will also depend on the random simulations;
#  results will be more stable with n3=3000 in sienaAlgorithmCreate)
# we can conclude that there is a tendency toward evidence
# for social influence.

# All estimation results are in the file CoEvol_results.out.
