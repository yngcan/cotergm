#===================================================================
# This file contains the 5 following MHP initializers, each
# prepended with 'InitMHP.'  All of these functions may also be
# found in the <InitMHP> file.
#      <formation>       <formationTNT>
#      <dissolution>
#===================================================================

InitMHP.formation <- function(arguments, nw, model) {
  MHproposal <- list(name = "Formation", inputs=NULL)
  MHproposal
}

InitMHP.formationTNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "FormationTNT", inputs=NULL)
  MHproposal
}

InitMHP.dissolution <- function(arguments, nw, model) {
  MHproposal <- list(name = "Dissolution", inputs=NULL)
  MHproposal
}

InitMHP.dissolutionTNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "DissolutionTNT", inputs=NULL)
  MHproposal
}
