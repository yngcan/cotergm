
#################################

InitMHP.formationPlusMLE <- function(arguments, nw){
	y0<-arguments$constraints$atleastnonminus$free.dyads()
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	y0.nw <- arguments$constraints$atleastnonminus$nw
	MHproposal <- list(name = "randomtoggleList_nodal", inputs=c(ergm.Cprepare.el(y0),sum(get.vertex.attribute(y0.nw,"status")=="2"),which(get.vertex.attribute(y0.nw,"status")=="2")), pkgname="ergm")
	MHproposal
	
}



InitMHP.formationMinusMLE <- function(arguments, nw){
	y0<-arguments$constraints$atleastnonplus$free.dyads()
	y0.nw <- arguments$constraints$atleastnonplus$nw
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	MHproposal <- list(name = "randomtoggleList_nodal", inputs=c(ergm.Cprepare.el(y0),sum(get.vertex.attribute(y0.nw,"status")==1),which(get.vertex.attribute(y0.nw,"status")==1)), pkgname="ergm")
	MHproposal
	
}



InitMHP.dissolutionPlusMLE <- function(arguments, nw){
	nw <- arguments$constraints$atmostnonminus$nw
	y0<-arguments$constraints$atmostnonminus$free.dyads()
	y0.nw <- arguments$constraints$atmostnonminus$nw
	
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	MHproposal <- list(name = "randomtoggleList_nodal", inputs=c(ergm.Cprepare.el(y0),sum(get.vertex.attribute(y0.nw,"status")==2),which(get.vertex.attribute(y0.nw,"status")==2)), pkgname="ergm")
	MHproposal
	
}



InitMHP.dissolutionMinusMLE <- function(arguments, nw){
	y0<-arguments$constraints$atmostnonplus$free.dyads()
	y0.nw <- arguments$constraints$atmostnonplus$nw
	
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	MHproposal <- list(name = "randomtoggleList_nodal", inputs=c(ergm.Cprepare.el(y0),sum(get.vertex.attribute(y0.nw,"status")==1),which(get.vertex.attribute(y0.nw,"status")==1)), pkgname="ergm")
	MHproposal
	
}

#################################
