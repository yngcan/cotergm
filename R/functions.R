# TODO: Add comment
# 
# Author: kirk
###############################################################################

setdiff.mat <- function(mat1,mat2){
	if(!is.matrix(mat1))
		mat1 <- t(as.matrix(mat1))
	if(!is.matrix(mat2))
		mat2 <- t(as.matrix(mat2))
	match(setdiff(apply(mat2,1,paste,collapse="_"), apply(mat1,1,paste,collapse="_")) ,apply(mat2,1,paste,collapse="_"))
	
}

setdiff(c(1,2,3),c(2,4))




FP.nw <- function(net1,net2,status){
	if (!all(net1 %v% status %in% c(2,1)))
	stop("net1 nodal status is not correctly specified")		
	index.mm <- which(outer(net1 %v% status , net1 %v% status,paste0)=="11",arr.ind = TRUE)

	# calculate FP dyads
	index <- setdiff.mat(mat1=as.edgelist(net1),mat2=as.edgelist(net2))
	formed.edge <- as.edgelist(net2)[index,,drop=FALSE]
	formed.edge.p <- formed.edge[setdiff.mat(index.mm,formed.edge),,drop=FALSE]
	net.tmp <- net1
	formed.net.p <- add.edges(net.tmp,formed.edge.p[,1],formed.edge.p[,2])
	
	# calculate FP nodal attributes
	formed.net.p %v% status <- net2 %v% status
	net1.al <- net1%v%status
	formed.net.p.al <- formed.net.p %v% status
	formed.net.p.al[which(net1.al==1)] <- 1
	formed.net.p %v% status <- formed.net.p.al 
	formed.net.p
}



FM.nw <- function(net1,net2,status){
	
	if (!all(net1 %v% status %in% c(2,1)))
		stop("net1 nodal status is not correctly specified")
	
		index.pp <- which(outer(net1 %v% status , net1 %v% status,paste0)=="22",arr.ind = TRUE)
	
	index <- setdiff.mat(mat1=as.edgelist(net1),mat2=as.edgelist(net2))
	formed.edge <- as.edgelist(net2)[index,,drop=FALSE]
	formed.edge.m <- formed.edge[setdiff.mat(index.pp,formed.edge),,drop=FALSE]
	net.tmp <- net1
	formed.net.m <- add.edges(net.tmp,formed.edge.m[,1],formed.edge.m[,2])
	
	
	# calculate FM nodal attributes
	formed.net.m %v% status <- net2 %v% status
	net1.al <- net1%v%status
	formed.net.m.al <- formed.net.m %v% status
	formed.net.m.al[which(net1.al==2)] <- 2
	formed.net.m %v% status <- formed.net.m.al 
	formed.net.m
	
}




DP.nw <- function(net1,net2,status){
	
	if (!all(net1 %v% status %in% c(2,1)))
		stop("net1 nodal status is not correctly specified")
	
		index.mm <- which(outer(net1 %v% status , net1 %v% status,paste0)=="11",arr.ind = TRUE)
	index <- setdiff.mat(mat1=as.edgelist(net2),mat2=as.edgelist(net1))
	dissolved.edge <- as.edgelist(net1)[index,,drop=FALSE]
	dissolved.edge.p <- dissolved.edge[setdiff.mat(index.mm,dissolved.edge),,drop=FALSE]
	#sort tail < head
	dissolved.edge.p <- t(apply(dissolved.edge.p, 1, sort))
	net.tmp <- net1
	if(length(net.tmp$mel)){
	el <- t(apply(t(sapply(net.tmp$mel,function(x)unlist(x)[1:2])), 1, sort))
	delete.edge.id.p <- match(apply(dissolved.edge.p,1,paste,collapse="_") , apply(el,1,paste,collapse="_"))
}
	dissolved.net.p <- net.tmp
	
	if(exists("delete.edge.id.p")){
	if(!any(is.na(delete.edge.id.p))){
		dissolved.net.p <- delete.edges(net.tmp,delete.edge.id.p)} else {warnings("no edge deleted")}
		dissolved.net.p <- net.tmp
	}
	
	# calculate DP nodal attributes
	dissolved.net.p %v% status <- net2 %v% status
	net1.al <- net1%v%status
	dissolved.net.p.al <- dissolved.net.p %v% status
	dissolved.net.p.al[which(net1.al==1)] <- 1
	dissolved.net.p %v% status <- dissolved.net.p.al 
	dissolved.net.p
}





DM.nw <- function(net1,net2,status){
	
	if (!all(net1 %v% status %in% c(2,1)))
		stop("net1 nodal status is not correctly specified")
	
		index.pp <- which(outer(net1 %v% status , net1 %v% status,paste0)=="22",arr.ind = TRUE)
	
	index <- setdiff.mat(mat1=as.edgelist(net2),mat2=as.edgelist(net1))
	dissolved.edge <- as.edgelist(net1)[index,,drop=FALSE]
	dissolved.edge.m <- dissolved.edge[setdiff.mat(index.pp,dissolved.edge),,drop=FALSE]
	#sort tail < head
	dissolved.edge.m <- t(apply(dissolved.edge.m, 1, sort))
	net.tmp <- net1
	if(length(net.tmp$mel)){
	el <- t(apply(t(sapply(net.tmp$mel,function(x)unlist(x)[1:2])), 1, sort))
	delete.edge.id.m <- match(apply(dissolved.edge.m,1,paste,collapse="_") , apply(el,1,paste,collapse="_"))
}
	dessolved.net.m <- net.tmp
	if(exists("delete.edge.id.m")){
	if(!any(is.na(delete.edge.id.m))){
	dissolved.net.m <- delete.edges(net.tmp,delete.edge.id.m)} else {warnings("no edge deleted")}}
	dissolved.net.m <- net.tmp
	

	# calculate DM nodal attributes
	dissolved.net.m %v% status <- net2 %v% status
	net1.al <- net1%v%status
	dissolved.net.m.al <- dissolved.net.m %v% status
	dissolved.net.m.al[which(net1.al==2)] <- 2
	dissolved.net.m %v% status <- dissolved.net.m.al 
	dissolved.net.m

}