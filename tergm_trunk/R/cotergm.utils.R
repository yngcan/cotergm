# TODO: Add comment
# 
# Author: kirk
###############################################################################
#
#
#nw <- network(10,directed=FALSE,density=0.5)
#set.vertex.attribute(nw,"status",sample(c(2,1),10,replace=TRUE))
#get.vertex.attribute(nw,"status")


get.MM.dyads <- function(nw,is.edgelist=FALSE,...){
	if(!all(unlist(lapply(nw$val,function(x)x$status)) %in% c(2,1))) stop("node status is not specifed correctly")
	node.status <- get.vertex.attribute(nw,"status")
	if(any(is.na(node.status))) stop("node status is missing")
	minus.edge.index <- which(node.status==1)
	nw.new <- network.initialize(network.size(nw))
	nw.new$gal$directed <- nw$gal$directed
	edgelist <- expand.grid(minus.edge.index,minus.edge.index)
	
	if(is.directed(nw)){
	edgelist.final <- edgelist[!apply(edgelist,1,function(x)(x[1]==x[2])),,drop=FALSE]
	} else{
	edgelist.temp <- unique(t(apply(edgelist,1,sort)))
	edgelist.final <- edgelist.temp[!apply(edgelist.temp,1,function(x)(x[1]==x[2])),,drop=FALSE]
	}
	
	add.edges(nw.new,tail=edgelist.final[,1],head=edgelist.final[,2])
	if(is.edgelist)
		as.edgelist(nw.new) else
		nw.new} 


get.PP.dyads <- function(nw,is.edgelist=FALSE,...){
	if(!all(unlist(lapply(nw$val,function(x)x$status)) %in% c(2,1))
			) stop("node status is not specifed correctly")
	node.status <- get.vertex.attribute(nw,"status")
	if(any(is.na(node.status))) stop("node status is missing")
	plus.edge.index <- which(node.status==2)
	nw.new <- network.initialize(network.size(nw))
	nw.new$gal$directed <- nw$gal$directed
	edgelist <- expand.grid(plus.edge.index,plus.edge.index)
	
	if(is.directed(nw)){
		edgelist.final <- edgelist[!apply(edgelist,1,function(x)(x[1]==x[2])),,drop=FALSE]
	} else{
		edgelist.temp <- unique(t(apply(edgelist,1,sort)))
		edgelist.final <- edgelist.temp[!apply(edgelist.temp,1,function(x)(x[1]==x[2])),,drop=FALSE]
	}
	
	add.edges(nw.new,tail=edgelist.final[,1],head=edgelist.final[,2])
	if(is.edgelist)
		as.edgelist(nw.new) else
		nw.new} 


