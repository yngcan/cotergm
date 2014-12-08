# FORTESS MONITORING 
# 
# Author: kirk
###############################################################################
library(xtable)
library(zoo)
library(forecast)
source("functions.R")
library("statnet.common",lib=cotergm)
library("network",lib=cotergm)
library("ergm",lib=cotergm)
library("tergm",lib=cotergm)
################################################



aaa <- function(x){
	vec <- NULL
	for(i1 in 1:2) for(j1 in 1:2) for(k1 in 0:1) 
				vec <-rbind(vec,c(i1,j1,k1))
	ll <- list()
	
	mat <- apply(vec,1,function(xx){
				show(c(x,xx))
			})	
	mat <- matrix(rowMeans(mat,na.rm=T),4,4,byrow=F)
	effstring <- c("edges", "plus", "nodematch.status.1", "nodematch.status.2")
	colnames(mat) <- c("FP,SS+","FM,SS-","DP,SI+","DM,SI-")
	rownames(mat) <- effstring
	signif(mat,2)
}

show <- function(x){ 
	net1 <- network.initialize(2,directed=FALSE)
	if(x[3])
		add.edges(net1,tail=c(1),head=c(2))
	net1%v%"status" <- c(x[1],x[2])
	net2 <- network.initialize(2,directed=FALSE)
	if(x[6])
		add.edges(net2,tail=c(1),head=c(2))
	net2%v%"status" <- c(x[4],x[5])
	y0 <- net1
	nw <- net2
	## ------------------------------------------------------------------------
# Formation Plus Network
	FP <- FP.nw(y0,nw,"status")
# Formation Minus Network
	FM <- FM.nw(y0,nw,"status")
# Dissolution Plus Network
	DP <- DP.nw(y0,nw,"status")
# Dissolution Minus Network
	DM <- DM.nw(y0,nw,"status")
	constraints=~atleastnonplus(y0)
	try(
			MHproposal <- MHproposal(constraints,nw=y0,arguments=list()),silent=T
	)
	constraints=~atleastnonminus(y0)
	try(
			MHproposal <- MHproposal(constraints,nw=y0,arguments=list()),silent=T
	)
	
	effstring <- c("edges", "plus", "nodematch.status.1", "nodematch.status.2")
	mat <- cbind(
			summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=1),MHproposal=MHproposal)[effstring],
			summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=2),MHproposal=MHproposal)[effstring],
	
			summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=1),MHproposal=MHproposal)[effstring],
			
			summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=2),MHproposal=MHproposal)[effstring])
	colnames(mat) <- c("FP,SS+","FM,SS-","DP,SI+","DM,SI-")
	rownames(mat) <- effstring
	signif(mat,2)
}


shinyServer(
		function(input, output,session){
			
			res <- reactive({
							x <- as.numeric(c(input$net1node1,input$net1node2,input$net1edge,input$net2node1,input$net2node2,input$net2edge))
							show(x)
					})
			
			res2 <- reactive({
						x <- as.numeric(c(input$net1node1,input$net1node2,input$net1edge,input$net2node1,input$net2node2,input$net2edge))
						aaa(x[c(1,2,3)])
					})
			
			output$plot <- renderPlot({
						
						x <- as.numeric(c(input$net1node1,input$net1node2,input$net1edge,input$net2node1,input$net2node2,input$net2edge))
						
						net1 <- network.initialize(2,directed=FALSE)
						if(x[3])
							add.edges(net1,tail=c(1),head=c(2))
						net1%v%"status" <- c(x[1],x[2])
						net2 <- network.initialize(2,directed=FALSE)
						if(x[6])
							add.edges(net2,tail=c(1),head=c(2))
						net2%v%"status" <- c(x[4],x[5])
						y0 <- net1
						nw <- net2
# Formation Plus Network
						FP <- FP.nw(y0,nw,"status")
# Formation Minus Network
						FM <- FM.nw(y0,nw,"status")
# Dissolution Plus Network
						DP <- DP.nw(y0,nw,"status")
# Dissolution Minus Network
						DM <- DM.nw(y0,nw,"status")
						
						net.base <- network.initialize(2)
						add.edges(net.base,c(1),c(2))
						p0 <- plot(net.base)

						layout_mat <- matrix(c(7,2,8,1,3,6,1,4,6,0,5,0),4,3,byrow=T)
						par(mar=c(2,2,3,2))
						layout(layout_mat)
						
						p1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=4,vertex.col="status",cex.lab=2,coord=p0,main="y0",cex.main=2,displaylabels=T,cex.lab=2)
						
						p1.1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=4,vertex.col="status",cex.lab=2,coord=p0,main="F+",cex.main=2,displaylabels=T,cex.lab=2,jitter=FALSE)
						plot(FP-y0,coord=p0,new=FALSE,vertex.cex=4,vertex.col=get.vertex.attribute(FP,"status"),cex.lab=2,edge.col=3,jitter=FALSE)
#
#
						p2 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=4,vertex.col="status",cex.lab=2,coord=p0,main="F-",cex.main=2,displaylabels=T,cex.lab=2,jitter=FALSE)
						plot(FM-y0,coord=p0,new=FALSE,vertex.cex=4,vertex.col=get.vertex.attribute(FM,"status"),cex.lab=2,edge.col=3,jitter=FALSE)
						
						p3 <- plot(y0-DP,edge.col="yellow",edge.lwd=5,vertex.cex=4,vertex.col=get.vertex.attribute(DP,"status"),cex.lab=2,coord=p0,main="D+",cex.main=2,displaylabels=T,cex.lab=2,jitter=FALSE)
						plot(DP,coord=p0,new=FALSE,vertex.cex=4,vertex.col="status",cex.lab=2,edge.col=1,jitter=FALSE)
						
						p4 <- plot(y0-DM,edge.col="yellow",edge.lwd=5,vertex.cex=4,vertex.col=get.vertex.attribute(DM,"status"),cex.lab=2,coord=p0,main="D-",cex.main=2,displaylabels=T,cex.lab=2,jitter=FALSE)
						plot(DM,coord=p0,new=FALSE,vertex.cex=4,vertex.col="status",cex.lab=2,edge.col=1,jitter=FALSE)
						
						p5 <- plot(nw,edge.col=1,edge.lwd=5,vertex.cex=4,vertex.col="status",cex.lab=2,coord=p0,main="y1",cex.main=2,displaylabels=T,cex.lab=2,vertex.sides=ifelse(y0 %v% "status" == nw %v% "status",50,3),jitter=FALSE)
#
						p7 <- plot.new()
						legend("topleft",legend=c("+","-"),pch=19,col=c(2,1),cex=4,bty="n")
						
						p8 <- plot.new()
						legend("topleft",legend=c("tie formed","tie dissolved","nodal status changed"),lty=c(1,1,NA),pch=c(NA,NA,24), col=c(3,"yellow",1),cex=2,bty="n")
						
					})
			output$test <- renderPrint({
						list(observed=res(),random=res2(),diff=res()-res2())
					})
			
		})