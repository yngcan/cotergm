# TODO: Add comment
# 
# Author: kirk
###############################################################################

# Nov 11, 2014

library("statnet.common",lib=cotergm)
library("network",lib=cotergm)
library("ergm",lib=cotergm)
library("tergm",lib=cotergm)
fls <-list.files("../../",recursive=T)
source(paste0("../../",fls[grep("kk_utils/utils.R",fls)]))
source(paste0("../../",fls[grep("cotergm/R/functions.R",fls)]))
fig.path = "figures/"


## ------------------------------------------------------------------------

net1 <- network.initialize(4,directed=FALSE)
add.edges(net1,tail=c(3,1),head=c(2,3))
net1%v%"status" <- c(2,2,1,2)

net2 <- network.initialize(4,directed=FALSE)
add.edges(net2,tail=c(2,4),head=c(3,1))
net2%v%"status" <- c(2,1,1,2)



y0 <- net1
nw <- net2
as.matrix(nw)
as.matrix(DP)
as.matrix(y0)
## ------------------------------------------------------------------------
# Formation Plus Network
FP <- FP.nw(y0,nw,"status")
# Formation Minus Network
FM <- FM.nw(y0,nw,"status")
# Dissolution Plus Network
DP <- DP.nw(y0,nw,"status")
# Dissolution Minus Network
DM <- DM.nw(y0,nw,"status")

as.matrix(DP)


net.base <- network.initialize(4)
add.edges(net.base,c(1,2,3,4),c(2,3,4,1))
p0 <- plot(net.base)

layout_mat <- matrix(c(7,2,8,1,3,6,1,4,6,0,5,0),4,3,byrow=T)
par(mar=c(2,2,3,2))
layout(layout_mat)


p1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=4,vertex.col="status",cex.lab=2,coord=p0,main="y0",cex.main=2,displaylabels=T,cex.lab=2)

p1.1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=4,vertex.col="status",cex.lab=2,coord=p0,main="F+",cex.main=2,displaylabels=T,cex.lab=2)
plot(FP-y0,coord=p0,new=FALSE,vertex.cex=4,vertex.col=get.vertex.attribute(FP,"status"),cex.lab=2,edge.col=3)
#
#
p2 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=4,vertex.col="status",cex.lab=2,coord=p0,main="F-",cex.main=2,displaylabels=T,cex.lab=2)
plot(FM-y0,coord=p0,new=FALSE,vertex.cex=4,vertex.col=get.vertex.attribute(FM,"status"),cex.lab=2,edge.col=3)

p3 <- plot(y0-DP,edge.col="yellow",edge.lwd=5,vertex.cex=4,vertex.col=get.vertex.attribute(DP,"status"),cex.lab=2,coord=p0,main="D+",cex.main=2,displaylabels=T,cex.lab=2)
plot(DP,coord=p0,new=FALSE,vertex.cex=4,vertex.col="status",cex.lab=2,edge.col=1)

p4 <- plot(y0-DM,edge.col="yellow",edge.lwd=5,vertex.cex=4,vertex.col=get.vertex.attribute(DM,"status"),cex.lab=2,coord=p0,main="D-",cex.main=2,displaylabels=T,cex.lab=2)
plot(DM,coord=p0,new=FALSE,vertex.cex=4,vertex.col="status",cex.lab=2,edge.col=1)

p5 <- plot(nw,edge.col=1,edge.lwd=5,vertex.cex=4,vertex.col="status",cex.lab=2,coord=p0,main="y1",cex.main=2,displaylabels=T,cex.lab=2,vertex.sides=ifelse(y0 %v% "status" == nw %v% "status",50,3))
#
p7 <- plot.new()
legend("topleft",legend=c("+","-"),pch=19,col=c(2,1),cex=2,bty="n")

p8 <- plot.new()
legend("topleft",legend=c("tie formed","tie dissolved","nodal status changed"),lty=c(1,1,NA),pch=c(NA,NA,24), col=c(3,"yellow",1),cex=2,bty="n")


constraints=~atleastnonminus(y0)

MHproposal <- MHproposal(constraints,nw=y0,arguments=list())

mat <- 
rbind(

# influence + 
summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=1),MHproposal=MHproposal)
,
# influence - 
summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=2),MHproposal=MHproposal)
,
# selection +
summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=2),MHproposal=MHproposal)
,
# selection -
summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=1),MHproposal=MHproposal)
)

rownames(mat) <- c("influence + to -","influence - to + ","selection +"," selection -")

colnames(mat) <- c("edges","plus","count of social effect")


mat