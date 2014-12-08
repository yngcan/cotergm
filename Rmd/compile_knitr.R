# TODO: Add comment
# 
# Author: kirk
###############################################################################



library(knitr)
setwd("Rmd")
#knit2html("cotergm.Rmd")
#purl("cotergm.Rmd")


knit2html("cotergm_school3.Rmd")


knit2html("cotergm_school2.Rmd")


knit2html("cotergm_school.Rmd")


purl("cotergm_school.Rmd")

#  scp /home/kirk/Dropbox/hiv_network/workspace/cotergm/Rmd/cotergm_school.html kirkli@madrid.stat.washington.edu:../../var/www/stat/people/kirkli/hivnetwork/cotergm_school.html













```{r}
# pdf(paste0(fig.path,"0922.pdf"),height=10,width=10)
		layout_mat <- matrix(c(7,2,8,1,3,6,1,4,6,0,5,0),4,3,byrow=T)
		par(mar=c(2,2,3,2))
		layout(layout_mat)
#
		p1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,main="y0",cex.main=2,displaylabels=T,cex.lab=1)
		
		p1.1 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="F+",cex.main=2,displaylabels=T,cex.lab=1)
		plot(FP,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=3)
#
#
		p2 <- plot(y0,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="F-",cex.main=2,displaylabels=T,cex.lab=1)
		plot(FM,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=3)
		
		p3 <- plot(y0,edge.col="yellow",edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="D+",cex.main=2,displaylabels=T,cex.lab=1)
		plot(DP,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=1)
		
		p4 <- plot(y0,edge.col="yellow",edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="D-",cex.main=2,displaylabels=T,cex.lab=1)
		plot(DM,coord=p1,new=FALSE,vertex.cex=2,vertex.col="status",cex.lab=2,edge.col=1)
		
		p5 <- plot(nw,edge.col=1,edge.lwd=5,vertex.cex=2,vertex.col="status",cex.lab=2,coord=p1,main="y1",cex.main=2,displaylabels=T,cex.lab=1,vertex.sides=ifelse(y0 %v% "status" == nw %v% "status",50,3))
#
		p7 <- plot.new()
		legend("topleft",legend=c("+","-"),pch=19,col=c(2,1),cex=2,bty="n")
		
		p8 <- plot.new()
		legend("topleft",legend=c("tie formed","tie dissolved","nodal status changed"),lty=c(1,1,NA),pch=c(NA,NA,24), col=c(3,"yellow",1),cex=2,bty="n")
#
#dev.off()
		```