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



show <- function(net1,net2){ 
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
MHproposal <- MHproposal(constraints,nw=y0,arguments=list())
)


constraints=~atleastnonminus(y0)
try(
		MHproposal <- MHproposal(constraints,nw=y0,arguments=list())
)

effstring <- c("edges", "plus", "nodematch.status.1", "nodematch.status.2")
mat <- cbind(
# selection +
summary(FP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=1),MHproposal=MHproposal)[effstring],
# selection -
summary(FM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0edge=2,keep=2),MHproposal=MHproposal)[effstring],
# influence - 
summary(DM~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=2),MHproposal=MHproposal)[effstring],
# influence + 
summary(DP~edges+plus+nodematch_cotergm("status",diff=TRUE,y0nodal=2,keep=1),MHproposal=MHproposal)[effstring])

colnames(mat) <- c("FP,SS+","FM,SS-","DP,SI+","DP,SI-")
rownames(mat) <- effstring
mat

}








vec <- NULL
for(i1 in 1:2) for(j1 in 1:2) for(k1 in 0:1) for(i2 in 1:2) for(j2 in 1:2) for(k2 in 0:1) 
						vec <-rbind(vec,c(i1,j1,k1,i2,j2,k2))

names(lis) <- apply(vec,1,paste,collapse=",",sep=",")
a <- data.frame(do.call(rbind,strsplit(names(lis),",")))
aa <- split(a[,4:6],a[,1:3])
namevec <- do.call(rbind,strsplit(names(aa),split="[:.:]"))
ll <- list()
for (i in 1:length(aa)){
	tt	<- aa[[i]]
	mat <- apply(tt,1,function(x){
show(c(as.numeric(namevec[i,]),as.numeric(x)))
				})	
		ll[[i]] <-matrix(rowMeans(mat,na.rm=T),4,4,byrow=F)
	}
ll			




aaa <- function(x){
	vec <- NULL
	for(i1 in 1:2) for(j1 in 1:2) for(k1 in 0:1) 
				vec <-rbind(vec,c(i1,j1,k1))
	ll <- list()
	
	mat <- apply(vec,1,function(xx){
				show(c(x,xx))
			})	
	matrix(rowMeans(mat,na.rm=T),4,4,byrow=F)
}






