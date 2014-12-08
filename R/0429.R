# TODO: Add comment
# 
# Author: KirkLi
###############################################################################
library(ergm,lib=performance_path)
#xvec1 <- c("S","S","L","L")
#xvec2 <- c("L","S","L","L")
########Jun 12, 2014################
n <- 6
net1 <- network(n,density=0.3,directed=FALSE)
net2 <- network(n,density=0.3,directed=FALSE)
al1 <- sample(0:1,n,replace=T)
al2 <- sample(0:1,n,replace=T)
net1 %v% "alcohol" <- al1
net2 %v% "alcohol" <- al2

## Formation +
index.mm <- which(outer(net1 %v% "alcohol" , net1 %v% "alcohol",paste0)=="00",arr.ind = TRUE)

setdiff.mat <- function(mat1,mat2){
    if(!is.matrix(mat1))
    mat1 <- t(as.matrix(mat1))
    if(!is.matrix(mat2))
    mat2 <- t(as.matrix(mat2))
    match(setdiff(apply(mat2,1,paste,collapse="_"), apply(mat1,1,paste,collapse="_")) ,apply(mat2,1,paste,collapse="_"))
    
}

index <- setdiff.mat(mat1=as.edgelist(net1),mat2=as.edgelist(net2))

formed.edge <- as.edgelist(net2)[index,]

formed.edge.p <- formed.edge[setdiff.mat(index.mm,formed.edge),,drop=FALSE]

net.tmp <- net1
formed.net.p <- add.edges(net.tmp,formed.edge.p[,1],formed.edge.p[,2])
formed.net.p %v% "alcohol" <- net1%v%"alcohol"
net2.al <- net2%v%"alcohol"
formed.net.p.al <- formed.net.p%v%"alcohol"
formed.net.p.al[which(formed.net.p.al==1)] <- net2.al[which(formed.net.p.al==1)]
formed.net.p %v% "alcohol" <- formed.net.p.al 


## Formation -
index.pp <- which(outer(net1 %v% "alcohol" , net1 %v% "alcohol",paste0)=="11",arr.ind = TRUE)

formed.edge.m <- formed.edge[setdiff.mat(index.pp,formed.edge),,drop=FALSE]
net.tmp <- net1
formed.net.m <- add.edges(net.tmp,formed.edge.m[,1],formed.edge.m[,2])
formed.net.m %v% "alcohol" <- net1%v%"alcohol"
net2.al <- net2%v%"alcohol"
formed.net.m.al <- formed.net.m%v%"alcohol"
formed.net.m.al[which(formed.net.m.al==0)] <- net2.al[which(formed.net.m.al==0)]
formed.net.m %v% "alcohol" <- formed.net.m.al 


net.tmp <- net1
formed.net <- add.edges(net.tmp,formed.edge[,1],formed.edge[,2])
col.vec <- rep(0,nrow(as.edgelist(formed.net)))
col.vec[setdiff.mat(as.edgelist(net1),as.edgelist(formed.net))] <- 1

#plot.network(net1,displaylabels=T,label.col=net1%v%"alcohol"+2)
layout_mat <- matrix(c(1,2,5,1,3,5,1,4,5),3,3,byrow=T)
par(mar=c(1,1,1,1))
layout(layout_mat)
p1 <- plot.network(net1,vertex.col=net1%v%"alcohol"+1,vertex.cex=4,displaylabels=T)
plot.network(formed.net.p,vertex.col=formed.net.p%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T )
plot.network(formed.net,vertex.col=net2%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T,edge.col=col.vec+1)
plot.network(formed.net.m,vertex.col=formed.net.m%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T )
plot.network(net2,vertex.col=net2%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T )
#legend("topright",legend=c("-","+"),col=c(2,3),pch=19,cex=4,bty="n")



## Dissolution +
index <- setdiff.mat(mat1=as.edgelist(net2),mat2=as.edgelist(net1))
dissolved.edge <- as.edgelist(net1)[index,]
dissolved.edge.p <- dissolved.edge[setdiff.mat(index.mm,dissolved.edge),,drop=FALSE]

net.tmp <- net1
delete.edge.id.p <- match(apply(dissolved.edge.p,1,paste,collapse="_") , apply(as.edgelist(net.tmp),1,paste,collapse="_")
)
dissolved.net.p <- delete.edges(net.tmp,delete.edge.id.p)
dissolved.net.p %v% "alcohol" <- net1%v%"alcohol"
net2.al <- net2%v%"alcohol"
dissolved.net.p.al <- dissolved.net.p%v%"alcohol"
dissolved.net.p.al[which(dissolved.net.p.al==1)] <- net2.al[which(dissolved.net.p.al==1)]
dissolved.net.p %v% "alcohol" <- dissolved.net.p.al 



# Dissolution -
index.pp <- which(outer(net1 %v% "alcohol" , net1 %v% "alcohol",paste0)=="11",arr.ind = TRUE)
dissolved.edge.m <- dissolved.edge[setdiff.mat(index.pp,dissolved.edge),,drop=FALSE]
net.tmp <- net1
delete.edge.id.m <- match(apply(dissolved.edge.m,1,paste,collapse="_") , apply(as.edgelist(net.tmp),1,paste,collapse="_")
)

dissolved.net.m <- delete.edges(net.tmp,delete.edge.id.m)
dissolved.net.m %v% "alcohol" <- net1%v%"alcohol"
net2.al <- net2%v%"alcohol"
dissolved.net.m.al <- dissolved.net.m%v%"alcohol"
dissolved.net.m.al[which(dissolved.net.m.al==0)] <- net2.al[which(dissolved.net.m.al==0)]
dissolved.net.m %v% "alcohol" <- dissolved.net.m.al 
#

net.tmp <- net1
delete.edge.id<- match(apply(dissolved.edge,1,paste,collapse="_") , apply(as.edgelist(net.tmp),1,paste,collapse="_")
)

dissolved.net <- delete.edges(net.tmp,delete.edge.id)
col.vec <- rep(0,nrow(as.edgelist(net1)))
col.vec[setdiff.mat(as.edgelist(dissolved.net), as.edgelist(net1))] <- 1




















#plot.network(net1,displaylabels=T,label.col=net1%v%"alcohol"+2)
layout_mat <- matrix(c(1,2,5,1,3,5,1,4,5),3,3,byrow=T)
par(mar=c(1,1,1,1))
layout(layout_mat)
p1 <- plot.network(net1,vertex.col=net1%v%"alcohol"+1,vertex.cex=4,displaylabels=T)
plot.network(dissolved.net.p,vertex.col=dissolved.net.p%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T )
plot.network(net1,vertex.col=net2%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T,edge.col=col.vec+1)
plot.network(dissolved.net.m,vertex.col=dissolved.net.m%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T )
plot.network(net2,vertex.col=net2%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T )
#legend("topright",legend=c("-","+"),col=c(2,3),pch=19,cex=4,bty="n")




#layout_mat <- matrix(c(1,2,6,1,3,6,1,4,6,1,5,6),4,3,byrow=T)
#par(oma=c(4,4,4,4),mar=c(5,5,5,5))
#layout(layout_mat)
pdf("fdplot.pdf")

p1 <- plot.network(net1,vertex.col=net1%v%"alcohol"+1,vertex.cex=4,displaylabels=T,label.cex=2)
plot.network(formed.net.p,vertex.col=formed.net.p%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T,label.cex=2)
#plot.network(formed.net,vertex.col=net2%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T,edge.col=col.vec+1)
plot.network(formed.net.m,vertex.col=formed.net.m%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T ,label.cex=2)
plot.network(dissolved.net.p,vertex.col=dissolved.net.p%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T ,label.cex=2)
#plot.network(net1,vertex.col=net2%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T,edge.col=col.vec+1)
plot.network(dissolved.net.m,vertex.col=dissolved.net.m%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T ,label.cex=2)
plot.network(net2,vertex.col=net2%v%"alcohol"+1,coord=p1,vertex.cex=4,displaylabels=T ,label.cex=2)

dev.off()








########Jun 12, 2014################
net1 <-as.network(read.table("s50_data/s50-network1.dat"))
net2 <-as.network(read.table("s50_data/s50-network2.dat"))
net3 <-as.network(read.table("s50_data/s50-network3.dat"))
alcohol <- read.table("s50_data/s50-alcohol.dat")
drugs <- read.table("s50_data/s50-drugs.dat")
net <- net1

d2ud <- function(net){
    for (i in 1:network.size(net)){
        for(j in 1:network.size(net)){
            net[i,j] = max(net[i,j],net[j,i])
        }}
    net$gal$directed=FALSE
    net
}

net1 <- d2ud(net1)
net2 <- d2ud(net2)
net3 <- d2ud(net3)



al1 <- as.numeric(alcohol[,1]>3)
al2 <- as.numeric(alcohol[,2]>3)
net1 %v% "alcohol" <- al1
net2 %v% "alcohol" <- al2

plot(net1)
plot(net2)













########Jun 9, 2014########randomized test########
# convert al t-1 and t to matrix
al1 <- alcohol[,2]
table(al1)
al1 <- as.numeric(al1>=3)
al1 <- outer(al1,al1,FUN="==")*1 # 1 is homophily and 0 is not

al2 <- alcohol[,3]
table(al2)
al2 <- as.numeric(al2>=3)
al2 <- outer(al2,al2,FUN="==")*1 # 1 is homophily and 0 is not

# 1
t1 <- matrix (c(sum(as.numeric(al1)==0 & as.numeric(al2)==1 & as.numeric(as.matrix(net1))==1),
sum(as.numeric(al1)==0 & as.numeric(al2)==1 & as.numeric(as.matrix(net1))==0),
sum(as.numeric(al1)==1 & as.numeric(al2)==0 & as.numeric(as.matrix(net1))==1),
sum(as.numeric(al1)==1 & as.numeric(al2)==0 & as.numeric(as.matrix(net1))==0)),2,2,byrow=T)

# 2
t2 <- matrix (c(sum(as.numeric(al1)==0 & as.numeric(al2)==0 & as.numeric(as.matrix(net1))==1),
sum(as.numeric(al1)==0 & as.numeric(al2)==0 & as.numeric(as.matrix(net1))==0),
sum(as.numeric(al1)==1 & as.numeric(al2)==1 & as.numeric(as.matrix(net1))==1),
sum(as.numeric(al1)==1 & as.numeric(al2)==1 & as.numeric(as.matrix(net1))==0)),2,2,byrow=T)


# 3
t3 <- matrix (c(sum(as.numeric(al1)==0 & as.numeric(al2)==1 & as.numeric(as.matrix(net1))==1),
                sum(as.numeric(al1)==0 & as.numeric(al2)==1 & as.numeric(as.matrix(net1))==0),
                sum(as.numeric(al1)==0 & as.numeric(al2)==0 & as.numeric(as.matrix(net1))==1),
                sum(as.numeric(al1)==0 & as.numeric(al2)==0 & as.numeric(as.matrix(net1))==0)),2,2,byrow=T)


# 4
(t4 <- matrix (c(sum(as.numeric(al1)==1 & as.numeric(al2)==0 & as.numeric(as.matrix(net1))==1),
                sum(as.numeric(al1)==1 & as.numeric(al2)==0 & as.numeric(as.matrix(net1))==0),
                sum(as.numeric(al1)==1 & as.numeric(al2)==1 & as.numeric(as.matrix(net1))==1),
                sum(as.numeric(al1)==1 & as.numeric(al2)==1 & as.numeric(as.matrix(net1))==0)),2,2,byrow=T))

library(sna)
# 5
(t5 <- matrix (c(sum(alcohol[,1]!=alcohol[,2] & degree(net1,gmode="graph")>0),
                 sum(alcohol[,1]!=alcohol[,2] & degree(net1,gmode="graph")==0),
                 sum(alcohol[,1]==alcohol[,2] & degree(net1,gmode="graph")>0),
                 sum(alcohol[,1]==alcohol[,2] & degree(net1,gmode="graph")==0)),2,2,byrow=T))


# 6
(t6 <- matrix (c(
    sum(
                    as.numeric(as.matrix(net1))==0 & 
                    as.numeric(as.matrix(net2))==1 &
                    as.numeric(al1)==1),

    sum(
                    as.numeric(as.matrix(net1))==0 & 
                    as.numeric(as.matrix(net2))==1 &
                    as.numeric(al1)==0),
    sum(
                    as.numeric(as.matrix(net1))==1 & 
                    as.numeric(as.matrix(net2))==0 &
                    as.numeric(al1)==1),
    sum(
                    as.numeric(as.matrix(net1))==1 & 
                    as.numeric(as.matrix(net2))==0 &
                    as.numeric(al1)==0)),2,2,byrow=T))
    
# 7
(t7 <- matrix (c(
                            sum(
                                    as.numeric(as.matrix(net1))==0 & 
                                            as.numeric(as.matrix(net2))==0 &
                                            as.numeric(al1)==1),
                            
                            sum(
                                    as.numeric(as.matrix(net1))==0 & 
                                            as.numeric(as.matrix(net2))==0 &
                                            as.numeric(al1)==0),
                            sum(
                                    as.numeric(as.matrix(net1))==1 & 
                                            as.numeric(as.matrix(net2))==1 &
                                            as.numeric(al1)==1),
                            sum(
                                    as.numeric(as.matrix(net1))==1 & 
                                            as.numeric(as.matrix(net2))==1 &
                                            as.numeric(al1)==0)),2,2,byrow=T))

# 8
(t8 <- matrix (c(
                            sum(
                                    as.numeric(as.matrix(net1))==0 & 
                                            as.numeric(as.matrix(net2))==1 &
                                            as.numeric(al1)==1),
                            
                            sum(
                                    as.numeric(as.matrix(net1))==0 & 
                                            as.numeric(as.matrix(net2))==1 &
                                            as.numeric(al1)==0),
                            sum(
                                    as.numeric(as.matrix(net1))==0 & 
                                            as.numeric(as.matrix(net2))==0 &
                                            as.numeric(al1)==1),
                            sum(
                                    as.numeric(as.matrix(net1))==0 & 
                                            as.numeric(as.matrix(net2))==0 &
                                            as.numeric(al1)==0)),2,2,byrow=T))


# 9
(t9 <- matrix (c(
                            sum(
                                    as.numeric(as.matrix(net1))==1 & 
                                            as.numeric(as.matrix(net2))==0 &
                                            as.numeric(al1)==1),
                            
                            sum(
                                    as.numeric(as.matrix(net1))==1 & 
                                            as.numeric(as.matrix(net2))==0 &
                                            as.numeric(al1)==0),
                            sum(
                                    as.numeric(as.matrix(net1))==1 & 
                                            as.numeric(as.matrix(net2))==1 &
                                            as.numeric(al1)==1),
                            sum(
                                    as.numeric(as.matrix(net1))==1 & 
                                            as.numeric(as.matrix(net2))==1 &
                                            as.numeric(al1)==0)),2,2,byrow=T))



# 10
(t10 <- matrix (c(
                            sum(
                                            as.numeric(al1)==0 &
                                            as.numeric(al2)==1 &
                                            as.numeric(as.matrix(net1))==0 & 
                                            as.numeric(as.matrix(net2))==1) ,
                            sum(
                                    as.numeric(al1)==0 &
                                            as.numeric(al2)==1 &
                                            as.numeric(as.matrix(net1))==1 & 
                                            as.numeric(as.matrix(net2))==0) ,
                            sum(
                                    as.numeric(al1)==1 &
                                            as.numeric(al2)==0 &
                                            as.numeric(as.matrix(net1))==0 & 
                                            as.numeric(as.matrix(net2))==1) ,
                            sum(
                                    as.numeric(al1)==1 &
                                            as.numeric(al2)==0 &
                                            as.numeric(as.matrix(net1))==1 & 
                                            as.numeric(as.matrix(net2))==0))              
                            ,2,2,byrow=T))
        
        
# 11
        (t11 <- matrix (c(
                                    sum(
                                            as.numeric(al1)==0 &
                                                    as.numeric(al2)==0 &
                                                    as.numeric(as.matrix(net1))==0 & 
                                                    as.numeric(as.matrix(net2))==0) ,
                                    sum(
                                            as.numeric(al1)==0 &
                                                    as.numeric(al2)==0 &
                                                    as.numeric(as.matrix(net1))==1 & 
                                                    as.numeric(as.matrix(net2))==1) ,
                                    sum(
                                            as.numeric(al1)==1 &
                                                    as.numeric(al2)==1 &
                                                    as.numeric(as.matrix(net1))==0 & 
                                                    as.numeric(as.matrix(net2))==0) ,
                                    sum(
                                            as.numeric(al1)==1 &
                                                    as.numeric(al2)==1 &
                                                    as.numeric(as.matrix(net1))==1 & 
                                                    as.numeric(as.matrix(net2))==1))              
                            ,2,2,byrow=T))

        
        
# 12
        (t12 <- matrix (c(
                                    sum(
                                            as.numeric(al1)==0 &
                                                    as.numeric(al2)==1 &
                                                    as.numeric(as.matrix(net1))==0 & 
                                                    as.numeric(as.matrix(net2))==1) ,
                                    sum(
                                            as.numeric(al1)==0 &
                                                    as.numeric(al2)==1 &
                                                    as.numeric(as.matrix(net1))==0 & 
                                                    as.numeric(as.matrix(net2))==0) ,
                                    sum(
                                            as.numeric(al1)==0 &
                                                    as.numeric(al2)==0 &
                                                    as.numeric(as.matrix(net1))==0 & 
                                                    as.numeric(as.matrix(net2))==1) ,
                                    sum(
                                            as.numeric(al1)==0 &
                                                    as.numeric(al2)==0 &
                                                    as.numeric(as.matrix(net1))==0 & 
                                                    as.numeric(as.matrix(net2))==0))              
                            ,2,2,byrow=T))
        
        
# 13
        (t13 <- matrix (c(
                                    sum(
                                            as.numeric(al1)==1 &
                                                    as.numeric(al2)==0 &
                                                    as.numeric(as.matrix(net1))==1 & 
                                                    as.numeric(as.matrix(net2))==0) ,
                                    sum(
                                            as.numeric(al1)==1 &
                                                    as.numeric(al2)==0 &
                                                    as.numeric(as.matrix(net1))==1 & 
                                                    as.numeric(as.matrix(net2))==1) ,
                                    sum(
                                            as.numeric(al1)==1 &
                                                    as.numeric(al2)==1 &
                                                    as.numeric(as.matrix(net1))==1 & 
                                                    as.numeric(as.matrix(net2))==0) ,
                                    sum(
                                            as.numeric(al1)==1 &
                                                    as.numeric(al2)==1 &
                                                    as.numeric(as.matrix(net1))==1 & 
                                                    as.numeric(as.matrix(net2))==1))              
                            ,2,2,byrow=T))


        

        
t1
t2
t3
t4
t5
t6
t7
t8
t9
t10
t11
t12
t13



chisq.test(t1)$p.value
chisq.test(t2)$p.value
chisq.test(t3)$p.value
chisq.test(t4)$p.value
chisq.test(t5)$p.value
chisq.test(t6)$p.value
chisq.test(t7)$p.value
chisq.test(t8)$p.value
chisq.test(t9)$p.value
chisq.test(t10)$p.value
chisq.test(t11)$p.value
chisq.test(t12)$p.value
chisq.test(t13)$p.value

########Jun 9, 2014########showing network auto-correlation using graphs########
library(plot3D)
scatter3D(rep(1:3,each=50),rep(1:50,3),c(drugs[,1],drugs[,2],drugs[,3]),xlab="wave",ylab="people",zlab="drugs",phi = 20, theta =20,colkey=FALSE,bty="b2")

for (i in 1:50)
	scatter3D(1:3,rep(i,3),c(drugs[i,1],drugs[i,2],drugs[i,3]),xlab="wave",ylab="people",zlab="drugs",phi = 20, theta =20,add=T	,type="l",colkey=FALSE,bty="b2")



scatter3D(rep(1:3,each=50),rep(1:50,3),c(alcohol[,1],alcohol[,2],alcohol[,3]),xlab="wave",ylab="people",zlab="alcohol",phi = 20, theta =20,colkey=FALSE,bty="b2")

for (i in 1:50)
scatter3D(1:3,rep(i,3),c(alcohol[i,1],alcohol[i,2],alcohol[i,3]),xlab="wave",ylab="people",zlab="alcohol",phi = 20, theta =20,add=T	,type="l",colkey=FALSE,bty="b2")




matplot(apply(drugs,1,jitter),type="l",ylab="durg level",xlab="wave")
matplot(apply(alcohol,1,jitter),type="l",ylab="alcohol level",xlab="wave")

pdf("figures/net1.pdf")
p1 <- plot(net1)
dev.off()

pdf("figures/net2.pdf")
p2 <- plot(net2,coord=p1)
dev.off()

pdf("figures/net3.pdf")
p3 <- plot(net3,coord=p1)
dev.off()


