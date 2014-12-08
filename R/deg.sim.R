degen <- function(network.size=5,stat=c("edges","triangle"),net0=network(5,density=0.4,directed=FALSE)){
	library(network)
	library(ergm)
	library(parallel)
	VEC <- strsplit(outer(1:network.size,1:network.size,paste)[upper.tri(outer(1:network.size,1:network.size,paste))]," ")
	VEC <- lapply(VEC,as.numeric)
	VEC <- do.call(rbind, VEC)

  # index of edges in original network
	edge0 <- intersect.mat(as.edgelist(net0),VEC)
	
  # all possible combinations of edge index
	index.list <- sapply(0:(2^choose(network.size,2)-1),function(x){ which((head(as.numeric(unlist(strsplit(paste(intToBits(x),collapse=""),""))),choose(network.size,2)*2+2)[seq(2,choose(network.size,2)*2+2,2)])==1)})
	

  attr.mat <- outer(as.numeric(net0 %v% "alcohol"),as.numeric(net0 %v% "alcohol"),FUN="+")
  # edgelist of ++
	pp.edge <- which(attr.mat==2,arr.ind = TRUE)  
  # edgelist of --
  mm.edge <- which(attr.mat==0,arr.ind = TRUE)  
  
  # edge index of ++
  edgepp <- intersect.mat(pp.edge,VEC) 
  
  # index of potential formed ++ edges
  edgepp.can.form <- setdiff(edgepp,edge0)
  
  # index of potential dissolved ++ edges 
  edgepp.can.diss <- intersect(edge0,edgepp)
  # edge index of --
  edgemm <- intersect.mat(mm.edge,VEC) 
  
  # index of potential formed -- edges
  edgemm.can.form <- setdiff(edgemm,edge0)
  
  # index of potential dissolved ++ edges 
  edgemm.can.diss <- intersect(edge0,edgemm)
  
   # whether is a candidate F+ network
   index.form.p <- unlist(lapply(index.list,function(x){all(edge0 %in% x, !(edgemm.can.form %in% x))}))
   
   # whether is a candidate F- network
   index.form.m <- unlist(lapply(index.list,function(x){all(edge0 %in% x, !(edgepp.can.form %in% x))}))
   
   
   index.diss.p <- unlist(lapply(index.list,function(x){all(x %in% edge0,edgemm.can.diss %in% x)}))
   
   index.diss.m <- unlist(lapply(index.list,function(x){all(x %in% edge0,edgepp.can.diss %in% x)}))
    
  
	parfn <- function(x,VEC, stat,network.size, net0){
		require(ergm)
		index <- x
		if(length(VEC[index,,drop=F])==0)
			net <- network.initialize(network.size,directed=FALSE) else 
			net <- as.network(VEC[index,,drop=F],directed=FALSE,matrix.type="edgelist")
  		net %v% "alcohol" <- net0  %v% "alcohol"
		summary.stat <- function(net,stat)
		{
			formula <- paste("net~", paste(stat,sep="",collapse = "+"),sep="")
			eval(parse(text=paste("summary(",formula,")")))
		}
		summary.stat(net,stat)
	}

  output.fn <- function(LIST,stat){
      POINT <- do.call(rbind,LIST)
      p.dat <- data.frame(POINT)
      unique.point.list <- split(x = p.dat, f = p.dat)
      size.list <- unlist(lapply(unique.point.list,nrow))
      TEMP <- cbind(do.call(rbind,strsplit(names(size.list),".",fixed=T)),size.list)  
      TEMP <- apply(TEMP, 2, as.numeric)
      TEMP <- TEMP[order(TEMP[,1],TEMP[,2]),]
      TEMP <- TEMP[TEMP[,3]!=0,]
      TEMP <- data.frame(TEMP)
      names(TEMP) <- c(stat,"count")
      TEMP
  }
  
	cl1 <- makeCluster(detectCores())
	LIST.form.p <- parLapply(cl1,index.list[index.form.p],parfn,VEC,stat,network.size, net0)
	#   LIST <- parLapply(cl1,0:10,parfn,VEC,stat,network.size)
	# ss <- lapply(LIST,function(x)x$stat)
  TEMP.form.p <- output.fn(LIST.form.p,stat)
  
  LIST.form.m <- parLapply(cl1,index.list[index.form.m],parfn,VEC,stat,network.size, net0)
  #   LIST <- parLapply(cl1,0:10,parfn,VEC,stat,network.size)
  # ss <- lapply(LIST,function(x)x$stat)
  TEMP.form.m <- output.fn(LIST.form.m,stat)
  
  LIST.diss.p <- parLapply(cl1,index.list[index.diss.p],parfn,VEC,stat,network.size, net0)
  #   LIST <- parLapply(cl1,0:10,parfn,VEC,stat,network.size)
  # ss <- lapply(LIST,function(x)x$stat)
  TEMP.diss.p <- output.fn(LIST.diss.p,stat)
  
  LIST.diss.m <- parLapply(cl1,index.list[index.diss.m],parfn,VEC,stat,network.size, net0)
  #   LIST <- parLapply(cl1,0:10,parfn,VEC,stat,network.size)
  # ss <- lapply(LIST,function(x)x$stat)
  TEMP.diss.m <- output.fn(LIST.diss.m,stat)
#  
#	save(LIST,file=paste(paste(stat,collapse="_"),"_LIST_form.RData",sep=""))
#	save(TEMP,file=paste(paste(stat,collapse="_"),"_TEMP_form.RData",sep=""))

	return(list(
                  TEMP.form.p=TEMP.form.p,
                  TEMP.form.m=TEMP.form.m,
                  TEMP.diss.p=TEMP.diss.p,
                  TEMP.diss.m=TEMP.diss.m))
  
}

library(statnet)

net0=network(5,density=0.4,directed=FALSE)
net0 %v% "alcohol" <- sample(c(0,1),size=network.size(net0),replace=TRUE)

summary(net0~edges+nodematch("alcohol"))


RES <- degen(network.size=5,stat=c("edges","nodematch('alcohol')"),net0)


load("edges_nodematch('alcohol')_TEMP_form.RData")
TEMP



source("C:/Dropbox/hiv_network/workspace/degeneracy/R/functions.R", echo=FALSE, encoding="GBK")


RES2 <- myfun(RES$TEMP.form.m,name=c("edge","hom"))


RES2$statmat






?nodematch