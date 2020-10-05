
######## Assessment of Realdata Clusters ##########

library(nloptr)
library(doMC)
registerDoMC(cores=40)
library(tools)
library(cluster)
library(clusterCrit)
library(reshape2)
library(factoextra)
library(dplyr) 
library(ggpubr)
library(reportRx)
library(xtable)


V<-PD
n.obs <- sum(V$class==2)
n.otu <- dim(V)[2]-41

test.set.idx <- 1:n.obs
gen.data.array <- V[,1:n.otu]
gen.data.array.cases <- gen.data.array[which(V$class==2),]
covfile <- V[,c((n.otu+1):ncol(V))]

clst.num = 10 # num of cluster to determine the optimal number based on SI
num.matrix = 6 # six methods to compare


#Use the distance matrices for clustering methods
l2.disct.pdf <- dist.mat.array.sum[,,1]+t(dist.mat.array.sum[,,1])
l2.disct.cdf <- dist.mat.array.sum[,,2]+t(dist.mat.array.sum[,,2])
l2.ctns.cdf <- dist.mat.array.sum[,,3]+t(dist.mat.array.sum[,,3])
dist.mat.array.eucl <- dist.mat.array.sum[,,4]+t(dist.mat.array.sum[,,4])
dist.mat.array.mht <- dist.mat.array.sum[,,5]+t(dist.mat.array.sum[,,5])
dist.mat.array.bc <- dist.mat.array.sum[,,6]+t(dist.mat.array.sum[,,6])

d.pool <- array(c(l2.disct.pdf,l2.disct.cdf,l2.ctns.cdf,dist.mat.array.eucl,
                  dist.mat.array.mht,dist.mat.array.bc), dim = c(dim(l2.disct.pdf),6) )






###############  internal indices  #################
SI <- array(0,dim=c((clst.num-1),num.matrix))
Dunn <- array(0,dim=c((clst.num-1),num.matrix))
Ball_Hall <- array(0,dim=c((clst.num-1),num.matrix))
Wemmert_Gancarski <- array(0,dim=c((clst.num-1),num.matrix)) 
Xie_Beni <- array(0,dim=c((clst.num-1),num.matrix)) #min best

Dunn.frame <- data.frame(cbind(Dunn,c(2:clst.num)))
colnames(Dunn.frame)[7]<-"Cluster.num"
Ball_Hall.frame <- data.frame(cbind(Ball_Hall,c(2:clst.num)))
colnames(Ball_Hall.frame)[7]<-"Cluster.num"
Wemmert_Gancarski.frame <- data.frame(cbind(Wemmert_Gancarski,c(2:clst.num)))
colnames(Wemmert_Gancarski.frame)[7]<-"Cluster.num"
Xie_Beni.frame <- data.frame(cbind(Xie_Beni,c(2:clst.num)))
colnames(Xie_Beni.frame)[7]<-"Cluster.num"

indices <- array(0,dim=c(num.matrix,5)) 
rownames(indices)<-c("L2.disct.pdf","L2.disct.cdf","L2.ctns.cdf","Manhattan","Euclidean","Bray-Curtis")
colnames(indices)<-c( "Dunn",
                      "SI",
                      "Ball_Hall",
                      "Wemmert_Gancarski",
                      "Xie_Beni")

for (i in 2:clst.num){
  
  print(paste0("Cluster number ", i))
  pr.disct.pdf<- pam(l2.disct.pdf, k=i,diss=TRUE,medoids = NULL,cluster.only = FALSE,
                     do.swap = TRUE)
  SI[i-1,1] <- pr.disct.pdf$silinfo$avg.width
  prclt.disct.pdf<-pr.disct.pdf$clustering
  Ball_Hall[i-1,1] <- as.numeric(intCriteria(gen.data.array.norm,prclt.disct.pdf,"Ball_Hall"))
  Wemmert_Gancarski[i-1,1] <- as.numeric(intCriteria(gen.data.array.norm,prclt.disct.pdf,"Wemmert_Gancarski"))
  Dunn[i-1,1] <- as.numeric(intCriteria(gen.data.array.norm,prclt.disct.pdf,"Dunn"))
  Xie_Beni[i-1,1] <- as.numeric(intCriteria(gen.data.array.norm,prclt.disct.pdf,"Xie_Beni"))
  
  
  pr.disct.cdf<- pam(l2.disct.cdf, k=i,diss=TRUE,medoids = NULL,cluster.only = FALSE,
                     do.swap = TRUE)
  SI[i-1,2] <- pr.disct.cdf$silinfo$avg.width
  prclt.disct.cdf<-pr.disct.cdf$clustering
  Ball_Hall[i-1,2] <- as.numeric(intCriteria(gen.data.array.norm,prclt.disct.cdf,"Ball_Hall"))
  Wemmert_Gancarski[i-1,2] <- as.numeric(intCriteria(gen.data.array.norm,prclt.disct.cdf,"Wemmert_Gancarski"))
  Dunn[i-1,2] <- as.numeric(intCriteria(gen.data.array.norm,prclt.disct.cdf,"Dunn"))
  Xie_Beni[i-1,2] <- as.numeric(intCriteria(gen.data.array.norm,prclt.disct.cdf,"Xie_Beni"))
  
  
  pr.ctns.cdf<- pam(l2.ctns.cdf, k=i,diss=T,medoids = NULL,cluster.only = FALSE,
                    do.swap = TRUE)
  SI[i-1,3] <- pr.ctns.cdf$silinfo$avg.width
  prclt.ctns.cdf<-pr.ctns.cdf$clustering
  Ball_Hall[i-1,3] <- as.numeric(intCriteria(gen.data.array.norm,prclt.ctns.cdf,"Ball_Hall"))
  Wemmert_Gancarski[i-1,3] <- as.numeric(intCriteria(gen.data.array.norm,prclt.ctns.cdf,"Wemmert_Gancarski"))
  Dunn[i-1,3] <- as.numeric(intCriteria(gen.data.array.norm,prclt.ctns.cdf,"Dunn"))
  Xie_Beni[i-1,3] <- as.numeric(intCriteria(gen.data.array.norm,prclt.ctns.cdf,"Xie_Beni"))
  
  
  pr.mht<- pam(dist.mat.array.mht, k=i,diss=T,medoids = NULL,cluster.only = FALSE,
               do.swap = TRUE)
  SI[i-1,4] <- pr.mht$silinfo$avg.width
  prclt.mht<-pr.mht$clustering
  Ball_Hall[i-1,4] <- as.numeric(intCriteria(gen.data.array.norm,prclt.mht,"Ball_Hall"))
  Wemmert_Gancarski[i-1,4] <- as.numeric(intCriteria(gen.data.array.norm,prclt.mht,"Wemmert_Gancarski"))
  Dunn[i-1,4] <- as.numeric(intCriteria(gen.data.array.norm,prclt.mht,"Dunn"))
  Xie_Beni[i-1,4] <- as.numeric(intCriteria(gen.data.array.norm,prclt.mht,"Xie_Beni"))
  
  pr.eucl<- pam(dist.mat.array.eucl, k=i,diss=T,medoids = NULL,cluster.only = FALSE,
                do.swap = TRUE)
  SI[i-1,5] <- pr.eucl$silinfo$avg.width
  prclt.eucl<-pr.eucl$clustering
  Ball_Hall[i-1,5] <- as.numeric(intCriteria(gen.data.trans,prclt.eucl,"Ball_Hall"))
  Wemmert_Gancarski[i-1,5] <- as.numeric(intCriteria(gen.data.trans,prclt.eucl,"Wemmert_Gancarski"))
  Dunn[i-1,5] <- as.numeric(intCriteria(gen.data.trans,prclt.eucl,"Dunn"))
  Xie_Beni[i-1,5] <- as.numeric(intCriteria(gen.data.trans,prclt.eucl,"Xie_Beni"))
  
  pr.bc<- pam(dist.mat.array.bc, k=i,diss=T,medoids = NULL,cluster.only = FALSE,
              do.swap = TRUE)
  SI[i-1,6] <- pr.bc$silinfo$avg.width
  prclt.bc<-pr.bc$clustering
  Ball_Hall[i-1,6] <- as.numeric(intCriteria(gen.data.trans,prclt.bc,"Ball_Hall"))
  Wemmert_Gancarski[i-1,6] <- as.numeric(intCriteria(gen.data.trans,prclt.bc,"Wemmert_Gancarski"))
  Dunn[i-1,6] <- as.numeric(intCriteria(gen.data.trans,prclt.bc,"Dunn"))
  Xie_Beni[i-1,6] <- as.numeric(intCriteria(gen.data.trans,prclt.bc,"Xie_Beni"))
}
for (a in 1:num.matrix){
  indices[a,] <- c(which.max(Dunn[,a]),
                   which.max(SI[,a]),
                   which.max(Ball_Hall[,a]),
                   which.max(Wemmert_Gancarski[,a]),
                   which.min(Xie_Beni[,a]))
}


