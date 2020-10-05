library(gtools)

accu.mat.c2 <- function(dist.mat, cls.label1, cls.label2, sample.size, pred.clst.num){
  confus.mat <- array(0,dim=c(pred.clst.num,2))
  lst = permutations(n = pred.clst.num, r = 2, repeats.allowed = F) #all combination with no repeats
  accu <- array(0,dim=c(nrow(lst),(ncol(lst)+3)))
  max.accu=0
  comb=c(0,0)
  accu[,1] = lst[,1]
  accu[,2] = lst[,2]
  for (p in 1:pred.clst.num){
    confus.mat[p,1] <- sum(dist.mat[1:cls.label1]==p)/sample.size
    confus.mat[p,2] <- sum(dist.mat[(cls.label1+1):cls.label2]==p)/sample.size
  }
  print(confus.mat) # confusion matrix
  for (a in 1:nrow(lst)){
    accu[a,3] = confus.mat[lst[a,1],1]
    accu[a,4] = confus.mat[lst[a,2],2]
    accu[a,5] = confus.mat[lst[a,1],1]+confus.mat[lst[a,2],2]
    #print(accu[a,])
    if (accu[a,5]>max.accu){
      max.accu = accu[a,5]
      comb = lst[a,]
      # print(a)
      # print(max.accu)
      # print(comb)
    }
  }
  accu.final = max(accu[,5])
  #print(accu)
  
  list(accu.mat=accu, max.accu=accu.final,comb.clst=comb,confus.mat=confus.mat)
}





accu.mat.c3 <- function(dist.mat, cls.label1, cls.label2, cls.label3, sample.size, pred.clst.num){
  if (pred.clst.num==2){
    pred.clst.num=3
  }
  confus.mat <- array(0,dim=c(pred.clst.num,3))
  lst = permutations(n = pred.clst.num, r = 3, repeats.allowed = F) #all combination with no repeats
  accu <- array(0,dim=c(nrow(lst),(ncol(lst)+4)))
  max.accu=0
  comb=c(0,0,0)
  accu[,1] = lst[,1]
  accu[,2] = lst[,2]
  accu[,3] = lst[,3]
  for (p in 1:pred.clst.num){
    confus.mat[p,1] <- sum(dist.mat[1:cls.label1]==p)/sample.size
    confus.mat[p,2] <- sum(dist.mat[(cls.label1+1):cls.label2]==p)/sample.size
    confus.mat[p,3] <- sum(dist.mat[(cls.label2+1):cls.label3]==p)/sample.size
  }
  print(confus.mat) # confusion matrix
  for (a in 1:nrow(lst)){
    accu[a,4] = confus.mat[lst[a,1],1]
    accu[a,5] = confus.mat[lst[a,2],2]
    accu[a,6] = confus.mat[lst[a,3],3]
    accu[a,7] = confus.mat[lst[a,1],1]+confus.mat[lst[a,2],2]+confus.mat[lst[a,3],3]
    #print(accu[a,])
    if (accu[a,7]>max.accu){
      max.accu = accu[a,7]
      comb = lst[a,]
      # print(a)
      # print(max.accu)
      # print(comb)
    }
  }
  accu.final = max(accu[,7])
  #print(accu)
  
  list(accu.mat=accu, max.accu=accu.final,comb.clst=comb,confus.mat=confus.mat)
}




########### Jaccard index caluclation
Jac_cal <- function(clust.ind,sample.size,true.cls){
  if (true.cls==2){
    lst = permutations(n = true.cls, r = 2, repeats.allowed = F)
    part1 <- array(0,dim=c(true.cls,sample.size))
    part1[1,] <- c(rep(lst[1,1],sample.size/true.cls),rep(lst[1,2],sample.size/true.cls))
    part1[2,] <- c(rep(lst[2,1],sample.size/true.cls),rep(lst[2,2],sample.size/true.cls))
    concord = max(length(which(part1[1,]==clust.ind)),length(which(part1[2,]==clust.ind)))
    concord/(sample.size*2-concord)
  }
  else{
    lst = permutations(n = true.cls, r = 3, repeats.allowed = F)
    part1 <- array(0,dim=c(true.cls,sample.size))
    part1[1,] <- c(rep(lst[1,1],sample.size/true.cls),rep(lst[1,2],sample.size/true.cls),rep(lst[1,3],sample.size/true.cls))
    part1[2,] <- c(rep(lst[2,1],sample.size/true.cls),rep(lst[2,2],sample.size/true.cls),rep(lst[2,3],sample.size/true.cls))
    part1[3,] <- c(rep(lst[3,1],sample.size/true.cls),rep(lst[3,2],sample.size/true.cls),rep(lst[3,3],sample.size/true.cls))
    concord = max(length(which(part1[1,]==clust.ind)),length(which(part1[2,]==clust.ind)),length(which(part1[3,]==clust.ind)))
    concord/(sample.size*2-concord)
    
  }

}
