##### calculate PD distance matrix #####

n.dir.rep <- 300 #number of boostrap replicates


library(nloptr)
library(tools)
###paralell tools
library(snow); library(iterators); library(foreach); library(doSNOW)
library(parallel)
library(cluster)


source(paste0(path, 'updated_code/scinet_part3_ex_fns_v1.r'))
source(paste0(path, 'updated_code/scinet_part2_new_function_v1.r'))
source(paste0(path, 'updated_code/scinet_estim_functions_v1.r'))

n.otu <- 280

##create gen.data.array to be the same as simulation
cls.sample <- as.numeric(as.factor(PD[, n.otu+1]))
table(cls.sample)

gen.data.array1 <-PD[which(cls.sample == 2),]
gen.data.array <- apply(gen.data.array1[, c(1:n.otu)], 2, as.numeric)
dim(gen.data.array)

ncount.list <- list()
ncount.list[[1]] <- table(cls.sample)[1]
ncount.list[[2]] <- table(cls.sample)[2]

## normalized the count data;
rds.min <- 10000; rds.max <- 20000
rds.sample <- lapply(as.list(1:length(ncount.list)), function(x){as.integer(runif(ncount.list[[x]], rds.min, rds.max))})
rds.sample.norm <- lapply(rds.sample, function(x){x/mean(x)})
cls.names <- 1:length(ncount.list)
gen.data.array.norm <- matrix(0, nrow=dim(gen.data.array)[1], ncol=dim(gen.data.array)[2])
gen.data.array.norm <- apply(gen.data.array, 2, function(x){x/rds.sample[[2]]})

gen.data.trans <- apply(gen.data.array.norm, 2, function(x) log(x+1))

###############################################################
#Define the parameters for the mixture distribution for the fit
###############################################################
lc.alp.cut <- 8
alp.max.q <- 0.85

subdiv.list <- c(0.05, 0.10, 0.30, 0.55) #percentage of data where overlap is in the ovlp list
ovlp.lvl.list <- c(0.65, 0.50, 0.35, 0.20) #cutoff for overlap to skip a distribution

## Define the model type to run ##
#Default - overspecify the model
alp.l <- c(1, 1, seq(2, lc.alp.cut-1)) #low count dist. parameters
bet.l <- c(2, 1, rep(1, lc.alp.cut-2))
ct.brk <- 20 #number of bins for the data

dist.opt <- 3 #the optimal distance to use - l2-discrete CDF
tot.mod.dist <- 3 #total distances to calculate

#precalculate the conditional matrices for all mixtures (select later for each dataset)
alg.indx <- 1 #algorithm to use (BFGS)
sig.level.round <- 2 #round the scaling factor for the conditional matrix
tmat.vals <- unique(round(seq(0.01, 8.6, by=0.01), sig.level.round))


init.param.list <- list(lc=lc.alp.cut, alpmax=alp.max.q, subdiv=subdiv.list, ovlp=ovlp.lvl.list, a=alp.l, b=bet.l, ct=ct.brk, dist=dist.opt, tdist=tot.mod.dist, alg=alg.indx, sig=sig.level.round, tmat=tmat.vals, ndir=n.dir.rep)

###################

n.obs <- dim(gen.data.array)[1]
n.gen <- dim(gen.data.array)[2]
## create training and test
test.set.idx <- 1:n.obs
train.set.idx <- 1:n.obs

for (jj in 1:n.otu) {
  print(c('Iteration:', jj))
  dat.iter <- cbind(gen.data.array[,jj], unlist(rds.sample.norm[[2]]))
  data.train <- dat.iter ## the training set only
  dat.run.dist <- dat.iter
  source(paste0(path, 'updated_code/scinet_part3_implementation_v1_ex_meth_run.r'))
}

####################
#Manhattan Distance:
t.t1 <- Sys.time()
dist.calc.pw.mat <- array(0, dim=c(dim(gen.data.array.norm)[1], dim(gen.data.array.norm)[1]))
mx.k <- ceiling(exp(log(dim(gen.data.array.norm)[2]) + log(dim(matrix.rc.list)[1]) - log(7e+7))) #about 1Gb blocks
iter.blk <- ceiling(dim(matrix.rc.list)[1]/mx.k)
for (k.indx in 1:mx.k) {
  if (k.indx == mx.k) {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- dim(matrix.rc.list)[1]; tot.idx <- idx.en - idx.st + 1
  } else {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- k.indx*iter.blk; tot.idx <- idx.en - idx.st + 1
  }
  comp.array <- array(0, dim=c(tot.idx, dim(gen.data.array.norm)[2], 2))
  comp.array[,,1] <- gen.data.array.norm[matrix.rc.list[idx.st:idx.en,1],]
  comp.array[,,2] <- gen.data.array.norm[matrix.rc.list[idx.st:idx.en,2],]
  dist.calc.pw.mat[matrix.rc.list[idx.st:idx.en,]] <- apply(abs(comp.array[,,1] - comp.array[,,2]), 1, sum)
}
dist.mat.array.mht <- dist.calc.pw.mat
print(Sys.time() - t.t1)

####################
#Euclidean Distance:
t.t1 <- Sys.time()
dist.calc.pw.mat <- array(0, dim=c(dim(gen.data.array.norm)[1], dim(gen.data.array.norm)[1]))
mx.k <- ceiling(exp(log(dim(gen.data.array.norm)[2]) + log(dim(matrix.rc.list)[1]) - log(7e+7))) #about 1Gb blocks
iter.blk <- ceiling(dim(matrix.rc.list)[1]/mx.k)
for (k.indx in 1:mx.k) {
  if (k.indx == mx.k) {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- dim(matrix.rc.list)[1]; tot.idx <- idx.en - idx.st + 1
  } else {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- k.indx*iter.blk; tot.idx <- idx.en - idx.st + 1
  }
  comp.array <- array(0, dim=c(tot.idx, dim(gen.data.array.norm)[2], 2))
  comp.array[,,1] <- log(gen.data.array.norm[matrix.rc.list[idx.st:idx.en,1],]+1)
  comp.array[,,2] <- log(gen.data.array.norm[matrix.rc.list[idx.st:idx.en,2],]+1)
  dist.calc.pw.mat[matrix.rc.list[idx.st:idx.en,]] <- sqrt(apply((comp.array[,,1] - comp.array[,,2])^2, 1, sum))
}
dist.mat.array.eucl <- dist.calc.pw.mat
print(Sys.time() - t.t1)

####################
#Bray-Curtis Distance:
t.t1 <- Sys.time()
dist.calc.pw.mat <- array(0, dim=c(dim(gen.data.array.norm)[1], dim(gen.data.array.norm)[1]))
mx.k <- ceiling(exp(log(dim(gen.data.array.norm)[2]) + log(dim(matrix.rc.list)[1]) - log(7e+7))) #about 1Gb blocks
iter.blk <- ceiling(dim(matrix.rc.list)[1]/mx.k)
for (k.indx in 1:mx.k) {
  if (k.indx == mx.k) {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- dim(matrix.rc.list)[1]; tot.idx <- idx.en - idx.st + 1
  } else {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- k.indx*iter.blk; tot.idx <- idx.en - idx.st + 1
  }
  comp.array <- array(0, dim=c(tot.idx, dim(gen.data.array.norm)[2], 2))
  comp.array[,,1] <- log(gen.data.array.norm[matrix.rc.list[idx.st:idx.en,1],]+1)
  comp.array[,,2] <- log(gen.data.array.norm[matrix.rc.list[idx.st:idx.en,2],]+1)
  dist.calc.pw.mat[matrix.rc.list[idx.st:idx.en,]] <- apply(abs(comp.array[,,1] - comp.array[,,2]),1, sum)/apply(comp.array[,,1] + comp.array[,,2]+0.0000001,1, sum)
}

dist.mat.array.bc <- dist.calc.pw.mat
print(Sys.time() - t.t1)


## Sum the distances over the different OTU
dist.mat.array.sum <- array(0, dim=c(length(train.set.idx), length(train.set.idx),6)) 

for (jj in 1:n.gen) {
  model.fit.init.fn <- paste0(f.dir, '/model_train_fit_otu_iteration', jj, '.RData')
  load(model.fit.init.fn)
  dist.mat.array.sum[,,1] <- dist.mat.array.sum[,,1] + dist.mat.array[,,1] #discrete pdf
  dist.mat.array.sum[,,2] <- dist.mat.array.sum[,,2] + dist.mat.array[,,2] #discrete cdf
  dist.mat.array.sum[,,3] <- dist.mat.array.sum[,,3] + dist.mat.array[,,3] #cont cdf
  
}


dist.mat.array.sum[,,4] <- dist.mat.array.mht
dist.mat.array.sum[,,5] <- dist.mat.array.eucl
dist.mat.array.sum[,,6] <- dist.mat.array.bc


